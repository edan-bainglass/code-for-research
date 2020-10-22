#!/bin/bash

#SBATCH -J Test_C
#SBATCH -o err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:10:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
source /work/01293/hudamn/Edan/genFunc

# K - KPOINTS ; E - ENCUT ; N - NBANDS ; S - SIGMA
tType=K
cycles=10
crit=0.001
# Test specific parameters
val=1
nBands=1
ismear=1
sig=1
step=1
even=
kInc=2

# CONVERGENCE TESTING ###############################################################################################
nElect=$(storePos e)
if ((nElect % 2 != 0)); then ((nElect++)); fi
if [[ $tType == S ]]; then setTag 'ISMEAR SIGMA' "$ismear $sig"; fi
if [[ $tType == K ]]; then
  k1=$(($(getKps | getNumber | sort -nr | head -1) + kInc))
  kPoints $k1 $even
fi
for i in $(seq $cycles); do
  ibrun vasp_std >Vasp.out
  getEnergy >>status
  cur=$(grep E0 status | tail -1 | cut -d ' ' -f 3 | getNumber)
  pre=$(grep E0 status | tail -2 | head -1 | cut -d ' ' -f 3 | getNumber)
  dE=$(abs $cur $pre | xargs printf %0.5f)
  if [[ $(echo "$dE <= $crit" | bc) == 1 ]]; then
    if [[ $tType == S ]]; then
      if [[ $(occCheck 1) == '' ]]; then
        e=$(grep 'energy  without entropy' OUTCAR | rmSpace | cut -d ' ' -f 4,7)
        i=$(storePos i)
        c=$(
          python3.7 - "$e" $i <<-END
							from sys import argv
							e = [float(i) for i in argv[1].split()]
							print(1) if abs(e[0] - e[1]) / float(argv[-1]) > 0.001 else print(0)
					END
        )
        if [[ $c == 1 ]]; then
          sed -i "$ s/$/ -> check entropy.../" status
          backUp C $ismear $(echo $sig | sed 's/^\./0\./') $(getKps | sed 's/ /-/g')
          break
        fi
      else
        sed -i "$ s/$/ -> check occupancies.../" status
        backUp C $ismear $(echo $sig | sed 's/^\./0\./') $(getKps | sed 's/ /-/g')
        break
      fi
    else
      sed -i "$ s/$/ -> dE = $dE/" status
      break
    fi
  fi
  case $tType in
  K)
    k1=$((k1 + kInc))
    kPoints $k1 $even
    ;;
  E)
    val=$((val + 25))
    setTag ENCUT $val
    ;;
  N)
    nBands=$((nBands + step))
    setTag NBANDS $nBands
    ;;
  S)
    backUp C $ismear $(echo $sig | sed 's/^\./0\./') $(getKps | sed 's/ /-/g')
    sig=$(echo $sig + $step | bc)
    setTag SIGMA $(echo $sig | sed 's/^\./0\./')
    ;;
  esac
done
if [[ $tType == S ]]; then cp status output/sigma/$ismear; fi
#####################################################################################################################
