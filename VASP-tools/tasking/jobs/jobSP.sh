#!/bin/bash

#SBATCH -J Test_SP
#SBATCH -o err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
source /work/01293/hudamn/Edan/genFunc

# INCAR setup

setTag LVHAR F
setTag LVTOT F
setTag LAECHG F
setTag LOPTICS F
enbTag ICORELEVEL

if [[ $(getTag LAECHG) == T ]]; then setTag LCHARG T; fi
if [[ -s CHGCAR ]]; then setTag ICHARG 1; elif [[ -s WAVECAR ]]; then setTag ICHARG 0; else setTag ICHARG 2; fi

if [[ $(getTag LOPTICS) == T ]]; then
  setTag NEDOS 2000
  setTag ISMEAR -5
  oldBands=$(getTag NBANDS)
  setTag NBANDS $(echo "$oldBands * 3" | bc)
  folder=optics
elif [[ $(getTag LHFCALC) == T ]]; then
  folder=hybrid/$(echo "$(getTag AEXX) * 100 / 1" | bc)
elif [[ $(getTag LVHAR) == T ]]; then
  folder=lv/har
elif [[ $(getTag LVTOT) == T ]]; then
  folder=lv/tot
else
  folder=totE
fi

ibrun vasp_std >Vasp.out

# output

getEnergy >>status
backUp SP $folder
cp KPOINTS SCF_KPOINTS

if [[ $(getTag LAECHG) == T ]]; then
  calcBader >baderLog
  rm AECCAR*
fi
if [[ $(occCheck 1) == '' ]]; then msg='Final energy value !!!'; else msg='Check occupancies !!!'; fi
sed -i "$ s/$/ -> $msg/" status
if [[ $(getTag LOPTICS) == T ]]; then setTag NBANDS $oldBands; fi
mv LOCPOT output/sp/$folder 2>/dev/null

# reset INCAR

setTag 'LWAVE LCHARG LAECHG LVTOT LVHAR LOPTICS NEDOS ISTART ICHARG ISMEAR' 'F F F F F F 1000 0 2 0'
disTag ICORELEVEL
