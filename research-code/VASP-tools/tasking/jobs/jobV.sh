#!/bin/bash

#SBATCH -J Test_V
#SBATCH -o out.err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:30:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
source /work/01293/hudamn/Edan/genFunc

wave=1
cycles=1
k1=3

# RELAX CELL ########################################################################################################
  setTag ISTART 0
  if [[ $wave == 1 ]]; then setTag LWAVE T; if [[ -s WAVECAR ]]; then setTag ISTART 1; fi; fi
  setTag 'NSW IBRION ISIF' '20 2 4'; kPoints $k1
  echo -e "\nEquilibrium Volume Search\n" >> status
  showTag IBRION >> status
  showTag ISIF >> status
  echo -n "" > EV
  scale=1.01
  step=0.001
  i=1
  while (( i <= cycles )); do
    if [ $i == 2 ]; then if [[ $wave == 1]]; then setTag ISTART 1; fi; fi
    sed -i "2s/.*/$scale/" POSCAR
    ibrun vasp_std > Vasp.out
    cp CONTCAR POSCAR; backUp V $i 
    tot=$(grep Elapsed OUTCAR | rev | cut -d ' ' -f 1 | rev); t=$(getTime $tot)
    Istep=$(tail -1 OSZICAR | grep -o '[0-9][0-9]* F' | cut -d ' ' -f 1)
    getEnergy $Istep $i $scale >> EV
    kPoints $k1 
    scale=$(echo "$scale - $step" | bc)
    ((i++)) 
  done
  echo -e "\n -> E/V data written to EV" >> status
  setTag 'LWAVE IBRION ISIF NSW' 'F -1 2 0'
#####################################################################################################################
