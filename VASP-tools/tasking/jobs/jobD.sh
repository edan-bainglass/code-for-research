#!/bin/bash

#SBATCH -J Test_D
#SBATCH -o err
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 80
#SBATCH -t 02:00:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
source /work/01293/hudamn/Edan/genFunc

# CALCULATE DOS #####################################################################################################
setTag ALGO normal
ismear=$(getTag ISMEAR)
setTag 'LCHARG ISTART ICHARG LORBIT ISMEAR' 'F 0 11 10 -5'
kTemp=$(getKps | getNumber | sort -r | head -1)
#  kPoints $(printf %0.f $(echo "$kTemp * 1.5" | bc))
kSca="$(head -4 KPOINTS | tail -1 | rmSpace | sed 's/ /x/g')"
echo -en "\nCalculating DOS... " >>status
ibrun vasp_std >Vasp.out
t=$(getTime)
kPoints $kTemp
backUp D $kSca
echo -n "generated on a ${kSca//x/ x } grid" >>status
if [[ $(occCheck 1) == '' ]]; then
  vbm=$(vbm i)
  echo "; Band Gap = $(getGap $vbm $((vbm + 1)) output/dos/$kSca); $t" >>status
else
  echo "; $t" '-> Check occupancies !!! ' >>status
fi
setTag 'ICHARG LORBIT ISMEAR' "1 0 $ismear"
#####################################################################################################################
