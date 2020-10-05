#!/bin/bash

#SBATCH -J Test_B
#SBATCH -o err
#SBATCH -p development
#SBATCH -N 4
#SBATCH -n 80
#SBATCH -t 02:00:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
source /work/01293/hudamn/Edan/genFunc

# CALCULATE BAND STRUCTURE ###########################################################################################################
	setTag ALGO normal
  setTag 'LCHARG ISTART ICHARG LORBIT' 'F 0 11 10'; enbTag ICORELEVEL
	echo -en "\nCalculating bands... " >> status
  ibrun vasp_std > Vasp.out
	t=$(getTime)
  x=$(grep Line_mode KPOINTS | wc -l)
  if (( x == 0 )) && [ ! -f status ]; then
		echo "Generated; $t" > status
  elif (( x == 1 )); then
		path="$(head -1 KPOINTS | sed 's/\\Gamma/Gm/g')"
    echo "generated along $path; $t" >> status
  elif (( x == 2 )); then
		path=full
		sed -i '1s/.*/Full_Spectrum/' KPOINTS
    echo "full band spectrum generated; $t" >> status
  fi
	if [[ $(occCheck 1) == '*' ]]; then sed -i '$ s/$/ -> Check occupancies !!!/' status; fi
  backUp B bands/$path; cp SCF_KPOINTS KPOINTS
  setTag 'LORBIT ICHARG' '0 1'; disTag ICORELEVEL
######################################################################################################################################
