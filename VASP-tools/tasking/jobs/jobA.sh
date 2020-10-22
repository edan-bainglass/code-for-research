#!/bin/bash

#SBATCH -J Test_A
#SBATCH -o out.err
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 01:00:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
source /work/01293/hudamn/Edan/genFunc

# set bands
a=$(grep "^#SBATCH -J" jobA | cut -d ' ' -f 3)
sed -i "s/[A-Z][a-z]$/$a/" POSCAR
storePos
setTag NBANDS $nElect

#####################################################################################################################
# update INCAR
setTag PREC accurate
setTag ISTART 0
setTag ICHARG 2
setTag NSW 0
setTag IBRION -1
setTag ISIF 2
setTag ISMEAR 2
setTag LWAVE F
setTag ISYM 1
# set KPOINTS
kPoints 1
# write header in status file
echo -e "Atomic Energies\n" >status
showTag ISMEAR >>status
echo "KPOINTS = $(getKps)" >>status
showTag ENCUT >>status
showTag NBANDS >>status
# loop through multiplicities
for i in 1 3 5; do
  setTag NUPDOWN $i
  mkdir -p output/$i
  echo -e "\n2S + 1 = $((i + 1))\n" >>status
  # execute VASP
  ibrun vasp_std >Vasp.out
  # get elapsed time
  tot=$(grep Elapsed OUTCAR | rev | cut -d ' ' -f 1 | rev)
  t=$(getTime $tot)
  # record data
  getEnergy >>status
  # save output
  mv CHG* EIGENVAL OUTCAR OSZICAR Vasp.out vasprun.xml output/$i
done
#####################################################################################################################
