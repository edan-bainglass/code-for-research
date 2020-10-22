#!/bin/bash

#SBATCH -J Test_IN
#SBATCH -o err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:30:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
#module load remora
source /work/01293/hudamn/Edan/genFunc

# parameters

lwave=F
lcharg=F
ediff=6

# setup

kPoints 1
sys="$(head -1 POSCAR | cut -d ' ' -f 1 | sed 's/ $//')"
setTag 'System ISTART ICHARG EDIFF LWAVE LCHARG' "$sys 0 2 1E-0$ediff $lwave $lcharg"
if [[ $(getTag LDAU) == T ]]; then u=$(grep DFT+U INCAR | getNumber); else u=""; fi
if [[ $u != "" ]]; then u="(Ueff = $u eV)"; fi
echo $sys $u >status
echo -e "\nInitial Values" >>status
realSpace >>status

#  echo -n "E0; EDIFF; K; ENCUT; NBANDS; NELM (ALGO); ISMEAR <(SIGMA)>; PREC; " >> status
#  echo -e "<R#ISIF#IBRION Trial(Cycle): NSW (EDIFFG)>; <V>; <(days)> h:m:s\n" >> status

#remora ibrun vasp_std > Vasp.out
ibrun vasp_std >Vasp.out

# output

getEnergy >>status
backUp IN $(getTag NPAR) $(getTag NBANDS)
