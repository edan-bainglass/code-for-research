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

setTag LWAVE F
setTag LCHARG F

ibrun vasp_std > Vasp.out

setTag 'LWAVE LCHARG' 'F F'
