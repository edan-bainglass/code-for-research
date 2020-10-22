#!/bin/bash

#SBATCH -J Test_R
#SBATCH -o err
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 16
#SBATCH -t 00:30:00
#SBATCH -A Photo-catalysts

module load vasp/5.4.4
source /work/01293/hudamn/Edan/genFunc

# parameters

wave=1
cycles=1
k1=1
eqVol=0
isif=3
ibrion=2
nsw=20
ediff=4
ediffg=01

# setup

disTag NELM
enbTag 'NELMDL NELMIN'
setTag 'ISTART ICHARG' '0 2'
setTag 'NELMIN NELMDL' '6 -12'
setTag 'ISIF IBRION NSW' "$isif $ibrion $nsw"
setTag EDIFF "1E-0$ediff"
setTag EDIFFG "-0.$ediffg"
setTag LCHARG F
kPoints $k1

if ((ediff == 6)); then setTag NELMIN 1; fi

if ((wave == 1)); then
	setTag LWAVE T
	if [[ -s WAVECAR ]]; then
		setTag ISTART 1
		disTag NELMDL
	fi
else
	setTag LWAVE F
fi

cp POSCAR tempPOS

if ((eqVol == 0)); then
	if [[ $(grep '; R' status) ]]; then
		n=$(grep -o '; R.*' status | tail -1 | cut -d ' ' -f 3 | cut -d '(' -f 1)
		((n++))
	else
		echo "" >>status
		n=1
	fi
else
	scalers=
	s=-$(((cycles - 1) / 2))
	echo -en "\nSearching for equilibrium volume..." >>status
	echo -n '' >EV
fi

# begin relaxation

for i in $(seq $cycles); do

	# start reading WAVECAR once available
	# NELMDL is only necessary for a bad guess with no WAVECAR
	if [[ -s WAVECAR ]]; then
		setTag 'ISTART ICHARG' '1 0'
		disTag NELMDL
		cp WAVECAR tempWave
	fi

	((eqVol == 1)) && scaleLattice $s "$scalers"

	ibrun vasp_std >Vasp.out

	# backup files for current run and recalculate k-point grid
	cp CONTCAR POSCAR
	backUp $(if ((eqVol == 0)); then echo R; else echo V; fi) $i
	kPoints $k1

	# adjust POTIM
	#  if (( ibrion == 2 )); then
	#    x=($(grep -o 'trialstep.*' Vasp.out | tail -1 | getNumber))
	#    s=$(echo "scale=5;${x[0]} * 10 ^ ${x[1]}" | bc)
	#    setTag POTIM $(echo "$(getTag POTIM) * $s" | bc | xargs printf %0.3f)
	#  fi

	# determine convergence and write output
	Istep=$(tail -1 OSZICAR | grep -o '[0-9][0-9]* F' | cut -d ' ' -f 1)
	if ((eqVol == 0)); then
		getEnergy $Istep $n $i $isif $ibrion >>status
		if [[ $(grep 'reached required accuracy' Vasp.out) ]]; then
			sed -i "$ s/$/ -> Done \!\!\!/" status
			if (($(grep E0 status | tail -1 | grep -o 1E... | cut -c 5) >= 6)); then
				realSpace >>status
				rm tempWave
			fi
			break
			#	elif [[ $(grep "can't locate minimum" Vasp.out) ]]; then
			#		setTag IBRION 2
			#		ibrion=2
			# TODO - if delE in meV range (close to minimum)
		fi
	else
		getEnergy $Istep $i $isif $ibrion >>EV
		((s++))
	fi
done

# final backup

if ((eqVol == 0)); then
	cp status output/relax
	mkdir -p output/relax/trial_$n
	mv output/relax/{1,2,3,4,5} output/relax/trial_$n 2>/dev/null
else
	cp EV output/relax/EV_data
	cp tempPOS POSCAR
fi

# reset INCAR

setTag 'LWAVE IBRION ISIF NSW EDIFF' 'F -1 2 0 1E-06'
disTag NELMDL NELMIN
