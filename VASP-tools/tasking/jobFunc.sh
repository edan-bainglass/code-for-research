#!/bin/bash

source /work/01293/hudamn/Edan/genFunc

FC() { # font color
    if [ $3 ]; then echo $(tput setaf $1)$2; else echo -n $(tput setaf $1)$2; fi
    tput sgr0
}

BC() { # background color
    if [ $3 ]; then echo $(tput setab $1)$2; else echo -n $(tput setaf $1)$2; fi
    tput sgr0
}

baderJob() {
    setTag "ISTART ICHARG LAECHG LCHARG" "0 2 T T"
    if [[ -s WAVECAR ]]; then setTag "ISTART ICHARG" "1 0"; fi
    if [[ -s CHGCAR ]]; then setTag ICHARG 1; fi
    echo -e "\ncalcBader > baderLog; rm AECCAR*" >>job
    echo -e "setTag \"LAECHG LCHARG\" \"F F\"" >>job
}

getFiles() { # copy work files
    cp INCAR oldINCAR 2>/dev/null
    cp /work/01293/hudamn/Edan/files/general/* ./
    if [ -f oldINCAR ]; then
        setTag 'NBANDS NPAR' "$(getTag NBANDS oldINCAR) $(getTag NPAR oldINCAR)"
        rm oldINCAR
    fi
    if [ -f POSCAR ]; then
        local func
        echo -ne "\nFunctional ( 0: LDA | 1: GGA | $(FC 2 2): PBE ) -> "
        read func
        SetUp $func
    else
        echo -e "\nNo POSCAR found. Get POSCAR and run SetUp!\n"
    fi
    # set AMIN if any lattice vector is > 50 angstrom
    if [[ $(getLat disp | awk '{if($1 > 50) print "1"}') ]]; then enbTag AMIN; fi
}

editJob() { # edit job template
    if [[ $1 ]] && [ $(inArray $1 "IN R C SP D B V") ]; then
        vi /work/01293/hudamn/Edan/files/jobs/job$1
    else
        vi /work/01293/hudamn/Edan/files/jobs/job
    fi
}

editINCAR() { # edit INCAR template
    vi /work/01293/hudamn/Edan/files/general/INCAR
}

name() { # set system name

    local name x term
    term=0

    if [ ! -f job ]; then
        echo -e "\nNo job file found\n"
        return
    fi

    if [ $1 ]; then x=$1; fi
    set ${x:=$(grep "\-J" job | cut -d '_' -f 2)}
    echo -ne "\nEnter system name -> "
    read name
    if [ ! $name ]; then
        if [ -f prevJob ]; then
            name=$(grep "^#SBATCH -J" prevJob | cut -d ' ' -f 3 | cut -d '_' -f 1)
        else
            echo -e "\nNo previous job found"
            name
            term=1
        fi
    fi
    if ((term == 0)); then
        sed -i "s/^#SBATCH -J.*/#SBATCH -J "$name"_"$x"/" job
        sed -i "1s/.*/$name/" POSCAR
    fi
}

queue() { # set queue
    local q term
    term=0
    if [ $1 ]; then
        q=$1
    else
        echo -n "Select queue ( 1-normal | 2-development ) -> "
        read q
        if [[ ! $q ]]; then
            if [ -f prevJob ]; then
                q=$(grep "^#SBATCH -p" prevJob | cut -d ' ' -f 3)
                if [[ $q == normal ]]; then q=1; else q=2; fi
            else
                echo -e "\nNo previous job found\n"
                queue
                term=1
            fi
        fi
    fi
    if ((q == 1)); then q=normal; else q=development; fi
    if ((term == 0)); then sed -i "s/^#SBATCH -p .*/#SBATCH -p $q/" job; fi
}

tLim() { # set time limit
    local tUnit t job term
    term=0
    if [ $1 ]; then tUnit=$1; else
        echo -n "Enter time unit ( 1-hours | 2-minutes ) -> "
        read tUnit
    fi
    if [ $3 ]; then file=$3; else file=job; fi
    if [ ! $tUnit ]; then
        if [ -f prevJob ]; then
            t=$(grep "^#SBATCH -t" prevJob | cut -d ' ' -f 3)
            sed -i "s/^#SBATCH -t .*/#SBATCH -t $t/" $file
        else
            echo -e "\nNo previous job found\n"
            tLim
            term=1
        fi
    else
        timeLimit() {
            local t
            if [ $1 ]; then t=$1; else
                echo -n "Enter time limit -> "
                read t
            fi
            if [ ! $t ] || [[ ! $t =~ ^[0-9]+$ ]]; then
                echo -e "\nInvalid time limit\n"
                timeLimit $1
            elif [ $tUnit == 1 ]; then
                sed -i "s/^#SBATCH -t .*/#SBATCH -t $(echo $t | xargs printf '%s:00:00')/" $file
            else
                if ((t < 10)); then t="0$t"; fi
                sed -i "s/^#SBATCH -t .*/#SBATCH -t $(echo $t | xargs printf '00:%s:00')/" $file
            fi
        }
        if ((term == 0)); then timeLimit $2; fi
    fi
}

cores() { # set number of cores
    local cores term
    term=0
    if [[ $1 ]]; then cores=$1; else
        echo -ne "Enter number of cores -> "
        read cores
    fi
    if [ ! $cores ]; then
        if [ -f prevJob ]; then
            cores=$(grep "^#SBATCH -n" prevJob | cut -d ' ' -f 3)
        else
            echo -e "\nNo previous job found\n"
            cores
            term=1
        fi
    fi
    if ((term == 0)); then sed -i "s/^\(\#SBATCH -n\).*/\1 $cores/1" job; fi
}

nodes() { # set number of nodes
    local nodes t
    term=0
    if [[ $1 ]]; then nodes=$1; else
        echo -ne "Enter number of nodes -> "
        read nodes
    fi
    if [ ! $nodes ]; then
        if [ -f prevJob ]; then
            nodes=$(grep "^#SBATCH -N" prevJob | cut -d ' ' -f 3)
        else
            echo -e "\nNo previous job found\n"
            nodes
            term=1
        fi
    fi
    if ((term == 0)); then sed -i "s/^\(\#SBATCH -N\).*/\1 $nodes/" job; fi
}

getJob() { # get job template

    cp INCAR oldINCAR 2>/dev/null

    local x nElect

    if [[ $1 ]] && [[ $(inArray $1 "IN R C SP D B") ]]; then x=$1; else x=''; fi
    if [ -f job ]; then cp job prevJob; fi
    cat /work/01293/hudamn/Edan/files/jobs/job$x >job # write job file
    echo -e "\n$(FC 2 '### Leave blank for previous job parameters ###')"
    case $x in
    IN)
        local ediff kpt lwave lcharg
        echo -en "\nEDIFF (1E-0$(FC 2 6)) -> "
        read ediff
        set ${ediff:=6}
        echo -n "K-Point Grid ($(FC 2 1)) -> "
        read kpt
        set ${kpt:=1}
        echo -n "Write WAVECAR? ($(FC 2 F)) -> "
        read lwave
        set ${lwave:=F}
        echo -n "Write CHGCAR? ($(FC 2 F)) -> "
        read lcharg
        set ${lcharg:=F}
        sed -i -e "s/ediff=.*/ediff=$ediff/" -e "s/lwave=.*/lwave=$lwave/" -e "s/lcharg=.*/lcharg=$lcharg/" job
        sed -i "s/kPoints 1/kPoints $kpt/" job
        ;;
    C)
        local tType cycles crit ismear sig step defEn cent val k1 nBands op max even kInc scheme
        echo -en "\nConvergence test ( $(FC 2 K): KPOINTS) | E: ENCUT | N: NBANDS | S: SIGMA ) -> "
        read tType
        set ${tType:=K}
        echo -n "Convergence criteria ($(FC 2 0.001)) -> "
        read crit
        set ${crit:=0.001}
        if [[ ! $tType == K ]]; then
            echo -n "k-point ($(FC 2 1)) -> "
            read k1
            set ${k1:=1}
        fi
        case $tType in
        K)
            echo -n "Allow even grids? ( 1: yes | $(FC 2 0): no ) -> "
            read even
            ((even == 1)) && even=even || even=
            echo -n "K-point increment ($(FC 2 2)) -> "
            read kInc
            set ${kInc:=2}
            sed -i -e "s/^even=.*/even=$even/" -e "s/^kInc=.*/kInc=$kInc/" job
            if [[ $even ]]; then scheme=Gamma; else scheme=MonkP; fi
            sed -i "3s/.*/$scheme/" KPOINTS
            #k1=$(grep E0 status | tail -1 | cut -d ';' -f 3 | getNumber | sort -r | head -1); kPoints $((k1+kInc)) $even
            cycles=10
            ;;
        E)
            defEn=$(defEn)
            cent=$(echo "$defEn / 100 * 100" | bc)
            if [[ $cent == $(echo "$defEn / 1" | bc) ]]; then
                val=$cent
            elif [[ $(echo "($defEn - $cent) <= 25" | bc) == 1 ]]; then
                val=$((cent + 25))
            elif [[ $(echo "($defEn - $cent) > 25 && ($defEn - $cent) <= 50" | bc) == 1 ]]; then
                val=$((cent + 50))
            elif [[ $(echo "($defEn - $cent) > 50 && ($defEn - $cent) <= 75" | bc) == 1 ]]; then
                val=$((cent + 75))
            else
                val=$((cent + 100))
            fi
            setTag ENCUT $val
            kPoints $k1
            sed -i "s/^val=.*/val=$val/" job
            ;;
        N)
            nElect=$(storePos e)
            if ((nElect % 2 != 0)); then ((nElect++)); fi
            step=$(getTag NPAR)
            nBands=$((nElect / 2 + step))
            sed -i -e "s/^nBands=.*/nBands=$nBands/" -e "s/^step=.*/step=$step/" job
            setTag 'NBANDS ICHARG' "$nBands 12"
            ;;
        S)
            kPoints $k1
            echo -n "Smearing type ( -1: fermi | $(FC 2 0): gauss | >0: MP ) -> "
            read ismear
            set ${ismear:=0}
            echo -n "Initial SIGMA ($(FC 2 INCAR))-> "
            read sig
            set ${sig:=$(getTag SIGMA)}
            echo -n "Step size ($(FC 2 0)) -> "
            read step
            set ${step:=0}
            echo -n "( 1: increment | 2: decrement ) -> "
            read op
            if [[ $(echo "$step != 0" | bc) == 1 ]]; then
                if [[ $op == 1 ]]; then
                    echo -n "Max value -> "
                    read max
                    cycles=$(echo "($max - $sig) / $step + 1" | bc)
                else
                    cycles=$(abs $(echo "$sig / $step - 1" | bc) 0) # check if -1 is necessary
                fi
            else
                cycles=1
            fi
            sed -i -e "s/^\(ismear=\).*/\1$ismear/" -e "s/^\(sig=\).*/\1$sig/" -e "s/^\(step=\).*/\1$step/" job
            ;;
        esac
        if [[ ! $tType =~ [SK] ]]; then
            echo -n "Number of iterations ($(FC 2 1))-> "
            read cycles
            set ${cycles:=1}
        fi
        sed -i -e "s/\(tType=\).*/\1$tType/" -e "s/\(cycles=\).*/\1$cycles/" -e "s/\(crit=\).*/\1$crit/" job
        ;;
    R)
        local wave cycles k1 isif ibrion nsw eqVol scale step ediff ediffg
        echo -en "\nWrite wavefunction ( $(FC 2 1): y | 0: n ) -> "
        read wave
        set ${wave:=1}
        echo -n "Number of iterations ($(FC 2 1)) -> "
        read cycles
        set ${cycles:=1}
        echo -n "k-point ($(FC 2 1)) -> "
        read k1
        set ${k1:=1}
        echo -n "EDIFF ($(FC 2 4)) -> "
        read ediff
        set ${ediff:=4}
        echo -n "EDIFFG (-0.$(FC 2 01)) -> "
        read ediffg
        set ${ediffg:=01}
        echo -n "ISIF ($(FC 2 3)) -> "
        read isif
        set ${isif:=3}
        echo -n "IBRION ($(FC 2 2)) -> "
        read ibrion
        set ${ibrion:=2}
        echo -n "NSW ($(FC 2 20)) -> "
        read nsw
        set ${nsw:=20}
        echo -n "Equilibrium volume search ( 1: y | $(FC 2 0): n ) -> "
        read eqVol
        set ${egVol:=0}
        if [[ $eqVol == 1 ]]; then
            echo -n "Initial scale for volume ($(FC 2 1.00)) -> "
            read scale
            set ${scale:=1.00}
            echo -n "Step size ($(FC 2 0.01)) -> "
            read step
            set ${step:=0.01}
            sed -i -e "s/\(eqVol=\).*/\1$eqVol/" -e "s/\(initScale=\).*/\1$scale/" -e "s/\(\tstep=\).*/\1$step/" job
        fi
        sed -i -e "s/^\(wave=\).*/\1$wave/" -e "s/^\(cycles=\).*/\1$cycles/" -e "s/^\(k1=\).*/\1$k1/" job
        sed -i -e "s/^\(isif=\).*/\1$isif/" -e "s/^\(ibrion=\).*/\1$ibrion/" -e "s/^\(nsw=\).*/\1$nsw/" job
        sed -i -e "s/^\(ediff=\).*/\1$ediff/" -e "s/^\(ediffg=\).*/\1$ediffg/" job
        ;;
    SP)
        local l xc o b
        echo -en "\nLocal potential ( 1: y | $(FC 2 0): n ) -> "
        read l
        set ${l:=0}
        if ((l == 1)); then
            l=T
            echo -n "Include XC? ( 1: y | $(FC 2 0): n ) -> "
            read xc
            set ${xc:=0}
            if ((xc == 1)); then xc=LVTOT; else xc=LVHAR; fi
        else
            l=F
        fi
        echo -n "Optical properties ( 1: y | $(FC 2 0): n ) -> "
        read o
        set ${o:=0}
        if ((o == 1)); then o=T; else o=F; fi
        echo -n "Bader charges ( 1: y | $(FC 2 0): n ) -> "
        read b
        set ${b:=0}
        if ((b == 1)); then b=T; else b=F; fi
        sed -i "s/\(setTag $xc\) [TF]/\1 $l/" job
        sed -i "s/\(setTag LOPTICS\) [TF]/\1 $o/" job
        sed -i "s/\(setTag LAECHG\) [TF]/\1 $b/" job
        ;;
    D) ;;
    B) ;;
    V) ;;
    *)
        local w c
        echo -n -e "\nWrite WAVECAR ( 1: y | $(FC 2 0): n ) -> "
        read w
        set ${w:=0}
        ((w == 1)) && w=T || w=F
        echo -n "Write CHGCAR ( 1: y | $(FC 2 0): n ) -> "
        read c
        set ${c:=0}
        ((c == 1)) && c=T || c=F
        sed -i -e "s/\(setTag LWAVE\) .*/\1 $w/" -e "s/\(setTag LCHARG\) .*/\1 $c/" job
        ;;
    esac
    name $x
    queue
    tLim
    nBands
    nPar
    kPar
    cores
    nodes

    rm oldINCAR

    echo ""
}

nBands() { # set number of bands
    local bands term
    term=0
    echo -n "Enter number of bands ($(storePos e) electrons; $(storePos i) ions) -> "
    read bands
    if [[ ! $bands ]]; then
        if [ -f oldINCAR ]; then
            bands=$(getTag NBANDS oldINCAR)
        else
            echo -e "\nPrevious INCAR not found\n"
            nBands
            term=1
        fi
    fi
    if ((term == 0)); then setTag NBANDS $bands; fi
}

nPar() { # set parallelization over bands
    local npar term
    term=0
    echo -ne "Enter number of bands to run in parallel (NPAR) -> "
    read npar
    if [[ ! $npar ]]; then
        if [ -f oldINCAR ]; then
            npar=$(getTag NPAR oldINCAR)
        else
            echo -e "\nPrevious INCAR not found\n"
            nPar
            term=1
        fi
    fi
    if ((term == 0)); then setTag NPAR $npar; fi
}

kPar() {
    local kpar term
    term=0
    if [[ $1 ]]; then
        kpar=$1
    else
        echo -n "Enter number of k-points to run in parallel (KPAR) -> "
        read kpar
        if [[ ! $kpar ]]; then
            if [ -f oldINCAR ]; then
                kpar=$(getTag KPAR oldINCAR)
            else
                echo -e "\nPrevious INCAR not found\n"
                kPar
                term=1
            fi
        fi
    fi
    if ((term == 0)); then setTag KPAR $kpar; fi
}

hybrid() { # turn hybrid on/off
    local hybrid hf
    echo -n "Hybrid funcational ( 1: On | $(FC 2 0): Off ) -> "
    read hybrid
    if [[ $hybrid != 1 ]]; then
        setTag LHFCALC F
        disTag HFSCREEN TIME PRECFOCK NKRED AEXX
    else
        echo -n "Enter HF enxchange contribution ($(FC 2 0.25)) -> "
        read hf
        set ${hf:=0.25}
        if [ ! -f status ]; then setTag LWAVE T; fi
        setTag 'LHFCALC AEXX' "T $hf"
        enbTag HFSCREEN TIME PRECFOCK NKRED AEXX
    fi
}

uParam() { # turn DFT+U on/off

    # TO-DO: Allow for simultaneous multiple U parameters?

    local c s i u site orbit ul uu uj l
    echo -n "DFT+U ( 1: On | $(FC 2 0): Off ) -> "
    read c
    if [[ $c != 1 ]]; then
        setTag LDAU F
        disTag LDAUTYPE LDAUL LDAUU LDAUJ LMAXMIX LASPH
        sed -i "s/DFT+U:.*/DFT+U:/" INCAR
    else
        s=($(storePos s))
        for i in ${!s[*]}; do s[$i]="$i-${s[$i]}"; done
        echo -n "Effective U -> "
        read u
        echo -n "On which site? ( $(echo "${s[*]}" | sed 's/ / | /g') ) -> "
        read site
        echo -n "On which orbital? ( 0-s | 1-p | 2-d | 3-f ) -> "
        read orbit
        setTag 'LDAU LASPH' 'T T'
        enbTag LDAUTYPE LDAUL LDAUU LDAUJ LMAXMIX LASPH

        ul=()
        uu=()
        uj=()
        for i in ${!s[*]}; do
            if [[ $(echo ${s[$i]} | cut -c 1) == $site ]]; then
                ul=(${ul[*]} $orbit)
                uu=(${uu[*]} $(printf %0.1f $(echo "$u + 1" | bc)))
                uj=(${uj[*]} 1)
            else
                ul=(${ul[*]} -1)
                uu=(${uu[*]} 0)
                uj=(${uj[*]} 0)
            fi
        done

        setTag LDAUL "${ul[*]}"
        setTag LDAUU "${uu[*]}"
        setTag LDAUJ "${uj[*]}"
        case $orbit in
        2) l=4 ;;
        3) l=6 ;;
        *)
            l=2
            setTag LASPH F
            ;;
        esac
        setTag LMAXMIX $l
        sed -i "s/DFT+U:.*/DFT+U: $u/" INCAR
    fi
}

VDW() { # turn Van Dar Waals forces on/off
    local vdw t
    echo -n "Van Dar Waals ( 1: On | $(FC 2 0): Off ) -> "
    read vdw
    if [[ $vdw != 1 ]]; then
        setTag IVDW 0
    else
        echo -n "Type ( 10: D2 | 11: D3 | 12: D3+ ) -> "
        read t
        setTag IVDW $t
    fi
}

symmetry() { # turn symmetry on/off
    local s
    echo -n "Symmetry ( 2: On | $(FC 2 0): Off ) -> "
    read s
    set ${s:=0}
    if ((s != 0)) && [[ $(getTag LHFCALC) == T ]]; then s=3; fi
    setTag ISYM $s
}

spinPolar() {
    local sp
    echo -n "Spin-polarized ( $(FC 2 1): On | 0: Off ) -> "
    read sp
    set ${sp:=1}
    setTag ISPIN $((sp + 1))
}

realProj() {
    local p
    echo -n "Projections in real space ( $(FC 2 1): On | 0: Off ) -> "
    read p
    set ${p:=1}
    if [[ $p == 1 ]]; then p=T; else p=F; fi
    setTag LREAL $p
}

magnetic() {
    local m
    echo -n "Magnetic moment expected? ( 0: Yes | $(FC 2 1): No ) -> "
    read m
    set ${m:=1}
    setTag NUPDOWN $((m - 1))
}

dipoleCorrect() {
    local d
    echo -n "Dipole corrections? ( 1: Yes | $(FC 2 0): No ) -> "
    read d
    set ${d:=0}
    if ((d == 1)); then
        setTag LDIPOL T
        setTag IDIPOL 3
    else
        setTag LDIPOL F
        setTag IDIPOL 0
    fi
}

SetUp() { # define initial job settings
    local s f
    s=$(head -6 POSCAR | tail -1 | rmSpace)
    if [ $1 ]; then f=$1; else f=2; fi
    potGen $f $s
    if [ -s POTCAR ]; then
        hybrid
        uParam
        VDW
        dipoleCorrect
        spinPolar
        symmetry
        realProj
        magnetic
        echo ""
    fi
}

potGen() { # generate POTCAR

    local func dir opts

    gen() {
        local elements i
        if [[ $1 ]]; then
            elements=$*
        else
            echo -ne "\nEnter system elements -> "
            read elements
        fi
        for i in $elements; do
            if [[ ! $(grep "^ $pat $i[_ ]" $dir) ]]; then
                echo -e "\nOne or more of the elements does not exist!"
                gen $*
            else
                sed -n "/^ $pat $i[_ ]/,/End/p" $dir >>POTCAR
            fi
        done
    }

    # define location of POTCAR database (LDA | GGA | PBE)
    if [ -f POTCAR ]; then rm POTCAR; fi
    if [[ ! $1 ]]; then
        echo -ne "\nFunctional ( 0: LDA | 1: GGA | $(FC 2 2): PBE ) -> "
        read func
        set ${func:=2}
    else func=$1; fi
    case $func in
    0)
        dir=/work/01293/hudamn/Edan/PAW_LDA
        pat=PAW
        ;;
    1)
        dir=/work/01293/hudamn/Edan/PAW_GGA
        pat=PAW
        ;;
    2)
        dir=/work/01293/hudamn/Edan/PAW_PBE
        pat=PAW_PBE
        ;;
    *)
        echo -e "\nUsage: potGen <0: LDA | 1: GGA | 2: PBE> <element list>\n"
        return
        ;;
    esac

    # generate POTCAR
    opts=($*)
    gen ${opts[*]:1}

    # verify POTCAR
    echo ""
    grep "^ $pat [A-Z]" POTCAR
    echo ""
}
