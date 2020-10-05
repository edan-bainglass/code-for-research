interfaceSeparationMod() {

    python3.7 - "$*" <<-END

    import numpy as np
    from sys import argv

    h, d = [float(i) for i in argv[1].split()]

    with open('POSCAR', 'r') as f, open('newPOS', 'w') as g:

        for _ in range(2):
            g.write(f.readline().lstrip())

        M = np.matrix([[float(col) for col in row.split()]
                    for row in [f.readline() for i in range(3)]])
        N = M.copy()

        N[2, 2] += 2 * d

        g.write('\n'.join('\t' + '\t'.join('%21.16f' % col for col in row)
                        for row in N.A) + '\n')

        for _ in range(4):
            g.write(f.readline().lstrip())

        h *= M[2, 2]
        for line in f.readlines():
            args = line.split()
            if len(args) > 0:
                z = float(args[2]) * M[2, 2]
                args[2] = str(round((z + d if z > h else z) / N[2, 2], 16))
                g.write('\t'.join(args) + '\n')
            else:
                break

END

}
scaleLattice() {

    python3.7 - $* <<-END

    from sys import argv
    import numpy as np

    sca = np.array([float(i) for i in argv[2:]])
    s = int(argv[1])

    with open('tempPOS', 'r') as f, open('scaPOS', 'w') as g:

        for _ in range(2):
            g.write(f.readline())

        M = np.matrix([[float(col) for col in row.split()]
                    for row in [f.readline() for i in range(3)]])
        T = np.matrix(np.identity(3) * (1 + s * sca))

        g.write('\n'.join('\t'.join('%0.9f' % col for col in row)
                        for row in (T.T * M).A) + '\n')

        for line in f.readlines():
            g.write(line)

END

    mv scaPOS POSCAR

}

rotateCell() {

    cp POSCAR origPOS

    python3.7 - "$*" <<-END

    from sys import argv
    import numpy as np

    vals = [float(i) for i in argv[1].split()]

    with open('POSCAR', 'r') as f, open('rotPOS', 'w') as g:

        for _ in range(2):
            g.write(f.readline())

        M = np.matrix([[float(col) for col in row.split()]
                    for row in [f.readline() for i in range(3)]])
        R = np.matrix([vals[:3], vals[3:6], vals[6:]])
        Mp = (M * R).T

        g.write('\n'.join('\t'.join('%0.9f' % col for col in row)
                        for row in Mp.A) + '\n')

        for _ in range(3):
            g.write(f.readline())

        newCoords = np.matrix([[float(col) for col in row.split()[:-1]]
                            for row in f.readlines()]) * R

        g.write('\n'.join('\t'.join('%0.9f' % (col + 1 if col < 0 else col)
                                    for col in row) for row in newCoords.A))

END

    mv rotPOS POSCAR
}

chgDiff() {
    /work/01293/hudamn/Edan/files/extra/chgdiff.py CHGCAR $*
}

sigCheck() {
    grep sigma OUTCAR | grep -o =.* | awk -v n=$(storePos i) '{print ($2 - $5)/n}'
}

compE() {
    local i
    for i in */status status; do
        echo -n "$i -> "
        grep E0 $i | tail -1 | getNumber | head -2 | tail -1
    done
}

getE() {
    local loc i
    if [[ $1 ]]; then loc=$*; else loc=.; fi
    for i in $loc; do grep E0 $i/status | tail -1 | getNumber | head -2 | tail -1; done | sort -nr
}

getAngle() {
    local mags mat a b c
    mags=($(getLat disp))
    mat=($(getLat pass))
    a=${mat[*]:0:3}
    b=${mat[*]:3:3}
    c=${mat[*]:6:3}

    case $1 in
    ab) trig acos $(echo "scale=16; $(dot "$a" "$b") / ${mags[0]} / ${mags[1]}" | bc) ;;
    cz) trig acos $(echo "scale=16; $(dot "$c" "0 0 1") / ${mags[2]}" | bc) ;;
    *) echo -e "\nmissing argument (ab or cz)\n" ;;
    esac
}

trig() {
    awk -v CONVFMT=%.16f -v t=$1 -v x=$2 'BEGIN \
  { \
    if (t == "acos") printf atan2((1.-x^2)^0.5,x); \
  }'
}

dot() {
    local a b
    a=($1)
    b=($2)
    echo "${a[0]} * ${b[0]} + ${a[1]} * ${b[1]} + ${a[2]} * ${b[2]}" | bc
}

originalSlab() {
    if [[ $1 ]]; then
        mkdir -p $(dirname $(find $1/output/init/ -name POSCAR | head -1))/unrelaxed && cp POSCAR $_
    else
        echo -e "\nmissing directory\n"
    fi
}

slabAnalysis() {

    local inDir currDir lat s n t a i j k z w st x o exp c

    if [[ ! -d output ]]; then
        echo -e "\nCheck location\n"
        return
    fi #verify current location

    # locate input POSCAR
    #if [[ $(find output/init -name POSCAR | wc -l) > 1 ]]; then echo -e "\nMultiple input POSCAR available\n"; return; fi
    inDir=$(dirname $(find output/init -name POSCAR | tail -1))

    # get coordinates from input and current POSCARs
    currDir=$(pwd)
    cd $inDir
    getCoords 1 >$currDir/inPOS
    cd $currDir
    getCoords 1 >currPOS

    # write atom list to file
    s=($(storePos s))
    n=($(storePos n))
    v=($(storePos v))
    echo -n "Formula unit (${s[*]}) -> "
    read -a a
    #a=(1 1 2 8) # formula unit
    echo -n "Sort by ( $(FC 2 0): layer | 1: element ) -> "
    read st
    set ${st:=0}
    echo -n "Limit sites numbers to unit cell? ( $(FC 2 0): no | 1: yes ) -> "
    read o
    set ${o:=0}
    echo -n "" >atomList
    w=0
    for i in ${!s[*]}; do
        t=$(((n[i] / a[i] - 1) / 2 + 1))
        x=1
        for j in $(seq 1 $(echo $((n[i] / a[i])))); do
            if ((j == t)); then z=BK; elif ((j > t)); then z=$((-j + 2 * t)); else z=-$j; fi
            for k in $(seq ${a[$i]}); do
                echo -n ${s[$i]} >>atomList
                #TODO (( o != 0 )) && ( (( a[i] > 1 )) && echo -n $(((x-1)%a[i]+1)) ) || echo -n $x >> atomList
                ((a[i] > 1)) && echo -n $(if ((o == 0)); then echo $x; else echo $(((x - 1) % a[i] + 1)); fi) >>atomList
                echo " $z $((j + w))" >>atomList
                ((x++))
            done
        done
        if ((st == 1)); then w=$((w + n[i])); fi
    done

    # write output
    title() {
        if [[ ! -f ACF.dat ]]; then
            echo " -----------------------------------------------------------"
            echo "     z         dx      dy      dz          dr        el   # "
            echo " -----------------------------------------------------------"
        else
            echo " ------------------------------------------------------------------------- "
            echo "     z         dx      dy      dz          dr        el   #         bader  "
            echo " ------------------------------------------------------------------------- "
        fi
    }

    paste inPOS currPOS atomList |
        awk -v z=$(dot "$(getLat pass | cut -d ' ' -f 7-9)" "0 0 1") -v mat="$(getLat pass)" -v t=.5 -v vac=20 \
    '{ \
        split(mat,m); h = $6*z; \
        da = $4-$1; db = $5-$2; dc = $6-$3; \
        da = da > t ? da - 1 : da < -t ? da + 1 : da; \
        db = db > t ? db - 1 : db < -t ? db + 1 : db; \
        dc = dc > t ? dc - 1 : dc < -t ? dc + 1 : dc; \
        h = h > z-vac/2 ? h - z : h; \
        dx = da*m[1]+db*m[4]+dc*m[7]; \
        dy = da*m[2]+db*m[5]+dc*m[8]; \
        dz = da*m[3]+db*m[6]+dc*m[9]; \
        dr = sqrt(dx*dx+dy*dy+dz*dz); \
        printf "%3s %7.3f || %7.3f %7.3f %7.3f  ->  %5.3f  || %4s  %2s\n", $9, h, dx, dy, dz, dr, $7, $8
    }' >slabAnalysis

    if [[ -f ACF.dat ]]; then
        # write valence charges to file
        for i in ${!n[*]}; do
            for j in $(seq 1 ${n[$i]}); do
                echo ${v[$i]}
            done
        done >vTemp

        # write bader charges to file
        sed -n '/--/,/--/p' ACF.dat | clip b 1 | awk '{print $5}' >bTemp

        paste bTemp vTemp | awk '{printf "    ||  %6.3f\n", $2-$1}' >baderAnalysis
        paste slabAnalysis baderAnalysis >sTemp
        mv sTemp slabAnalysis
    fi

    # define block separator
    if ((st == 0)); then
        #TODO get number of layers and atoms per layer automatically in general
        exp=$(
            c=0
            for i in $(seq 1 ${n[0]}); do
                c=$((i * 12))
                echo -n $c\G\;
            done
        )
    else
        exp=$(
            c=0
            for i in ${n[*]}; do
                let c+=$i
                echo -n $c\G\;
            done
        )
    fi

    title
    cat slabAnalysis | sort -nr -k1,1 -k2,2 | sed 's/^  *\w*//' | sed "s/\t//" | tac | eval "sed '$exp'" | clip e 1 #| awk '{print $9, $1, $7, $12}'
    title

    rm inPOS currPOS atomList vTemp bTemp sTemp baderAnalysis slabAnalysis 2>/dev/null
}

diskSpace() {
    /usr/local/etc/taccinfo
}

surfTerm() {

    local s n ox l i j line el num ind nc

    if [[ ! -f POSCAR ]]; then
        echo -e "\nMissing POSCAR\n"
        return
    fi

    s=($(storePos s))
    n=($(storePos n))
    nc=(${n[*]})

    #echo -en "\nExpected oxidation states (${s[*]}) -> "; read -arr ox; echo ""

    ox=(1 3 6 -2)
    echo ""

    echo -n "" >surfPos

    l=9
    for i in ${!s[*]}; do
        for j in $(seq 1 ${n[$i]}); do
            head -$l POSCAR | tail -1 | rmSpace | sed "s/^/${s[$i]} $j /" | cut -d ' ' -f 1,2,5 >>surfPos
            ((l++))
        done
    done

    cat surfPos | sort -rk3 | while read line; do
        el=$(echo $line | cut -d ' ' -f 1)
        num=$(echo $line | cut -d ' ' -f 2)
        ind=$(index $el "${s[*]}")
        echo -n "$el$num "
        if [[ $1 == show ]]; then echo -n "${n[*]} "; fi
        echo "$(math sum $(for i in ${!s[*]}; do echo $((n[i] * ox[i])); done)) / 1" | bc
        n[$ind]=$((n[$ind] - 1))
        if [[ $1 == show ]]; then
            for i in ${n[*]}; do
                if ((i == 0)); then
                    for j in ${!s[*]}; do
                        n[$j]=$((n[j] + nc[j]))
                    done
                    break
                fi
            done
        fi
    done

    echo ""
    rm surfPos

}

getCoords() {
    local n i
    if [ -f POTCAR ]; then
        n=$(storePos i)
    else
        n=0
        for i in $(storePos n); do
            let n+=$i
        done
    fi
    grep -EA $n '((d|D)irect|(c|C)artesian)' POSCAR | tail -n +2 | rmSpace |
        if [[ $1 ]]; then cut -d ' ' -f 1-3; else cut -d ' ' -f 1-; fi
}

scaleVac() {

    local s n i co cn x inc

    s=($(storePos s))
    n=($(storePos n))

    # copy POSCAR header to new POSCAR
    head -2 POSCAR >scaledPos

    if [[ ! $1 ]]; then inc=0; else inc=$1; fi

    python3.7 - $inc <<-END

    import numpy as np
    from sys import argv
    from pymatgen.core import structure, Lattice
    from pymatgen.io.vasp.inputs import Poscar

    vac = float(argv[1])

    slab_in = structure.Structure.from_file('POSCAR')
    lat = slab_in.as_dict().get('lattice')

    lat = Lattice.from_parameters(lat.get('a'),
                                  lat.get('b'),
                                  lat.get('c'),
                                  lat.get('alpha'),
                                  lat.get('beta'),
                                  lat.get('gamma'),
                                  vesta=True)

    lat = Lattice.from_parameters(lat.a,
                                  lat.b,
                                  lat.c + vac / lat.matrix[2, 2] * lat.c,
                                  lat.alpha,
                                  lat.beta,
                                  lat.gamma,
                                  vesta=True)

    with open('scaledPos', 'a') as f:
        f.write(np.array2string(lat.matrix) + "\n")

END

    # format matrix and finish writing header
    sed -i -e 's/[[]//g' -e 's/[]]//g' scaledPos
    head -9 POSCAR | tail -4 >>scaledPos

    # determine scaling ratios in long axis for slab
    co=$(getLat disp | tail -1)
    cn=$(getLat disp scaledPos | tail -1)

    # split POSCAR by element
    # scale z-axis coordinates to fit scaled cell
    x=0
    for i in ${!s[*]}; do
        getCoords | head -$((n[i] + x)) | tail -${n[$i]} |
            awk -v co=$(printf %0.3f $co) -v cn=$(printf %0.3f $cn) \
        -v h=$(echo "scale=16; $co / $cn" | bc) -v o=$1 -v vac=20 -v CONVFMT=%.16f \
        '{$3 * co > co - vac/2 ? $3=$3*h + (cn - co)/cn : $3=$3*h}1' >${s[$i]}_temp
        let x+=${n[$i]}
    done

    # combine elements
    for i in ${!s[*]}; do
        cat ${s[$i]}_temp >>scaledPos
    done

    sed -i 's/ 0 / 0\. /g' scaledPos

    rm *_temp
    mv scaledPos POSCAR

}

scaleSlab() {

    local x s n i a co cn h lo ln l j sc sb t hi

    s=($(storePos s))
    n=($(storePos n))

    # get user input
    echo -en "\nLayer height -> "
    read h
    echo -n "Current number of layers -> "
    read lo
    echo -n "Target number of layers -> "
    read ln
    echo -n "Number of atoms per element in a single layer (${s[*]}) -> "
    read -a a

    if [[ -z $a ]]; then
        for i in ${n[*]}; do
            a+=($((i / lo)))
        done
    fi

    # copy POSCAR header to new POSCAR
    head -2 POSCAR >scaledPos

    # scale current slab
    python3.7 - $h $lo $ln <<-END

    import numpy as np
    from sys import argv
    from pymatgen.core import structure, Lattice
    from pymatgen.io.vasp.inputs import Poscar

    h = float(argv[1])
    lo, ln = [int(i) for i in argv[2:]]

    slab_in = structure.Structure.from_file('POSCAR')
    lattice = slab_in.as_dict().get('lattice')

    old_lattice = Lattice.from_parameters(lattice.get('a'),
                                          lattice.get('b'),
                                          lattice.get('c'),
                                          lattice.get('alpha'),
                                          lattice.get('beta'),
                                          lattice.get('gamma'),
                                          vesta=True)

    new_lattice = Lattice.from_parameters(old_lattice.a,
                                          old_lattice.b,
                                          old_lattice.c + (h * (ln - lo)),
                                          old_lattice.alpha,
                                          old_lattice.beta,
                                          old_lattice.gamma,
                                          vesta=True)

    with open('scaledPos', 'a') as f:
        f.write(np.array2string(new_lattice.matrix) + "\n")

END

    # format matrix and finish writing header
    sed -i -e 's/[[]//g' -e 's/[]]//g' scaledPos
    head -9 POSCAR | tail -4 >>scaledPos

    # modify number of each element
    for i in ${!s[*]}; do
        sed -i "7s/[0-9][0-9]*/$((n[i] + a[i] * (ln - lo)))/$((i + 1))" scaledPos
    done

    # determine scaling ratios in long axis for slab and bulk
    co=$(getLat disp | tail -1)
    cn=$(getLat disp scaledPos | tail -1)
    sc=$(echo "scale=16; $co / $cn" | bc)
    sb=$(echo "scale=16; $h / $cn" | bc)

    # split POSCAR by element
    # scale z-axis coordinates to fit scaled cell
    x=0
    for i in ${!s[*]}; do
        getCoords | head -$((n[i] + x)) | tail -${n[$i]} | awk -v c=$sc -v CONVFMT=%.16f '{$3=$3*c}1' >${s[$i]}_temp
        let x+=${n[$i]}
    done

    # insert new layers
    # use scaled layer height to increment z-axis coordinates
    # modify selective dynamics
    for i in ${!s[*]}; do
        l=${a[$i]}
        x=$(grep -n 'F F F' ${s[$i]}_temp | tail -1 | cut -d ':' -f 1)
        set ${x:=$((n[i] / 2))}
        for j in $(seq $((ln - lo))); do
            # insert layer
            t="$(grep -B $((l - 1)) -- "$(sed -n "$x"p ${s[$i]}_temp)" ${s[$i]}_temp)"
            hi=$(echo "$t" | tail -1 | cut -d ' ' -f 3) # increment z for all z > hi
            awk -v t="$t" "NR==$((x + 1)){print t}1" ${s[$i]}_temp >$i
            mv $i ${s[$i]}_temp
            # increment z-axis coordinate for new layers and above
            awk -v l=$x -v s=$sb -v hi=$hi -v CONVFMT=%.16f '{if (NR > l || $3 > hi) $3=$3+s}1' ${s[$i]}_temp >$i
            mv $i ${s[$i]}_temp
            sed -i 's/ 0 / 0\. /g' ${s[$i]}_temp
            let x+=$l
        done
        # modify selective dynamics
        m=$(cat ${s[$i]}_temp | wc -l)
        ((m % 2 != 0)) && ((m++)) || let m+=2
        awk -v m=$m -v l=$l 'NR < m/2 - l/2 || NR >= m/2 + l/2 {sub("F F F","T T T")}1' ${s[$i]}_temp >${s[$i]}_temp1
        #cat ${s[$i]}_temp1 | sort -k3 > ${s[$i]}_temp; rm ${s[$i]}_temp1
        mv ${s[$i]}_temp1 ${s[$i]}_temp
    done

    # combine elements
    for i in ${!s[*]}; do
        cat ${s[$i]}_temp >>scaledPos
    done

    rm *_temp
    mv scaledPos POSCAR
    echo ""

}

totRelaxTime() {
    if [ ! -f status ]; then
        echo -e "\nMissing status file\n"
        return
    fi
    local runTimes tot i
    runTimes=$(grep E0 status | grep ' R' | cut -d ' ' -f 20 | xargs)
    tot=0
    for t in ${runTimes[*]}; do
        let tot+=$(date -u -d "jan 1 1970 $t" +%s)
    done
    if ((tot > 86400)); then echo -n "($((tot / 86400))) "; fi
    date -u -d @${tot} +%T
}

SD() {

    if [[ ! $1 ]]; then return; fi

    local ax long low high

    long=$(getLat disp | sort -n | tail -1)
    ax=$(index $long "$(getLat disp)")
    low1=$1
    high1=$2
    low2=$3
    high2=$4

    python3.7 - $ax $long $low1 $high1 $low2 $high2 <<-END

    from sys import argv

    axis = int(argv[1])
    l, low1, high1 = [float(i) for i in argv[2:5]]
    low2 = high2 = 0
    if len(argv) == 7: low2, high2 = [float(i) for i in argv[5:]]

    coordType = 'Direct'

    l = 1  # set l to 1 for direct coordinates

    with open('POSCAR', 'r') as f, open('SD', 'w') as g:

        for line in range(5):
            g.write(f.readline())

        atoms = f.readline()
        g.write(atoms)
        atoms = atoms.split()
        nums = f.readline()
        g.write(nums)
        nums = nums.split()
        g.write('SD\n')
        g.write(f.readline())

        for s in range(len(nums)):

            lines = [[float(value) for value in f.readline().split()]
                    for line in range(int(nums[s]))]

            for line in sorted(lines, key=lambda line: line[axis]):
                SD = ' F F F' if low1 * l < line[
                    axis] < high1 * l or low2 * l < line[
                        axis] < high2 * l else ' T T T'
                g.write('%18.9f %18.9f %18.9f' % tuple(line) + SD + ' \n')

END

    mv SD POSCAR

}

clip() { # remove lines from file edges
    if [[ $# < 2 ]]; then
        echo -e "\nUSAGE: clip <s(start) | e(end) | b(both)> <# of lines> *<file name>\n"
    else
        if [[ $3 ]]; then file=$3; else file=''; fi
        case $1 in
        s) tail -n +$(($2 + 1)) $file ;;
        e) head -n -$2 $3 ;;
        b) tail -n +$(($2 + 1)) $file | head -n -$2 ;;
        esac
    fi
}

bandDecomp() {

    local name path rType range c

    run() {
        if [ $1 == 1 ]; then
            echo -ne "\nEnter range name -> "
            read name
            path=charge/$name
            mkdir -p $path
            if [ ! -d charge ]; then
                echo -e "\nFailed to create necessary folders\n"
                return
            fi
            cp POSCAR POTCAR KPOINTS WAVECAR job $path
            sed -i "s/\(\-J\).*/\1 $name/" $path/job
            echo -e "mv PARCHG ../$name.vasp; rm WAVECAR" >>$path/job
            echo -n "Enter range type ( 0: energy | 1: bands ) -> "
            read rType
            set ${rType:=0}
            echo -n "Enter range -> "
            read range
            if [[ $rType == 0 ]]; then
                rType=EINT
                range=$(float $range)
            else
                rType=IBAND
            fi

            # Consider fixing disTag and enbTag to act on file in other locations

            disTag IBAND EINT
            enbTag $rType
            cp INCAR $path
            setTag "$rType LPARD ISTART ICHARG NPAR KPAR" "$(joinStr , $range) T 1 0 1 1" $path
            sbatch -D $path $path/job
            echo -ne "\nProcess another range? ( 1: yes | 0: no ) -> "
            read c
            set ${c:=0}
            run $c
        else
            echo ""
            cp prevJob job
        fi
    }
    sed -i -e "s/\(setTag LWAVE\) .*/\1 F/" -e "s/\(setTag LCHARG\) .*/\1 F/" job
    queue 1
    nodes 1
    cores 24
    tLim 2 5
    run 1
}

archTest() {

    local a n c

    run() {
        if [ $1 == 1 ]; then
            a=(${a[*]} $n)
            getJob IN
            mkdir -p $n
            cp POSCAR INCAR POTCAR KPOINTS job $n
            echo -ne "Process another architecture? ( 1: yes | 0: no ) -> "
            read c
            set ${c:=0}
            ((n++))
            run $c
        else
            for i in ${a[*]}; do
                sbatch -D $i $i/job
            done
            disTag NELM
            echo ""
            return
        fi
    }

    a=()
    n=$(find [0-9] -maxdepth 0 -type d 2>/dev/null | tail -1)
    if [[ $n ]]; then ((n++)); else n=1; fi
    enbTag NELM
    setTag "LCHARG NELM" "F 1"
    run 1
}

float() {
    local a
    for i in $*; do
        a=(${a[*]} $(echo "$i * 1.0" | bc))
    done
    echo "${a[*]}"
}

plot() {
    if [[ ! $1 ]]; then
        echo -e "\nUsage: plot <type> <loc>\n"
    elif [[ $1 == edit ]]; then
        vi /work/01293/hudamn/Edan/files/plot.py
    else
        python3.7 /work/01293/hudamn/Edan/files/plot.py ${@:1}
        echo ""
    fi
}

waveFind() {

    local d

    if [[ $1 == d ]]; then d='-delete'; fi
    find -name WAVECAR -not -empty $d
}

getChgDen() {

    local n nIon i paw j pos

    mkdir charge 2>/dev/null
    if [ ! -d charge ]; then
        echo -e "\nFailed creating folder\n"
    elif [[ ! $(grep "\-J" job | grep SP) ]]; then
        echo -e "\nMissing SP job file\n"
    else
        s=($(storePos s))
        n=($(storePos n))
        nIon=$(storePos i)
        if [[ $(grep PAW_PBE POTCAR) ]]; then paw=2; else paw=0; fi
        queue 1
        cp INCAR charge
        cp SCF_KPOINTS charge/KPOINTS
        cp job charge/tempJOB
        cat POSCAR >charge/tempPOS
        cd charge
        if [[ $(head -8 tempPOS | grep '^[Ss]') ]]; then sed -i 8d tempPOS; fi
        head -$((nIon + 8)) tempPOS >POSCAR
        setTag 'ISTART ICHARG' '0 2'

        j=1
        for i in ${!s[*]}; do
            mkdir -p ${s[$i]}
            pos=$((8 + j))
            sed -n "$pos,$((pos + ${n[$i]} - 1))p;1,8p" POSCAR >${s[$i]}/POSCAR
            sed -i -e "6s/.*/   ${s[$i]}/" -e "7s/.*/   ${n[$i]}/" ${s[$i]}/POSCAR
            cp tempJOB ${s[$i]}/job
            echo -e "\ne=-3" >>${s[$i]}/job
            echo -en "c=(\$(grep : OSZICAR | rmSpace | tail -1 | cut -d ' ' -f 4,5 | " >>${s[$i]}/job
            echo "grep -o E... | sed 's/0//' | getNumber | xargs))" >>${s[$i]}/job
            echo -e "if (( \${c[0]} <= \$e )) && (( \${c[1]} <= \$e )); then" >>${s[$i]}/job
            echo -e "  mv CHGCAR ../${s[$i]}.vasp; rm CHG*" >>${s[$i]}/job
            echo -e "fi" >>${s[$i]}/job
            sed -i -e "s/\(setTag LWAVE\) .*/\1 F/" -e "s/\(setTag LCHARG\) .*/\1 T/" ${s[$i]}/job
            potGen $paw ${s[$i]}
            cp INCAR POTCAR KPOINTS ${s[$i]}
            if [[ $1 ]]; then sbatch -D ${s[$i]} ${s[$i]}/job; fi
            j=$((j + ${n[$i]}))
        done

        rm INCAR POTCAR KPOINTS tempJOB tempPOS POSCAR
        cd ../
    fi
}

vbm() {

    local nElect vbmI vbm cbm t

    nElect=$(head -6 EIGENVAL | tail -1 | getNumber | head -1)
    if ((nElect % 2 != 0)); then ((nElect++)); fi
    if [[ $(head -9 EIGENVAL | tail -1 | rmSpace | wc -w) == 3 ]]; then col='2'; else col='2,3'; fi

    vbmI=$((nElect / 2))
    vbm=($(tail -n +8 EIGENVAL | grep " $vbmI " | rmSpace | cut -d ' ' -f $col | sort -n | tail -1))
    cbm=($(tail -n +8 EIGENVAL | grep " $((vbmI + 1)) " | rmSpace | cut -d ' ' -f $col | sort -nr | tail -1))

    if [ $2 ]; then t=$2; else t='1-'; fi

    case $1 in
    i) echo $vbmI ;;
    eV) echo ${vbm[*]} | cut -d ' ' -f $t ;;
    eC) echo ${cbm[*]} | cut -d ' ' -f $t ;;
    *) echo -e "\nUsage: vbm <i: index | eV: vbm energies | eC: cbm energies>\n" ;;
    esac

}

occCheck() {

    local vbm cbm patV patC n nElect

    if [ ! -f EIGENVAL ]; then
        echo -e '\nNo EIGENVAL file found\n'
        return
    fi

    nElect=$(head -6 EIGENVAL | tail -1 | getNumber | head -1)
    vbm=$(vbm i)
    if [[ $1 ]]; then
        if [[ $(head -9 EIGENVAL | tail -1 | rmSpace | wc -w) == 3 ]]; then
            if ((nElect % 2 != 0)); then patV=0.500000; else patV=1.000000; fi
            patC=0.000000
        else
            if ((nElect % 2 != 0)); then patV='1.000000 0.000000'; else patV='1.000000 1.000000'; fi
            patC='0.000000 0.000000'
        fi
        n=$(head -6 EIGENVAL | tail -1 | rmSpace | cut -d ' ' -f 2)
        if [[ $(grep " $vbm " EIGENVAL | rmSpace | grep "$patV" | wc -l) == $n ]] &&
            [[ $(grep " $((vbm + 1)) " EIGENVAL | rmSpace | grep "$patC" | wc -l) == $n ]]; then echo ''; else echo '*'; fi
    else
        grep -C 1 " $vbm " EIGENVAL
    fi
}

backUp() {
    if [ ! $1 ]; then
        echo -e "\nUsage: backUp <jobType> <args>\n"
    else
        case $1 in
        IN)
            mkdir -p output/init/$2/$3
            cp INCAR KPOINTS POSCAR OUTCAR EIGENVAL Vasp.out status output/init/$2/$3
            ;;
        R)
            mkdir -p output/relax/$2
            cp INCAR POSCAR OUTCAR EIGENVAL XDATCAR Vasp.out output/relax/$2
            ;;
        V)
            mkdir -p output/relax/EV/$2
            cp INCAR POSCAR OUTCAR EIGENVAL XDATCAR Vasp.out output/relax/EV/$2
            ;;
        C)
            mkdir -p output/sigma/$2/$3/$4
            cp EIGENVAL OUTCAR Vasp.out output/sigma/$2/$3/$4
            ;;
        SP)
            mkdir -p output/sp/$2
            cp INCAR EIGENVAL OUTCAR Vasp.out vasprun.xml output/sp/$2
            ;;
        D)
            mkdir -p output/dos/$2
            cp INCAR POTCAR OUTCAR DOSCAR EIGENVAL Vasp.out output/dos/$2
            cp vasprun.xml output/dos/$2/dos.xml
            ;;
        B)
            mkdir -p output/$2
            cp INCAR EIGENVAL KPOINTS OUTCAR output/$2
            cp vasprun.xml output/$2/bands.xml
            ;;
        esac
    fi
}

math() {

    if [ ! $1 ]; then
        echo -e "\nUsage: math <avg sum min max ceil floor triProd> \"<args>\"\n"
    elif [[ ! $2 ]]; then
        echo -e "\nEmpty argument list\n"
    else
        python3.7 - $1 "${@:2}" <<-END

        from sys import argv
        import numpy as np

        if argv[1] == 'avg':
            print(round(np.average([float(i) for i in argv[2:]]), 3))
        elif argv[1] == 'sum':
            print(round(np.sum([float(i) for i in argv[2:]]), 3))
        elif argv[1] == 'min':
            print(round(np.min([float(i) for i in argv[2:]]), 3))
        elif argv[1] == 'max':
            print(round(np.max([float(i) for i in argv[2:]]), 3))
        elif argv[1] == 'ceil':
            print(int(np.ceil([float(i) for i in argv[2:]])[0]))
        elif argv[1] == 'floor':
            print(int(np.floor([float(i) for i in argv[2:]])[0]))
        elif argv[1] == 'triProd':
            v = [np.array([float(j) for j in argv[2+i:5+i]]) for i in (0, 3, 6)]
            if len(v) == 3:
                print(np.dot(v[0], np.cross(v[1], v[2])))
END
    fi

}

oxyCore() {
    local c a path
    if [ $1 ]; then
        path=$1
    elif [ -d OUTCAR ] && [ -f vasprun.xml ] || [ -f bands.xml ]; then
        path=./
    else
        echo -e "\n -> Missing files...\n"
        return
    fi
    file=$(find $1/*.xml)
    c=$(grep ' <rc><c>O </c><c> ' $file | wc -l)
    a=($(grep 1s $path/OUTCAR | tail -$c | rmSpace | cut -d ' ' -f 3))
    if [ ${#a[*]} != 0 ]; then
        echo -e "\n Avg O 1s = $(math avg ${a[*]}) eV\n"
    else
        echo -e "\n -> No 1s states available...\n"
    fi
}

goto() {
    local dir
    #echo -ne "\nSelect directory: \n\n"
    #echo -ne "\n-> "
    #read dir

    home=/work/01293/hudamn/stampede2/EDAN

    case $dir in
    *) cd $home ;;
    esac
    #echo ''
}

makeSuper() {
    if (($# != 3)); then
        echo -e "\nUsage: makeSuper n m l (n m l are multiples of lattice parameters A B C)\n"
    else
        cp POSCAR unitCell
        python3.7 - $* <<-END

        from sys import argv
        from pymatgen.core.structure import Structure
        from pymatgen.io.vasp.inputs import Poscar

        st = Structure.from_file('POSCAR')
        st.make_supercell(argv[1:])

        Poscar(st).write_file('POSCAR')

END
    fi
}

highSymK() {

    local args gen

    args=(0.01 0.0 0.0)
    gen=False

    if [[ $1 ]]; then
        if [[ $1 == edit ]]; then
            vi /work/01293/hudamn/Edan/files/makeKpoints.py
            return
        elif [[ $1 == disp ]] || [[ $1 == gen ]]; then
            if [[ $1 == gen ]]; then gen=True; fi
            if [[ $2 ]]; then args=($2); fi
        fi
    else
        echo -e "\nUSAGE: HighSymK <edit | disp | gen> [symprec, angTol, AbsTol]>\n"
        return
    fi

    if [ ! -f SCF_KPOINTS ] && [ -f KPOINTS ]; then cp KPOINTS SCF_KPOINTS; fi
    python3.7 /work/01293/hudamn/Edan/files/makeKpoints.py "${args[*]}" $gen
    echo ""
}

dispFunc() { # display all available user-defined functions
    echo ""
    grep '^  [A-Za-z][A-Za-z]* ()' /work/01293/hudamn/Edan/genFunc
    echo ""
}

rmSpace() { # removes extra whitespace including tabs
    sed -e "s/\t/ /g" -e 's/  */ /g' -e 's/^ //'
}

joinStr() { # Usage: joinStr <delimiter> <args>
    IFS=$1
    local args
    args=$*
    echo "${args[*]:2}"
}

getNumber() {
    grep -o "\-*[0-9][0-9]*\.*[0-9]*"
}

range() { # Usage: range <upper bound>
    local rng ind
    rng=()
    if [ $2 ]; then ind=$2; else ind=0; fi
    while (($ind < $1)); do
        rng+=($ind)
        ((ind++))
    done
    echo ${rng[*]}
}

index() { # return index of element in array
    local a i
    a=($2)
    for i in ${!a[*]}; do if [[ ${a[$i]} == $1 ]]; then echo $i; fi; done
}

inArray() { # Usage: inArray <var/str> <"array">
    local i
    for i in $2; do if [[ $i == $1 ]]; then echo 1; fi; done
}

getPrimitive() {
    if [ ! -f POSCAR ]; then
        echo -e "\nNo POSCAR found in current directory\n"
        return
    else
        if [ -f fullPOS ]; then cp fullPOS POSCAR; fi
        cp POSCAR fullPOS

        python3.7 <<-END

        from pymatgen.core.structure import Structure
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer as sga
        from pymatgen.io.vasp.inputs import Poscar

        st = Structure.from_file('POSCAR')

        Poscar(sga(st).get_primitive_standard_structure(international_monoclinic=False)).write_file('POSCAR')

END
        echo ""
        cat POSCAR
        echo -e "\n -> written to POSCAR...\n"
    fi
}

genMovie() { # Usage: genMovie <space-separated list of relaxation trial numbers>
    local dir str files i
    #    dir="output/relax"
    dir=./
    str=$(joinStr , $*)
    files=($(find $dir/trial_[$str]/*/XDATCAR))
    cat ${files[0]} >xDat
    for i in ${files[*]:1}; do
        tail -n +8 $i >>xDat
    done
    cat $(find $dir/trial_[$str]/*/OUTCAR) | grep "energy  without entropy" >xOut
    cat $(find $dir/trial_[$str]/*/OUTCAR) | grep FORCES >>xOut
    /work/01293/hudamn/Edan/files/extra/xdat2xyz.pl
    rm xDat xOut
}

duplicateLines() {
    local n prev
    n=0
    cat $1 | while read line; do
        if ((n > 0)); then
            if [ ! "$line" == "$prev" ]; then echo $line; fi
        else
            echo $line
        fi
        prev="$line"
        ((n++))
    done >$1_temp
    mv $1_temp $1
}

kList() { # Usage: kList <dir>
    local folder file
    if [ ! $1 ] || [ ! -d $1 ]; then
        echo -e "\nMissing path!\n"
        return
    else folder=$1; fi
    file=$(find *.xml)
    sed -n '/kpointlist/,/weights/p' $folder/$file | grep '<v>' |
        rmSpace | cut -d ' ' -f 2-4 | xargs printf "%0.6f %0.6f %0.6f\n" >kList
    duplicateLines kList
}

bandEdges() { # obtain band edges for selected spin (provide any argument to print band edges)

    local s col vbmI cbmI vbm cbm eV kV kVn eC kC kCn

    kList ./

    vbmI=$(vbm i)
    cbmI=$((vbmI + 1))
    vbm=($(vbm eV))
    cbm=($(vbm eC))

    echo -ne "\nSpin ( 1-up | 2-down ) -> "
    read s

    if [[ $s =~ ^[12]$ ]]; then
        eV=${vbm[$((s - 1))]}
        kV=$(tail -n +8 EIGENVAL | rmSpace | grep -B $vbmI "$vbmI $(if ((s == 2)); then echo ".* "; fi)$eV" |
            head -1 | cut -d ' ' -f 1-3 | xargs printf "%0.6f ")
        kV=$(for i in $(echo $kV | cut -d ' ' -f 1-); do if [[ $i == -0.000000 ]]; then
            echo -n '0.000000 '
        else echo -n "$i "; fi; done)
        kVn=$(grep -n -- "^$(echo $kV)" kList | cut -d ':' -f 1 | head -1)
        eC=${cbm[$((s - 1))]}
        kC=$(tail -n +8 EIGENVAL | rmSpace | grep -B $cbmI "$cbmI $(if ((s == 2)); then echo ".* "; fi)$eC" |
            head -1 | cut -d ' ' -f 1-3 | xargs printf "%0.6f ")
        kCn=$(grep -n -- "^$(echo $kC)" kList | cut -d ':' -f 1 | head -1)
        if ((s == 1)); then echo -e "Spin Up\n"; else echo -e "Spin Down\n"; fi
        printf "  VBM: %0.3f eV @ %d [ %s]\n" $eV $kVn "$kV"
        printf "  CBM: %0.3f eV @ %d [ %s]\n" $eC $kCn "$kC"
        printf "  Gap: %0.3f eV\n" $(echo $eC - $eV | bc)
        echo ""
    else
        echo -e "\nBad spin value\n"
    fi

    if [[ ! $1 ]]; then rm kList; fi
}

emc() {
    local k band spin st step
    if [[ $1 == 1 ]]; then
        echo -ne "\nK-point -> "
        read k
        if [[ ! $k ]]; then
            if [ -f KPOINTS ]; then k="$(head -4 KPOINTS | tail -1 | rmSpace | cut -d ' ' -f 1-3)"; else return; fi
        fi
        echo -n "Band number -> "
        read band
        if [[ ! $band ]]; then
            if [ -f KPOINTS ]; then band=$(head -1 KPOINTS | cut -d ' ' -f 1); else return; fi
        fi
        echo -n "Spin ( 1-up | 2-down ) -> "
        read spin
        if [[ ! $spin ]]; then
            if [ -f KPOINTS ]; then spin=$(head -1 KPOINTS | cut -d ' ' -f 2); else return; fi
        fi
        echo -n "Stencil for central difference ( 3 | 5 ) -> "
        read st
        if [[ ! $st ]]; then
            if [ -f KPOINTS ]; then st=$(head -1 KPOINTS | cut -d ' ' -f 4); else return; fi
        fi
        echo -n "Step size -> "
        read step
        if [[ ! $step ]]; then
            if [ -f KPOINTS ]; then step=$(head -1 KPOINTS | cut -d ' ' -f 3); else return; fi
        fi
        python3.7 /work/01293/hudamn/Edan/files/emc.py 'gen' $band $spin $step $st "$k"
        echo ""
    elif [[ $1 == 2 ]]; then
        python3.7 /work/01293/hudamn/Edan/files/emc.py 'calc'
        highSymK
    elif [[ $1 == 3 ]]; then
        vi /work/01293/hudamn/Edan/files/emc.py
        return
    else
        echo -e "\nUsage: emc <option> -> 1: KPOINTS generation; 2: Effective mass calulation; 3: edit EMC\n"
        return
    fi
}

genK() {
    local c k kc kb kf step n
    echo -ne "\nHole (h) or Electron (e) -> "
    read c
    if [[ ! $c =~ ^[he]$ ]]; then
        echo -e "\nInvalid selection\n"
        genK
    fi
    bandEdges 1
    echo -n "k-point index -> "
    read k
    if [[ $k =~ ^[0-9]+$ ]]; then
        if ((k <= $(cat kList | wc -l))); then
            kc=($(head -$k kList | tail -1 | cut -d ' ' -f 1-3))
            kb=($(head -$((k - 1)) kList | tail -1 | cut -d ' ' -f 1-3))
            if [[ ! $kb ]]; then kb=0; fi
            kf=($(head -$((k + 1)) kList | tail -1 | cut -d ' ' -f 1-3))
            if [ $(head -$((k + 1)) kList | wc -l) == $(head -$k kList | wc -l) ]; then kf=0; fi
            echo -n "Step size (2pi/A) -> "
            read step
            if [[ $step =~ ^[0-9]+\.*[0-9]*$ ]]; then
                echo -n "Number of steps away from k-point -> "
                read n
                if [[ ! $n =~ ^[0-9]+$ ]]; then
                    echo -e "\nValue must be a positive integer"
                    genK
                fi
            else
                echo -e "\nValue must be a positive integer or float"
                genK
            fi
        else
            echo -e "\nk-point out of range"
            genK
        fi
    else
        echo -e "\nValue must be a positive integer"
        genK
    fi

    mkdir -p $c
    cp ../../../{INCAR,POTCAR,CHGCAR,POSCAR,job} $c

    python3.7 - "${kb[*]}" "${kc[*]}" "${kf[*]}" $step $n $c <<-END

    import numpy as np
    from numpy import pi, dot, array
    from numpy.linalg import norm, inv
    from sys import argv
    from pymatgen.core.structure import Structure as st

    rec = st.from_file(argv[6] + '/POSCAR').lattice.reciprocal_lattice
    kb, kc, kf = [
        rec.get_cartesian_coords(i) if len(i) > 1 else "0" for i in [
            array([float(j) for j in k.split()]) if len(k) > 1 else "0"
            for k in argv[1:4]
        ]
    ]

    # define unit vectors
    if len(kb) > 1: ub = (kb - kc) / norm(kb - kc)
    if len(kf) > 1: uf = (kf - kc) / norm(kf - kc)

    k = []

    if len(kb) > 1:
        for i in range(int(argv[5])):
            k.append(rec.get_fractional_coords(kc + ub * (i + 1) * float(argv[4])))
        k.reverse()

    k.append(rec.get_fractional_coords(kc))

    if len(kf) > 1:
        for i in range(int(argv[5])):
            k.append(rec.get_fractional_coords(kc + uf * (i + 1) * float(argv[4])))

    with open('kTest.txt', 'w') as f:
        f.write("Step size = {0:0.3f} 2pi/A\n{1:d}\nrec\n".format(
            float(argv[4]), len(k)))
        for i in k:
            f.write(" {0:11.8f} {1:11.8f} {2:11.8f}  0.01\n".format(
                i[0], i[1], i[2]))

    print("\n -> k-point grid generated!")

END

    mv kTest.txt $c/KPOINTS
    echo ""
}

effM() { # FIX COMPATABILITY WITH ISPIN=1

    # conversion factor -> (6.582e-16)^2[eVs]^2 * dk^2[1/A]^2 / d2E[eV] * {(3.00e18[A/s])^2 / 511[keV/c^2]} =~ 7.63

    local vbmI cbmI col vbm cbm s bn e k delK kb kf st step kpts

    kList ./

    vbmI=$(vbm i)
    cbmI=$((vbmI + 1))
    vbm=($(vbm eV))
    cbm=($(vbm eC))

    echo -ne "\nSpin ( 1-up | 2-down ) -> "
    read s
    if [[ ! $s =~ ^[12]$ ]]; then s=1; fi
    echo ''

    echo -n "Band index -> "
    read bn # get band index
    if [[ $bn =~ ^[0-9]+$ ]]; then
        if (($bn <= $(grep -m 1 NBANDS *.xml | getNumber))); then
            tail -n +8 EIGENVAL | grep " $bn " | rmSpace | cut -d ' ' -f 2,3 >eList
            duplicateLines eList
        else
            echo -e "\nBand number out of range"
            effM
        fi
    else
        echo -e "\nBand value must be a positive integer"
        effM
    fi

    echo ""
    if [ $bn == $vbmI ]; then
        e=${vbm[$((s - 1))]}
    elif [ $bn == $cbmI ]; then
        e=${cbm[$((s - 1))]}
    elif ((bn < vbmI)); then
        e=$(cat eList | cut -d ' ' -f $s | sort | tail -1)
    else
        e=$(cat eList | cut -d ' ' -f $s | sort -r | tail -1)
    fi
    grep -m 1 -n -A 20 -B 20 $e eList | sed 's/[:-]/ /' | cut -d ' ' -f 1,$((s + 1)) | sed "s/\($e\)/\1 <- \*/"
    echo ""

    # get k-point indices
    echo -n "k-point index -> "
    read k
    if [[ $k =~ ^[0-9]+$ ]]; then
        if ((k <= $(cat kList | wc -l))); then
            echo -n "Step size (points away from k) -> "
            read delK
            if [[ $delK =~ ^[0-9]+$ ]]; then
                kb=$((k - delK))
                kf=$((k + delK))
            else
                echo -e "\nDelta k must be a positive integer"
                effM
            fi
            echo -n "Stencil ( 3 | 5 ) -> "
            read st
            if [ $st == 3 ] || [ $st == 5 ]; then step=$(((st - 1) / 2)); else
                echo -e "\nStencil must be 3 or 5"
                return
            fi
            if ((k + step * (kf - k) > $(cat kList | wc -l))) || ((k - step * (kf - k) < 0)); then
                echo -e "\nk-points out of range"
                effM
            else
                if [[ $step == 2 ]]; then
                    kpts=($((2 * kb - k)) $kb $k $kf $((2 * kf - k)))
                else
                    kpts=($kb $k $kf)
                fi
            fi
        else
            echo -e "\nk-point out of range"
            effM
        fi
    else
        echo -e "\nk-point index must be a positive integer"
        effM
    fi

    # convert kpoints to cartesian coordinates and calculate effective mass
    python3.7 - "${kpts[*]}" $s <<-END
    
    from numpy import array, abs
    from numpy.linalg import norm
    from sys import argv
    from pymatgen.core.structure import Structure as st

    rec = st.from_file('POSCAR').lattice.reciprocal_lattice

    indices = [int(i) - 1 for i in argv[1].split(' ')]

    kf = []  # fractional coordinates
    with open('kList', 'r') as f:
        for line in enumerate(f.readlines()):
            if line[0] in indices:
                kf.append(
                    array([float(i) for i in line[1].split('\n')[0].split()[0:3]]))

    k = [rec.get_cartesian_coords(i)
        for i in kf]  # converted to cartesian coordinates
    kn = [i / max(abs(i)) if norm(i) != 0 else i for i in k]  # normalize

    e = []  # energies in eV
    with open('eList', 'r') as f:
        for line in enumerate(f.readlines()):
            if line[0] in indices:
                e.append(float(line[1].split('\n')[0].split(' ')[int(argv[2]) -
                                                                1]))

    print("")
    print(e)

    if len(k) == 3:
        M = 7.63050245928016 / ((e[0] - 2.0 * e[1] + e[2]) / norm(k[2] - k[1])**2)
    else:
        M = 7.63050245928016 / (
            (-e[0] + 16.0 * e[1] - 30.0 * e[2] + 16.0 * e[3] - e[4]) /
            (12.0 * norm(k[3] - k[2])**2))

    print("")
    for i, kpt in enumerate(zip(kf, kn)):
        print(
            "k ( frac | cart ) = [%7.4f %7.4f %7.4f ] | [%7.4f %7.4f %7.4f ] ; e = %9.6f"
            % (kpt[0][0], kpt[0][1], kpt[0][2], kpt[1][0], kpt[1][1], kpt[1][2],
            e[i]))

    print("")
    index = int((len(k) + 1) / 2 - 1)
    print("backward = [ {0:0.3f} eV ; {1:0.3f} 2pi/A ]".format(
        abs(e[index] - e[index - 1]), norm(k[index] - k[index - 1])))
    print("forward  = [ {0:0.3f} eV ; {1:0.3f} 2pi/A ]".format(
        abs(e[index + 1] - e[index]), norm(k[index + 1] - k[index])))

    print("\nM_eff = %0.3f m_e\n" % M)

END

    #rm eList kList

}

calcBader() { # calculate bader charges
    if [ -e AECCAR0 ] && [ -e AECCAR2 ]; then
        /work/01293/hudamn/Edan/files/extra/chgsum.pl AECCAR0 AECCAR2
    else
        echo -e "\nCore charges are missing!\n"
        return
    fi
    if [ -e CHGCAR_sum ]; then
        /work/01293/hudamn/Edan/files/extra/bader CHGCAR -ref CHGCAR_sum
        rm CHGCAR_sum
    else
        echo -e "\nMissing sum over core charges!\n"
        return
    fi
}

baderCharges() { # Usage: getBaderCharges ...

    local s n k j x i a c z

    s=($(storePos s))
    n=($(storePos n))
    v=($(storePos v))

    if [[ ! $1 ]]; then
        echo -e "\nAverage net charges for (${s[*]}):\n"
        k=(${n[*]})
        j=1
        x=0
        for i in ${!s[*]}; do
            a=0
            while [ $j -le ${n[$i]} ]; do
                c=$(grep " $j " ACF.dat | rmSpace | cut -d ' ' -f 5)
                printf "%s_%d = %6.3f e\n" ${s[$i]} $((j - x)) $(echo "${v[$i]} - $c" | bc)
                a=$(echo $a + $c | bc)
                ((j++))
            done
            z=$(printf "%0.3f e\n" $(echo "scale=16;${v[$i]} - $a / ${k[$i]}" | bc))
            echo -e "\n${s[$i]}_avg = $z\n"
            n[$i + 1]=$((n[$i + 1] + n[$i]))
            x=$((j - 1))
        done
    else
        a=()
        a=$*
        for i in ${a[*]}; do
            s=()
            n=$(echo $i | cut -d '_' -f 1)
            if [ $(echo $i | grep _) ]; then s=($(echo $i | cut -d '_' -f 2 | sed 's/,/ /g')); fi
            c=($(grep " $n " ACF.dat | sed 's/  */ /g' | cut -d ' ' -f 3,4,5,6))
            if [ ${s[0]} ]; then
                x=($(head -3 POSCAR | tail -1))
                c[0]=$(echo "${c[0]} + ${x[0]} * ${s[0]}" | bc)
                c[1]=$(echo "${c[1]} + ${x[1]} * ${s[0]}" | bc)
                c[2]=$(echo "${c[2]} + ${x[2]} * ${s[0]}" | bc)
            fi
            if [ ${s[1]} ]; then
                x=($(head -4 POSCAR | tail -1))
                c[0]=$(echo "${c[0]} + ${x[0]} * ${s[1]}" | bc)
                c[1]=$(echo "${c[1]} + ${x[1]} * ${s[1]}" | bc)
                c[2]=$(echo "${c[2]} + ${x[2]} * ${s[1]}" | bc)
            fi
            if [ ${s[2]} ]; then
                x=($(head -5 POSCAR | tail -1))
                c[0]=$(echo "${c[0]} + ${x[0]} * ${s[2]}" | bc)
                c[1]=$(echo "${c[1]} + ${x[1]} * ${s[2]}" | bc)
                c[2]=$(echo "${c[2]} + ${x[2]} * ${s[2]}" | bc)
            fi
            printf "%0.3f %0.3f %0.3f %0.3f\n" ${c[3]} ${c[0]} ${c[1]} ${c[2]}
        done
    fi
}

# get direct gaps
getDirectGap() {
    local vbm up dn pat s
    if [ ! $1 ]; then
        echo -e "\nUsage: getDirectGap vbm\n"
        return
    fi
    vbm=$1
    if [ $(getTag ISPIN) == 1 ]; then s=2; else s=2,3; fi
    tail -n +8 EIGENVAL | grep " $vbm " | rmSpace | cut -d ' ' -f $s >vbm
    tail -n +8 EIGENVAL | grep " $((vbm + 1)) " | rmSpace | cut -d ' ' -f $s >cbm
    if [ $(getTag ISPIN) == 1 ]; then
        paste vbm cbm | awk '{printf "%0.3f\n", $2-$1}' >directGap
    else
        paste vbm cbm | awk '{printf "%0.3f %0.3f\n", $3-$1, $4-$2}' >directGap
    fi
    rm vbm cbm
    up=$(cat directGap | cut -d ' ' -f 1 | sort -nr | tail -1)
    printf "\nSpin up   -> %0.3f eV" $up
    if [ $(getTag ISPIN) == 2 ]; then
        dn=$(cat directGap | cut -d ' ' -f 2 | sort -nr | tail -1)
        printf "\nSpin down -> %0.3f eV\n" $dn
    fi
    rm directGap
    echo ""
}

# generate DOS plot
parseD() {

    local path vbm num fermi s n v OFS j m k maxO i ions ranges range rMin rMax \
    a cS cP cD cSu cSd cPu cPd cDu cDd l col xS xP xD xSu xSd xPu xPd xDu xDd \
    p pos1 pos2 species c1 o x gap

    if [ -d PDOS ]; then rm -fr PDOS; fi
    mkdir PDOS/

    if [ $1 ]; then
        path=$1
    else
        path=output/dos/$(grep DOS status | tail -1 | cut -d ' ' -f 6-10 | sed 's/ //g')
    fi

    vbm=$(math floor $(math min $(vbm eV)))
    num=$(getTag NEDOS)

    echo -e "\n$(grep '\.... ' $path/DOSCAR | rmSpace | head -$num | grep "^$vbm\.... ") \n"
    echo -n "Enter Fermi level -> "
    read fermi

    s=($(storePos s))
    n=($(storePos n))

    grep '\.... ' $path/DOSCAR | rmSpace >dosData

    head -$num dosData | tail -$num >PDOS/total
    if [[ $(getTag ISPIN) == 2 ]]; then
        awk '{print $1,$2*1,$3*-1,$4,$5}' OFS=" " PDOS/total | sed 's/e/E/g' >PDOS/temp
        mv PDOS/temp PDOS/total
    fi

    l=2
    for j in ${!s[*]}; do
        for k in $(seq ${n[$j]}); do
            m=$((num * l))
            ((l++))
            head -$m dosData | tail -$num >PDOS/${s[$j]}$k
        done
    done
    rm dosData

    # sum over requested ions
    maxO=$(tail -1 $path/DOSCAR | grep -o '[0-9]*\.[0-9]*' | wc -l)
    for j in ${!s[*]}; do

        # get list of sites for current element
        echo -n "Enter ${s[$j]} ions to sum over (1-${n[$j]}) -> "
        read ions
        echo $ions >list
        ranges=($(grep -o '[0-9][0-9]*-[0-9][0-9]*' list))
        for i in ${ranges[*]}; do
            range=()
            rMin=$(echo $i | cut -d '-' -f 1)
            rMax=$(echo $i | cut -d '-' -f 2)
            while ((rMin <= rMax)); do
                range=(${range[*]} $rMin)
                ((rMin++))
            done
            sed -i "s/$i/${range[*]}/" list
        done

        # collect ionic contributions for requested list
        a=()
        if [ $(getTag ISPIN) == 1 ]; then
            cS=('$2 ')
            cP=('$3 ')
            cD=('$4 ')
        else
            cSu=('$2 ')
            cSd=('$3 ')
            cPu=('$4 ')
            cPd=('$5 ')
            cDu=('$6 ')
            cDd=('$7 ')
        fi
        l=1
        for k in $(grep '' list); do
            a=(${a[*]} "PDOS/${s[$j]}$k")
            col=$(echo "$maxO * $l" | bc)
            if [ $(getTag ISPIN) == 1 ]; then
                xS=$((2 + col))
                xP=$((3 + col))
                xD=$((4 + col))
                cS=(${cS[*]} '+ $'"$xS"'')
                cP=(${cD[*]} '+ $'"$xP"'')
                cD=(${cP[*]} '+ $'"$xD"'')
            else
                xSu=$((2 + col))
                xSd=$((3 + col))
                xPu=$((4 + col))
                xPd=$((5 + col))
                xDu=$((6 + col))
                xDd=$((7 + col))
                cSu=(${cSu[*]} '+ $'"$xSu"'')
                cSd=(${cSd[*]} '+ $'"$xSd"'')
                cPu=(${cPu[*]} '+ $'"$xPu"'')
                cPd=(${cPd[*]} '+ $'"$xPd"'')
                cDu=(${cDu[*]} '+ $'"$xDu"'')
                cDd=(${cDd[*]} '+ $'"$xDd"'')
            fi
            ((l++))
        done
        if [ $(getTag ISPIN) == 1 ]; then
            c=("${cS[*]}" "${cP[*]}" "${cD[*]}")
        else
            c=("${cSu[*]}" "${cSd[*]}" "${cPu[*]}" "${cPd[*]}" "${cDu[*]}" "${cDd[*]}")
        fi
        pos1=("%0.5f")
        pos2=()
        for i in $(seq 0 $((maxO - 1))); do
            pos1=(${pos1[*]} " %0.5f")
            pos2=(${pos2[*]} "${c[$i]}")
        done
        col=$(echo ${pos2[*]} | sed -e 's/ + /+/g' -e 's/ /,/g')
        paste $(echo ${a[*]}) | awk '{ printf "'"${pos1[*]}"'\n", $1,'"$col"' }' >PDOS/${s[$j]}
        if [ $(getTag ISPIN) == 2 ]; then
            awk '{print $1,$2*1,$3*-1,$4*1,$5*-1,$6*1,$7*-1}' OFS=" " PDOS/${s[$j]} >PDOS/temp
            mv PDOS/temp PDOS/${s[$j]}
        fi
    done
    rm list

    # combine results for plotting
    vbm=$(vbm i)
    gap=$(getGap $vbm $((vbm + 1)) $path 1)
    if [ $(getTag ISPIN) == 1 ]; then
        echo -n "" "Total" >PDOS/DOS_$gap.dos
    else
        echo -n "" "Total_u Total_d" >PDOS/DOS_$gap.dos
    fi
    species=()
    if [ $(getTag ISPIN) == 1 ]; then
        o=(s p d)
        c1=("%0.5f" "%0.5f")
        x=5
    else
        o=(s_u s_d p_u p_d d_u d_d)
        c1=("%0.5f" "%0.5f" "%0.5f")
        x=7
    fi
    for i in ${!s[*]}; do
        species=(${species[*]} "PDOS/${s[$i]}")
        for j in $(seq 0 $((maxO - 1))); do
            if [[ ! ${s[$i]}_${o[$j]} =~ O_d_[ud] ]]; then
                echo -n " ${s[$i]}_${o[$j]}" >>PDOS/DOS_$gap.dos
                c1=(${c1[*]} "%0.5f")
                echo -n '$' >>col
                echo -n "$x " >>col
            fi
            ((x++))
        done
        ((x++))
    done
    if [ $(getTag ISPIN) == 2 ]; then echo '$3' $(grep '' col) >col; fi
    col=$(grep '' col | sed -e 's/ /,/g' -e 's/,$//')
    rm col
    echo "" >>PDOS/DOS_$gap.dos
    perl -lane 'BEGIN{$"=" "}$F[0]=$F[0]-'"$fermi"';print "@F"' PDOS/total >PDOS/totalFermi
    paste PDOS/totalFermi ${species[*]} | awk '{ printf "'"${c1[*]}"'\n", $1,$2,'"$col"'}' >>PDOS/DOS_$gap.dos
    echo ""
}

intDOS() { # integrate total DOS

    local min max l u lNum uNum Iu Id a b I

    min=$(echo "$(grep '' PDOS/totalFermi | head -1 | cut -d ' ' -f 1) / 1" | bc)
    max=$(echo "$(grep '' PDOS/totalFermi | tail -1 | cut -d ' ' -f 1) / 1" | bc)
    echo -e "\nRange: $min eV - $max eV"
    echo -ne "\nEnter lower bound -> "
    read l
    echo -ne "Enter upper bound -> "
    read u
    lNum=$(grep -n "^$l.[0-9]*" PDOS/totalFermi | tail -1 | cut -d ':' -f 1)
    uNum=$(grep -n "^$u.[0-9]*" PDOS/totalFermi | tail -1 | cut -d ':' -f 1)
    echo $lNum $uNum
    Iu=0
    Id=0
    while [ $lNum -lt $uNum ]; do
        a=($(grep -n '' PDOS/totalFermi | grep "^$lNum:" | cut -d ':' -f 2 | cut -d ' ' -f 1-3))
        b=($(grep -n '' PDOS/totalFermi | grep "^$((lNum + 1)):" | cut -d ':' -f 2 | cut -d ' ' -f 1-3))
        Iu=$(echo "$Iu + $(echo "(${b[0]} - ${a[0]}) * (${a[1]} + ${b[1]}) * 0.5" | bc)" | bc)
        Id=$(echo "$Id + $(echo "(${b[0]} - ${a[0]}) * (${a[2]} + ${b[2]}) * -0.5" | bc)" | bc)
        ((lNum++))
    done
    I=$(echo "$Iu + $Id" | bc)
    printf "\nTotal -> %0.3f \n" $I
    printf "Spin up -> %0.3f\n" $Iu
    printf "Spin down -> %0.3f\n\n" $Id
}

# calculate percent error to experimental
perErr() {
    local err
    err=$(echo "($1 / $2 - 1) * 100" | bc -l | xargs printf %0.1f)
    if (($(echo "$err < 0" | bc -l))); then
        echo "$err%"
    else
        echo "+$err%"
    fi
}

getLat() { # get lattice parameters (any argument to print)
    if [[ $2 ]]; then file=$2; else file=POSCAR; fi
    local i j l a1 a2 a3
    l=3
    for i in a1 a2 a3; do
        for j in 0 1 2; do
            eval "$i[$j]=$(head -$l $file | tail -1 | rmSpace | cut -d ' ' -f $((j + 1)))"
        done
        ((l++))
    done
    mags=($(
        python3.7 - "${a1[*]}" "${a2[*]}" "${a3[*]}" <<-END

            from sys import argv
            import numpy as np

            for i in range(3): print(np.linalg.norm([float(j) for j in argv[i+1].split()]))
END
    ))
    if [[ $1 == disp ]]; then
        printf "%s\n" ${mags[*]}
    elif [[ $1 == pass ]]; then
        echo "${a1[*]}" "${a2[*]}" "${a3[*]}"
    elif [[ $1 == area ]]; then
        echo "${mags[0]} * ${mags[1]}" | bc
    fi
}

# get current POSCAR parameters
storePos() {

    if [ ! -f POSCAR ]; then
        echo -e "\nMissing POSCAR\n"
        return
    fi

    local i s n v elec ion

    s=($(head -6 POSCAR | tail -1))
    n=($(head -7 POSCAR | tail -1))

    if [ $1 ]; then
        case $1 in
        s)
            echo "${s[*]}"
            return
            ;;
        n)
            echo "${n[*]}"
            return
            ;;
        *)
            if [ ! -f POTCAR ]; then
                echo -e "\nMissing POTCAR\n"
                return
            fi
            v=()
            elec=0
            ion=0
            for i in ${!s[*]}; do
                if [[ $(grep PAW_PBE POTCAR) ]]; then
                    #TODO - fix C/Cu
                    v=(${v[*]} $(echo "$(grep -A 1 "^ PAW_PBE ${s[$i]}.*" POTCAR | tail -1 | getNumber) / 1" | bc))
                else
                    v=(${v[*]} $(echo "$(grep -A 1 "^ PAW ${s[$i]}.*" POTCAR | tail -1 | getNumber) / 1" | bc))
                fi
                elec=$(echo "$elec + (${n[$i]} * ${v[$i]})" | bc)
                ion=$(echo "$ion + ${n[$i]}" | bc)
            done
            ;;
        esac
    fi
    case $1 in
    v) echo "${v[*]}" ;;
    e) echo $elec ;;
    i) echo $ion ;;
    esac
}

# get ionic forces
forces() {

    if [ ! -f OUTCAR ]; then
        echo -e "\n -> No OUTCAR found"!"\n"
        return
    fi

    local n ax

    n=$(grep NIONS OUTCAR | rev | cut -d ' ' -f 1 | rev)
    #               sed -n '/TOTAL-FORCE/,/total drift/p' OUTCAR | tail -$((n + 4)) # print Forces as they appear in OUTCAR
    #               return

    if [[ $1 ]]; then ax=$1; else ax=$(index $(getLat disp | sort -n | tail -1) "$(getLat disp)" | head -1); fi
    if [[ ! $(grep TOTAL-FORCE OUTCAR) ]]; then
        echo -e "\nNo forces written yet\n"
        return
    fi
    sed -n '/TOTAL-FORCE/,/total drift/p' OUTCAR | tail -$((n + 4)) | clip b 2 >temp1
    grep -o ' [TF] .*' POSCAR >temp2
    paste temp1 temp2 | rmSpace | sort -k$((ax + 1)) -n | cut -d ' ' -f 1- >temp3 # sorted by longest lattice parameter

    # TODO: allow for mixed T/F cases

    forcesHeader() {

        python3.7 <<-END

        print(' ' + '-' * 64)
        print('%6s %10s %10s %10s %10s %10s' % ('x', 'y', 'z', 'Fx', 'Fy', 'Fz'))
        print(' ' + '-' * 64)

END

    }

    forcesHeader

    python3.7 - $(getTag EDIFFG | cut -d '-' -f 2) <<-END | sort -nr -k3

    from sys import argv
    from colorama import Fore


    def check(range):
        for n in range:
            if abs(float(args[n])) <= cutoff:
                args.append(eval('Fore.GREEN + "1"'))
            else:
                args.append(eval('Fore.RED + "0"'))
        args.append(eval('Fore.WHITE'))


    cutoff = float(argv[1])

    with open('temp3', 'r') as f:
        for line in f.readlines():
            args = line.split()
            l = len(args) - 3
            if 'T' in args or l == 3:
                args.append('--->')
                check(range(3, 6))

            vars = tuple([float(i) for i in args[:6]] + args[6:])
            if len(args) == 11:
                print('%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %s %s %s %s %s' %
                    vars)
            elif len(args) == 9:
                print('%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f  <%s %s %s>' %
                    vars)
            else:
                print(
                    '%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f  <%s %s %s> %s %s %s %s %s'
                    % vars)

END

    forcesHeader

    rm temp[123]
}

# get nearest neighbor table
nearest() {
    local n
    n=$(grep NIONS OUTCAR | rev | cut -d ' ' -f 1 | rev)
    sed -n "/nearest neighbor/,/ $n  /p" OUTCAR | sed '/^$/d'
}

enbTag() { # enable INCAR tag
    local i
    for i in $*; do sed -i "/#$i /s/#//" INCAR; done
}

disTag() { # disable INCAR tag
    local i
    for i in $*; do sed -i "/^$i /s/^ */#/" INCAR; done
}

setTag() { # update an INCAR tag

    local tag tags vals row value path

    if [[ $1 ]]; then
        tags=($1)
        vals=($2)
        path=$3
        set ${path:=./}
        if [[ ${#tags[*]} != ${#vals[*]} ]]; then vals=${vals[*]}; fi
        for tag in ${!tags[*]}; do sed -i "/^#*${tags[$tag]} /s/=.*/= ${vals[$tag]/,/ }/" $path/INCAR; done
    else
        echo ""
        grep -n = INCAR | cut -d '=' -f 1 | sed 's/#//'
        echo ""
        echo -n "Which tag? (enter numeric value) -> "
        read row
        echo -n "Enter new value -> "
        read value
        tag=$(head -$row INCAR | tail -1 | cut -d '=' -f 1)
        setTag $tag $value
        echo ""
    fi
}

getTag() { # get INCAR tag value
    local file
    if [ $2 ]; then file=$2; else file=INCAR; fi
    echo -n $(cat $file | rmSpace | grep "^#*$1 " | cut -d ' ' -f 3)
}

showTag() { # write INCAR tag to file
    if [ $1 == ISMEAR ] && [ $(getTag $1) != -5 ]; then
        echo -n "$1 = $(getTag $1) ($(getTag SIGMA))"
    else echo -n "$1 = $(getTag $1)"; fi
    if [[ ! $2 ]]; then echo ""; else echo -n " $2 "; fi
}

defEn() { # get default ENCUT from POTCAR
    grep ENMAX POTCAR | cut -d ' ' -f 8 | cut -c 1-7 | sort -nr | head -1 | xargs printf '%0.f'
}

kPoints() { # set KPOINTS
    if [ ! "$*" ]; then
        echo "Usage: kPoints <k1> <k2> <k3> -> k1 only to calculate k2 & k3 according to lattice"
        return
    elif [ ! $2 ] || [ ! $3 ]; then

        local a1 a2 a3 k

        a1=($(getLat pass | cut -d ' ' -f 1-3))
        a2=($(getLat pass | cut -d ' ' -f 4-6))
        a3=($(getLat pass | cut -d ' ' -f 7-9))

        if [[ $2 == even ]]; then odd=False; else odd=True; fi

        k=($(

            python3.7 - $1 "${a1[*]}" "${a2[*]}" "${a3[*]}" $odd <<-END

            from sys import argv
            import numpy as np

            odd = argv[-1]

            k1 = int(argv[1])
            k1 = k1 + 1 if odd == 'True' and k1 % 2 == 0 else k1

            # define real lattice
            a = [np.array([float(j) for j in i.split()]) for i in argv[2:5]]
            realVecs = [np.linalg.norm(r) for r in a]

            # define reciprocal lattice
            b = [(np.dot(a[0], np.cross(a[1], a[2])))**-1 * np.cross(i, j) for i, j in [a[1:3], [a[2], a[0]], a[0:2]]]
            recVecs = [np.linalg.norm(r) for r in b]

            # define ratios with respect to shortest vector
            longIndex = recVecs.index(max(recVecs))
            r21 = recVecs[(longIndex + 1) % len(recVecs)] / recVecs[longIndex]
            r31 = recVecs[(longIndex + 2) % len(recVecs)] / recVecs[longIndex]

            k2 = int(k1 * r21)
            k3 = int(k1 * r31)

            k2, k3 = [i + 1 if odd == 'True' and i % 2 == 0 else i for i in [k2, k3]]

            k = []
            for i in range(3):
            if i == longIndex:
                k.append(k1)
            elif i == (longIndex + 1) % len(recVecs):
                k.append(k2)
            else:
                k.append(k3)

            print(k[0], k[1], k[2])

END

        ))

    else
        k=($*)
    fi
    if [[ $(inArray 0 "${k[*]}") ]]; then
        echo -e "\n -> Null index detected. Please use larger k value.\n"
    else
        sed -i "4s/.*/${k[0]}  ${k[1]}  ${k[2]}/" KPOINTS
    fi
}

getKps() { # get KPOINTS
    echo "$(head -4 KPOINTS | tail -1 | rmSpace | cut -d ' ' -f 1-3)"
}

setKps() {
    local k
    echo -n "Kpoints (enter kx only to extrapolate) -> "
    read k
    if [[ ! $k =~ ^[0-9]+$ ]]; then k=1; fi
    kPoints $k
}

getTime() { # get elapsed time of last VASP run

    local raw t_days days t_hr hr r_min min t_sec sec

    raw=$(grep Elapsed OUTCAR | rev | cut -d ' ' -f 1 | rev)
    t_days=$(echo "$raw / 86400" | bc -l)
    days=$(echo "$t_days / 1" | bc)
    t_hr=$(echo "($t_days - $days) * 24" | bc)
    hr=$(echo "$t_hr / 1" | bc)
    t_min=$(echo "($t_hr - $hr) * 60" | bc)
    min=$(echo "$t_min / 1" | bc)
    t_sec=$(echo "($t_min - $min) * 60" | bc)
    sec=$(printf %0.f "$t_sec")
    if ((sec == 60)); then
        sec=0
        ((min++))
    fi
    if ((min == 60)); then
        min=0
        ((hr++))
    fi
    if ((hr == 24)); then
        hr=0
        ((days++))
    fi
    if [ "$(echo "$sec < 10" | bc)" == 1 ]; then sec="0$sec"; fi
    if [ "$(echo "$min < 10" | bc)" == 1 ]; then min="0$min"; fi
    if [ "$(echo "$hr < 10" | bc)" == 1 ]; then hr="0$hr"; fi
    echo "$(if ((days != 0)); then echo -n "($days) "; fi)$hr:$min:$sec"
}

getVolume() { # get current cell volume

    local vol

    if [ -f OUTCAR ]; then
        vol=$(grep volume OUTCAR | tail -1 | getNumber)
    else
        vol=$(math triProd $(getLat pass) | xargs printf %0.2f)
    fi

    printf "VOLUME = %0.2f" $vol

}

getEnergy() { # display energy (eV) with calculation parameters

    local algo cycl smear volume energy occ ediff

    occ=$(occCheck 1)
    file=OSZICAR
    algo=$(grep -B 1 ' F= ' $file | tail -2 | head -1 | cut -c 1-3)
    cycl=$(grep -B 1 ' F= ' $file | tail -2 | head -1 | rmSpace | cut -d ' ' -f 2)
    smear=$(getTag ISMEAR)
    volume=$(getVolume | getNumber)
    energy=$(grep "free  energy" OUTCAR | tail -1 | getNumber | xargs printf %0.5f)
    ediff=$(grep -B 1 E0 $file | tail -2 | head -1 | rmSpace | cut -d ' ' -f 4 | grep -o E.*)
    echo -ne "E0 = $energy; 1$ediff; ($(getKps)); $(getTag ENCUT); $(getTag NBANDS); $cycl ($algo)"
    echo -n "; $smear$(if [ $smear != -5 ]; then echo -n " ($(getTag SIGMA | xargs printf %0.3f))"; fi)"$occ"; $(getTag PREC)"
    echo -n "$(if [[ $1 ]]; then echo "; R$4$5 $2($3): $1 ($(getTag EDIFFG)); $volume"; fi)"
    echo "; $(getTime)"
}

realSpace() { # write real space matrix
    echo ""
    local a1 a2 a3
    a1=($(getLat pass | cut -d ' ' -f 1-3))
    a2=($(getLat pass | cut -d ' ' -f 4-6))
    a3=($(getLat pass | cut -d ' ' -f 7-9))
    printf "( %0.3f  %0.3f  %0.3f )\n\n" ${a1[0]} ${a1[1]} ${a1[2]}
    printf "( %0.3f  %0.3f  %0.3f )\n\n" ${a2[0]} ${a2[1]} ${a2[2]}
    printf "( %0.3f  %0.3f  %0.3f ) -> $(getVolume)\n\n" ${a3[0]} ${a3[1]} ${a3[2]}
}

abs() { # Usage: abs <value 1> <value 2>
    if [ "$(echo "$1 - $2 >= 0" | bc)" == 1 ]; then
        echo "($1 - $2) * 1" | bc
    else
        echo "($1 - $2) * -1" | bc
    fi
}

getGap() { # Usage: getGap <vbm> <cbm> <path> -> default path = current dir
    local vbm cbm path vbMax cbMin s
    if [ $1 ]; then vbm=$1; else vbm=$(vbm i); fi
    cbm=$((vbm + 1))
    if [ $3 ]; then path=$3; else path=.; fi
    if [ $(getTag ISPIN) == 1 ]; then s=2; else s=3; fi
    vbMax=$(sed '1,7d' $path/EIGENVAL | grep " $vbm " | rmSpace | cut -d ' ' -f $s | sort -n | tail -1)
    cbMin=$(sed '1,7d' $path/EIGENVAL | grep " $cbm " | rmSpace | cut -d ' ' -f $s | sort -rn | tail -1)
    printf %0.3f "$(echo "$cbMin - $vbMax" | bc)"
    if [ ! $4 ]; then echo ""; fi
}

# modify POSCAR #####################################################################################################
posMod() {

    local dir def vac s n s1 i el loc defect x k l ind c

    getElement() {
        s1=()
        for i in ${!s[*]}; do s1+=("$i-${s[$i]}"); done
        echo -n "Select element to modify ( $(echo "${s1[*]}" | sed 's/ / | /g') ) -> "
        read el
        if [ ! $(inArray $el "${!s[*]}") ]; then
            echo 'Element does not exist in system!'
            getElement
        fi
    }

    getLocation() {
        echo -n "Select location (1-${n[$el]}) -> "
        read -a loc
        for i in ${loc[*]}; do
            if [ ! $(inArray $i "$(range $((${n[$el]} + 1)) 1)") ]; then
                echo 'Not a valid location!'
                getLocation
            fi
        done
    }

    getReplace() {
        echo -n "Choose dopent element (leave blank for vacancy) -> "
        read defect
        if [[ $defect ]] && [[ ! $(grep "^ PAW_PBE $defect[_ ]" $dir) ]]; then
            echo 'Not a valid defect!'
            getReplace
        fi
    }

    runMod() {

        s=($(storePos s))
        n=($(storePos n))
        getElement
        getLocation
        getReplace

        for k in ${loc[*]}; do
            l=$((k + 8 - x))
            for i in $(seq 0 $el); do l=$((l + ${n[$i]})); done # determine line of selected defect
            if [ ! $defect ]; then
                sed -i -e "7s/[0-9][0-9]*/$((${n[$el]} - 1))/$((el + 1))" POSCAR # -e "$l"d POSCAR
            else
                if [ ! $(inArray $defect "${s[*]}") ]; then
                    sed -i -e "6s/.*/&   $defect/" -e "7s/.*/&    1/" POSCAR
                else
                    ind=$(index $defect "${s[*]}")
                    sed -i "7s/[0-9][0-9]*/$((${n[$ind]} + 1))/$((ind + 1))" POSCAR
                fi
            fi
        done

        head -$l POSCAR | tail -1 >>POSCAR # move coordinate to end of file
        sed -i "${l}d" POSCAR              # delete old coordinate
        return
        sed -i "6s/\(.*\)/$(tail -n +6 POSCAR | head -1 | xargs printf %4s%4s%4s%4s)/" POSCAR
        sed -i "7s/\(.*\)/$(tail -n +7 POSCAR | head -1 | xargs printf %4s%4s%4s%4s)/" POSCAR
    }

    # execute modifications
    if [ ! -d tempPos ]; then cp POSCAR tempPos; fi
    sed -i -e '1s/  */ /g' -e '1s/ $//' POSCAR # reformat POSCAR title
    blank=$(grep -n '^[\t ]*$' POSCAR | head -1 | cut -d ':' -f 1)
    if [[ $blank ]]; then sed -i "$blank,\$d" POSCAR; fi
    dir=/work/01293/hudamn/Edan/PAW_PBE

    echo ""
    runMod
    echo ""
    head -7 POSCAR
    echo ""

    # prompt user for subsequent modification
    c=1
    while [ $c != 0 ]; do
        echo -n "Enter 1 for additional modification; 0 to exit -> "
        read c
        if [ $c != 1 ] && [ $c != 0 ]; then echo 'Not a valid selection!'; fi
        if [ $c == 1 ]; then
            runMod
            echo ""
            head -7 POSCAR
            echo ""
        else
            # getFiles
            # getJob IN
            echo ""
        fi
    done
}
#####################################################################################################################

parseB() {

    # set fermi level
    vbm=$(vbm i)
    fermi=$(vbm eV)

    # get paths
    echo -n "Enter k-path -> "
    read path
    if [ $path ]; then
        paths=($path)
    else
        paths=($(grep '' paths))
    fi
    # get number of bands
    b=$(head -6 EIGENVAL | tail -1 | sed 's/  */ /g' | cut -d ' ' -f 4)

    # parse EIGENVAL ####################################################################################################
    n=0
    while [ $n -lt ${#paths[*]} ]; do
        # set k-points path
        p=($(echo ${paths[$n]} | sed 's/-/ /g'))
        # reformat path name
        path_name=$(echo ${paths[$n]})
        if [ ! -d bands/$path_name ]; then mkdir bands/$path_name; fi
        # initialize record
        if [ $(getTag ISPIN) == 1 ]; then
            echo -n "" >bands/$path_name/BANDS
        else
            echo -n "" >bands/$path_name/BANDS_u
            echo -n "" >bands/$path_name/BANDS_d
        fi
        # collect eigenvalues from all path segments for current path
        k1=0
        k2=1
        m=1
        while [ $k1 -lt $(echo "${#p[*]} - 1" | bc) ]; do
            eigVal=$(find bands/*/${p[$k1]}-${p[$k2]}/ -name EIGENVAL)
            pattern="^ *[0-9][0-9]*  *[-0-9][0-9]*\.[0-9]"
            kps=$(head -2 $(find bands/*/${p[$k1]}-${p[$k2]}/ -name KPOINTS) | tail -1)
            if [ $(getTag ISPIN) == 1 ]; then
                # record band energies for all k-points
                i=1
                x=1
                while [ $i -le $kps ]; do
                    a="$(printf "%0.4f " $(grep "$pattern" $eigVal | head -$((x * $b)) | tail -$b | sed 's/  */ /g' | cut -d ' ' -f 3))"
                    echo $m $i $a "" >>bands/$path_name/BANDS
                    ((x++))
                    ((i++))
                    ((m++))
                done
                ((k1++))
                ((k2++))
            else
                # record band energies for all k-points
                i=1
                x=1
                while [ $i -le $kps ]; do
                    a1="$(printf "%0.4f " $(grep "$pattern" $eigVal | head -$((x * $b)) | tail -$b | sed 's/  */ /g' | cut -d ' ' -f 3))"
                    a2="$(printf "%0.4f " $(grep "$pattern" $eigVal | head -$((x * $b)) | tail -$b | sed 's/  */ /g' | cut -d ' ' -f 4))"
                    echo $m $i $a1 >>bands/$path_name/BANDS_u
                    echo $m $i $a2 >>bands/$path_name/BANDS_d
                    ((x++))
                    ((i++))
                    ((m++))
                done
                ((k1++))
                ((k2++))
            fi
        done
        if [ $(getTag ISPIN) == 1 ]; then
            echo -n "$path_name 100" >>bands/$path_name/BANDS
        else
            echo -n "$path_name 100" >>bands/$path_name/BANDS_u
            echo -n "$path_name 100" >>bands/$path_name/BANDS_d
        fi
        ((n++))
    done

    # combine results for plotting ######################################################################################
    # determine spacing
    i=0
    tot_kps=0
    while [ $i -lt ${#paths[*]} ]; do
        if [ $(getTag ISPIN) == 1 ]; then
            tot_kps=$(echo "$tot_kps + $(grep '^[0-9]' bands/${paths[$i]}/BANDS | wc -l)" | bc)
        else
            tot_kps=$(echo "$tot_kps + $(grep '^[0-9]' bands/${paths[$i]}/BANDS_u | wc -l)" | bc)
        fi
        ((i++))
    done
    ratio=$(printf %0.f $(echo "(0.05 * $tot_kps)" | bc))
    i=1
    while [ $i -lt $ratio ]; do
        echo -n "\\n" >>file
        ((i++))
    done
    spacing=$(grep '' file)
    rm file
    # concatenate bands
    i=0
    if [ $(getTag ISPIN) == 1 ]; then
        while [ $i -lt ${#paths[*]} ]; do
            grep '' bands/${paths[$i]}/BANDS >>bands/BANDS
            if [ $i -lt $(echo "${#paths[*]} - 1" | bc) ]; then echo -e "$spacing" >>bands/BANDS; fi
            ((i++))
        done
        # set Fermi level at 0 eV
        perl -pe 's/([-0-9][0-9]*\.[0-9][0-9]*)/$1-'"$fermi"'/eg' bands/BANDS >bands/Fermi
        sed -i 's/\( [-0-9][0-9]*\) /\1\.0000 /g' bands/Fermi
        sed -i 's/\( [-0-9][0-9]*\.[0-9]\) /\1000 /g' bands/Fermi
        sed -i 's/\( [-0-9][0-9]*\.[0-9][0-9]\) /\100 /g' bands/Fermi
        # calculate band gap
        cbm=$((vbm + 1))
        cbMin=$(printf %0.4f $(grep "^[0-9]" bands/Fermi | cut -d ' ' -f $((cbm + 2)) | sort -nr | tail -1))
        gap=$(printf %0.3f $cbMin)
        # determine total number of rows in BANDS
        rows=$(wc -l bands/BANDS | cut -d ' ' -f 1)
        # record VBM/CBM locations
        kVBM=($(grep -n ' 0\.0000 ' bands/Fermi | cut -d ':' -f 1))
        kCBM=($(grep -n " $cbMin " bands/Fermi | cut -d ':' -f 1))
        for i in ${!kVBM[*]}; do
            if [ ${kVBM[$i]} -gt $rows ]; then kVBM[$i]="$(echo ${kVBM[$i]} - $rows | bc)(d)"; fi
        done
        for i in ${!kCBM[*]}; do
            if [ ${kCBM[$i]} -gt $rows ]; then kCBM[$i]="$(echo ${kCBM[$i]} - $rows | bc)(d)"; fi
        done
        echo -e "\nVBM at ${kVBM[*]}\n\nCBM at ${kCBM[*]}" >>status
        # rename BANDS
        mv bands/Fermi bands/Fermi_"$rows"_"$gap".bands
    else
        while [ $i -lt ${#paths[*]} ]; do
            grep '' bands/${paths[$i]}/BANDS_u >>bands/BANDS_u
            grep '' bands/${paths[$i]}/BANDS_d >>bands/BANDS_d
            if [ $i -lt $(echo "${#paths[*]} - 1" | bc) ]; then
                echo -e "$spacing" >>bands/BANDS_u
                echo -e "$spacing" >>bands/BANDS_d
            fi
            ((i++))
        done
        # set Fermi level at 0 eV
        perl -pe 's/([-0-9][0-9]*\.[0-9][0-9]*)/$1-'"$fermi"'/eg' bands/BANDS_u >bands/Fermi_u
        perl -pe 's/([-0-9][0-9]*\.[0-9][0-9]*)/$1-'"$fermi"'/eg' bands/BANDS_d >bands/Fermi_d
        sed -i 's/\( [-0-9][0-9]*\) /\1\.0000 /g' bands/Fermi*
        sed -i 's/\( [-0-9][0-9]*\.[0-9]\) /\1000 /g' bands/Fermi*
        sed -i 's/\( [-0-9][0-9]*\.[0-9][0-9]\) /\100 /g' bands/Fermi*
        # calculate band gap
        cbm=$((vbm + 1))
        cbMin=$(printf %0.4f $(cat bands/Fermi_u bands/Fermi_d | grep "^[0-9]" | cut -d ' ' -f $((cbm + 2)) | sort -nr | tail -1))
        gap=$(printf %0.3f $cbMin)
        # determine total number of rows in BANDS
        rows_u=$(wc -l bands/Fermi_u | cut -d ' ' -f 1)
        rows_d=$(wc -l bands/Fermi_d | cut -d ' ' -f 1)
        # record VBM/CBM locations
        kVBM=($(cat bands/Fermi_u bands/Fermi_d | grep -n ' 0\.0000 ' | cut -d ':' -f 1))
        kCBM=($(cat bands/Fermi_u bands/Fermi_d | grep -n " $cbMin " | cut -d ':' -f 1))
        for i in ${!kVBM[*]}; do
            if [ ${kVBM[$i]} -gt $rows_u ]; then kVBM[$i]="$(echo ${kVBM[$i]} - $rows_u | bc)(d)"; fi
        done
        for i in ${!kCBM[*]}; do
            if [ ${kCBM[$i]} -gt $rows_u ]; then kCBM[$i]="$(echo ${kCBM[$i]} - $rows_u | bc)(d)"; fi
        done
        echo -e "\nVBM at ${kVBM[*]}\n\nCBM at ${kCBM[*]}" >>status
        # rename BANDS
        mv bands/Fermi_u bands/Fermi_u_"$rows_u"_"$gap".bands
        mv bands/Fermi_d bands/Fermi_d_"$rows_d"_"$gap".bands
    fi

}
