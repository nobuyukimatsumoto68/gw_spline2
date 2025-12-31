#/!/bin/bash
set -e


Ns=32
npts=40
# declare -a nbetas=(21 9 8 8)
# declare -a masses=(0p1000 0p2000 0p3000 0p4000)
declare -a nbetas=(21)
declare -a masses=(0p1000)
# declare -a nbetas=(9 8 8)
# declare -a masses=(0p2000 0p3000 0p4000)

########################

# Ns=24
# npts=32
# declare -a nbetas=(13 10 9 8)
# declare -a masses=(0p1000 0p2000 0p3000 0p4000)

########################

# Ns=16
# npts=24
# declare -a nbetas=(14 12 6 6)
# declare -a masses=(0.1 0.2 0.3 0.4)


ibetamin=0
ibetamax=1000

for (( i=0; i<"${#masses[@]}"; i++))
do
    mass=${masses[$i]}
    nbeta=${nbetas[$i]}
    ./get_potential_jk.o $mass $ibetamin $ibetamax $Ns $npts $nbeta
done

