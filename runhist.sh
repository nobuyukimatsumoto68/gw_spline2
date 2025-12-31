#!/bin/bash
set -e


Ns=32
npts=40

# mass=0p4000
# mass=0p3000
# mass=0p2000
# mass=0p1000
declare -a masses=(0p1000 0p2000 0p3000 0p4000)

########################

# Ns=24
# npts=32

# # mass=0p4000
# # mass=0p3000
# # mass=0p2000
# # mass=0p1000
# declare -a masses=(0p1000 0p2000 0p3000 0p4000)

########################

# Ns=16
# npts=24

# declare -a masses=(0.1 0.2 0.3 0.4)


ibetamin=0
ibetamax=1000

for mass in "${masses[@]}"
do
    ./get_potential.o $mass $ibetamin $ibetamax $Ns $npts
done

