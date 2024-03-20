#!/usr/bin/env bash

mlist=$(python -c "import numpy as np; print(' '.join(map(str, np.geomspace(1, 2, 50))))")
echo "" > data/Stc_power_m.csv
for m in $mlist; do
    echo $m
    ./critical-stokes.py $m -f power --maxstep 1e-2 --tmax=1e2 >> data/Stc_power_m.csv
done

mlist=$(python -c "import numpy as np; print(' '.join(map(str, np.geomspace(1, 2, 50))))")
echo "" > data/Stc_power_m_x2.csv
for m in $mlist; do
    echo $m
    ./critical-stokes.py $m -f power --maxstep 1e-2 --tmax=1e2 -x 2 >> data/Stc_power_m_x2.csv
done

mlist=$(python -c "import numpy as np; print(' '.join(map(str, np.geomspace(1, 2, 50))))")
echo "" > data/Stc_power_m_x05.csv
for m in $mlist; do
    echo $m
    ./critical-stokes.py $m -f power --maxstep 1e-2 --tmax=1e2 -x 0.5 >> data/Stc_power_m_x05.csv
done

Relist=$(python -c "import numpy as np; print(' '.join(map(str, np.geomspace(1, 1e8, 50))))")
echo "" > data/Stc_hiemenz.csv
for Re in $Relist; do
    echo $Re
    ./critical-stokes.py $Re -f hiemenz >> data/Stc_hiemenz.csv
done

Relist=$(python -c "import numpy as np; print(' '.join(map(str, np.geomspace(1, 1e8, 50))))")
echo "" > data/Stc_hiemenz_x2.csv
for Re in $Relist; do
    echo $Re
    ./critical-stokes.py $Re -f hiemenz -x 2 >> data/Stc_hiemenz_x2.csv
done

Relist=$(python -c "import numpy as np; print(' '.join(map(str, np.geomspace(1, 1e8, 50))))")
echo "" > data/Stc_hiemenz_x05.csv
for Re in $Relist; do
    echo $Re
    ./critical-stokes.py $Re -f hiemenz -x 0.5 >> data/Stc_hiemenz_x05.csv
done

for m in 1.05 1.10 1.15 1.20 1.25 1.30 1.50 2.00; do
    ./efficiency.py $m -v -f power -s 1e-7 1 100 --maxstep 1e-2 --tmax 1e2 > data/efficiency_power_m=$m.csv &
done
wait

for Re in 1 2 5 10 50 100 500 1000 2000 5000 10000 100000 1000000 10000000 100000000; do
    ./efficiency.py $Re -v -f hiemenz -s 1e-5 1 100 > data/efficiency_hiemenz_Re=$Re.csv &
done
wait

./efficiency.py 0.15 -v -f kuwabara -s 1e-5 1 100 > data/efficiency_kuwabara_alpha=0.15_2.csv &

./efficiency.py -f shm -s 1e-3 1 100 --maxstep 1e-2 > data/efficiency_shm.csv &

for x0 in $(seq -0.1 -0.1 -1); do
    ./efficiency.py 1 -f chord --yguess 1e-4 -x $x0 --maxstep 1e-2 --tmax 1e2 --stokes 1e-5 1 50 > data/efficiency_chord_x0=$x0.csv &
done

path = data/Stc_chord_x0.csv
rm -f $path
for x0 in $(seq -0.1 -0.1 -2); do
    read dummy Stc <<< $(./critical-stokes.py 1 -f chord --maxstep 1e-2 --tmax=1e2 -x $x0)
    echo $x0 $Stc >> $path
done

wait