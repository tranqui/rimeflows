for m in 1.05 1.10 1.15 1.20 1.25 1.30 1.50 2.00; do
    ./efficiency.py $m -v -f power -s 1e-7 1 100 --maxstep 1e-2 --tmax 1e2 > data/efficiency_power_m=$m.csv &
done

for Re in 1 2 5 10 50 100 500 2000 5000 10000 100000 1000000 10000000 100000000; do
    ./efficiency.py $Re -f hiemenz -s 1e-5 1 100 > data/efficiency_hiemenz_Re=$Re.csv &
done

./efficiency.py -f shm -s 1e-3 1 100 --maxstep 1e-2 > data/efficiency_shm.csv &