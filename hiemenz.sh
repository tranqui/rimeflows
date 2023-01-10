for Re in 1 2 5 10 50 100 500 2000 5000 10000 100000 1000000 10000000 100000000; do
    ./efficiency.py $Re -f hiemenz -s 1e-5 1 100 > data/efficiency_hiemenz_Re=$Re.csv &
done

./efficiency.py -f shm -s 1e-3 1 100 --maxstep 1e-2 > data/efficiency_shm.csv &
