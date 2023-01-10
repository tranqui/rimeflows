for Re in 1 2 5 10 50 100 500 2000 5000 10000 100000 1000000 10000000 100000000; do
    ReOut=$(printf "%04d" $Re)
    ./efficiency.py $Re -f hiemenz -s 1e-5 1 100 --niters 100 --tmax 1e6 > data/efficiency_hiemenz_Re=$ReOut.csv &
done
