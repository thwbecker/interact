for i in src/*.cpp; do
    printf "src/"
    g++ -I. -MM $i
done > make.depends
