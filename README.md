# project-4
Compile and link as g++ -O3 4b.cpp -o 4b.exe -larmadillo -std=c++11 ./executable [Outputfile] [initial temp] [final temp] [tempstep] [#spins] {#MC cycles] ./4b.exe output 1 1 1 2 10000

Compile and link parallel version as mpic++ -O3 4e.cpp -o 4e.exe -larmadillo -std=c++11

mpirun -np [Number of processes] [Program] [Outputfile] [initial temp] [final temp] [tempstep] [#spins] {#MC cycles]

mpirun -np 4 ./4e.exe Poutput 2.4 2.4 0.1 20 1000000
