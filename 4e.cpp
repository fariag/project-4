/* 
Program to solve the two-dimensional Ising model 
The coupling constant J = 1
Boltzmann's constant = 1
*/

#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <random>
#include <armadillo>
#include <string>
#include <mpi.h>

using namespace  std;
using namespace arma;
// output file
ofstream outfile;

// inline function for Periodic boundary conditions
inline int Periodic(int i, int limit, int add) { 
    return (i+limit+add) % (limit);
}

// function to initialize energy, spin matrix and magnetization
void Initialize(int n_spin, mat &spin_mat,  double& energy, double& magmoment) {
    // setup spin matrix and initial magnetization

    // Ordered spin orientation
    for(int x =0; x < n_spin; x++) {
        for (int y= 0; y < n_spin; y++){
            spin_mat(x,y) = 1.0; // spin orientation for the ground state
            magmoment +=  1.0;
        }
    }

    // Random spin orientation
    /*random_device get_rand;
    mt19937_64 gen(get_rand());
    uniform_real_distribution<double> dist(0.0, 1.0);
    for(int x =0; x < n_spin; x++) {
        for (int y= 0; y < n_spin; y++){
            if (dist(gen) < 0.5) {
                spin_mat(x,y) = 1.0; // spin orientation for the ground state
                magmoment +=  1.0;
            }
        }
    }*/

    // setup initial energy
    for(int x =0; x < n_spin; x++) {
        for (int y= 0; y < n_spin; y++){
            energy -=  (double) spin_mat(x,y)*
                       (spin_mat(Periodic(x,n_spin,-1),y) +
                        spin_mat(x,Periodic(y,n_spin,-1)));
        }
    }
}

void Metropolis(int n_spin, int mcs, double temp, vec &expect, int my_rank, int num_procs) {

    // Initialize needed variables 
    random_device get_rand;
    mt19937_64 gen(get_rand());

    uniform_real_distribution<double> dist(0.0, 1.0);

    double energy = 0.;
    double magmoment = 0.;
    int accept_count = 0;

    mat spin_mat = zeros<mat>(n_spin,n_spin);

    // Find all the different energy values
    vec energydiff = zeros<mat>(17);
    for( int de =-8; de <= 8; de+=4) {
        energydiff(de+8) = exp(-de/temp);
    }

    Initialize(n_spin, spin_mat, energy, magmoment);

    MPI_Status status;

    int my_start = (((double) mcs) / ((double) num_procs))*my_rank + 1;
    int my_finish = (((double) mcs) / ((double) num_procs))*(my_rank + 1);

    for (int i = my_start; i <= my_finish; i++){
        // accept_count = 0;

        // Looping over all the spins in the lattice
        for(int x =0; x < n_spin; x++) {
            for (int y= 0; y < n_spin; y++){
                int ix = (int) (dist(gen)*(double)n_spin);
                int iy = (int) (dist(gen)*(double)n_spin);
                
                int delta = 2*spin_mat(ix,iy)*
                            (spin_mat(ix,Periodic(iy,n_spin,-1))
                            + spin_mat(Periodic(ix,n_spin,-1),iy)
                            + spin_mat(ix,Periodic(iy,n_spin,1)) 
                            + spin_mat(Periodic(ix,n_spin,1),iy));
                
                if ( dist(gen) <= energydiff(delta+8) ) {
                    spin_mat(ix,iy) *= -1.0;                  // Randomly flip on of the spins
                    magmoment += (double) 2*spin_mat(ix,iy);
                    energy += (double) delta;
                    accept_count++;
                }
            }
        }

        expect(0) += energy;
        expect(1) += energy*energy;
        expect(2) += magmoment;    
        expect(3) += magmoment*magmoment; 
        expect(4) += fabs(magmoment);
    }

    // Record expectation values
    if (my_rank == 0) { 
        /* Assemble the whole picture in process 0 */
        double temp = 0;
        double temp1 = 0;
        double temp2 = 0;
        double temp3 = 0;
        double temp4 = 0;
        double temp5 = 0;

        for (int i=1; i<num_procs; i++) {
            MPI_Recv(&temp, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
            expect(0) += temp;

            MPI_Recv(&temp1, 1, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
            expect(1) += temp1;
            
            MPI_Recv(&temp2, 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
            expect(2) += temp2;

            MPI_Recv(&temp3, 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &status);
            expect(3) += temp3;

            MPI_Recv(&temp4, 1, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &status);
            expect(4) += temp4;
        }
    }

    if (my_rank != 0) {
        /* Each process except 0 sends each piece of the whole image to process 0 */
        MPI_Send(&expect(0), 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
        MPI_Send(&expect(1), 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        MPI_Send(&expect(2), 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
        MPI_Send(&expect(3), 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
        MPI_Send(&expect(4), 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
    }
}

void output(int n_spin, int mcs, double temperature, vec expect) {
  // Prints all the results to a text file 
    double norm = 1.0/((double) (mcs));
    double E_expect = expect(0)*norm;
    double E2_expect = expect(1)*norm;
    double M_expect = expect(2)*norm;
    double M2_expect = expect(3)*norm;
    double Mabs_expect = expect(4)*norm;

    double Evariance = (E2_expect- E_expect*E_expect)/n_spin/n_spin;          // Find the energy variance
    double Mvariance = (M2_expect - Mabs_expect*Mabs_expect)/n_spin/n_spin;   // Find the magnetic variance
    double spec_heat = E_expect/(temperature*temperature);
    outfile << setiosflags(ios::showpoint | ios::uppercase);
    outfile << setw(20) << "temp";
    outfile << setw(20) << "Energy Expectation";
    outfile << setw(20) << "Specific Heat";
    outfile << setw(20) << "Mag Expectation";
    outfile << setw(20) << "Mag Variance";  
    outfile << setw(20) << "|Mag Expectation|";
    outfile << endl;
    outfile << setw(20) << setprecision(8) << temperature;
    outfile << setw(20) << setprecision(8) << E_expect/n_spin/n_spin;
    outfile << setw(20) << setprecision(8) << Evariance/temperature/temperature;
    outfile << setw(20) << setprecision(8) << M_expect/n_spin/n_spin;
    outfile << setw(20) << setprecision(8) << Mvariance/temperature;
    outfile << setw(20) << setprecision(8) << Mabs_expect/n_spin/n_spin << endl;
}

int main(int argc, char* argv[]) {

    if (argc <= 5) {
        cout << "error: " << argv[0] << 
          " Require name of outputfile, amount of spins, Number of MC cycles, intial and final temperature, and Steps between temperatures" << endl;
        exit(1);
    }

    // Take arguments from the command line
    string filename=argv[1];
    double temp0 = atof(argv[2]);
    double tempf = atof(argv[3]);
    double temp = atof(argv[4]);
    int n_spin = atoi(argv[5]);
    int mcs = atoi(argv[6]);    

    int my_rank, num_procs;
    
    // Create/Replace text file for recording results
    string fileout = filename;
    string argument = to_string(n_spin);
    fileout.append("_Spin=" +argument);
    fileout.append(".txt");
    outfile.open(fileout);

    /*outfile << setw(20) << setprecision(8) << mcs;
    outfile << endl;*/

    // Begin simulation
    MPI_Init(&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size (MPI_COMM_WORLD, &num_procs);
    for (double t = temp0; t <= tempf; t+=temp){
        vec expect = zeros<mat>(5);
        Metropolis(n_spin, mcs, t, expect, my_rank, num_procs);
        /*output(n_spin, mcs, t, expect);
        outfile << endl;*/
    }
    MPI_Finalize();
    outfile.close();
    return 0;
}
