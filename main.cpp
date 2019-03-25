#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <time.h>
#include <fstream>

using namespace std;

int main ( int argc, char *argv[] );
pair<int, int> random_walks ( int a, int b, int x, double p, int N );

//****************************************************************************80

int main ( int argc, char *argv[] )
{
    int a = 0;
    int b = 1000;
    int x = 300;
    pair<int, int> hits;
    int hit_total = 0;
    int avg_live = 0;
    int ierr;
    double p = 0.6;
    const int master = 0;
    double pdf_estimate;
    int process_num;
    int process_rank;
    int N = 100000;
    int trial_total;
    int seed;

    clock_t tStart = clock();

    ierr = MPI_Init(&argc, &argv);

    if (ierr != 0) {
        exit(1);
    }

    ierr = MPI_Comm_size(MPI_COMM_WORLD, &process_num);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);

    seed = 123456789 + process_rank * 100;
    srand(seed);

    if (process_rank == 0) {
        cout << "The number of processes is " << process_num << endl;
        cout << endl;
        cout << "A = " << a << endl;
        cout << "B = " << b << endl;
        cout << "Start position X = " << x << endl;
        cout << "Probability to go right P = " << p << endl;
    }

        N = 1000;

        hits = random_walks(a, b, x, p, N);

        ierr = MPI_Reduce(&hits.first, &hit_total, 1, MPI_INT, MPI_SUM, master,
                          MPI_COMM_WORLD);
        ierr = MPI_Reduce(&hits.second, &avg_live, 1, MPI_INT, MPI_SUM, master,
                          MPI_COMM_WORLD);

        if (process_rank == 0) {
            trial_total = N * process_num;
            avg_live /= trial_total;
            pdf_estimate = (double) (hit_total) / (double) (trial_total);


            cout << "  " << setw(8) << trial_total
                 << "  " << setw(8) << hit_total
                 << "  " << setw(8) << pdf_estimate
                 << "  " << setw(8) << avg_live;
        }

        if (process_rank == 0) {
            cout << endl;
            printf("Time taken: %.2fs\n", (double) (clock() - tStart) / CLOCKS_PER_SEC);

            ofstream myfile;
            myfile.open("output.txt");
            myfile << (double) (hit_total) / (double) (N * process_num) << ' ' << avg_live << endl;
            myfile.close();

            myfile.open("stat.txt");
            myfile << a << ' ' << b << endl;
            myfile << x << endl;
            myfile << p << endl;
            myfile << N * process_num << endl;
            myfile << (double) (clock() - tStart) / CLOCKS_PER_SEC << endl;
            myfile << process_num << endl;
            myfile.close();

            myfile.open("time_from_p_n=1000.txt", std::ios_base::app);
            myfile << process_num << ' ' << (double) (clock() - tStart) / CLOCKS_PER_SEC << endl;
            myfile.close();

        }
    MPI_Finalize();
    return 0;
}
//****************************************************************************80
pair<int, int> random_walks(int a, int b, int x, double p, int N)
{
    pair<int,int> hits = make_pair(0, 0);
    int i;
    for (i = 0; i < N; ++i) {
        int cur_x = x;
        while ((cur_x > a) and (cur_x < b)) {
            double cur_rand = ((double) std::rand() / (RAND_MAX));
            if (cur_rand < p) {  //to right
                cur_x += 1;
            } else {
                cur_x -= 1;
            }
            hits.second += 1;
        }
        if (cur_x == b) {
            hits.first += 1;
        }
    }
    return hits;
}

