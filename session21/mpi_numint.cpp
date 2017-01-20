#include <mpi.h>
#include <printf.hpp>

double f(double x) {
    return 4.0/(1+x*x);
}

double Integral(double (*functionPtr)(double),
                double lowerBound,
                double upperBound,
                int nof_subintervals)
{
    double h = (upperBound-lowerBound)/nof_subintervals;
    double integral = ((*functionPtr)(lowerBound)
                    + (*functionPtr)(upperBound))/2.0
                    + 2.0*(*functionPtr)((2*upperBound - h)/2);
    for(int k=1; k<nof_subintervals; ++k) {
        integral += ((*functionPtr)(lowerBound + k*h)
                     + 2.0 * (*functionPtr)((2*lowerBound+(2*k-1)*h)/2));
    }
    integral *= (h/3);

    return integral;
}


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int nof_processes; MPI_Comm_size(MPI_COMM_WORLD, &nof_processes);

    double (*func)(double)=&f;
    double lowerBound = 0;
    double upperBound = 1;
    double h = (upperBound - lowerBound)/nof_processes;

    /* if nof_subintervals not devisible by nof_processes:
       nof_subintervals will be reduced to
       nof_subintervals/nofprocesses * nof_processes
    */
    int nof_subintervals = 200;

    double pi = 0;
    for (int k=0; k<nof_processes; ++k) {
        double subSol;
        if(rank == k) {
            subSol = Integral(func, 0.0+k*h, 0.0+(k+1)*h,
                              nof_subintervals/nof_processes);
            MPI_Send(&subSol, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }
    if (!rank) {
        for (int k=0; k<nof_processes; ++k) {
            MPI_Status status;
            double erg;
            MPI_Recv(&erg, 1, MPI_DOUBLE, MPI_ANY_SOURCE,
                     0, MPI_COMM_WORLD, &status);
            pi += erg;
        }
        fmt::printf("Integral: %1.20lf\n", pi);
    }

MPI_Finalize();
}
