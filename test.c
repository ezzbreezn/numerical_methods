#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <math.h>

double 
f(double x)
{
    //return - (x * x) + 2 * x - 2 + 12 * exp(-x);
    //return -3 * x * x - 6 * x + 15 * exp(x - 1) - 5;
    //return exp(x) + 0.5 * sin(x) - 0.5 * cos(x);
    //return 4 * exp(x) / 3 * (2 - 3 * x) - 2 * exp(-2 * x) / 3; 
    //return exp(x) * (5 - 3 * x) / 3 - 2 * exp(-2 * x) / 3;
    //return sqrt(x);
    //return x + exp(-x) - exp(-1);
    //return exp(x) - 2;
    return (2 * exp(2) * x - 2 * x + exp(1-x)-exp(x+1))/(1-exp(2));
}



int
main(int argc, char *argv[])
{
    double start, end;
    int n;
    sscanf(argv[2], "%lf", &start);
    sscanf(argv[3], "%lf", &end);
    sscanf(argv[4], "%d", &n);
    double step = (end - start) / n;
    int file = open(argv[1], O_WRONLY | O_CREAT | O_TRUNC, 0666);
    dup2(file, 1);
    close(file);
    printf("x y\n");
    for (int i = 0; i <= n; ++i) {
        printf("%.6lf %.6lf\n", start, f(start));
        start += step;
    }
    return 0;
}
