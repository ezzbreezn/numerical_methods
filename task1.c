#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <math.h>

int make_system(int argc, char *argv[], double ***pa, double **pf, int *n);
int read_file(FILE *fin, double ***pa, double **pf, int n);
void free_matrix(double **a, int n); 
void fill_matrix(double **a, int n); 
double **allocate_matrix(int n); 
void fill_vector(double *f, int n);

void 
swap_double(double *x, double *y)
{
    double r = *x;
    *x = *y;
    *y = r;
}

void 
swap_int(int *x, int *y)
{
    int r = *x;
    *x = *y;
    *y = r;
}

//Check whether the string contains nothing but the digits' symbols
int 
isdigits(char *s)
{
    while (*s) {
        if (*s < '0' || *s > '9') {
            return 0;
        }
        ++s;
    }
    return 1;
}

//Prints the help info
void 
help(void)
{
    printf("GAUSS'S METHOD, ZEIDEL'S METHOD AND UPPER RELAXATION METHOD\n");
    printf("Usage format:\n");
    printf("./program <ARG1> [<ARG2>] [<ARG3>]\n");
    printf("ARG1:\n");
    printf("\t-help (-h) to print help info\n");
    printf("\t-input (-i) to read matrix from input file\n");
    printf("\t-formula (-f) to specify elements using the formula\n");
    printf("ARG2:\n");
    printf("\tmatrix size as integer if ARG1 is set to -formula (-f)\n");
    printf("\tpath to input file if if ARG1 is set to -input (-i)\n");
    printf("ARG3:\n");
    printf("\tx parameter for calculating the system if ARG1 is set to -formula (-f)\n");
}

void 
memory_error(void)
{
    printf("Error: out of memory\n");
    exit(1);
}

void 
pid_error(void) 
{
    printf("Error: can't create process\n");
    exit(1);
}

//Parsing program arguments and reading or calculating the matrix and the vector according to 
//chosen arguments
//Return value: 0 if the system of linear equations successfully created, 1 in the case of fail
int 
make_system(int argc, char *argv[], double ***pa, double **pf, int *pn)
{
    int n;
    double **a, *f;
    if (argc < 2) {
        printf("Wrong input: first parameter missing\n");
        return 1;
    }
    if (strcmp(argv[1], "-help") == 0 || strcmp(argv[1], "-h") == 0) {
        help();
        exit(0);
    } else if (strcmp(argv[1], "-input") == 0 || strcmp(argv[1], "-i") == 0) {
        if (argc < 3) {
            printf("Wrong input: input file path missing\n");
            return 1;
        }
        FILE *fin = fopen(argv[2], "r");
        if (fin == NULL) {
            printf("Error: can't read the file\n");
            return 1;
        }
        if (fscanf(fin, "%d", &n) != 1) {
            printf("Wrong input: wrong input data format\n");
            fclose(fin);
            return 1;
        }
        int res = read_file(fin, &a, &f, n);
        fclose(fin);
        if (res != 0) {
            if (res == 2) {
                memory_error();
            } else {
                printf("Error: can't read data\n");
            }
            return 1;
        }
    } else if (strcmp(argv[1], "-formula") == 0 || strcmp(argv[1], "-f") == 0) {
        if (argc < 3 || !isdigits(argv[2])) {
            printf("Wrong input: matrix size is not specified\n");
            return 1;
        }
        sscanf(argv[2], "%d", &n);
        f = realloc(NULL, n * sizeof(double));
        a = allocate_matrix(n);
        if (f == NULL || a == NULL) {
            if (f) {
                free(f);
            }
            if (a) {
                free_matrix(a, n);
            }
            printf("Error: can't allocate memory\n");
            return 1;
        }
        fill_matrix(a, n);
        fill_vector(f, n);
    } else {
        printf("Error: wrong first parameter\n");
        return 1;
    }
    *pa = a;
    *pf = f;
    *pn = n;
    return 0;
}

double **
allocate_matrix(int n)
{
    double **a = realloc(NULL, n * sizeof(*a));
    if (a == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        a[i] = realloc(NULL, n * sizeof(double));
        if (a[i] == NULL) {
            for (int j = 0; j < i; ++j) {
                free(a[j]);
            }
            free(a);
            return NULL;
        }
    }
    return a;
}

void 
print_system(double **a, double *f, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%12.6f", a[i][j]);
        }
        printf(" | %12.6f\n", f[i]);
    }
}

double **
copy_matrix(double **a, int n)
{
    double **res = realloc(NULL, n * sizeof(*res));
    for (int i = 0; i < n; ++i) {
        res[i] = realloc(NULL, n * sizeof(double));
        for (int j = 0; j < n; ++j) {
            res[i][j] = a[i][j];
        }
    }
    return res;
}

double **
matrix_x_matrix(double **a, double **b, int n)
{
    double **res = realloc(NULL, n * sizeof(*res));
    if (res == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        res[i] = realloc(NULL, n * sizeof(double));
        if (res[i] == NULL) {
            return NULL;
        }
        for (int j = 0; j < n; ++j) {
            for (int c = 0; c < n; ++c) {
                res[i][j] += a[i][c] * b[c][j];
            }
        }
    }
    return res;
}

double *
matrix_x_vector(double **a, double *f, int n)
{
    double *res = realloc(NULL, n * sizeof(double));
    if (res == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i] += a[i][j] * f[j];
        }
    }
    return res;
}

//spoils the matrix
void 
transpose(double **a, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            swap_double(&a[i][j], &a[j][i]);
        }
    }
}

void 
free_matrix(double **a, int n)
{
    for (int i = 0; i < n; ++i) {
        free(a[i]);
    }
    free(a);
}

void 
print_matrix(double **a, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%12.6f", a[i][j]);
        }
        printf("\n");
    }
}

//Create the system matrix according to the formula
void 
fill_matrix(double **a, int n)
{
    double m;
    printf("Type m parameter for filling the matrix:\n");
    scanf("%lf", &m);
    m = 1.001 - 2 * m / 1000;
    double pow1 = m * m, pow2 = (m - 1) * (m - 1);
    int maxpow = 2 * n, i;
    for (int pow = 2; pow <= maxpow; ++pow) {
        for (int j = 1; j < pow; ++j) {
            i = pow - j;
            if (i <= n && j <= n) {
                a[i - 1][j - 1] = i == j ? pow2 : pow1 + 0.1 * (j - i);
            }
        }
        pow1 *= m;
        pow2 *= m - 1;
    }
}

//static void fill_matrix1(double **a, int n)
//{
//    for (int i = 0; i < n; ++i)
//        for (int j = 0; j < n; ++j)
//            a[i][j] = rand() % (i * j + 2) * 1.0 / (rand() % (j + 1) + 1) ;
//}

//Create the system vector according to the formula
void 
fill_vector(double *f, int n)
{
    double x;
    printf("Type x parameter for filling the vector:\n");
    scanf("%lf", &x);
    for (int i = 1; i <= n; ++i) {
        f[i - 1] = x * exp(x / i) * cos(x / i); 
    }
}

void 
print_vector(double *f, int n)
{
    for (int i = 0; i < n; ++i) {
        printf(" x%d = %8.6f\n", i + 1, f[i]);
    }
}

double *
copy_vector(double *f, int n)
{
    double *copy = realloc(NULL, n * sizeof(double));
    if (copy == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        copy[i] = f[i];
    }
    return copy;
}

//Allocate memory for system matrix and vector, read their
//values from the input file and write it in the according pointers
//Return value: 0 if data succsefully read and written to pointers,
//1 if data is incorrect in the input file, 2 if memory allocation failed
int 
read_file(FILE *in, double ***pa, double **pf, int n)
{
    double *f = realloc(NULL, n * sizeof(double));
    double **a = allocate_matrix(n);
    if (f == NULL || a == NULL) {
        if (f) {
            free(f);
        }
        if (a) {
            free_matrix(a, n);
        }
        return 2;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (fscanf(in, "%lf", &a[i][j]) != 1) {
                free_matrix(a, n);
                free(f);
                return 1;
            }
        }
        if (fscanf(in, "%lf", &f[i]) != 1) {
            free_matrix(a, n);
            free(f);
            return 1;
        }
    }
    *pa = a;
    *pf = f;
    return 0;
}

double 
matrix_norm(double **a, int n)
{
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            norm += a[i][j] * a[i][j];
        }
    }
    return sqrt(norm);
}

double 
vector_norm(double *f, int n)
{
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        norm += f[i] * f[i];
    }
    return sqrt(norm);
}

//Calculate matrix determinant
double 
determinant(double **a, int n)
{
    int sign = 1, flag = 1, i, j, c;
    double first;
    for (i = 0; i < n - 1; ++i) {
        if (a[i][i] == 0) {
            for (j = 0; j < n; ++j) {
                if (a[j][i] != 0) {
                    flag = 0;
                    break;
                }
            }
            if (flag++) {
                return 0; //no non-zero elements in the column
            }
            for (c = 0; c < n; ++c) {
                swap_double(&a[i][c], &a[j][c]);
            }
            sign = -sign;
        }
        for (j = i + 1; j < n; ++j) {
            first = a[j][i];
            for (c = i; c < n; ++c) {
                a[j][c] -= first * a[i][c] / a[i][i];
            }
        }
    }
    double res = 1;
    for (int i = 0; i < n; ++i) {
        res *= a[i][i];
    }
    return res * sign;
}

//Gauss's method implementation
//converting the system to the triangular form
//computing the solution and returning it
double *
gauss_method(double **a, double *f, int n)
{
    double div;
    for (int i = 0; i < n; ++i) {
        if (a[i][i] == 0) {
            for (int j = i + 1; j < n; ++j) {
                if (a[j][i] != 0) {
                    for (int c = i; c < n; ++c) {
                        swap_double(&a[j][c], &a[i][c]);
                    }
                    swap_double(&f[j], &f[i]);
                }
            }
        }
        div = a[i][i];
        for (int j = i; j < n; ++j) {
            a[i][j] /= div;
        }
        f[i] /= div;
        for (int j = i + 1; j < n; ++j) {
            for (int c = i + 1; c < n; ++c) {
                a[j][c] -= a[j][i] * a[i][c];
            }
            f[j] -= a[j][i] * f[i];
            a[j][i] = 0;
        }
    }
    double *res = realloc(NULL, n * sizeof(double));
    if (res == NULL) {
        return NULL;
    }
    res[n - 1] = f[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        res[i] = f[i];
        for (int j = n - 1; j > i; --j) {
            res[i] -= res[j] * a[i][j];
        }
    }
    return res;
}

//Gauss's method with the choose of main element implementation
double *
gauss_with_choice(double **a, double *f, int n)
{
    double div, max;
    int colmax;
    int *pos = realloc(NULL, n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        pos[i] = i;
    }
    for (int i = 0; i < n; ++i) {
        max = fabs(a[i][i]);
        colmax = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(a[i][j]) > max) {
                max = fabs(a[i][j]);
                colmax = j;
            }
        }
        if (i != colmax) {
            for (int c = 0; c < n; ++c) {
                swap_double(&a[c][colmax], &a[c][i]);
            }
            swap_int(&pos[i], &pos[colmax]);
        }
        div = a[i][i];
        for (int j = i; j < n; ++j) {
            a[i][j] /= div;
        }
        f[i] /= div;
        for (int j = i + 1; j < n; ++j) {
            for (int c = i + 1; c < n; ++c) {
                a[j][c] -= a[j][i] * a[i][c];
            }
            f[j] -= a[j][i] * f[i];
            a[j][i] = 0;
        }
    }
    double *res = realloc(NULL, n * sizeof(double));
    if (res == NULL) {
        return NULL;
    }
    for (int i = n - 2; i >= 0; --i) {
        for (int j = n - 1; j > i; --j) {
            f[i] -= f[j] * a[i][j];
        }
    }
    for (int i = 0; i < n; ++i) {
        res[pos[i]] = f[i];
    }
    free(pos);
    return res;
}

//Seidel's method implementation
void 
seidel_method(double **a, double *f, int n)
{
    //the matrix A^T * A is always self-adjoint and positive definite
    //Ax=f => A^T*Ax = A^T*f
    double eps;
    printf("Enter eps: ");
    scanf("%lf", &eps);
    eps /= 1000;
    double **trans = copy_matrix(a, n);
    double *res = f;
    transpose(trans, n);
    a = matrix_x_matrix(trans, a, n);
    f = matrix_x_vector(trans, f, n);
    double *prev = realloc(NULL, n * sizeof(double));
    double *next = realloc(NULL, n * sizeof(double));
    double *diff = realloc(NULL, n * sizeof(double));
    double s;
    do {
        for (int i = 0; i < n; ++i) {
            s = 0;
            for (int j = 0; j < i; ++j) {
                s += next[j] * a[i][j];
            }
            for (int j = i + 1; j < n; ++j) {
                s += prev[j] * a[i][j];
            }
            next[i] = (f[i] - s) / a[i][i];
        }
        for (int i = 0; i < n; ++i) {
            diff[i] = prev[i] - next[i];
            prev[i] = next[i];
        }

    } while (vector_norm(diff, n) > eps);

    for (int i = 0; i < n; ++i) res[i] = next[i];
    free(prev);
    free(next);
    free(f);
    free_matrix(a, n);
    free_matrix(trans, n);
}

//Upper relaxation method implementation
void 
upper_relaxation(double **a, double *f, double w, int n)
{
    //the matrix A^T * A is always self-adjoint and positive definite
    //Ax=f => A^T*Ax = A^T*f
    double eps;
    printf("Enter eps: ");
    scanf("%lf", &eps);
    eps /= 1000;
    int iter = 1;
    double **trans = copy_matrix(a, n);
    double *res = f;
    transpose(trans, n);
    a = matrix_x_matrix(trans, a, n);
    f = matrix_x_vector(trans, f, n);
    double *prev = realloc(NULL, n * sizeof(double));
    double *next = realloc(NULL, n * sizeof(double));
    double *diff = realloc(NULL, n * sizeof(double));
    double s;
    do {
        for (int i = 0; i < n; ++i) {
            s = 0;
            for (int j = 0; j < i; ++j) {
                s += next[j] * a[i][j];
            }
            for (int j = i; j < n; ++j) {
                s += prev[j] * a[i][j];
            }
            next[i] = prev[i] + w * (f[i] - s) / a[i][i];
        }
        for (int i = 0; i < n; ++i) {
            diff[i] = prev[i] - next[i];
            prev[i] = next[i];
        }
        iter++;
    } while (vector_norm(diff, n) > eps);
    printf("\n\n%d iterations, w = %8.3f\n\n", iter, w);
    for (int i = 0; i < n; ++i) {
        res[i] = next[i];
    }
    free(prev);
    free(next);
    free(f);
    free_matrix(a, n);
    free_matrix(trans, n);
}

double 
residual_norm(double **a, double *f, double *x, int n)
{
    double *tmp = matrix_x_vector(a, x, n);
    for (int i = 0; i < n; ++i) {
        tmp[i] -= f[i];
    }
    double norm = vector_norm(tmp, n);
    free(tmp);
    return norm;
}

//Calculate inverse matrix for A matrix
double **
inverse_matrix(double **a, int n)
{
    double **inv = realloc(NULL, n * sizeof(double));
    if (inv == NULL) {
        memory_error();
    }
    for (int i = 0; i < n; ++i) {
        inv[i] = realloc(NULL, n * sizeof(double));
        if (inv[i] == NULL) {
            memory_error();
        }
        inv[i][i] = 1;
    }
    double div, max;
    int strmax;
    for (int i = 0; i < n; ++i) {
        max = fabs(a[i][i]);
        strmax = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(a[j][i]) > max) {
                max = fabs(a[j][i]);
                strmax = j;
            }
        }
        for (int c = 0; c < n; ++c) {
            swap_double(&a[strmax][c], &a[i][c]);
            swap_double(&inv[strmax][c], &inv[i][c]);
        }
        div = a[i][i];
        for (int j = 0; j < n; ++j) {
            a[i][j] /= div;
            inv[i][j] /= div;
        }
        //printf("\n");
        //print_matrix(a,n);
        for (int j = i + 1; j < n; ++j) {
            for (int c = i + 1; c < n; ++c) {
                a[j][c] -= a[j][i] * a[i][c];
            }
            for (int c = 0; c < n; ++c) {
                inv[j][c] -= a[j][i] * inv[i][c];
            }
            a[j][i] = 0;
        }
    }
    for (int col = 0; col < n; ++col) {
        for (int i = 0; i < col; ++i) {
            div = a[i][col];
            for (int j = 0; j < n; ++j) {
                a[i][j] -= a[col][j] * div;
                inv[i][j] -= inv[col][j] * div;
            }
        }
    }
    return inv;
}

int 
main(int argc, char *argv[])
{
    int n;
    double **a, *f;
    if (make_system(argc, argv, &a, &f, &n)) {
        return 1;
    }
    int fd = open("out.txt", O_WRONLY | O_APPEND);
    dup2(fd, STDOUT_FILENO);
    printf("The system of linear equations:\n\n");
    print_system(a, f, n);
    pid_t pid = fork();
    if (pid == -1) {
        pid_error();
    }
    if (!pid) {
        double **acp = copy_matrix(a, n);
        double *fcp = copy_vector(f, n);
        double *res = gauss_method(a, f, n);
        if (res == NULL || acp == NULL || fcp == NULL) {
            memory_error();
        }
        printf("\nGauss method:\nSolution:\n");
        print_vector(res, n);
        printf("\nResidual norm: %30.20f\n", residual_norm(acp, fcp, res, n));
        free_matrix(a, n);
        free(f);
        free_matrix(acp, n);
        free(fcp);
        free(res);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        double **acp = copy_matrix(a, n);
        double *fcp = copy_vector(f, n);
        double *res = gauss_with_choice(a, f, n);
        if (res == NULL || acp == NULL || fcp == NULL) {
            memory_error();
        }
        //printf("\n\nGauss method (choosing maximum value in the string):\nThe solution vector is:\n");
        printf("\n\nGauss method with choosing the main element:\nSolution:\n");
        print_vector(res, n);
        printf("\nResidual norm: %30.20f\n", residual_norm(acp, fcp, res, n));
        free_matrix(a, n);
        free(f);
        free_matrix(acp, n);
        free(fcp);
        free(res);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        free(f);
        printf("\n\ndet A = %20.15f\n", determinant(a, n));
        free_matrix(a, n);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        free(f);
        double matnorm = matrix_norm(a, n);
        double **inv = inverse_matrix(a, n);
        if (inv == NULL) {
            memory_error();
        }
        double invnorm = matrix_norm(inv, n);
        printf("\nInverse matrix:\n");
        print_matrix(inv, n);
        printf("\nCondition number = %8.3f\n\n", matnorm * invnorm);
        free_matrix(a, n);
        free_matrix(inv, n);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        double **acp = copy_matrix(a, n);
        double *fcp = copy_vector(f, n);
        if (acp == NULL || fcp == NULL) {
            memory_error();
        }
        printf("Seidel method:\n");
        seidel_method(a, f, n);
        printf("Solution:\n");
        print_vector(f, n);
        printf("\nResidual norm: %8.6f\n", residual_norm(acp, fcp, f, n));
        free_matrix(a, n);
        free(f);
        free_matrix(acp, n);
        free(fcp);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        double **acp = copy_matrix(a, n);
        double *fcp = copy_vector(f, n);
        if (acp == NULL || fcp == NULL) {
            memory_error();
        }
        double w;
        printf("\nUpper relaxation method:\n");
        printf("Enter w parameter:\n");
        scanf("%lf", &w);
        upper_relaxation(a, f, w, n);
        printf("Solution:\n");
        print_vector(f, n);
        printf("\nResidual norm: %8.6f\n", residual_norm(acp, fcp, f, n));
        free_matrix(a, n);
        free(f);
        free_matrix(acp, n);
        free(fcp);
        return 0;
    }
    wait(NULL);
    free_matrix(a, n);
    free(f);
    return 0;
}
