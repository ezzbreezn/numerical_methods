#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <math.h>

int make_system(int argc, char *argv[], double ***pA, double **pf, int *n);
int read_file(FILE *fin, double ***pA, double **pf, int n);
void free_matrix(double **A, int n); 
void fill_matrix(double **A, int n); 
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

//Check if the string contains only digits
int 
digit_check(char *s)
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
    printf("./program <ARG1> [<ARG2>]\n");
    printf("ARG1:\n");
    printf("\t-help (-h) to print help info\n");
    printf("\t-input (-i) to read matrix from input file\n");
    printf("\t-formula (-f) to specify elements using the formula\n");
    printf("ARG2:\n");
    printf("\tmatrix size as integer if ARG1 is set to -formula (-f)\n");
    printf("\tpath to input file if if ARG1 is set to -input (-i)\n");
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
make_system(int argc, char *argv[], double ***pA, double **pf, int *pn)
{
    int n;
    double **A, *f;
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
        int res = read_file(fin, &A, &f, n);
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
        if (argc < 3 || !digit_check(argv[2])) {
            printf("Wrong input: matrix size is not specified\n");
            return 1;
        }
        sscanf(argv[2], "%d", &n);
        f = realloc(NULL, n * sizeof(double));
        A = allocate_matrix(n);
        if (f == NULL || A == NULL) {
            if (f) {
                free(f);
            }
            if (A) {
                free_matrix(A, n);
            }
            printf("Error: can't allocate memory\n");
            return 1;
        }
        fill_matrix(A, n);
        fill_vector(f, n);
    } else {
        printf("Error: wrong first parameter\n");
        return 1;
    }
    *pA = A;
    *pf = f;
    *pn = n;
    return 0;
}

double **
allocate_matrix(int n)
{
    double **A = realloc(NULL, n * sizeof(*A));
    if (A == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        A[i] = realloc(NULL, n * sizeof(double));
        if (A[i] == NULL) {
            for (int j = 0; j < i; ++j) {
                free(A[j]);
            }
            free(A);
            return NULL;
        }
    }
    return A;
}

void 
print_system(double **A, double *f, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%12.6f", A[i][j]);
        }
        printf(" | %12.6f\n", f[i]);
    }
}

double **
copy_matrix(double **A, int n)
{
    double **res = realloc(NULL, n * sizeof(*res));
    for (int i = 0; i < n; ++i) {
        res[i] = realloc(NULL, n * sizeof(double));
        for (int j = 0; j < n; ++j) {
            res[i][j] = A[i][j];
        }
    }
    return res;
}

double **
matrix_x_matrix(double **A, double **B, int n)
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
                res[i][j] += A[i][c] * B[c][j];
            }
        }
    }
    return res;
}

double *
matrix_x_vector(double **A, double *f, int n)
{
    double *res = realloc(NULL, n * sizeof(double));
    if (res == NULL) {
        return NULL;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            res[i] += A[i][j] * f[j];
        }
    }
    return res;
}

//Inverse matrix calculation
void 
transpose(double **A, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            swap_double(&A[i][j], &A[j][i]);
        }
    }
}

void 
free_matrix(double **A, int n)
{
    for (int i = 0; i < n; ++i) {
        free(A[i]);
    }
    free(A);
}

void 
print_matrix(double **A, int n)
{
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            printf("%12.6f", A[i][j]);
        }
        printf("\n");
    }
}

//Create the system matrix according to the formula
void 
fill_matrix(double **A, int n)
{
    double m;
    printf("Enter M for the formula:\n");
    scanf("%lf", &m);
    m = 1.001 - 2 * m / 1000;
    double pow1 = m * m, pow2 = (m - 1) * (m - 1);
    int maxpow = 2 * n;
    for (int pow = 2; pow <= maxpow; ++pow) {
        for (int j = 1; j < pow; ++j) {
            int i = pow - j;
            if (i <= n && j <= n) {
                A[i - 1][j - 1] = i == j ? pow2 : pow1 + 0.1 * (j - i);
            }
        }
        pow1 *= m;
        pow2 *= (m - 1);
    }
}


//Create the system vector according to the formula
void 
fill_vector(double *f, int n)
{
    double x;
    printf("Enter x for the formula:\n");
    scanf("%lf", &x);
    for (int i = 1; i <= n; ++i) {
        f[i - 1] = x * exp(x / i) * cos(x / i); 
    }
}

void 
print_vector(double *f, int n)
{
    for (int i = 0; i < n; ++i) {
        printf(" x%d = %12.6f\n", i + 1, f[i]);
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
read_file(FILE *in, double ***pA, double **pf, int n)
{
    double *f = realloc(NULL, n * sizeof(double));
    double **A = allocate_matrix(n);
    if (f == NULL || A == NULL) {
        if (f) {
            free(f);
        }
        if (A) {
            free_matrix(A, n);
        }
        return 2;
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (fscanf(in, "%lf", &A[i][j]) != 1) {
                free_matrix(A, n);
                free(f);
                return 1;
            }
        }
        if (fscanf(in, "%lf", &f[i]) != 1) {
            free_matrix(A, n);
            free(f);
            return 1;
        }
    }
    *pA = A;
    *pf = f;
    return 0;
}

double 
matrix_norm(double **A, int n)
{
    double norm = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            norm += A[i][j] * A[i][j];
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
determinant(double **A, int n)
{
    int sign = 1, flag = 1, i, j, c;
    double first;
    for (i = 0; i < n - 1; ++i) {
        flag = 1;
        if (A[i][i] == 0) {
            for (j = 0; j < n; ++j) {
                if (A[j][i] != 0) {
                    flag = 0;
                    break;
                }
            }
            if (flag) {
                return 0;
            }
            for (c = 0; c < n; ++c) {
                swap_double(&A[i][c], &A[j][c]);
            }
            sign *= -1;
        }
        for (j = i + 1; j < n; ++j) {
            first = A[j][i];
            for (c = i; c < n; ++c) {
                A[j][c] -= first * A[i][c] / A[i][i];
            }
        }
    }
    double res = 1;
    for (int i = 0; i < n; ++i) {
        res *= A[i][i];
    }
    return res * sign;
}

//Gauss's method implementation
//converting the system to the triangular form
//computing the solution and returning it
double *
gauss_method(double **A, double *f, int n)
{
    double div;
    for (int i = 0; i < n; ++i) {
        if (A[i][i] == 0) {
            for (int j = i + 1; j < n; ++j) {
                if (A[j][i] != 0) {
                    for (int c = i; c < n; ++c) {
                        swap_double(&A[j][c], &A[i][c]);
                    }
                    swap_double(&f[j], &f[i]);
                }
                break;
            }
        }
        div = A[i][i];
        for (int j = i; j < n; ++j) {
            A[i][j] /= div;
        }
        f[i] /= div;
        for (int j = i + 1; j < n; ++j) {
            for (int c = i + 1; c < n; ++c) {
                A[j][c] -= A[j][i] * A[i][c];
            }
            f[j] -= A[j][i] * f[i];
            A[j][i] = 0;
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
            res[i] -= res[j] * A[i][j];
        }
    }
    return res;
}

//Gauss's method with the choose of main element implementation
double *
gauss_with_choice(double **A, double *f, int n)
{
    double div, max;
    int colmax;
    int *pos = realloc(NULL, n * sizeof(int));
    for (int i = 0; i < n; ++i) {
        pos[i] = i;
    }
    for (int i = 0; i < n; ++i) {
        max = fabs(A[i][i]);
        colmax = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(A[i][j]) > max) {
                max = fabs(A[i][j]);
                colmax = j;
            }
        }
        if (i != colmax) {
            for (int c = 0; c < n; ++c) {
                swap_double(&A[c][colmax], &A[c][i]);
            }
            swap_int(&pos[i], &pos[colmax]);
        }
        div = A[i][i];
        for (int j = i; j < n; ++j) {
            A[i][j] /= div;
        }
        f[i] /= div;
        for (int j = i + 1; j < n; ++j) {
            for (int c = i + 1; c < n; ++c) {
                A[j][c] -= A[j][i] * A[i][c];
            }
            f[j] -= A[j][i] * f[i];
            A[j][i] = 0;
        }
    }
    double *res = realloc(NULL, n * sizeof(double));
    if (res == NULL) {
        return NULL;
    }
    for (int i = n - 2; i >= 0; --i) {
        for (int j = n - 1; j > i; --j) {
            f[i] -= f[j] * A[i][j];
        }
    }
    for (int i = 0; i < n; ++i) {
        res[pos[i]] = f[i];
    }
    free(pos);
    return res;
}

//Calculate inverse matrix for A matrix
double **
inverse_matrix(double **A, int n)
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
        max = fabs(A[i][i]);
        strmax = i;
        for (int j = i + 1; j < n; ++j) {
            if (fabs(A[j][i]) > max) {
                max = fabs(A[j][i]);
                strmax = j;
            }
        }
        for (int c = 0; c < n; ++c) {
            swap_double(&A[strmax][c], &A[i][c]);
            swap_double(&inv[strmax][c], &inv[i][c]);
        }
        div = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= div;
            inv[i][j] /= div;
        }
        for (int j = i + 1; j < n; ++j) {
            for (int c = i + 1; c < n; ++c) {
                A[j][c] -= A[j][i] * A[i][c];
            }
            for (int c = 0; c < n; ++c) {
                inv[j][c] -= A[j][i] * inv[i][c];
            }
            A[j][i] = 0;
        }
    }
    for (int col = 0; col < n; ++col) {
        for (int i = 0; i < col; ++i) {
            div = A[i][col];
            for (int j = 0; j < n; ++j) {
                A[i][j] -= A[col][j] * div;
                inv[i][j] -= inv[col][j] * div;
            }
        }
    }
    return inv;
}

//Seidel's method implementation
void 
seidel_method(double **A, double *f, int n)
{
    //To guarantee the applicability for Ax=f consider the equivalent 
    //system A^T*Ax = A^T*f
    double epsilon;
    printf("Enter the solution accuracy: ");
    scanf("%lf", &epsilon);
    epsilon /= 1000.0;
    double **trans = copy_matrix(A, n);
    double *res = f;
    transpose(trans, n);
    A = matrix_x_matrix(trans, A, n);
    f = matrix_x_vector(trans, f, n);
    double *prev = realloc(NULL, n * sizeof(double));
    double *next = realloc(NULL, n * sizeof(double));
    double *diff = realloc(NULL, n * sizeof(double));
    double s;
    do {
        for (int i = 0; i < n; ++i) {
            s = 0;
            for (int j = 0; j < i; ++j) {
                s += next[j] * A[i][j];
            }
            for (int j = i + 1; j < n; ++j) {
                s += prev[j] * A[i][j];
            }
            next[i] = (f[i] - s) / A[i][i];
        }
        for (int i = 0; i < n; ++i) {
            diff[i] = prev[i] - next[i];
            prev[i] = next[i];
        }

    } while (vector_norm(diff, n) > epsilon);

    for (int i = 0; i < n; ++i) {
        res[i] = next[i];
    }
    free(prev);
    free(next);
    free_matrix(trans, n);
}

//Upper relaxation method implementation
void 
upper_relaxation(double **A, double *f, double w, int n)
{
    //To guarantee the applicability for Ax=f consider the equivalent 
    //system A^T*Ax = A^T*f
    double epsilon;
    printf("Enter the solution accuracy: ");
    scanf("%lf", &epsilon);
    epsilon /= 1000;
    int iter = 1;
    double **trans = copy_matrix(A, n);
    double *res = f;
    transpose(trans, n);
    A = matrix_x_matrix(trans, A, n);
    f = matrix_x_vector(trans, f, n);
    double *prev = realloc(NULL, n * sizeof(double));
    double *next = realloc(NULL, n * sizeof(double));
    double *diff = realloc(NULL, n * sizeof(double));
    double s;
    do {
        for (int i = 0; i < n; ++i) {
            s = 0;
            for (int j = 0; j < i; ++j) {
                s += next[j] * A[i][j];
            }
            for (int j = i; j < n; ++j) {
                s += prev[j] * A[i][j];
            }
            next[i] = prev[i] + w * (f[i] - s) / A[i][i];
        }
        for (int i = 0; i < n; ++i) {
            diff[i] = prev[i] - next[i];
            prev[i] = next[i];
        }
        iter++;
    } while (vector_norm(diff, n) > epsilon);
    printf("\n\nRequired %d iterations, w = %.6f\n\n", iter, w);
    for (int i = 0; i < n; ++i) {
        res[i] = next[i];
    }
    free(prev);
    free(next);
    free_matrix(trans, n);
}

double 
residual_norm(double **A, double *f, double *x, int n)
{
    double *tmp = matrix_x_vector(A, x, n);
    for (int i = 0; i < n; ++i) {
        tmp[i] -= f[i];
    }
    double norm = vector_norm(tmp, n);
    free(tmp);
    return norm;
}


int 
main(int argc, char *argv[])
{
    int n;
    double **A, *f;
    if (make_system(argc, argv, &A, &f, &n)) {
        return 1;
    }
    printf("The system of linear equations:\n\n");
    print_system(A, f, n);
    pid_t pid;

    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        free(f);
        printf("\ndet A = %.6f\n", determinant(A, n));
        free_matrix(A, n);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        free(f);
        double matnorm = matrix_norm(A, n);
        double **inv = inverse_matrix(A, n);
        if (inv == NULL) {
            memory_error();
        }
        double invnorm = matrix_norm(inv, n);
        printf("\nInverse matrix:\n");
        print_matrix(inv, n);
        printf("\nCondition number Ma = %.6f\n", matnorm * invnorm);
        free_matrix(A, n);
        free_matrix(inv, n);
        return 0;
    }
    wait(NULL);


    pid = fork();
    if (pid == -1) {
        pid_error();
    }
    if (!pid) {
        double **A_cpy = copy_matrix(A, n);
        double *f_cpy = copy_vector(f, n);
        double *res = gauss_method(A, f, n);
        if (res == NULL || A_cpy == NULL || f_cpy == NULL) {
            memory_error();
        }
        printf("\nGauss method:\nSolution:\n");
        print_vector(res, n);
        printf("\nResidual norm: %.12f\n", residual_norm(A_cpy, f_cpy, res, n));
        free_matrix(A, n);
        free(f);
        free_matrix(A_cpy, n);
        free(f_cpy);
        free(res);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        double **A_cpy = copy_matrix(A, n);
        double *f_cpy = copy_vector(f, n);
        double *res = gauss_with_choice(A, f, n);
        if (res == NULL || A_cpy == NULL || f_cpy == NULL) {
            memory_error();
        }
        printf("\n\nGauss method with choosing the main element:\nSolution:\n");
        print_vector(res, n);
        printf("\nResidual norm: %.12f\n\n", residual_norm(A_cpy, f_cpy, res, n));
        free_matrix(A, n);
        free(f);
        free_matrix(A_cpy, n);
        free(f_cpy);
        free(res);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        double **A_cpy = copy_matrix(A, n);
        double *f_cpy = copy_vector(f, n);
        if (A_cpy == NULL || f_cpy == NULL) {
            memory_error();
        }
        printf("Seidel method:\n");
        seidel_method(A, f, n);
        printf("Solution:\n");
        print_vector(f, n);
        printf("\nResidual norm: %.12f\n", residual_norm(A_cpy, f_cpy, f, n));
        free_matrix(A, n);
        free(f);
        free_matrix(A_cpy, n);
        free(f_cpy);
        return 0;
    }
    wait(NULL);
    if ((pid = fork()) == -1) {
        pid_error();
    }
    if (!pid) {
        double **A_cpy = copy_matrix(A, n);
        double *f_cpy = copy_vector(f, n);
        if (A_cpy == NULL || f_cpy == NULL) {
            memory_error();
        }
        double w;
        printf("\nUpper relaxation method:\n");
        printf("Enter w parameter:\n");
        scanf("%lf", &w);
        upper_relaxation(A, f, w, n);
        printf("Solution:\n");
        print_vector(f, n);
        printf("\nResidual norm: %.12f\n", residual_norm(A_cpy, f_cpy, f, n));
        free_matrix(A, n);
        free(f);
        free_matrix(A_cpy, n);
        free(f_cpy);
        return 0;
    }
    wait(NULL);
    free_matrix(A, n);
    free(f);
    return 0;
}
