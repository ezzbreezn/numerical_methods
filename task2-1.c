#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>

void 
memory_err(void)
{
    printf("Error: memory error\n");
    _exit(1);
}

//Functions coefficients examples

double 
p1(double x)
{
    return -x;
}

double 
q1(double x)
{
    return 2;
}

double 
f1(double x)
{
    return -x + 1;
}

double 
p2(double x)
{
    return 1;
}

double 
q2(double x)
{
    return 0;
}

double 
f2(double x)
{
    return -1;
}

double 
p3(double x)
{
    return 1;
}

double 
q3(double x)
{
    return -x;
}

double 
f3(double x)
{
    return (4 * x * x * x - 2 * x + 1) / (4 * x * sqrt(x));
}

double 
p4(double x)
{
    return -1;
}

double 
q4(double x)
{
    return 0;
}

double 
f4(double x)
{
    return 0;
}

double
p5(double x)
{
    return 0;
}

double
q5(double x)
{
    return -1;
}

double
f5(double x)
{
    return -2 * x;
}


void 
print(double *f, double x, double h, int n)
{
    for (int i = 0; i < n; ++i) {
        printf("%15.6f %15.6f\n", x, f[i]);
        x += h;
    }
}

//Reduction method
double 
*reduction_method(int n, double last, double *alpha, double *beta)
{
    double *res = realloc(NULL, (n + 1) * sizeof(*res));
    if (res == NULL) {
        return NULL;
    }
    res[n] = last;
    for (int i = n - 1; i >= 0; --i) {
        res[i] = alpha[i] * res[i + 1] + beta[i];
    }
    return res;
}

//According differential equation solution
double 
*ODE_solution(int n, double a, double b, double s1, double g1, double d1,
        double s2, double g2, double d2, double (*p)(double), double (*q)(double), double (*f)(double))
{
    double h = (b - a) / n, x = a + h, A, B, C, tmp;
    double *alpha = realloc(NULL, n * sizeof(*alpha));
    if (alpha == NULL) {
        return NULL;
    }
    double *beta = realloc(NULL, n * sizeof(*beta));
    if (beta == NULL) {
        free(alpha);
        return NULL;
    }

    tmp = s1 * h - g1;
    alpha[0] = -g1 / tmp;
    beta[0] = (d1 * h) / tmp;
    for (int i = 0; i < n - 1; i++) {
        tmp = p(x) / (2 * h);
        A = 1.0 / (h * h);
        B = A;
        C = -2 * A;
        A -= tmp;
        B += tmp;
        C += q(x);

        tmp = C + A * alpha[i];
        alpha[i + 1] = -B / tmp;
        beta[i + 1] = (-f(x) - beta[i] * A) / tmp;
        x += h;
    }

    tmp = 1 / (g2 + h * s2);

//p_n
    A = g2 / tmp;
//q_n
    tmp *= h * d2;
//y_n
    tmp = (A * beta[n - 1] + tmp) / (1 - A * alpha[n - 1]); 

    return reduction_method(n, tmp, alpha, beta);
}
int 
main(int argc, char *argv[])
{
    int n, num;
    double a, b, s1, s2, g2, g1, d2, d1;
    printf("The available differential equations:\n");
    printf("General form:\n");
    printf("y'' + p(x) * y' + q(x) * y = -f(x)\n");
    printf("1) y'' - xy' + 2y = x - 1\n");
    printf("2: y'' + y' = 1\n");
    printf("3: y'' + y' - xy = (-4x^3+2x-1) / ( 4x*sqrt(x) )\n");
    printf("4: y'' - y' = 0\n");
    printf("5: y'' - y = 2x\n");
    printf("Enter the equation number: ");
    scanf("%d", &num);
    printf("Enter the left border: ");
    scanf("%lf", &a);
    printf("Enter the right border: ");
    scanf("%lf", &b);
    printf("Enter the border conditions:\ns1*y(a)+g1*y'(a)=d1\ns2*y(b)+g2*y'(b)=d2\n");
    printf("Enter the s1: ");
    scanf("%lf", &s1);
    printf("Enter the g1: ");
    scanf("%lf", &g1);
    printf("Enter d1: ");
    scanf("%lf", &d1);
    printf("Enter s2: ");
    scanf("%lf", &s2);
    printf("Enter g2: ");
    scanf("%lf", &g2);
    printf("Enter d2: ");
    scanf("%lf", &d2);
    printf("Enter the number of grid nodes: ");
    scanf("%d", &n);

    void *f, *p, *q;
    switch (num) {
    case 1:
        f = f1;
        p = p1;
        q = q1;
        break;
    case 2:
        f = f2;
        p = p2;
        q = q2;
        break;
    case 3:
        f = f3;
        p = p3;
        q = q3;
        break;
    case 4:
        f = f4;
        p = p4;
        q = q4;
        break; 
    default:
        f = f5;
        p = p5;
        q = q5;
        break;
    }

    double *res = ODE_solution(n, a, b, s1, g1, d1, s2, g2, d2, p, q, f);
    if (res == NULL) {
        memory_err();
    }
    printf("\n\nSolution y(x):\n");
    print(res, a, (b-a) / n, n + 1);
    printf("\n");
    free(res);
    return 0;
}
