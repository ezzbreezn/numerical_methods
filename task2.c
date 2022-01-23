#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <string.h>

void 
mem_error(void)
{
    printf("Error: emory error\n");
    _exit(1);
}


double
f1(double x, double y)
{
    return -y - (x * x);
}

double
f2(double x, double y)
{
    return y + 3 * x * x - 1;
}

double
f3(double x, double y)
{
    return cos(x) + y;
}



double 
f4(double x, double y)
{
    return exp(x) - 3 * y;
}

double 
f5(double x, double y)
{
    return (x - x * x) * y;
}


double 
sf11(double x, double u, double v)
{
    return sin(x + u) - 1.1 * v;
}

double 
sf12(double x, double u, double v)
{
    return 2.5 * u - (x + v) * (x + v);
}

double 
sf21(double x, double u, double v)
{
    return (u - v) / x;
}

double 
sf22(double x, double u, double v)
{
    return 2 * v / x + u * x;
}

double 
sf31(double x, double u, double v)
{
    return cos(u + 1.3 * v) - 2.1;
}

double 
sf32(double x, double u, double v)
{
    return 1.3 / (x + 2.3 * u * u) + x + 1;
}

double 
sf41(double x, double u, double v)
{
    return 2 * u + v + exp(x);
}

double 
sf42(double x, double u, double v)
{
    return -2 * x + 2 * x;
}


double 
sf51(double x, double u, double v)
{
    return 2 * u - 4 * v;
}

double 
sf52(double x, double u, double v)
{
    return u - 3 * v + 3 * exp(x);
}

void 
output(double *f, double x, double h, int n)
{
    for (int i = 0; i < n; ++i) {
        printf("%15.6f %15.6f\n", x, f[i]);
        x += h;
    }
}

//Runge-Kutta second-order method implementation for the differential equation
//Returns an array of n values of the approximating function in the grid nodes
double *
Runge_Kutta2(double a, double b, double yy0, int n, double alpha, double (*f)(double, double))
{
    double *y = realloc(NULL, (n + 1) * sizeof(*y));
    if (y == NULL) { 
        return NULL;
    }

    y[0] = yy0;
    double h = (b - a) / n, x = a, ff;

    for (int i = 0; i < n; ++i) {
        //y_i+1 = y_i + ((1 - a)f(x_i, y_i) + a*f( x_i + h/2a, y_i + h/2a * f(x_i, y_i) ) )h
        //a - alpha
        ff = f(x, y[i]);
        //(1 - alpha) * ff + alpha * f(x + h / (2 * alpha), y[i] + h * ff / (2 * alpha) ) );
        y[i + 1] = y[i] + ( (1 - alpha) * ff + alpha * f(x + h / (2 * alpha), y[i] + h * ff / (2 * alpha) ) ) * h;
        x += h;
    }
    return y;
}

//Runge-Kutta fourth-order method for the differential equation
//Returns an array of n values of the approximating function in the grid nodes
double *
Runge_Kutta4(double a, double b, double yy0, int n, double (*f)(double, double))
{
    double *y = realloc(NULL, (n + 1) * sizeof(*y));
    if (y == NULL) {
        return NULL;
    }

    y[0] = yy0;
    double h = (b - a) / n, x = a, k1, k2, k3, k4;

    for (int i = 0; i < n; ++i) {
        //y_i+1 = y_i + h/6(k1 + 2k2 + 2k3 + k4)
        //k1 = f(x_i, y_i)
        //k2 = f(x_i + h/2, y_i + k1*h/2)
        //k3 = f(x_i + h/2, y_i + k2*h/2)
        //k4 = f(x_i + h, y_i + hk3)

        k1 = f(x, y[i]);
        k2 = f(x + h / 2, y[i] + k1 * h / 2);
        k3 = f(x + h / 2, y[i] + k2 * h / 2);
        k4 = f(x + h, y[i] + k3 * h);

        y[i + 1] = y[i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        x += h;
    }

    return y;
}

//Runge-Kutta second-order method implementation for the system of differential equations
//Returns two arrays of n values of the approximating functions in the grid nodes
double **
Runge_Kutta_sys2(double a, double b, double u0, double v0, int n, double alpha, 
double (*f1)(double, double, double), double (*f2)(double, double, double))
{
    double **uv = realloc(NULL, 2 * sizeof(*uv));
    if (uv == NULL) {
        return NULL;
    }
    //u
    uv[0] = realloc(NULL, (n + 1) * sizeof(*uv));
    if (uv[0] == NULL) {
        free(uv);
        return NULL;
    }
    //v
    uv[1] = realloc(NULL, (n + 1) * sizeof(*uv));
    if (uv[1] == NULL) {
        free(uv[0]);
        free(uv);
        return NULL;
    }

    uv[0][0] = u0;
    uv[1][0] = v0;
    double h = (b - a) / n, x = a, ff1, ff2;

    for (int i = 0; i < n; ++i) {
        //u_i+1 = y_i + ((1 - alpha)f1(x_i, u_i, v_i) + alpha*f1(x_i + h/2alpha, u_i + h/2alpha * f1(x_i, u_i, v_i), v_i + h/2alpha * f2(x_i, y_i, v_i) ) )h
        //v_i+1 = v_i + ((1 - alpha)f2(x_i, u_i, v_i) + alpha*f2(x_i + h/2alpha, u_i + h/2alpha * f1(x_i, u_i, v_i), v_i + h/2alpha * f2(x_i, y_i, v_i) ) )h
        ff1 = f1(x, uv[0][i], uv[1][i]);
        ff2 = f2(x, uv[0][i], uv[1][i]);
        uv[0][i + 1] = uv[0][i] + ( (1 - alpha) * ff1 + alpha * f1(x + h / (2 * alpha),
                    uv[0][i] + h * ff1 / (2 * alpha), uv[1][i] + h * ff2 / (2 * alpha) ) ) * h;
        uv[1][i + 1] = uv[1][i] + ( (1 - alpha) * ff2 + alpha * f2(x + h / (2 * alpha),
                    uv[0][i] + h * ff1 / (2 * alpha), uv[1][i] + h * ff2 / (2 * alpha) ) ) * h;
        x += h;
    }

    return uv;
}

//Runge-Kutta fourth-order method implementation for the system of differential equations
//Returns an array of n values of the approximating functions in the grid nodes 
double **
Runge_Kutta_sys4(double a, double b, double u0, double v0, int n, double (*f1)(double, double, double), double (*f2)(double, double, double))
{
    double **uv = realloc(NULL, 2 * sizeof(*uv));
    if (uv == NULL) {
        return NULL;
    }
    //u
    uv[0] = realloc(NULL, (n + 1) * sizeof(*uv));
    if (uv[0] == NULL) {
        free(uv);
        return NULL;
    }
    //v
    uv[1] = realloc(NULL, (n + 1) * sizeof(*uv));
    if (uv[1] == NULL) {
        free(uv[0]);
        free(uv);
        return NULL;
    }

    uv[0][0] = u0;
    uv[1][0] = v0;
    double h = (b - a) / n, x = a, k1, k2, k3, k4, l1, l2, l3, l4;

    for (int i = 0; i < n; ++i) {
        //u_i+1 = u_i + h/6(k1 + 2k2 + 2k3 + k4)
        //v_i+1 = v_i + h/6(l1 + 2l2 + 2l3 + l4)

        //k1 = f1(x_i, u_i, v_i)
        //l1 = f2(x_i, u_i, v_i)
        k1 = f1(x, uv[0][i], uv[1][i]);
        l1 = f2(x, uv[0][i], uv[1][i]);

        //k2 = f1(x_i + h/2, u_i + k1*h/2, v_i + l1*h/2)
        //l2 = f2(x_i + h/2, u_i + k1*h/2, v_i + l1*h/2)
        k2 = f1(x + h / 2, uv[0][i] + k1 * h / 2, uv[1][i] + l1 * h / 2);
        l2 = f2(x + h / 2, uv[0][i] + k1 * h / 2, uv[1][i] + l1 * h / 2);

        //k3 = f1(x_i + h/2, u_i + k2*h/2, v_i + l2*h/2)
        //l3 = f2(x_i + h/2, u_i + k2*h/2, v_i + l2*h/2)
        k3 = f1(x + h / 2, uv[0][i] + k2 * h / 2, uv[1][i] + l2 * h / 2);
        l3 = f2(x + h / 2, uv[0][i] + k2 * h / 2, uv[1][i] + l2 * h / 2);

        //k4 = f1(x_i + h, u_i + k3*h, v_i + l3*h)
        //l4 = f2(x_i + h, u_i + k3*h, v_i + l3*h)
        k4 = f1(x + h, uv[0][i] + k3 * h, uv[1][i] + l3 * h);
        l4 = f2(x + h, uv[0][i] + k3 * h, uv[1][i] + l3 * h);

        uv[0][i + 1] = uv[0][i] + h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
        uv[1][i + 1] = uv[1][i] + h / 6 * (l1 + 2 * l2 + 2 * l3 + l4);

        x += h;
    }

    return uv;
}

int main(int argc, char *argv[])
{
    int n, num;
    double a, b, yy0, u0, v0, alpha;
    if (!strcmp(argv[1], "-e") || !strcmp(argv[1], "-equation")) {
        printf("The available differential equations:\n");
        printf("\n1) y' = -y + x ^ 2\n");
        printf("\n2) y' = y + 3x ^ 2 - 1\n");
        printf("\n3) y' = cos(x) + y\n");
        printf("\n4) y' = e ^ x - 3y\n");
        printf("\n5) (x - x ^ 2)y\n");
        printf("Choose the equation: ");
        scanf("%d", &num);

        printf("Enter the left border x0: ");
        scanf("%lf", &a);
        printf("Enter the  right border: ");
        scanf("%lf", &b);
        printf("Enter the initial condition y0 = y(x0): ");
        scanf("%lf", &yy0);
        printf("Enter the number of grid nodes (except x0): ");
        scanf("%d", &n);
        printf("Enter the Runge-Kutta 2nd-order method parameter (alpha): ");
        scanf("%lf", &alpha);

        void *function;
        switch (num) {
        case 1:
            function = f1;
            break;
        case 2:
            function = f2;
            break;
        case 3:
            function = f3;
            break;
        case 4:
            function = f4;
            break;
        default:
            function = f5;
        }

        if (!fork()) {
            double *rk2 = Runge_Kutta2(a, b, yy0, n, alpha, function);
            if (rk2 == NULL) {
                mem_error();
            }
            printf("Runge-Kutta2, Numerical solution:\n\n");
            printf("         x_i           y(x_i)\n");
            output(rk2, a, (b - a) / n, n + 1);
            printf("\n");
            free(rk2);
            return 0;
        }
        wait(NULL);

        if (!fork()) {
            double *rk4 = Runge_Kutta4(a, b, yy0, n, function);
            if (rk4 == NULL) {
                mem_error();
            }
            printf("\nRunge-Kutta4, Numerical solution:\n");
            printf("         x_i           y(x_i)\n");
            output(rk4, a, (b - a) / n, n + 1);
            printf("\n");
            free(rk4);
            return 0;
        }
        wait(NULL);
    } else if (!strcmp(argv[1], "-s") || !strcmp(argv[1], "-system")) {
        printf("The available systems of differential equations:\n");
        printf("\n1)\nu' = sin(x + u) - 1.1v\nv' = 2.5u - (x + v) ^ 2\n");
        printf("\n2)\nu' = (u - v)/x\nv' = 2v/x + ux\n");
        printf("\n3)\nu' = cos(u + 1.3v) - 2.1\nv' = 1.3/(x + 2.3u ^ 2) + x + 1\n");
        printf("\n4)\nu' = 2u + v + e ^ x\nv' = -2u + 2x\n");
        printf("\n5)\nu' = 2u - 4v\nv' = u - 3v + 3e ^ x\n");
        printf("Choose the system: ");
        scanf("%d", &num);
        printf("Enter the left border (x0): ");
        scanf("%lf", &a);
        printf("Enter the right border: ");
        scanf("%lf", &b);
        printf("Enter the initial condition u0 = u(x0): ");
        scanf("%lf", &u0);
        printf("Enter the initial condition v0 = v(x0): ");
        scanf("%lf", &v0);
        printf("Enter the number of grid nodes (except x0): ");
        scanf("%d", &n);
        printf("Enter the Runge-Kutta 2nd-order method parameter (alpha): ");
        scanf("%lf", &alpha);

        void *function1, *function2;
        switch (num) {
        case 1:
            function1 = sf11;
            function2 = sf12;
            break;
        case 2:
            function1 = sf21;
            function2 = sf22;
            break;
        case 3:
            function1 = sf31;
            function2 = sf32;
            break;
        case 4:
            function1 = sf41;
            function2 = sf42;
            break;
        default:
            function1 = sf51;
            function2 = sf52;
        }

        if (!fork()) {
            double **rk2 = Runge_Kutta_sys2(a, b, u0, v0, n, alpha, function1, function2);
            if (rk2 == NULL) {
                mem_error();
            }
            printf("Runge-Kutta2, Numerical solution:\nu(x):\n");
            printf("         x_i           u(x_i)\n");
            output(rk2[0], a, (b - a) / n, n + 1);
            printf("\nv(x):\n");
            printf("         x_i           v(x_i)\n");
            output(rk2[1], a, (b - a) / n, n + 1);
            printf("\n");
            free(rk2);
            return 0;
        }
        wait(NULL);

        if (!fork()) {

            double **rk4 = Runge_Kutta_sys4(a, b, u0, v0, n, function1, function2);
            if (rk4 == NULL) {
                mem_error();
            }
            printf("Runge-Kutta4, Numerical solution:\nu(x):\n");
            printf("         x_i           u(x_i)\n");
            output(rk4[0], a, (b - a) / n, n + 1);
            printf("\nv(x):\n");
            printf("         x_i           v(x_i)\n");
            output(rk4[1], a, (b - a) / n, n + 1);
            printf("\n");
            free(rk4);
            return 0;
        }
        wait(NULL);
    }
    return 0;
}
