/*
Numerical Methods
Pratical 6 variant 17 (Danilo Nascimento)

Description:
Solve the system of nonlinear equations (separate and refine a positive root)
using the fixed-point iteration method and the Newton method with an accuracy of ε=0.0001

Author: Danilo Nascimento, M6O-311Bki-21
Date: 30.10.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8
*/

#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

// main Functions ------------------------------------
// f1(x1,x2) 
double f(double x1, double x2)
{
    return 3 * x1 - cos(x2);
}

// f2(x1,x2)
double g(double x1, double x2)
{
    return 3 * x2 - exp(x1);
}

// Fixed Point Iteration f1
double phi1(double x1, double x2)
{
    return x1 - 0.1 * f(x1, x2);
}
// Fixed Point Iteration f2
double phi2(double x1, double x2)
{
    return x2 - 0.1 * g(x1, x2);
}
// END main Functions ------------------------------------

// partial derivatives ----------------------------
double f1(double x1, double x2)
{
    return 3;
}

double f2(double x1, double x2)
{
    return sin(x2);
}

double g1(double x1, double x2)
{
    return -exp(x1);
}
double g2(double x1, double x2)
{
    return 3;
}
// END partial derivatives ----------------------------

void FPI_seidel_solver(double a, double b, double c, double d, double e, bool seidelMethod)
{
    // double dphi = max(abs(f1(b, d)) + abs(f2(b, d)), abs(g1(b, d)) + abs(g2(b, d))); // this is > 1

    double q = 0.99; // manually setting the q to 0.99 gives a smallest # at sanity check // abs(dphi);

    if (q >= 1 || q <= 0)
    {
        return;
    }

    double x = (a + b) / 2;
    double y = (c + d) / 2;
    double r = 2 * e;
    int k = 0;

    double u, v, s;

    while (r > e)
    {
        u = x;
        v = y;

        if (!seidelMethod)
        {
            // for FPI Method
            x = phi1(u, v);
            y = phi2(u, v);
        }
        else
        {
            // for seidel Method
            x = phi1(u, v);
            y = phi2(x, v);
        }
        if (x < a || x > b || y < c || y > d)
        {
            return;
        }
        r = (x > u) ? x - u : u - x;
        s = (y > v) ? y - v : v - y;
        if (s > r)
        {
            r = s;
        }
        r *= q / (1 / q);
        k++;
    }

    // print results
    if (!seidelMethod)
    {
        cout << "Fixed-point Iteration Method" << endl;
    }
    else
    {
        cout << "Seidel Method" << endl;
    }
    cout << "'x': " << x << endl;
    cout << "'y': " << y << endl;
    cout << "Number of Iterations 'k': " << k << endl;
    cout << "Sanity Check f1(x1,x2) = " << f(x, y) << endl;
    cout << "Sanity Check f2(x1,x2) = " << g(x, y) << endl;
}

void newton_solver(double a, double b, double c, double d, double e)
{
    double x, y, r;
    int k;

    x = (a + b) / 2;
    y = (c + d) / 2;
    r = 2 * e;

    k = 0;

    double u, v, s;

    while (r > e)
    {
        u = x;
        v = y;

        s = f1(u, v) * g2(u, v) - f2(u, v) * g1(u, v);

        if (s == 0)
        {
            return;
        }

        x = u - (f(u, v) * g2(u, v) - g(u, v) * f2(u, v)) / s;
        y = v - (f1(u, v) * g(u, v) - f(u, v) * g1(u, v)) / s;

        r = (x > u) ? x - u : u - x;

        s = (y > v) ? y - v : v - y;

        if (s > r)
        {
            r = s;
        }

        k++;
    }

    cout << "Newton Method" << endl;
    cout << "'x': " << x << endl;
    cout << "'y': " << y << endl;
    cout << "Number of Iterations 'k': " << k << endl;
    cout << "Sanity Check f1(x1,x2) = " << f(x, y) << endl;
    cout << "Sanity Check f2(x1,x2) = " << g(x, y) << endl;
}

int main()
{

    const double e = 0.0001;
    cout << fixed << setprecision(4);
    double a, b, c, d;

    a = 0.2;
    b = 0.5;
    c = 0.2;
    d = 0.5;

    // FPI Method
    FPI_seidel_solver(a, b, c, d, e, false);
    cout << endl;
    // Seidel Method
    FPI_seidel_solver(a, b, c, d, e, true);
    cout << endl;
    // Newton Method
    newton_solver(a, b, c, d, e);

    return 0;
}