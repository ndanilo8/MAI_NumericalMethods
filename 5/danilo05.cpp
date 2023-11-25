/*
Numerical Methods
Pratical 5 variant 17 (Danilo Nascimento)

Description:
Refining the value of the root of the equation (4^x)-5x-2
1. Interval halving method
2. Fixed-point iteration method
3. Newton method
4. Secant Method

Author: Danilo Nascimento, M6O-311Bki-21
Date: 23.10.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8
*/

#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>
#include <limits>
#include <functional>

using namespace std;

using func = function<double(double)>;

// ----------------------------------------------------------------------
// variant 17

// function f
const double f(double x)
{
    return pow(4, x) - 5 * x - 2; // // variant 17 function
}

// function phi
const double phi(double x)
{
    return log(5 * x + 2) / log(4);
}
// ----------------------------------------------------------------------

// function to find the derivative of any function of type 'func'
double derivative(func f, double x)
{
    double h = sqrt(numeric_limits<double>::epsilon()); // small step
    return (f(x + h) - f(x - h)) / (2.0 * h);
}


// Function to calculate NL equation with interval halving method
void interval_halving_solver(double a, double b, double e)
{
    if (f(a) * f(b) >= 0)
    {
        return;
    }

    double x = (a + b) / 2;
    int k = 0;

    while ((b - a) / 2 > e)
    {
        if (f(x) == 0)
        {
            break;
        }
        if (f(a) * f(x) < 0)
        {
            b = x;
        }
        if (f(x) * f(b) < 0)
        {
            a = x;
        }
        x = (a + b) / 2;
        k++;
    }

    cout << "'x': " << x << endl;
    cout << "Number of Interations 'k':" << k << endl;
    cout << "Sanity Check f(x) = " << f(x) << endl;
}

// Function to calculate NL equation with Fixed-point iteration method
void fixed_point_iteration_solver(double a, double b, double e)
{
    // calculate q = |max(phi'(x))| , in this case x = b (upper limit)
    double dphi = derivative(phi, b);
    double q = abs(dphi);

    double y, x, r = 0;
    int k = 0;

    if (q >= 1 || q <= 0)
    {
        return;
    }

    x = (a + b) / 2;
    r = 2 * e;
    k = 0;

    while (r > e)
    {
        y = x;
        x = phi(y);

        if (x < a || x > b)
        {
            return;
        }

        r = (x > y) ? x - y : y - x;

        r *= q / (1 - q);
        k++;
    }

    cout << "'x': " << x << endl;
    cout << "Number of Interations 'k':" << k << endl;
    cout << "Sanity Check f(x) = " << f(x) << endl;
}

// First Derivative of f(x)
double g(double x)
{
    return pow(4, x) * log(4) - 5;
}

// Second Derivative of f(x)
double h(double x)
{
    return pow(4, x) * log2(4);
}

// Function to calculate NL equation with newton method
void newton_solver(double a, double b, double e)
{
    double x, y;

    if (f(a) * h(a) > 0) // h(a) = second_derivative(f, a)
    {
        x = a;
    }
    else if (f(b) * h(b) > 0) // h(b) = second_derivative(f, b)
    {
        x = b;
    }
    else
    {
        return;
    }

    double r = 2 * e;
    int k = 0;

    while (r > e)
    {
        y = x;

        x = y - f(y) / g(y); // g(y) = derivative(f, y)

        r = (x > y) ? x - y : y - x;

        k++;
    }

    cout << "'x': " << x << endl;
    cout << "Number of Interations 'k':" << k << endl;
    cout << "Sanity Check f(x) = " << f(x) << endl;
}

// Function to calculate NL equation with secant method
void secant_solver(double c, double d, double e)
{
    double y, x, r = 0;
    double z = 0;

    y = c;
    x = d;
    r = 2 * e;

    int k = 1;

    while (r > e)
    {
        z = y;
        y = x;
        x = y - f(y) * (y - z) / (f(y) - f(z));
        r = (x > y) ? x - y : y - x;
        k++;
    }
    cout << "'x': " << x << endl;
    cout << "Number of Interations 'k':" << k << endl;
    cout << "Sanity Check f(x) = " << f(x) << endl;
}

int main()
{

    double a = 0;
    double b = 2;

    const double e = 0.001;

    cout << fixed << setprecision(4);

    cout << "INTERVAL HALVING METHOD" << endl;
    interval_halving_solver(a, b, e);

    cout << endl;

    cout << "FIXED POINT ITERATION METHOD" << endl;
    fixed_point_iteration_solver(a, b, e);

    cout << endl;

    cout << "NEWTON METHOD" << endl;
    newton_solver(a, b, e);

    cout << endl;

    cout << "SECANT METHOD" << endl;
    a = 2;
    b = 1.9;
    secant_solver(a, b, e);

    return 0;
}