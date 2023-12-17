/*
Numerical Methods
Pratical 9 variant 17 (Danilo Nascimento)

Description:
Using a table of values of the function y(x)
at points xi, i=0,1,2,3,4 construct a Newton interpolation polynomial
to calculate numerically the first and second derivatives of the function at point x*,
calculate the value of the interpolation error at point x* for each of the derivatives:

Author: Danilo Nascimento, M6O-311Bki-21
Date: 07.11.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8
*/

#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

// comment to run example from PDF "Lecture 09" instead
//#define VARIANT_17

using namespace std;

/*
@brief Calculates the values of the first and second derivatives of a function at a given point using a Newton interpolation polynomial on the entire table of function values.

@param n size of x values
@param x 1D vector for the points
@param y 1D vector for priori-calculated values of x with the function y(x)
@param u point x* 

@return Print (v), (w) - 1st, 2nd Derivate AND Interpolation errors at point x* respectively
*/
void numericalDifferentiation(int n, vector<double> &x, vector<double> &y, double u);


/*
@brief divides 2 integers
@param a first number 
@param b second number

@return result of division
*/
int DIV(int a, int b)
{
    return a / b;
}
/*
@brief calculates the remainder of 2 integers
@param a first number 
@param b second number

@return remainder
*/
int REMAINDER(int a, int b)
{
    return a % b;
}


double firstDerivative(double u);
double secondDerivative(double u);

int main()
{
    cout << fixed << setprecision(4);

// Variant 17 Danilo Nascimento -----------------
#ifdef VARIANT_17
    int n = 5;

    vector<double> x(n, 0.0);
    x = {-2, -1, 0, 1, 2};

    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        y[i] = exp(x[i]) + x[i];
    }
    double u = -2;
    

#else
    // Example PDF ----------------

    int n = 5;

    vector<double> x(n, 0.0);
    x = {0.1, 0.5, 0.9, 1.3, 1.7};

    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        y[i] = x[i] * log(x[i]);
    }

    double u = 0.9;

#endif

    numericalDifferentiation(n, x, y, u);

    return 0;
}

void numericalDifferentiation(int n, vector<double> &x, vector<double> &y, double u)
{

    vector<vector<double>> d(n, vector<double>(n, 0.0));
    double f = 0, g = 0, v = 0, w = 0;

    // table of divided differences
    // rows
    for (int i = 0; i < n; i++)
    {
        d[i][0] = y[i];
    }
    // columns
    for (int j = 1; j < n; j++)
    {
        // rows
        for (int i = 0; i < n - j; i++)
        {
            d[i][j] = (d[i][j - 1] - d[i + 1][j - 1]) / (x[i] - x[i + j]);
        }
    }

    // FIRST DERIVATIVE
    v = 0;
    // summands
    for (int j = 1; j < n; j++)
    {
        f = 0;
        // summands
        for (int k = 0; k < j; k++)
        {
            g = d[0][j];

            // multipliers
            for (int i = 0; i < j; i++)
            {
                if (i != k)
                {
                    g *= (u - x[i]);
                }
            }
            f += g;
        }
        v += f;
    }

    // SECOND DERIVATIVE
    w = 0;

    for (int j = 3; j <= n; j++) // summands
    {
        f = 0;
        for (int k = 1; k <= (j - 1) * (j - 2); k++) // summands
        {
            g = d[0][j - 1];
            for (int i = 1; i < j; i++) // multipliers
            {
                if ((i != 1 + DIV(k - 1, j - 2)) && (i != 1 + REMAINDER(k + DIV(k - 1, j - 1), j - 1)))
                {
                    g *= (u - x[i - 1]);
                    //  cout << j << " " << k << " " << i << endl; // DEBUG
                }
            }
            f += g;
        }
        w += f;
    }

// PRINT OUTPUTS
    // print output v and w
    cout << "v = " << v << endl;
    cout << "w = " << w << endl; 
    cout << endl;

    // Value of interpolation error at point x*
    double epsilon1 = abs(v - firstDerivative(u));
    cout << "1st Derivative:" << endl;
    cout << "Absolute Error at " << u << " = " << epsilon1 << endl;

    double epsilon2 = abs(w - secondDerivative(u));
    cout << "2nd Derivative:" << endl;
    cout << "Absolute Error at " << u << " = " << epsilon2 << endl;
}

double firstDerivative(double u)
{
// Variant 17
#ifdef VARIANT_17
    return exp(u) + 1;

#else
    // Example PDF
    return log(u) + 1;
#endif
}

double secondDerivative(double u)
{
// Variant 17
#ifdef VARIANT_17
    return exp(u);

#else
    // Example PDF
    return 1 / u;
#endif
}