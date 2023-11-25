/*
Numerical Methods
Pratical 7 variant 17 (Danilo Nascimento)

Description:
Using a table of values of the function y(x) at points xi, i=0,1,2,3
construct Lagrange and Newton interpolation polynomials,
Calculate the value of the interpolation error at point x*

Author: Danilo Nascimento, M6O-311Bki-21
Date: 14.11.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8
*/

#include <iostream>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>

using namespace std;

template <typename T>
void printVector(const vector<T> &vec)
{
    for (const T &element : vec)
    {
        cout << element << " ";
    }
    cout << endl;
}

double cot(double a)
{
    return cos(a) / sin(a);
}

// Lagrange Polynomial solver
void computeLagrangePolynomial(int n, vector<double> &x, vector<double> &y, double u)
{
    double v, f;
    v = 0;

    vector<double> w(n, 0.0);

    // Summands
    for (int j = 0; j < n; j++)
    {
        w[j] = y[j];
        f = y[j];

        // Multipliers
        for (int i = 0; i < n; i++)
        {
            if (i != j)
            {
                w[j] /= (x[j] - x[i]);
                f *= (u - x[i]) / (x[j] - x[i]);
            }
        }
        v += f;
    }

    // PRINT OUTPUT
    cout << "w = ";
    printVector(w);
    cout << "v = " << v << endl;
    // Value of interpolation error at point x*
    double epsilon = abs((v - (cot(u) + u)));
    cout << "Absolute Error at " << u << " = " << epsilon << endl;

    // sanity check (Computing the Vandermonde determinant)
    double det = 1;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            det *= (x[j] - x[i]);
        }
    }
    cout << "Sanity Check: " << det << endl;
}

// Newton Polynomial solver
void computeNewtonPolynomial(int n, vector<double> &x, vector<double> &y, double u)
{
    vector<vector<double>> d(n, vector<double>(n, 0.0));
    vector<double> w(n, 0.0);
    double f, v;

    // table of divided differences
    for (int i = 0; i < n; i++) // rows
    {
        d[i][0] = y[i];
    }
    for (int j = 1; j < n; j++) // columns
    {
        for (int i = 0; i < n - j; i++) // rows
        {
            d[i][j] = (d[i][j - 1] - d[i + 1][j - 1]) / (x[i] - x[i + j]);
        }
    }
    v = 0;
    for (int j = 0; j < n; j++) // summands
    {
        w[j] = d[0][j];
        f = d[0][j];

        for (int i = 0; i < j; i++) // multipliers
        {
            f *= (u - x[i]);
        }
        v += f;
    }

    // PRINT OUTPUT
    cout << "w = ";
    printVector(w);
    cout << "v = " << v << endl;
    // Value of interpolation error at point x*
    double epsilon = abs((v - (cot(u) + u)));
    cout << "Absolute Error at " << u << " = " << epsilon << endl;

    // sanity check (Computing the Vandermonde determinant)
    double det = 1;
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            det *= (x[j] - x[i]);
        }
    }
    cout << "Sanity Check: " << det << endl;
}

int main()
{
    cout << fixed << setprecision(4);

    int n = 4;

    vector<double> x(n, 0.0);
    x = {M_PI / 8, 2 * M_PI / 8, 3 * M_PI / 8, 4 * M_PI / 8};

    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        y[i] = cot(x[i]) + x[i];
    }
    double u = M_PI / 3;

    /*
      // Example from pdf
      int n = 4;

      vector<double> x(n, 0.0);
      x = {0.1, 0.5, 0.9, 1.3};

      vector<double> y(n, 0.0);
      for (int i = 0; i < n; i++)
      {
          y[i] = x[i] * log(x[i]);
      }

      double u = 0.7;
   */
    cout << "Lagrange Polynomial" << endl;
    computeLagrangePolynomial(n, x, y, u);

    cout << endl;

    cout << "Newton Polynomial" << endl;
    computeNewtonPolynomial(n, x, y, u);

    return 0;
}