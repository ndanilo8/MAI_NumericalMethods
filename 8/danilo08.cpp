/*
Numerical Methods
Pratical 8 variant 17 (Danilo Nascimento)

Description:
Using a table of values of the function y(x) at points xi, i=0,1,2,3,4,5
1- construct a cubic spline interpolation polynomial,
2- calculate the value of the interpolation error at point x*;

3- solve the normal system of the method of least squares to
construct approximating polynomials of the 1st and 2nd degree,
4- calculate the sum of squared errors for each of the approximating polynomials

Author: Danilo Nascimento, M6O-311Bki-21
Date: 20.11.2023
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

template <typename T>
void printVector(const vector<vector<T>> &vec2D)
{
    for (const vector<T> &row : vec2D)
    {
        for (const T &element : row)
        {
            cout << element << " ";
        }
        cout << endl;
    }
}

double cot(double a)
{
    return cos(a) / sin(a);
}

void computeCubicSpline(int n, vector<double> &x, vector<double> &y, double u)
{
    double copy_u = u;
    vector<double> a(n - 1, 0.0);
    vector<double> b(n - 1, 0.0);
    vector<double> c(n - 1, 0.0);
    vector<double> d(n - 1, 0.0);
    vector<double> p(n - 2, 0.0);
    vector<double> q(n - 2, 0.0);

    int j = 0;
    double v = 0;

    // "Coefficients a" loop
    for (int i = 0; i < n - 1; i++)
    {
        a[i] = y[i];
    }

    c[0] = 0;

    // Sweep method for coefficients c
    p[0] = -(x[2] - x[1]) / (2 * x[2] - 2 * x[0]);
    q[0] = 3 * ((y[2] - y[1]) / (x[2] - x[1]) - (y[1] - y[0]) / (x[1] - x[0])) / (2 * x[2] - 2 * x[0]);

    // "Forward Path" loop
    for (int i = 1; i < n - 2; i++)
    {
        p[i] = -(x[i + 2] - x[i + 1]) / ((2 * x[i + 2] - 2 * x[i]) + (x[i + 1] - x[i]) * p[i - 1]);
        q[i] = (3 * ((y[i + 2] - y[i + 1]) / (x[i + 2] - x[i + 1]) - (y[i + 1] - y[i]) / (x[i + 1] - x[i])) - (x[i + 1] - x[i]) * q[i - 1]) / ((2 * x[i + 2] - 2 * x[i]) + (x[i + 1] - x[i]) * p[i - 1]);
    }

    c[n - 2] = q[n - 3];

    // "Backward Path" loop
    for (int i = n - 3; i >= 1; i--)
    {

        c[i] = p[i - 1] * c[i + 1] + q[i - 1];
    }

    // "Coefficients b and d" loop
    for (int i = 0; i < n - 2; i++)
    {

        b[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (x[i + 1] - x[i]) * (c[i + 1] + 2 * c[i]) / 3;
        d[i] = (c[i + 1] - c[i]) / (x[i + 1] - x[i]) / 3;
    }

    b[n - 2] = (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]) - (x[n - 1] - x[n - 2]) * 2 * c[n - 2] / 3;
    d[n - 2] = -c[n - 2] / (x[n - 1] - x[n - 2]) / 3;

    j = 0;

    // "Interval" loop
    for (int i = 0; i < n - 2; i++)
    {
        if (u >= x[i] && u <= x[i + 1])
        {
            j = i;
        }
    }

    if (j == 0)
    {
        return; // End the program if j is still 0 after the loop
    }

    u -= x[j];
    v = a[j] + b[j] * u + c[j] * u * u + d[j] * u * u * u;

    cout << "a = " ;//<< endl;
    printVector(a);
    cout << "b = " ;//<< endl;
    printVector(b);
    cout << "c = " ;//<< endl;
    printVector(c);
    cout << "d = " ;//<< endl;
    printVector(d);
    cout << "v = " << v << endl;

    // Value of interpolation error at point x*
    double epsilon = abs(v - (cot(copy_u) + copy_u));
    cout << "Absolute Error at " << copy_u << " = " << epsilon << endl;
}

void computeLeastSquares(int n, vector<double> &x, vector<double> &y, int m)
{
    vector<vector<double>> a(m, vector<double>(m, 0.0));
    vector<double> b(m, 0.0);
    vector<double> z(m, 0.0);
    double g = 0;
    double f = 0;

    // normal system of the method of least squares
    // rows
    for (int i = 0; i < m; i++)
    {
        // cout << i << endl;
        //  columns
        for (int j = 0; j < m; j++)
        {
            // cout << "\t" << j << endl;
            a[i][j] = 0;

            // summands
            for (int k = 0; k < n; k++)
            {
                // cout << "\t\t" << k << endl;
                a[i][j] += pow(x[k], i + j);
            }
        }

        b[i] = 0;
        // summands
        for (int k = 0; k < n; k++)
        {
            b[i] += +y[k] * pow(x[k], i);
        }
    }
    // cramer's method for the normal system
    if (m == 1)
    {
        z[0] = b[0] / a[0][0];
    }
    else if (m == 2)
    {
        g = a[0][0] * a[1][1] - a[1][0] * a[0][1];
        z[0] = (b[0] * a[1][1] - b[1] * a[1][0]) / g;
        z[1] = (a[0][0] * b[1] - a[1][0] * b[0]) / g;
    }
    else if (m == 3)
    {
        g = a[0][0] * a[1][1] * a[2][2] + a[1][0] * a[2][1] * a[0][2] + a[2][0] * a[0][1] * a[1][2] - a[2][0] * a[1][1] * a[0][2] - a[0][0] * a[2][1] * a[1][2] - a[1][0] * a[0][1] * a[2][2];
        z[0] = (b[0] * a[1][1] * a[2][2] + b[1] * a[2][1] * a[0][2] + b[2] * a[0][1] * a[1][2] - b[2] * a[1][1] * a[0][2] - b[0] * a[2][1] * a[1][2] - b[1] * a[0][1] * a[2][2]) / g;
        z[1] = (a[0][0] * b[1] * a[2][2] + a[1][0] * b[2] * a[0][2] + a[2][0] * b[0] * a[1][2] - a[2][0] * b[1] * a[0][2] - a[0][0] * b[2] * a[1][2] - a[1][0] * b[0] * a[2][2]) / g;
        z[2] = (a[0][0] * a[1][1] * b[2] + a[1][0] * a[2][1] * b[0] + a[2][0] * a[0][1] * b[1] - a[2][0] * a[1][1] * b[0] - a[0][0] * a[2][1] * b[1] - a[1][0] * a[0][1] * b[2]) / g;
    }

    vector<double> aproxPolynomial(n, 0.0);

    // sum of squared errors
    f = 0;
    // points
    for (int k = 0; k < n; k++)
    {
        g = 0;
        // summands
        for (int i = 0; i < m; i++)
        {
            g += z[i] * pow(x[k], i);
        }
        f += (g - y[k]) * (g - y[k]); // result of sum of squared errors
        aproxPolynomial[k] = g; // approx polynomial at the interval 
    }

    cout << fixed << setprecision(2); // just to make the matrix with less decimals places
    cout << "a = " << endl;           // (before coeffs) result of the system of the method of least squares
    printVector(a);

    cout << endl;
    cout << fixed << setprecision(4);

    cout << "b = "; // ( after coeffs ) result of the system of the method of least squares
    printVector(b);

    cout << "z = "; // coefficients
    printVector(z);

    cout << "f = " << f << endl; // sum of squared errors phi1,2

    cout << endl;

    cout << "P"<< m-1 << "(x)= " ;
    printVector(aproxPolynomial); // this prints the approximation polynomial at the given intervcl
}

int main()
{
    cout << fixed << setprecision(4);

      /* Example from PDF
    int n = 6;

    vector<double> x(n, 0.0);
    x = {0.1, 0.5, 0.9, 1.3, 1.7, 2.1};

    vector<double> y(n, 0.0);

    for (int i = 0; i < n; i++)
    {
        y[i] = x[i] * log(x[i]);
    }

    double u = 0.7;
    int m = 1;
*/
    // Variant 17 Danilo Nascimento
  
    int n = 6;

    vector<double> x(n, 0.0);
    x = {M_PI / 8, 2 * M_PI / 8, 3 * M_PI / 8, 4 * M_PI / 8, 5 * M_PI / 8, 6 * M_PI / 8};

    vector<double> y(n, 0.0);
    for (int i = 0; i < n; i++)
    {
        y[i] = cot(x[i]) + x[i];
    }
    double u = M_PI / 3;
    int m = 1;
    

    cout << "-------------Cubic Spline-------------" << endl;
    computeCubicSpline(n, x, y, u);
    cout << endl;
    cout << "-------------Least Squares 1st Degree-------------" << endl;
    computeLeastSquares(n, x, y, m + 1);
    cout << endl;
    cout << "-------------Least Squares 2nd Degree-------------" << endl;
    computeLeastSquares(n, x, y, m + 2);

    return 0;
}