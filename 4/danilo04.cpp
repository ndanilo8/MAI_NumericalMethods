/*
Numerical Methods
Pratical 4 variant 17 (Danilo Nascimento)

Description:
algorithm for finding the eigenvalues and eigenvectors
of a symmetric matrix using the Jacobi rotation method.

Author: Danilo Nascimento, M6O-311bki-21
Date: 10.10.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8

*/
#include <iostream>
#include <vector>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

// template function with overloading that can print 1d and 2d vectors
// this will allow to use the same printing function even of different data types for next practical works

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

void jacobiRotationMethodSolver(int n, vector<vector<double>> &a, double e)
{
    vector<vector<double>> A(n, vector<double>(n, 0.0));
    A = a; // copy of a matrix

    vector<vector<double>> v(n, vector<double>(n, 0.0));

    // ROWS
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            v[i][j] = (i == j) ? 1 : 0;
        }
    }
    double f = 0.0;

    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            f += a[i][j] * a[i][j];
        }
    }
    f = sqrt(f);
    int k = 0;

    vector<vector<double>> u(n, vector<double>(n, 0.0));
    vector<vector<double>> b(n, vector<double>(n, 0.0));
    double h = 0.0;

    // iterations
    while (f > e)
    {
        double g = 0.0;
        int l = 1;
        int m = 2;

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (a[i][j] > g)
                {
                    g = a[i][j];
                    l = i;
                    m = j;
                }
                if (-a[i][j] > g)
                {
                    g = -a[i][j];
                    l = i;
                    m = j;
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                u[i][j] = (i == j) ? 1 : 0;
            }
        }

        if (a[l][l] == a[m][m])
        {
            h = M_PI / 4;
        }
        else
        {
            h = atan(2 * a[l][m] / (a[l][l] - a[m][m])) / 2.0;
        }

        u[l][l] = cos(h);
        u[l][m] = -sin(h);
        u[m][l] = sin(h);
        u[m][m] = cos(h);

        // ROWS
        for (int i = 0; i < n; i++)
        {
            // COLUMNS
            for (int j = 0; j < n; j++)
            {
                if (i == l || i == m || j == l || j == m)
                {
                    b[i][j] = 0;
                    // ELEMENT
                    for (int p = 0; p < n; p++)
                    {
                        b[i][j] += u[p][i] * a[p][j];
                    }
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == l || i == m || j == l || j == m)
                {
                    a[i][j] = 0;
                    // ELEMENT
                    for (int p = 0; p < n; p++)
                    {
                        a[i][j] += b[i][p] * u[p][j];
                    }
                }
            }
        }

        // ROWS
        for (int i = 0; i < n; i++)
        {
            // COLUMNS
            for (int j = 0; j < n; j++)
            {
                if (i == l || i == m || j == l || j == m)
                {
                    b[i][j] = 0;
                    // ELEMENT
                    for (int p = 0; p < n; p++)
                    {
                        b[i][j] += v[i][p] * u[p][j];
                    }
                }
            }
        }

        // ROWS
        for (int i = 0; i < n; i++)
        {
            // COLUMNS
            for (int j = 0; j < n; j++)
            {
                if (i == l || i == m || j == l || j == m)
                {
                    v[i][j] = b[i][j];
                }
            }
        }

        f = 0.0;

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                f += a[i][j] * a[i][j];
            }
        }

        f = sqrt(f);
        k++;
    }

    vector<double> w(n);

    // ROWS
    for (int i = 0; i < n; i++)
    {
        w[i] = a[i][i];
    }

    // PRINT RESULTS
    // print w, v, k
    cout << "Eigenvalues 'w':" << endl;
    printVector(w);
    cout << endl;

    cout << "Eigenvectors  'v':" << endl;
    printVector(v);
    cout << endl;

    cout << "Number of Iterations: " << k << ", given the accuracy e = " << e << endl;
    cout << endl;

    // sanity check A * v - w * v = 0
    vector<vector<double>> lhs(n, vector<double>(n, 0.0)); // Left Hand Side A * v
    vector<vector<double>> rhs(n, vector<double>(n, 0.0)); // Right Hand Side w * v
    vector<vector<double>> result(n, vector<double>(n, 0.0));

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            rhs[i][j] += w[j] * v[i][j]; // Perform the matrix-vector multiplication w_1x3 * V_3x3 (w * v)

            for (int k = 0; k < 3; k++)
            {
                lhs[i][j] += A[i][k] * v[k][j]; // Perform the matrix-vector multiplication A_3x3 * V_3x3 (A * v)
            }
        }
    }
    // perform lhs_3x3 - rhs_3x3 (A * v - w * v)
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] += lhs[i][j] - rhs[i][j];
        }
    }

    cout << "Sanity Check:" << endl;
    cout << "pair 1 "
         << "pair 2 "
         << "pair 3" << endl;
    printVector(result);
    cout << endl;
}

void maxModulusMatrixEigenvalue_eigenvectorSolver(int n, vector<vector<double>> &a, double e)
{
    vector<double> y(n);
    for (int i = 0; i < n; i++)
    {
        y[i] = 1;
    }
    double f = 2 * e;
    double w = 0.0;
    int k = 0;

    vector<double> z(n);
    // ITERATION
    while (f > e)
    {
        double r = 0;
        // ROWS
        for (int i = 0; i < n; i++)
        {
            z[i] = 0.0;
            // ELEMENT
            for (int j = 0; j < n; j++)
            {
                z[i] += a[i][j] * y[j];
            }
            if (z[i] > r)
            {
                r = z[i];
            }
            if (-z[i] > r)
            {
                r = -z[i];
            }
        }

        f = (k > 0) ? w : 2 * e + z[0] / y[0];

        w = z[0] / y[0];

        for (int i = 0; i < n; i++)
        {
            y[i] = z[i] / r;
        }

        f = w - f;
        if (f < 0)
        {
            f = -f;
        }
        k++;
    }

    // PRINT RESULTS
    // print w, y, k
    cout << "'w': " << w << endl;
    cout << endl;

    cout << "'y': ";
    printVector(y);
    cout << endl;

    cout << "Number of Iterations: " << k << ", given the accuracy e = " << e << endl;
    cout << endl;
}

int main()
{
    cout << fixed << setprecision(3);
    int n = 3;
    vector<vector<double>> a(n, vector<double>(n, 0.0));

    // variant 17
    a = {{5, -3, -4},
         {-3, -3, 4},
         {-4, 4, 0}};

    double e = 0.01;

    vector<vector<double>> A(n, vector<double>(n, 0.0));
    A = a; // copy of a matrix

    cout << "JACOBI ROTATION METHOD" << endl
         << endl;

    jacobiRotationMethodSolver(n, a, e);

    cout << "---------------------------------------" << endl
         << endl;

    cout << "POWER ITERATION METHOD" << endl
         << endl;

    maxModulusMatrixEigenvalue_eigenvectorSolver(n, A, e);

    return 0;
}