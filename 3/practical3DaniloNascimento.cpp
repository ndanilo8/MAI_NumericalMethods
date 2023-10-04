/*
Numerical Methods
Pratical 3 variant 17 (Danilo Nascimento)

Description:
Solving SLAE using the fixed-point iteration method or the Seidel method.

Author: Danilo Nascimento, M6O-311bki-21
Date: 02.10.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8

*/

#include <iostream>
#include <vector>
#include <cmath> 
#include <chrono>
using namespace std;

// flag SeidelMethod to select solving mode FPI or Seidel (true, false)
void fixedPointIterationMethod_SeidelMethod(int n, vector<vector<double>> &c, vector<double> &d, double e, bool SeidelMethod)
{
    // copy of c matrix
    vector<vector<double>> A(n, vector<double>(n, 0.0)); // copy of c matrix
    vector<double> B(n);                                 // copy of d matrix
    A = c;                                               // assign the c matrix to A for the sanity check
    B = d;                                               // assign the d matrix to B for the sanity check

    // clock for timing
    auto start_time = chrono::high_resolution_clock::now();

    vector<vector<double>> a(n, vector<double>(n, 0.0));
    vector<double> b(n, 0.0);

    // Reduction to an equivalent form using the Jacobi method
    for (int i = 0; i < n; i++) // rows
    {
        if (c[i][i] == 0)
        {
            return;
        }
        for (int j = 0; j < n; j++) // columns
        {
            a[i][j] = (i != j) ? -c[i][j] / c[i][i] : 0;
        }
        b[i] = d[i] / c[i][i];
    }
    // Calculation of the norm of the matrix
    double g = 0.0;
    for (int i = 0; i < n; i++) // rows
    {
        double h = 0.0;
        for (int j = 0; j < n; j++)
        {
            h = (a[i][j] > 0) ? h + a[i][j] : h - a[i][j];
        }
        if (h > g)
            g = h;
    }

    if (g >= 1)
    {
        cout << "The Convergence condition is not fulfilled!" << endl;
        return;
    }
    double f = 0.0;
    vector<double> x(n, 0.0);

    for (int i = 0; i < n; i++) // rows
    {
        x[i] = b[i];

        if (x[i] > f)
        {
            f = x[i];
        }
        if (-x[i] > f)
        {
            f = -x[i];
        }
    }

    f = f * g / (1 - g);

    int k = 0;
    vector<double> y(n, 0.0);

    while (f > e) // Iteration
    {
        for (int i = 0; i < n; i++)
        {
            y[i] = x[i];
        }

        for (int i = 0; i < n; i++) // rows
        {
            x[i] = b[i];

            if (!SeidelMethod)
            {
                // For the fixed-point iteration method
                for (int j = 0; j < n; j++) // Element
                {
                    x[i] += a[i][j] * y[j];
                }
            }
            else
            {
                // Alternatively for the Seidel method
                for (int j = 0; j < i; j++) // element
                {
                    x[i] += a[i][j] * x[j];
                }
                for (int j = i; j < n; j++) // element
                {
                    x[i] += a[i][j] * y[j];
                }
            }
        }
        f = 0;

        for (int i = 0; i < n; i++)
        {
            if (x[i] - y[i] > f)
            {
                f = x[i] - y[i];
            }

            if (y[i] - x[i] > f)
            {
                f = y[i] - x[i];
            }
        }
        f = f * g / (1 - g);
        k++;
    }

    // Measure the elapsed time
    // deltaT = (final - initial)
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(end_time - start_time);

    if (!SeidelMethod)
    {
        cout << "Fixed Point Iteration Method took " << duration.count() << " microseconds." << endl;
    }
    else
    {
        cout << "Seidel Method took " << duration.count() << " microseconds." << endl;
    }

    // Output the results
    cout << "Solution vector 'x':" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "Number of Iterations: " << k << ", given the accuracy e = " << e << endl;

    // Sanity Check
    // Verify all calculations are correct by computing A*x-b
    vector<double> result(n);

    for (int i = 0; i < n; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            result[i] += A[i][j] * x[j];
        }
        result[i] -= B[i];
    }

    cout << "\nSanity Check (A * x - b):" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << result[i] << " ";
    }
    
    double tolerance = e;
    cout << "\nRounded Sanity check to: " << tolerance << " of tolerance" << endl;
       for (int i = 0; i < n; i++)
    {
        double roundedResult = round(result[i] / tolerance)* tolerance;
        cout << roundedResult << " ";
    }
    cout << endl;
}

int main()
{

    int n = 4;
    vector<vector<double>> c(n, vector<double>(n, 0.0));

    // variant 17
    c = {{-19, 2, -1, -8},
         {2, 14, 0, -4},
         {6, -5, -20, -6},
         {-6, 4, -2, 15}};

    vector<double> d(n, 0.0);
    d = {38, 20, 52, 43};
    double error = 0.01;

    // solve with Fixed-point iteration method
    cout << endl;
    fixedPointIterationMethod_SeidelMethod(n, c, d, error, false);
    // solve with Seidel method
    cout << endl;
    cout << "-----------------------------\n"
         << endl;
    fixedPointIterationMethod_SeidelMethod(n, c, d, error, true);
    cout << endl;

    return 0;
}