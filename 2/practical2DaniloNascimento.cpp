/*
Numerical Methods
Pratical 2 variant 17 (Danilo Nascimento)
"Gaussian method with pivot element selection"

Description:
Solve the system of linear algebraic equations 'x'
using the Gaussian method with pivot element selection,
calculate the determinant 'y' and the inverse matrix 'z'

Author: Danilo Nascimento, M6O-311bki-21
Date: 19.09.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8

*/

#include <iostream>
#include <vector>
#include <iomanip> // Include the <iomanip> header for the printer formatting of the matrices 'setprecision'
using namespace std;

// Function to print a matrix with a fixed number of decimal places
void printMatrix(vector<vector<double>> &mat)
{
    int n = mat.size();

    // Set the output format for floating-point numbers to fixed with 7 decimal places
    // ive used this because the format when printing the matrix wasnt corrent with only the \t
    // thus i defined a fixed size for all printing elements
    cout << fixed << setprecision(7);

    // Loop through rows
    for (int i = 0; i < n; i++)
    {
        // Loop through columns
        for (int j = 0; j < n; j++)
        {
            // Print the matrix element with fixed formatting
            cout << mat[i][j] << "\t";
        }
        // Print a newline character to move to the next row
        cout << endl;
    }
}

// Function to perform Gaussian elimination with pivot element selection
void gaussianPivotElementSelection(int n, vector<vector<double>> &a, vector<double> &b)
{
    vector<vector<double>> z(n, vector<double>(n, 0.0)); // Inverse Matrix
    vector<double> x(n);                                 // Solution vector
    double y = 0.0;                                      // determinant

    // Rows
    for (int i = 0; i < n; i++)
    {
        // Columns
        for (int j = 0; j < n; j++)
        {
            z[i][j] = (i == j) ? 1 : 0;
        }
    }

    int p = 0;

    double f = 0.0;
    double m = 0.0;

    // forward path
    for (int k = 0; k < n - 1; k++)
    {
        // For the Gaussian method with pivot element selection
        double g = (a[k][k] > 0) ? a[k][k] : -a[k][k];
        m = k;

        for (int i = k + 1; i < n; i++)
        {
            double f = (a[i][k] > 0) ? a[i][k] : -a[i][k];
            if (f > g)
            {
                g = f;
                m = i;
            }
        }
        if (g == 0)
        {
            cout << "Solution is not unique or not available!" << endl;
            return;
        }
        if (m != k)
        {
            for (int j = k; j < n; j++)
            {
                // this section bellow can be entirely replaced by
                // swap(a[k][j], a[m][j]);
                g = a[k][j];
                a[k][j] = a[m][j];
                a[m][j] = g;
            }
            // this section bellow can be entirely replaced by
            // swap(b[k], b[m]);
            g = b[k];
            b[k] = b[m];
            b[m] = g;

            for (int j = 0; j < n; j++)
            {
                // this section bellow can be entirely replaced by
                //  swap(z[k][j], z[m][j]);
                g = z[k][j];
                z[k][j] = z[m][j];
                z[m][j] = g;
            }
            p++;
        }

        // General part of the Gaussian method

        // rows
        for (int i = k + 1; i < n; i++)
        {
            double h = (-a[i][k]) / a[k][k];

            // columns
            for (int j = k; j < n; j++)
            {
                a[i][j] += h * a[k][j];
            }
            b[i] += h * b[k];
            // columns
            for (int j = 0; j < n; j++)
            {
                z[i][j] += h * z[k][j];
            }
        }
    }

    y = 1 - 2 * (p % 2);

    for (int k = n - 1; k >= 0; k--)
    {
        y *= a[k][k];
        double h = 0.0;

        for (int i = k + 1; i < n; i++)
        {
            h += a[k][i] * x[i];
        }
        x[k] = (b[k] - h) / a[k][k];

        for (int j = 0; j < n; j++)
        {
            h = 0.0;

            for (int i = k + 1; i < n; i++)
            {
                h += a[k][i] * z[i][j];
            }

            z[k][j] = (z[k][j] - h) / a[k][k];
        }
    }

    // Output the results
    cout << "\nSolution vector 'x':" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "\nDeterminant 'y': " << y << endl;

    cout << "\nInverse Matrix 'z':" << endl;
    printMatrix(z);

    // Sanity Check
    // Verify all calculations are correct by computing A*x-b
    // and if the result is 0 then its correct
    vector<double> result(n);
    for (int i = 0; i < n; i++)
    {
        result[i] = 0.0;
        for (int j = 0; j < n; j++)
        {
            result[i] += a[i][j] * x[j];
        }
        result[i] -= b[i];
    }

    cout << "\nSanity Check (A * x - b):" << endl;
    for (int i = 0; i < n; i++)
    {
        cout << result[i] << " ";
    }
    cout << endl;
}

int main()
{

    int n = 4;
    vector<vector<double>> a(n, vector<double>(n, 0.0));
    // Variant 17
    a = {
        {8, 8, -5, -8},
        {8, -5, 9, -8},
        {5, -4, -6, -2},
        {8, 3, 6, 6}

    };
    vector<double> b = {13, 38, 14, -95};

    gaussianPivotElementSelection(n, a, b);

    return 0;
}