/*
Numerical Methods
Pratical 1 variant 17 (Danilo Nascimento)

Description:
Solve the system of linear algebraic equations
with a tridiagonal matrix using the sweep method,
calculate the determinant

Author: Danilo Nascimento, M6O-311bki-21
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8
*/

#include <iostream>
#include <vector>
using namespace std;

/*
pass arguments to sweepMethod function as "const" to avoid data modification.
pass as references to avoid creating a copy of the arrays
thus improving efficiency, memory usage and overall performance

the indexes where changed by -1 compared to the given pseudocode,
this was because of the natural behavior of arrays, vectors, etc in C/C++
where they start from index 0 and not from 1.
*/
void SweepMethod(int n, const vector<double> &a, const vector<double> &b, const vector<double> &c, const vector<double> &d)
{
    vector<double> p(n);
    vector<double> q(n);
    double y = 0.0;

    p[0] = -c[0] / b[0];
    q[0] = d[0] / b[0];
    y = b[0];

    // Forward Path
    for (int i = 1; i < n; ++i)
    {
        /*
        Since this section gets repeated
        I stored it into a "temp" to reduce unecessary extra computations
        */
        double temp = b[i] + a[i] * p[i - 1];

        p[i] = -c[i] / temp;
        q[i] = (d[i] - a[i] * q[i - 1]) / temp;
        y *= temp;
    }

    vector<double> x(n);

    x[n - 1] = q[n - 1];

    // Backward Path
    for (int i = n - 2; i >= 0; --i)
    {
        x[i] = p[i] * x[i + 1] + q[i];
    }

    // Print the results
    cout << "x: ";
    for (int i = 0; i < n; ++i)
    {
        cout << x[i] << " ";
    }
    cout << endl;

    cout << "y: " << y << endl;
}

int main()
{
    int n = 5; // array size

    // input Matrix values
     vector<double> a = {0, -1, -9, -1, 9};
     vector<double> b = {-6, 13, -15, -7, -18};
     vector<double> c = {5, 6, -4, 1, 0};
     vector<double> d = {51, 100, -12, 47, -90};

    /*
    //     debug Raj results
    vector<double> a = {0, -8, 6, -8, 5};
    vector<double> b = {10, 16, -16, 16, -13};
    vector<double> c = {-1, 1, 6, -5, 0};
    vector<double> d = {16, -110, 24, -3, 87};
*/
    SweepMethod(n, a, b, c, d);

    return 0;
}
