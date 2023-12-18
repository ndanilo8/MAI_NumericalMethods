/*
Numerical Methods
Pratical 10 variant 17 (Danilo Nascimento)

Description:
1- “Numerical integration”
program algorithm for the numerical
calculation of the value of a definite integral of a function on an interval with the
left-point and right-point rectangle rule, midpoint quadrature rule, trapezoid rule,
Simpson's rule.

2 - “Pyramid of refinements”
program algorithm for the numerical
calculation of the value of a definite integral of a function on an interval with the
trapezoid rule refined with the Runge-Romberg-Richardson rule.

Author: Danilo Nascimento, M6O-311Bki-21
Date: 14.12.2023
E-mail: ndanilo8@gmail.com
github: https://github.com/ndanilo8
*/

#include <iostream>
#include <vector>
#include <iomanip>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

#define VARIANT_17 // comment to run example from PDF

#define ARR_SIZE n

double f(double x)
{
#ifdef VARIANT_17
    return pow(M_E, x) + x;
#else
    return x * log(x);
#endif
}

void numericalIntegration(int n, vector<double> &x);

/*
@brief Numerical Calculation of the value of a definite integral of a function on an interval with the
trapezoid rule refined with the Runge-Romberg-Richardson rule.
@param n size of array
@param x 1D vector for the points

@return result the definite integral with trapezoid rule refined with the Runge-Romberg-Richardson rule.
*/
void pyramidOfRefinements(int n, vector<double> &x);

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

int main()
{
    cout << fixed << setprecision(4);

#ifdef VARIANT_17
    const int n = 5;

    vector<double> x(ARR_SIZE, 0.0);
    x = {-2, -1, 0, 1, 2};

#else
    const int n = 5;

    vector<double> x(ARR_SIZE, 0.0);
    x = {0.1, 0.5, 0.9, 1.3, 1.7};

#endif

    cout << "NUMERICAL INTEGRATION" << endl;
    numericalIntegration(n, x);
    cout << endl;
    cout << "PYRAMID OF REFINEMENTS" << endl;
    pyramidOfRefinements(n, x);

    return 0;
}

void numericalIntegration(int n, vector<double> &x)
{
    vector<double> y(n, 0.0);

    // values
    for (int i = 0; i < n; i++)
    {
        y[i] = f(x[i]);
    }

    double s = 0, t = 0, u = 0, v = 0, w = 0;

    // All nodes are equally spaced thus h = x_i - x_i-1
    const double h = x[1] - x[0];

    // Summands
    for (int i = 0; i < n - 1; i++)
    {
        /*
         s += y[2 * i - 1] * (x[2 * i + 1] - x[2 * i - 1]); // left point rectangle rule
         t += y[2 * i + 1] * (x[2 * i + 1] - x[2 * i - 1]); // right point rectangle rule
         u += y[2 * i] * (x[2 * i + 1] - x[2 * i - 1]); // midpoint quadrature rule
         v += (y[2 * i + 1] + y[2 * i - 1]) / 2 * (x[2 * i + 1] - x[2 * i - 1]);                // trapezoid rule
         w += (y[2 * i + 1] + 4 * y[2 * i] + y[2 * i - 1]) / 6 * (x[2 * i + 1] - x[2 * i - 1]); // simpsons rule
 */

        s += y[i];                                  // left point rectangle rule
        t += y[i + 1];                              // right point rectangle rule
        u += f(x[i] + h / 2);                       // midpoint quadrature rule
        v += y[i + 1] + y[i];                       // trapezoid rule
        w += y[i + 1] + 4 * f(x[i] + h / 2) + y[i]; // simpsons rule
    }

    s *= h;
    t *= h;
    u *= h;
    v *= h / 2;
    w *= h / 6;

#ifdef VARIANT_17
    const double I = (pow(x[n - 1], 2) / 2) + pow(M_E, x[n - 1]) - (pow(x[n - n], 2) / 2) + pow(M_E, x[n - n]);
#else
    const double I = (pow(x[n - 1], 2) / 4) * (2 * log(x[n - 1]) - 1) - (pow(x[n - n], 2) / 4) * (2 * log(x[n - n]) - 1);
#endif
    double s_error = abs(s - I);
    double t_error = abs(t - I);
    double u_error = abs(u - I);
    double v_error = abs(v - I);
    double w_error = abs(w - I);

    // print Output
    cout << "exact value of the definte integral I = " << I << endl;
    cout << endl;

    cout << "left-point rectangular rule" << endl;
    cout << "s = " << s << " | Abs. Interp. error = " << s_error << endl;
    cout << endl;

    cout << "right-point rectangular rule" << endl;
    cout << "t = " << t << " | Abs. Interp. error = " << t_error << endl;
    cout << endl;

    cout << "midpoint rule" << endl;
    cout << "u = " << u << "| Abs. Interp. error = " << u_error << endl;
    cout << endl;

    cout << "trapezoid rule" << endl;
    cout << "v = " << v << "| Abs. Interp. error = " << v_error << endl;
    cout << endl;

    cout << "Simposon's rule rule" << endl;
    cout << "w = " << w << "| Abs. Interp. error = " << w_error << endl;
}

void pyramidOfRefinements(int n, vector<double> &x)
{

    vector<double> y(n, 0.0);

    // values
    for (int i = 0; i < n; i++)
    {
        y[i] = f(x[i]);
    }

    int m = 2, p = 2;
    int q = p + log(n - 1) / log(m);
    int k = (q - p) + 1;

    if (k > q - p + 1)
    {
        return;
    }

    vector<vector<double>> d(k, vector<double>(k, 0.0));

    // rows
    for (int j = 1; j <= k; j++)
    {
        d[j - 1][0] = 0;
        q = pow(m, j - 1);

        // summands
        for (int i = 1; i <= (n - 1) / q; i++)
        {
            d[j - 1][0] += (y[q * i] + y[q * (i - 1)]) / 2 * (x[q * i] - x[q * (i - 1)]);
            //  cout << "(y[ " <<  q * i<< " ] + y[ " <<  q * (i - 1)<< " ]) / 2*(x[ " <<  q * i<< " ] - x[ " <<  q * (i - 1) << "])" << endl;
        }
    }

    // columns
    for (int j = 1; j < k; j++)
    {
        q = pow(m, p + j - 1);

        for (int i = 0; i < k - j; i++)
        {
            d[i][j] = d[i][j - 1] + (d[i][j - 1] - d[i + 1][j - 1]) / (q - 1);
        }
    }

    // print output
    cout << "d = " << endl;
    printVector(d);
}

// different method
/*
void pyramidOfRefinements(int n, vector<double> &x )
{
    vector<double> y(ARR_SIZE, 0.0);

    // values
    for (int i = 0; i < ARR_SIZE; i++)
    {
        y[i] = f(x[i]);
    }

    double v_h = 0;
    double v_2h = 0;
    double v_4h = 0;

    // All nodes are equally spaced thus h = x_i - x_i-1
    const double h = x[1] - x[0];

    const int m = 2; // partition ratio coefficient
    const int p = 2; // order of accuracy of the trapezoid rule
    const double q = p + log(n - 1) / log(m);
    const int k = q - p + 1;

#ifdef VARIANT_17
    const double I = (pow(x[n - 1], 2) / 2) + pow(M_E, x[n - 1]) - (pow(x[n - n], 2) / 2) + pow(M_E, x[n - n]);
#else
    const double I = (pow(x[n - 1], 2) / 4) * (2 * log(x[n - 1]) - 1) - (pow(x[n - n], 2) / 4) * (2 * log(x[n - n]) - 1);
#endif

    vector<vector<double>> d(k, vector<double>(k, 0.0));

    // Summands
    for (int i = 0; i < n - 1; i++)
    {
        v_h += y[i + 1] + y[i];
    }

    v_h *= h / 2;
    double v_h_error = abs(v_h - I);

    cout << "trapezoid rule" << endl;
    cout << "v = " << v_h << "| Abs. Interp. error = " << v_h_error << endl;
    cout << endl;

    v_2h = y[0] + 2 * y[2] + y[4];
    v_2h *= 2 * h / 2;
    double v_2h_error = abs(v_2h - I);

    cout << "trapezoid rule" << endl;
    cout << "v = " << v_2h << " | Abs. Interp. error = " << v_2h_error << endl;
    cout << endl;

    v_4h = y[0] + y[4];
    v_4h *= 4 * h / 2;
    double v_4h_error = abs(v_4h - I);

    cout << "trapezoid rule" << endl;
    cout << "v = " << v_4h << " | Abs. Interp. error = " << v_4h_error << endl;
    cout << endl;

    double l = (k * (k - 1)) / 2;

    double I_tilda_h_t = v_h + ((v_h - v_2h) / (pow(m, p) - 1));
    double I_tilda_h_t_error = abs(I_tilda_h_t - I);
    cout << "~I_h_t = " << I_tilda_h_t << " | Abs. Interp. error = " << I_tilda_h_t_error << endl;
    cout << endl;

    double I_tilda_2h_t = v_2h + ((v_2h - v_4h) / (pow(m, p) - 1));
    double I_tilda_2h_t_error = abs(I_tilda_2h_t - I);
    cout << "~I_2h_t = " << I_tilda_2h_t << " | Abs. Interp. error = " << I_tilda_2h_t_error << endl;
    cout << endl;

    double I_tilda_tilda_h_t = I_tilda_h_t + ((I_tilda_h_t - I_tilda_2h_t) / (pow(m, p+1) - 1));
    double I_tilda_tilda_h_t_error = abs(I_tilda_tilda_h_t - I);
    cout << "~~I_h_t = " << I_tilda_tilda_h_t << " | Abs. Interp. error = " << I_tilda_tilda_h_t_error << endl;
    cout << endl;

}
 */