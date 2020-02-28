#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

class RegrPoli {
public:
    int N;      // Numar puncte
    int *x;  // Vector abscise
    double *y;  // Vector ordonate (rata)
    double *yN;
    int g;      // Grad polinom
    double *c;  // Vector coeficienti

    RegrPoli(int grad=0, int Nel=0, int xin[]=nullptr, double yin[]=nullptr)
    {
        N = Nel;
        g = grad;
        x = xin;//new int[N];
        y = yin;//new double[N];
        c = new double[g+1];
    }
    void print()
    {
        for(int i=0; i<N; i++)
            cout << x[i] << ' ' << y[i] << endl;
    }
    void printX()
    {
        for(int i=0; i<N; i++)
            cout << x[i] << ' ';
        cout << endl;
    }
    void printY()
    {
        for(int i=0; i<N; i++)
            cout << y[i] << ' ';
        cout << endl;
    }
    void printC()
    {
        for(int i=g; i>=0; i--)
            printf("a%d= %.8f\t", i, c[i]);
        cout << endl;
    }
    double predict(int year)
    {
        double prediction=0;
        for(int i=0; i<=g; i++)
        {
            prediction += c[i] * pow(year, i);
        }
        //cout << "An: " << year << '\t' << "Sinucideri: " << prediction << endl;
        return prediction;
    }

    void Fit()
{
    // Calculare valori xi^2n
    double X[2 * g + 1];
    for(int i = 0; i < 2 * g + 1; i++)
    {
        X[i] = 0;
        for (int j = 0; j < N; j++)
        {
            X[i] += pow(x[j], i);
        }
    }
    // Calculare matrice normala augmentata
    double A[g + 1][g + 2];
    for (int i = 0; i <= g; i++)
    {
        for (int j = 0; j <= g; j++)
        {
            A[i][j] = X[i + j];
        }
    }
    // Calculare valori xi^n * yi
    double Y[g + 1];
    for (int i = 0; i < g + 1; i++)
    {
        Y[i] = 0;
        for (int j = 0; j < N; j++)
        {
            Y[i] += pow(x[j], i) * y[j];
        }
    }
    // Memorare valori pe ultima coloana a matricii A
    for (int i = 0; i <= g; i++)
    {
        A[i][g + 1] = Y[i];
    }
    g++;

    // Pivotare matrice A
    for (int i = 0; i < g; i++)
    {
        for (int k = i + 1; k < g; k++)
            if (A[i][i] < A[k][i])
                for (int j = 0; j <= g; j++) {
                    double temp = A[i][j];
                    A[i][j] = A[k][j];
                    A[k][j] = temp;
                }
    }
    // Eliminarea Gaussiana (creare zerouri sub pivoti)
    for (int i = 0; i < g - 1; i++)
    {
        for (int k = i + 1; k < g; k++)
        {
            double t = A[k][i] / A[i][i];
            for (int j = 0; j <= g; j++)
            {
                A[k][j] -= t * A[i][j];
            }
        }
    }
    // UTRIS (Rezolvare Ax = b)
    for (int i = g - 1; i >= 0; i--)
    {
        c[i] = A[i][g];
        for (int j = i + 1; j < g; j++)
        {
            c[i] -= A[i][j] * c[j];
        }
        c[i] /= A[i][i];
    }
    g--;
}
};
