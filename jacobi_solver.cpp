/*
This program solves a system of linear equations (A \cdot x = b)
using the iterative Jacobi method. It reads a square matrix (A), 
a free vector (b), and the desired accuracy (\varepsilon). 
If the method is not applicable (e.g., (A) is not diagonally dominant), 
it prints an error message. Otherwise, it computes the decomposition 
into (\alpha) and (\beta), iteratively approximates the solution vectors (x_i), 
prints each iteration with its current accuracy (e), and finally outputs 
the approximate solution (\tilde{x}). All numbers are formatted to four decimal places.
*/

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;
class Matrix
{
public:
    int n;
    vector<vector<double>> array;
    Matrix(int n)
    {
        this->n = n;
        this->array = vector<vector<double>>(n, vector<double>(n));

        // filling the vector with "empty" values (NaN)
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (j == 0)
                {
                    array[i][j] = 1;
                } else
                {
                    this->array[i][j] = nan("");
                }
            }
        }
    }

    virtual void input()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cin >> array[i][j];
            }
        }
    }

    void output()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n - 1; j++)
            {
                if (!isnan(this->array[i][j]))
                {
                    cout << fixed << setprecision(4) << array[i][j] << " ";
                }
            }
            if (!isnan(this->array[i][n - 1]))
            {
                cout << fixed << setprecision(4) << array[i][n - 1];
            }
            cout << '\n';
        }
    }
};

class ColumnVector : public Matrix
{
private:
    vector<double> vectorArr;
public:
    ColumnVector(int n) : Matrix(n)
    {
        this->vectorArr = vector<double>(n);
    }

    void input() override
    {
        for (int i = 0; i < n; i++)
        {
            cin >> vectorArr[i];
        }
    }

    void outputVector()
    {
        for (int j = 0; j < n - 1; j++)
        {
            cout << fixed << setprecision(4) << vectorArr[j] << endl;
        }
        cout << fixed << setprecision(4) << vectorArr[n - 1];
        cout << '\n';
    }

    ColumnVector operator-(const ColumnVector& vector) {
        ColumnVector result(n);
        for (int i = 0; i < n; i++) {
            result.vectorArr[i] = (double)vectorArr[i]-vector.vectorArr[i];
        }
        return result;
    }

    ColumnVector operator+(const ColumnVector& vector) {
        ColumnVector result(n);
        for (int i = 0; i < n; i++) result.vectorArr[i] = vectorArr[i] + vector.vectorArr[i];
        return result;
    }

    ColumnVector operator*(const Matrix& matrix)
    {
        ColumnVector result(n);
        for (int i = 0; i < matrix.n; i++)
        {
            for (int j = 0; j < matrix.n; j++)
            {
                result.vectorArr[i] += vectorArr[j] * matrix.array[i][j];
            }
        }
        return result;
    }

    double condition(const ColumnVector& vector) {
        ColumnVector tempVector(n);
        tempVector = this->operator-(vector);
        double result = 0;
        for (int i = 0; i < n; i++) {
            result += pow(tempVector.vectorArr[i], 2);
        }
        return sqrt(result);
    }

    void JacobiMethod(const Matrix& A, double e) {
        int size = A.n;
        vector<vector<double>> array;
        vector<vector<double>> a = A.array;
        vector<double> b = vectorArr;
        vector<double> beta = vector<double>(size, 0);
        // Let's express from all the rows the xᵢ:
        for (int i = 0; i < size; i++) {
            beta[i] = b[i] / a[i][i];
            double sum = 0;
            for (int j = 0; j < size; j++) {
                if (i != j) sum += abs(a[i][j]);
            }
            if (sum >= abs(a[i][i])) {
                cout << "The method is not applicable";
                exit(0);
            }
        }
        // αᵢⱼ = -аᵢⱼ / аᵢᵢ
        vector<vector<double>> alpha = vector<vector<double>>(size, vector<double>(size, 0));
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if (i != j) alpha[i][j] = -a[i][j] / a[i][i];
            }
        }
        // Next, we start an iterative process:
        ColumnVector xARRAY1(n);
        ColumnVector xARRAY2(n);
        Matrix MatrixALPHA(n);
        MatrixALPHA.array = alpha;
        cout << "alpha:" << endl;
        MatrixALPHA.output();
        ColumnVector VectorBETA(n);
        VectorBETA.vectorArr = beta;
        cout << "beta:" << endl;
        VectorBETA.outputVector();

        xARRAY1.vectorArr = beta;
        int k = 1;
        xARRAY2 = xARRAY1*MatrixALPHA+VectorBETA;
        ColumnVector Answer(n);
        Answer.vectorArr = xARRAY2.vectorArr;
        double ex = xARRAY2.condition(xARRAY1);
        while (ex > e) {
            cout << "x(" << k << ")" << endl;
            k += 1;
            xARRAY2.outputVector();
            ex = xARRAY2.condition(xARRAY1);
            cout << "e: " << fixed << setprecision(4) << ex << endl;
            Answer.vectorArr = xARRAY2.vectorArr;
            xARRAY1.vectorArr = xARRAY2.vectorArr;
            xARRAY2 = xARRAY1*MatrixALPHA+VectorBETA;
        }
        cout << "x~:" << endl;
        Answer.outputVector();
    }
};

int main()
{
    int n;
    cin >> n;
    Matrix A(n);
    A.input();
    int m;
    cin >> m;
    ColumnVector b(m);
    b.input();
    double e;
    cin >> e;
    b.JacobiMethod(A, e);
    return 0;
}
