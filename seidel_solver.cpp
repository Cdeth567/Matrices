/*
This program solves a system of linear equations (A \cdot x = b) 
using the iterative Seidel (Gauss–Seidel) method. It reads a square matrix (A),
a free vector (b), and the approximation accuracy. 
If the method is not applicable (e.g., matrix does not allow convergence), 
it prints an error message. Otherwise, it computes matrices alpha and beta, 
the lower-triangular (B) and upper-triangular (C) decompositions, (I-B), 
and its inverse ((I-B)^{-1}). The program iteratively computes approximation vectors (x_i), 
prints each step with its current accuracy (e), and finally outputs the approximate solution. 
All numbers are formatted to four decimal places.
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

    void output(vector<vector<double>> tempArr)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                cout << fixed << setprecision(4) << array[i][j] << " ";
            }
            for (int j = 0; j < n - 1; j++)
            {
                cout << fixed << setprecision(4) << tempArr[i][j] << " ";
            }
            cout << fixed << setprecision(4) << tempArr[i][n - 1];
            cout << '\n';
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

    void operator=(const Matrix& ma)
    {
        array = ma.array;
        n = ma.n;
    }

    Matrix operator-(const Matrix& ma)
    {
        vector<vector<double>> temp = vector<vector<double>>(n, vector<double>(n));
        vector<vector<double>> a1 = ma.array;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                temp[i][j] = array[i][j] - a1[i][j];
            }
        }
        Matrix Temporary(n);
        Temporary.array = temp;
        return Temporary;
    }

    vector<vector<double>> makePermutation(int q, int w, vector<vector<double>> tempArr)
    {
        if (q != w)
        {
            swap(tempArr[q], tempArr[w]);
            swap(array[q], array[w]);
        }
        return tempArr;
    }

    vector<vector<double>> makeElimination(int col, vector<vector<double>> tempArr, bool flag)
    {
        if (flag)
        {
            for (int i = col+1; i < n; i++) // walk through the rows
            {
                double k = array[i][col] / array[col][col];
                if (array[i][col] != 0)
                {
                    for (int j = 0; j < n; j++)
                    {
                        tempArr[i][j] -= k * tempArr[col][j];
                        array[i][j] -= k * array[col][j];
                    }
                    step++;
                }
            }
        } else {
            for (int g = col - 1; g >= 0; g--) // walk through the rows
            {
                double k = array[g][col] / array[col][col];
                if (array[g][col] != 0) {
                    for (int j = n - 1; j >= 0; j--) {
                        tempArr[g][j] -= k * tempArr[col][j];
                        array[g][j] -= k * array[col][j];
                    }
                    step++;
                }
            }
        }
        return tempArr;
    }

    // array[i][j] - k * array[i-1][j]
    // array[i][j] - k * array[i-1][j] = 0 => k = array[i][j] / array[i-1][j]
    int step = 1;
    vector<vector<double>> GaussianProcess(vector<vector<double>> tempArray)
    {
        int column = 0;
        for (int count = 0; count < n; count++) { // counter for column detection
            double mx = 0;
            int temp = 0;
            for (int i = count; i < n; i++) {
                if (abs(array[i][column]) > mx) {
                    mx = abs(array[i][column]);
                    temp = i;
                }
            }
            column += 1;
            if (temp > 0 && abs(mx) != abs(array[count][count])) {
                tempArray = makePermutation(temp, count, tempArray);
                cout << "step #" << step << ": permutation" << endl;
                output(tempArray);
                step++;
            }
            tempArray = makeElimination(count, tempArray, true);
        }
        // reversed
        for (int count = n - 1; count > 0; count--) { // counter for column detection
            tempArray = makeElimination(count, tempArray, false);
        }
        return tempArray;
    }

    vector<vector<double>> diagonalNormalization(vector<vector<double>> tempArray)
    {
        for (int i = 0; i < n; i++)
        {
            double k = array[i][i];
            for (int j = 0; j < n; j++)
            {
                array[i][j] /= k;
                tempArray[i][j] /= k;
            }
        }
        array = tempArray;
        return tempArray;
    }

    void lowerTriangular()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                array[i][j] = 0;
            }
        }
        output();
    }

    void upperTriangular()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                array[i][j] = 0;
            }
        }
        output();
    }

    Matrix operator*(const Matrix& ma)
    {
        if (n == ma.n)
        {
            vector<vector<double>> temp = vector<vector<double>>(n, vector<double>(ma.n, 0));
            vector<vector<double>> a1 = ma.array;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < ma.n; j++)
                {
                    for (int c = 0; c < n; c++)
                    {
                        temp[i][j] += array[i][c] * a1[c][j];
                    }
                }
            }
            Matrix Temporary(n);
            Temporary.array = temp;
            return Temporary;
        } else
        {
            cout << "Error: the dimensional problem occurred" << endl;
            vector<vector<double>> temp = vector<vector<double>>(0, vector<double>(0));
            Matrix Temporary(0);
            Temporary.array = temp;
            return Temporary;
        }
    }
};

class IdentityMatrix : public Matrix
{
public:
    IdentityMatrix(int n) : Matrix(n)
    {
        this->n = n;
        this->array = vector<vector<double>>(n, vector<double>(n, 0));
    }
    vector<vector<double>> change()
    {
        for (int i = 0; i < n; i++)
        {
            array[i][i]=1;
        }
        return array;
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

    void SeidelMethod(const Matrix& A, double e) {
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

        Matrix B(n);
        B = MatrixALPHA;
        cout << "B:" << endl;
        B.lowerTriangular();
        Matrix C(n);
        C = MatrixALPHA;
        cout << "C:" << endl;
        C.upperTriangular();

        IdentityMatrix I(n);
        I.change();
        cout << "I-B:" << endl;
        Matrix IB = I - B;
        IB.output();
        cout << "(I-B)_-1:" << endl;
        IB.diagonalNormalization(IB.GaussianProcess(I.array));
        IB.output();

        xARRAY1.vectorArr = beta;
        int k = 1;
        Matrix IBC = IB * C;
        xARRAY2 = xARRAY1 * IBC + VectorBETA * IB;
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
            IBC = IB * C;
            xARRAY2 = xARRAY1 * IBC + VectorBETA * IB;
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
    b.SeidelMethod(A, e);
    return 0;
}
