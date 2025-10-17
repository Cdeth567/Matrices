/*
This program solves a system of linear equations (A \cdot x = b)
using Gaussian elimination with pivoting. It implements a `ColumnVector` 
class with operations for input/output, addition, multiplication, and norm calculation. 
During computation, it prints each elimination step and the transformed matrix
and vector under the section “Gaussian process:”, 
followed by diagonal normalization in “Diagonal normalization:” 
and the final solution in “[Result:](Result:)”. 
If the matrix is singular, it prints “Error: matrix A is singular.”
*/

#include <iostream>
#include <vector>
#include <iomanip>

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

    void output(vector<double> tempArr)
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n - 1; j++)
            {
                cout << fixed << setprecision(2) << array[i][j] << " ";
            }
            cout << fixed << setprecision(2) << array[i][n - 1];
            cout << '\n';
        }
        for (int j = 0; j < n - 1; j++)
        {
            cout << fixed << setprecision(2) << tempArr[j] << " ";
        }
        cout << fixed << setprecision(2) << tempArr[n - 1];
        cout << '\n';
    }

    void makePermutation1(int q, int w)
    {
        if (q != w)
        {
            swap(array[q], array[w]);
        }
    }

    void makeElimination1(int col)
    {
        for (int i = col+1; i < n; i++) // walk through the rows
        {
            double k = array[i][col] / array[col][col];
            if (array[i][col] != 0)
            {
                for (int j = 0; j < n; j++)
                {
                    array[i][j] -= k * array[col][j];
                }
            }
        }
    }

    void upperTriangular()
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
                makePermutation1(temp, count);
            }
            makeElimination1(count);
        }
    }

    vector<double> makePermutation(int q, int w, vector<double> tempArr)
    {
        if (q != w)
        {
            swap(tempArr[q], tempArr[w]);
            swap(array[q], array[w]);
        }
        return tempArr;
    }

    // array[i][j] - k * array[i-1][j]
    // array[i][j] - k * array[i-1][j] = 0 => k = array[i][j] / array[i-1][j]
    vector<double> makeElimination(int col, vector<double> tempArr, bool flag)
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
                        array[i][j] -= k * array[col][j];
                    }
                    tempArr[i] -= k * tempArr[col];
                    cout << "step #" << step << ": elimination" << endl;
                    output(tempArr);
                    step++;
                }
            }
        } else {
            for (int g = col - 1; g >= 0; g--) // walk through the rows
            {
                double k = array[g][col] / array[col][col];
                if (array[g][col] != 0) {
                    for (int j = n - 1; j >= 0; j--) {
                        array[g][j] -= k * array[col][j];
                    }
                    tempArr[g] -= k * tempArr[col];
                    cout << "step #" << step << ": elimination" << endl;
                    output(tempArr);
                    step++;
                }
            }
        }
        return tempArr;
    }

    int step = 1;
    vector<double> GaussianProcess(vector<double> tempArray)
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

    void findDeterminant()
    {
        double det = 1;
        for (int i = 0; i < n; i++)
        {
            det *= array[i][i];
        }
        if (det == 0)
        {
            cout << "Error: matrix A is singular";
            exit(0);
        }
    }

    vector<double> diagonalNormalization(vector<double> tempArray)
    {
        for (int i = 0; i < n; i++)
        {
            double k = array[i][i];
            for (int j = 0; j < n; j++)
            {
                if (array[i][j] != 0.00)
                {
                    array[i][j] = array[i][j] / k;
                }
            }
            if (tempArray[i] != 0.00)
            {
                tempArray[i] = tempArray[i] / k;
            }
        }
        cout << "Diagonal normalization:" << endl;
        output(tempArray);
        return tempArray;
    }

    void operator=(const Matrix& ma)
    {
        array = ma.array;
        n = ma.n;
    }
};

class IdentityMatrix : public Matrix
{
public:
    IdentityMatrix(int n) : Matrix(n)
    {
        this->n = n;
        for (int i = 0; i < n; i++)
        {
            array[i][i]=1;
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
        this->vectorArr = ::vector<double>(n);
    }

    void input() override
    {
        for (int i = 0; i < n; i++)
        {
            cin >> vectorArr[i];
        }
    }

    vector<double> getArrayVector()
    {
        return vectorArr;
    }

    void setArrayVector(vector<double> t)
    {
        vectorArr = std::move(t);
    }

    void outputVector()
    {
        for (int j = 0; j < n - 1; j++)
        {
            cout << fixed << setprecision(2) << vectorArr[j] << " ";
        }
        cout << fixed << setprecision(2) << vectorArr[n - 1];
        cout << '\n';
    }
};

int main()
{
    int n1, n2;
    cin >> n1;
    Matrix A(n1); // matrix A
    A.input();
    cin >> n2;

    ColumnVector vector(n2);
    vector.input();

    Matrix TempMatrix(n1);
    TempMatrix = A;
    TempMatrix.upperTriangular();
    TempMatrix.findDeterminant();

    cout << "Gaussian process:" << endl;
    vector.setArrayVector(A.diagonalNormalization(A.GaussianProcess(vector.getArrayVector())));

    cout << "Result:" << endl;
    vector.outputVector();
    return 0;
}
