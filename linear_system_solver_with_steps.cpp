/*
This program solves a system of linear equations (A \cdot x = b) 
using Gaussian elimination with full pivoting and prints all intermediate steps. 
It defines a `Matrix` class for square matrices with methods for row permutations, 
elimination, upper-triangular conversion, Gaussian process, and diagonal normalization. 
A `ColumnVector` class represents the free vector (b) with input/output and setters/getters. 
The program prints each elimination or permutation step with the updated matrix and vector, 
performs diagonal normalization, and outputs the final solution. 
If the matrix is singular, it prints an error message.
*/

#include <iostream>
#include <utility>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;
class Matrix
{
public:
    int n;
    int m;
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

    void setM(int m1)
    {
        this->m = m1;
    }

    int getM()
    {
        return m;
    }

    vector<double> inputArrayVector()
    {
        vector<double> vector = ::vector<double>(n, 0);
        for (int i = 0; i < n; i++)
        {
            cin >> array[i][1];
            cin >> vector[i];
        }
        return vector;
    }

    void fillingInTheCells(int m1)
    {
        m1 += 1;
        this-> m = m1;
        if (m1 >= 2)
        {
            for (int i = 0; i < n; i++)
            {
                for (int j = 2; j < m1; j++)
                {
                    array[i][j] = pow(array[i][1], j);
                }
            }
        }
    }

    void transposed()
    {
        vector<vector<double>> temp = vector<vector<double>>(m, vector<double>(n));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                temp[j][i] = array[i][j];
            }
        }
        this->array = vector<vector<double>>(m, vector<double>(n));
        int t = m;
        this->m = n;
        this->n = t;
        this->array = std::move(temp);
    }

    void output()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m - 1; j++)
            {
                if (!isnan(this->array[i][j]))
                {
                    cout << fixed << setprecision(4) << array[i][j] << " ";
                }
            }
            if (!isnan(this->array[i][m - 1]))
            {
                cout << fixed << setprecision(4) << array[i][m - 1];
            }
            cout << '\n';
        }
    }

    void clear()
    {
        int m1 = 0;
        vector<vector<double>> new_array = vector<vector<double>>(n, vector<double>(n));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                if (!isnan(this->array[i][j]))
                {
                    m1 += 1;
                    new_array[i][j] = array[i][j];
                }
            }
        }
        this->array = std::move(new_array);
    }

    vector<vector<double>> multiplication(int n2, int m2, vector<vector<double>> arr2)
    {
        vector<vector<double>> temp = vector<vector<double>>(n, vector<double>(m2, 0));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m2; j++)
            {
                for (int c = 0; c < m; c++)
                {
                    temp[i][j] += array[i][c] * arr2[c][j];
                }
            }
        }
        return temp;
    }

    void setArray(vector<vector<double>> t)
    {
        array = std::move(t);
    }

    vector<vector<double>> getArray()
    {
        return array;
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
        return tempArray;
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
        this->vectorArr = ::vector<double>(n);
    }

    void input()
    {
        for (int i = 0; i < n; i++)
        {
            cin >> vectorArr[i];
        }
    }

    void setArrayVector(vector<double> t)
    {
        vectorArr = std::move(t);
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

    void clearVector()
    {
        int m1 = 0;
        vector<double> new_vector = vector<double>(n);
        for (int i = 0; i < n; i++)
        {
            if (!isnan(this->vectorArr[i]))
            {
                m1 += 1;
                new_vector[i] = vectorArr[i];
            }
        }
        this->vectorArr = std::move(new_vector);
    }

    vector<double> vectorMultiplication(int n2, int m2, vector<vector<double>> arr2)
    {
        vector<double> tempVector = vector<double>(n2, 0);
        for (int i = 0; i < n2; i++)
        {
            for (int j = 0; j < m2; j++)
            {
                tempVector[i] += vectorArr[j] * arr2[i][j];
            }
        }
        return tempVector;
    }
};

int main()
{
    int m, n;
    cin >> m;
    Matrix matrix(m);
    ColumnVector Vector_b(m);
    Vector_b.setArrayVector(matrix.inputArrayVector());
    cin >> n;
    matrix.fillingInTheCells(n);
    cout << "A:" << endl;
    Matrix A(m);
    A.setM(matrix.getM());
    A.setArray(matrix.getArray());
    A.output();
    A.clear();
    matrix.clear();

    // transposed
    Matrix A_T(m);
    A_T.setM(matrix.getM());
    A_T.setArray(matrix.getArray());
    A_T.transposed();

    // multiplication
    cout << "A_T*A:" << endl;
    n += 1;
    Matrix MatrixA_TA(n);
    Matrix A_TCopy(n);
    A_TCopy.setArray(A_T.getArray());
    A_TCopy.setM(A_T.getM());
    MatrixA_TA.setM(n);
    MatrixA_TA.setArray(A_T.multiplication(m, n, matrix.getArray()));
    if (!MatrixA_TA.getArray().empty())
        MatrixA_TA.output();

    cout << "(A_T*A)^-1:" << endl;
    Matrix A_Inv(n);
    IdentityMatrix identityArray(n);
    identityArray.change();
    vector<vector<double>> tempArr = MatrixA_TA.GaussianProcess(identityArray.getArray()); // start for gaussian process
    vector<vector<double>> res = MatrixA_TA.diagonalNormalization(tempArr);
    MatrixA_TA.setArray(res);
    A_Inv.setArray(MatrixA_TA.getArray());
    A_Inv.setM(MatrixA_TA.getM());
    A_Inv.output();

    cout << "A_T*b:" << endl;
    Matrix MatrixA_Tb(m);
    MatrixA_Tb.setM(1);
    Vector_b.clearVector();
    ColumnVector Result(n);
    Result.setArrayVector(Vector_b.vectorMultiplication(n, m, A_TCopy.getArray()));
    Result.outputVector();

    cout << "x~:" << endl;
    ColumnVector Answer(n);
    Answer.setArrayVector(Result.vectorMultiplication(n, n, A_Inv.getArray()));
    Answer.outputVector();
    return 0;
}
