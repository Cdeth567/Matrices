/*
This program extends the matrix framework to compute 
the inverse of a square matrix using Gaussian elimination
with pivoting by the maximum absolute element. It first prints
the augmented matrix, then performs and displays each elimination
and normalization step with labels. If the matrix is singular, 
it prints “Error: matrix A is singular.” 
All results are printed with two decimal places.
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

    void input()
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
                cout << fixed << setprecision(2) << array[i][j] << " ";
            }
            for (int j = 0; j < n - 1; j++)
            {
                cout << fixed << setprecision(2) << tempArr[i][j] << " ";
            }
            cout << fixed << setprecision(2) << tempArr[i][n - 1];
            cout << '\n';
        }
    }

    void output()
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
                        tempArr[g][j] -= k * tempArr[col][j];
                        array[g][j] -= k * array[col][j];
                    }
                    cout << "step #" << step << ": elimination" << endl;
                    output(tempArr);
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

int main()
{
    int n;
    cin >> n;
    Matrix A(n); // temporary matrix
    A.input();

    Matrix TempMatrix(n);
    TempMatrix = A;
    TempMatrix.upperTriangular();
    TempMatrix.findDeterminant();
    IdentityMatrix identityArray(n);

    cout << "Augmented matrix:" << endl;
    A.output(identityArray.getArray());

    cout << "Gaussian process:" << endl;
    A.setArray(A.diagonalNormalization(A.GaussianProcess(identityArray.getArray())));

    cout << "Result:" << endl;
    A.output();
    return 0;
}
