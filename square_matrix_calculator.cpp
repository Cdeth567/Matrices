/*
This program works with square matrices (n×n) 
using inheritance and operator overloading. It supports input/output, 
addition, subtraction, multiplication, and transposition using upcasting 
and downcasting with pointers. If matrix sizes are incompatible, 
it prints “Error: the dimensional problem occurred.”
*/

#include <iostream>
#include <vector>

using namespace std;
class Matrix
{
private:
    int n;
    vector<vector<int>> array;
public:
    Matrix(int n)
    {
        this->n = n;
        this->array = vector<vector<int>>(n, vector<int>(n));
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

    Matrix transposed()
    {
        vector<vector<int>> temp = vector<vector<int>>(n, vector<int>(n));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                temp[j][i] = array[i][j];
            }
        }
        Matrix ma(n);
        ma.array = temp;
        return ma;
    }

    void output()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n - 1; j++)
            {
                cout << array[i][j] << " ";
            }
            cout << array[i][n - 1];
            cout << '\n';
        }
    }

    Matrix operator+(const Matrix& ma)
    {
        if (ma.n != n)
        {
            cout << "Error: the dimensional problem occurred" << endl;
            vector<vector<int>> temp = vector<vector<int>>(0, vector<int>(0));
            Matrix Temporary(0);
            Temporary.array = temp;
            return Temporary;
        } else
        {
            vector<vector<int>> temp = vector<vector<int>>(n, vector<int>(n));
            vector<vector<int>> a1 = ma.array;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    temp[i][j] = array[i][j] + a1[i][j];
                }
            }
            Matrix Temporary(n);
            Temporary.array = temp;
            return Temporary;
        }
    }

    Matrix operator-(const Matrix& ma)
    {
        if (ma.n != n)
        {
            cout << "Error: the dimensional problem occurred" << endl;
            vector<vector<int>> temp = vector<vector<int>>(0, vector<int>(0));
            Matrix Temporary(0);
            Temporary.array = temp;
            return Temporary;
        } else
        {
            vector<vector<int>> temp = vector<vector<int>>(n, vector<int>(n));
            vector<vector<int>> a1 = ma.array;
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
    }

    Matrix operator*(const Matrix& ma)
    {
        if (n == ma.n)
        {
            vector<vector<int>> temp = vector<vector<int>>(n, vector<int>(ma.n, 0));
            vector<vector<int>> a1 = ma.array;
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
            vector<vector<int>> temp = vector<vector<int>>(0, vector<int>(0));
            Matrix Temporary(0);
            Temporary.array = temp;
            return Temporary;
        }
    }

    void operator=(const Matrix& ma)
    {
        array = ma.array;
        n = ma.n;
    }
};

int main()
{
    int n1, n2, n3;
    cin >> n1;
    Matrix A(n1);
    A.input();
    cin >> n2;
    Matrix B(n2);
    B.input();
    cin >> n3;
    Matrix C(n3);
    C.input();

    // D = A + B
    Matrix D = A + B;
    D.output();

    // E = B - A
    Matrix E = B - A;
    E.output();

    // F = C * A
    Matrix F = C * A;
    F.output();

    // G = A^T
    Matrix G = A.transposed();
    G.output();
    return 0;
}
