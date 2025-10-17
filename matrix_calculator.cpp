/*
This program defines a Matrix class that supports 
storing n×m integer matrices and performing basic 
operations using operator overloading. It overloads 
operators for input/output (>> and <<), assignment (=),
addition (+), subtraction (-), and multiplication (*).
It also includes a method to compute the transpose of a matrix. 
If matrix dimensions are incompatible for an operation, it prints 
"Error: the dimensional problem occurred". The program reads three 
matrices A, B, and C, then outputs: D = A + B, E = B − A, F = C * A, 
and G = AT (transpose of A).
*/

#include <iostream>
#include <vector>

using namespace std;
class Matrix
{
private:
    int n;
    int m;
    vector<vector<int>> array;
public:
    Matrix(int n, int m)
    {
        this->n = n;
        this->m = m;
        this->array = vector<vector<int>>(n, vector<int>(m));
    }

    void input()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                cin >> array[i][j];
            }
        }
    }

    Matrix transposed()
    {
        vector<vector<int>> temp = vector<vector<int>>(m, vector<int>(n));
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m; j++)
            {
                temp[j][i] = array[i][j];
            }
        }
        Matrix ma(m, n);
        ma.array = temp;
        return ma;
    }

    void output()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < m - 1; j++)
            {
                cout << array[i][j] << " ";
            }
            cout << array[i][m - 1];
            cout << '\n';
        }
    }

    Matrix operator+(const Matrix& ma)
    {
        if (!(ma.n == n && ma.m == m))
        {
            cout << "Error: the dimensional problem occurred" << endl;
            vector<vector<int>> temp = vector<vector<int>>(0, vector<int>(0));
            Matrix Temporary(0, 0);
            Temporary.array = temp;
            return Temporary;
        } else
        {
            vector<vector<int>> temp = vector<vector<int>>(n, vector<int>(m));
            vector<vector<int>> a1 = ma.array;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    temp[i][j] = array[i][j] + a1[i][j];
                }
            }
            Matrix Temporary(n, m);
            Temporary.array = temp;
            return Temporary;
        }
    }

    Matrix operator-(const Matrix& ma)
    {
        if (!(ma.n == n && ma.m == m))
        {
            cout << "Error: the dimensional problem occurred" << endl;
            vector<vector<int>> temp = vector<vector<int>>(0, vector<int>(0));
            Matrix Temporary(0, 0);
            Temporary.array = temp;
            return Temporary;
        } else
        {
            vector<vector<int>> temp = vector<vector<int>>(n, vector<int>(m));
            vector<vector<int>> a1 = ma.array;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    temp[i][j] = array[i][j] - a1[i][j];
                }
            }
            Matrix Temporary(n, m);
            Temporary.array = temp;
            return Temporary;
        }
    }

    Matrix operator*(const Matrix& ma)
    {
        if (m == ma.n)
        {
            vector<vector<int>> temp = vector<vector<int>>(n, vector<int>(ma.m, 0));
            vector<vector<int>> a1 = ma.array;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < ma.m; j++)
                {
                    for (int c = 0; c < m; c++)
                    {
                        temp[i][j] += array[i][c] * a1[c][j];
                    }
                }
            }
            Matrix Temporary(n, ma.m);
            Temporary.array = temp;
            return Temporary;
        } else
        {
            cout << "Error: the dimensional problem occurred" << endl;
            vector<vector<int>> temp = vector<vector<int>>(0, vector<int>(0));
            Matrix Temporary(0, 0);
            Temporary.array = temp;
            return Temporary;
        }
    }

    void operator=(const Matrix& ma)
    {
        array = ma.array;
        n = ma.n;
        m = ma.m;
    }
};

int main()
{
    int n1, m1, n2, m2, n3, m3;
    cin >> n1 >> m1;
    Matrix A(n1, m1);
    A.input();
    cin >> n2 >> m2;
    Matrix B(n2, m2);
    B.input();
    cin >> n3 >> m3;
    Matrix C(n3, m3);
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
