/*
This program extends the matrix hierarchy to include identity,
elimination, and permutation matrices using inheritance.
It supports creating and multiplying these special square matrices 
(e.g., B = E21 * A, C = P21 * A). Identity matrices represent I(n*n),
elimination matrices nullify specific elements, and permutation matrices swap rows.
All matrices use integer elements.
*/

#include <iostream>
#include <vector>

using namespace std;
class Matrix
{
public:
    int n;
    vector<vector<int>> array;
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

    void setArray(vector<vector<int>> t)
    {
        array = std::move(t);
    }

    vector<vector<int>> getArray()
    {
        return array;
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
        this->array = vector<vector<int>>(n, vector<int>(n, 0));
        for (int i = 0; i < n; i++)
        {
            array[i][i]=1;
        }
    }
};

class EliminationMatrix : public Matrix
{
public:
    EliminationMatrix(int n) : Matrix(n)
    {
        this->n = n;
        this->array = vector<vector<int>>(n, vector<int>(n));
    }

    void makeElimination(int q, int w, int el)
    {
        array = vector<vector<int>>(n, vector<int>(n, 0));
        q -= 1;
        w -= 1;
        for (int i = 0; i < n; i++)
        {
            array[i][i]=1;
        }
        array[q][w] = el * (-1);
    }

    int change(int q, int w)
    {
        q -= 1;
        w -= 1;
        if (q < w) // above the main diagonal
        {
            return array[q][w] / array[q][q];
        } else // under the main diagonal
        {
            return array[q][w] / array[w][w];
        }
    }
};

class PermutationMatrix : public Matrix
{
public:
    PermutationMatrix(int n) : Matrix(n)
    {
        this->n = n;
        this->array = vector<vector<int>>(n, vector<int>(n));
    }

    void makePermutation(int q, int w)
    {
        q -= 1;
        w -= 1;
        swap(array[q], array[w]);
    }
};

int main()
{
    int n;
    cin >> n;
    Matrix A(n); // temporary matrix
    A.input();

    IdentityMatrix I(3); // I
    I.output();

    EliminationMatrix E21(n); // E_21
    E21.setArray(A.getArray());
    E21.makeElimination(2, 1, E21.change(2, 1));
    E21.output();

    Matrix B = E21 * A; // B = E_21 * A
    B.output();

    PermutationMatrix P21(n); // P_21
    IdentityMatrix IP(n);
    P21.setArray(IP.getArray());
    P21.makePermutation(2, 1);
    P21.output();

    Matrix C = P21 * A;
    C.output();
    return 0;
}

