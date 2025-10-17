/*
This program extends the previous matrix system 
to calculate the determinant of a square matrix using 
Gaussian elimination with pivoting by the maximum absolute element. 
After each permutation or elimination step, it prints the updated 
matrix and the step type. Finally, it outputs the determinant formatted to two decimal places.
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

    void output()
    {
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n - 1; j++)
            {
                cout << fixed << setprecision(2) <<array[i][j] << " ";
            }
            cout << fixed << setprecision(2) << array[i][n - 1];
            cout << '\n';
        }
    }

    virtual void makePermutation(int q, int w)
    {
        if (q != w)
        {
            swap(array[q], array[w]);
        }
    }

    // array[i][j] - k * array[i-1][j]
    // array[i][j] - k * array[i-1][j] = 0 => k = array[i][j] / array[i-1][j]
    int makeElimination(int col, int st)
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
                cout << "step #" << st << ": elimination" << endl;
                output();
                st++;
            }
        }
        return st;
    }

    int upperTriangular()
    {
        int column = 0;
        int step = 1;
        int st = 0;
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

                makePermutation(temp, count);
                st += 1;
                cout << "step #" << step << ": permutation" << endl;
                output();
                step ++;
            }
            step = makeElimination(count, step);
        }
        return st;
    }

    void findDeterminant(bool flag)
    {
        double det = 1;
        for (int i = 0; i < n; i++)
        {
            det *= array[i][i];
        }
        cout << "result:" << endl;
        if (det == 0)
        {
            cout << fixed << setprecision(2) << abs(det);
        } else
        {
            if (flag)
            {
                cout << fixed << setprecision(2) << det;
            } else
            {
                cout << fixed << setprecision(2) << det * (-1);
            }
        }
    }
};

int main()
{
    int n;
    cin >> n;
    Matrix A(n); // temporary matrix
    A.input();
    int st = A.upperTriangular();
    if (st % 2 != 0)
    {
        A.findDeterminant(false);
    } else
    {
        A.findDeterminant(true);
    }
    return 0;
}
