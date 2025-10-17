/*
This program models a predator-prey system using the Lotka-Volterra equations: 
(v'(t)=\alpha_1 v(t)-\beta_1 v(t) k(t)) and (k'(t)=-\alpha_2 k(t)+\beta_2 v(t) k(t)). 
It takes initial populations (v_0) and (k_0), coefficients alpha_1, beta_1, alpha_2, beta_2, 
a time limit (T), and the number of approximation points (N). 
Using analytical solutions derived from linear algebra methods,
it computes the populations of victims (v(t_i)) and predators (k(t_i)) 
at each time step (t_i). The output prints arrays of time moments, victim counts, 
and predator counts, all formatted with 2 decimal precision and space-separated.
*/

#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>

using namespace std;
int main()
{
    int numberOfVictims;
    cin >> numberOfVictims;
    int numberOfKillers;
    cin >> numberOfKillers;
    double alpha1, beta1, alpha2, beta2;
    cin >> alpha1 >> beta1 >> alpha2 >> beta2;
    double timeLimit;
    cin >> timeLimit;
    int numberOfThePoints;
    cin >> numberOfThePoints;
    double v0 = numberOfVictims - alpha2 / beta2;
    double k0 = numberOfKillers - alpha1 / beta1;
    double t = 0;
    vector<double> allT = vector<double>(numberOfThePoints+1, 0);
    vector<double> allV = vector<double>(numberOfThePoints+1, 0);
    vector<double> allK = vector<double>(numberOfThePoints+1, 0);
    int i = 0;
    while (t <= timeLimit)
    {
        allT[i] = t;
        double tempV = v0*cos(sqrt(alpha1*alpha2)*t) - k0*((sqrt(alpha2)*beta1)/(beta2*sqrt(alpha1)))*sin(sqrt(alpha1*alpha2)*t);
        double tempK = v0*((sqrt(alpha1)*beta2)/(beta1*sqrt(alpha2)))*sin(sqrt(alpha1*alpha2)*t)+k0*cos(sqrt(alpha1*alpha2)*t);
        allV[i] = tempV + (alpha2/beta2);
        allK[i] = tempK + (alpha1/beta1);
        i++;
        t += (double) timeLimit/numberOfThePoints;
    }
    for (int in = 0; in < 3; in++)
    {
        if (in == 0)
        {
            cout << "t:" << endl;
        } else if (in == 1)
        {
            cout << endl;
            cout << "v:" << endl;
        } else
        {
            cout << endl;
            cout << "k:" << endl;
        }
        for (int j = 0; j <= numberOfThePoints; j++)
        {
            if (in == 0)
            {
                cout << fixed << setprecision(2) << allT[j] << " ";
            } else if (in == 1)
            {
                cout << fixed << setprecision(2) << allV[j] << " ";
            } else
            {
                cout << fixed << setprecision(2) << allK[j] << " ";
            }
        }
    }
    return 0;
}
