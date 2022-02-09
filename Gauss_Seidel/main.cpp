#include <iostream>
#include <fstream> //read and write file
#include <iomanip> //format output to the screen
#include <cmath> //for absolute value functions
/*
The solution is obtained iteratively via
x_(k+1) = D^-1(b - (L+U)x_k)
Do this formula until convergence is reached

*/
using namespace std;
int main()
{
	//---------- Set up ----------//

	//Variables
	int n = 0;
	int iteration = 0;
	double actualValue, checkValue = 0.0;
	double error = 0.0;
	double tolerableError = 0.00001;
	double sumOfb = 0.0; //used to figure out the check value

	//Matrices
	double **A; //LHS
	double **L_U; //L + U
	double **D; //Diagonal

	//Vectors
	double *b; //RHS
	double *newB; //used to calculate actual value
	double *x; //solution

	//File streams
	ifstream fin;
	ofstream fout;

	//Open file containing the matrix
	fin.open("q13.txt");

	//Find dimensions
	fin >> n;


	//Finish constructing matrices
	A = new double *[n]; //LHS
	for (int i = 0; i < n; i++)
		A[i] = new double[n];

    D = new double *[n]; //LHS
	for (int i = 0; i < n; i++)
		D[i] = new double[n];

    L_U = new double *[n];
	for (int i = 0; i < n; i++)
		L_U[i] = new double[n];

	//Fill A, D, L_U
	for(int row = 0; row < n; row++)
	{
		for (int col = 0; col < n; col++)
		{
			fin >> A[row][col];
			if(row == col)
			{
			    if(A[row][col] != 0)
                    D[row][col] = 1/A[row][col]; //D gets initialized inverted
                else
                    D[row][col] = 0.0;

			    L_U[row][col] = 0.0;
			}
			else
			{
			    L_U[row][col] = A[row][col];
			    D[row][col] = 0.0;
			}
		}
	}

	//Finish constructing vectors
	b = new double[n]; //RHS
	newB = new double[n];
	x = new double[n]; //solution vector

	//Fill RHS and solution vector
	for(int row = 0; row < n; row++)
	{
		fin >> b[row];
		x[row] = 0.0; //make first guess here
		newB[row] = 0.0;
	}
	fin.close();

	//Initialize check value
	for(int i = 0; i < n; i++)
    {
        sumOfb += b[i] * b[i];
    }
    checkValue = sqrt(sumOfb);

    /*
    cout << "--- DEBUG ---" << endl << "---   A   ---" << endl;
    for(int i = 0; i < n; i++){ for(int j = 0; j < n; j++){ cout << A[i][j] << " "; }cout << endl; }
    cout << "-------------" << endl;

    cout << "---  L_U  ---" << endl;
    for(int i = 0; i < n; i++){ for(int j = 0; j < n; j++){ cout << L_U[i][j] << " "; }cout << endl; }
    cout << "-------------" << endl;

    cout << "---   D   ---" << endl;
    for(int i = 0; i < n; i++){ for(int j = 0; j < n; j++){ cout << D[i][j] << " "; } cout << endl; }
    cout << "-------------" << endl;

    cout << "---   b   ---" << endl;
    for(int i = 0; i < n; i++){ cout << b[i] << " "; } cout << endl;
    cout << "-------------" << endl;
    */

	//---------- Gauss-Seidel Iteration ----------//
	do
    {
        for(int i = 0; i < n; i++)
        {
            //(L+U)x
            double sumOfThis = 0.0;
            for(int j = 0; j < n; j++)
            {
                sumOfThis += (L_U[i][j] * x[j]);
            }
            double newX = sumOfThis;
            //b-((L+U)x)
            newX = b[i] - newX;
            //D^-1(b-((L+U)x))
            x[i] = D[i][i] * newX;
        }

        //calculate newB vector
        for(int i = 0; i < n; i++)
        {
            double sum = 0.0;
            for(int j = 0; j < n; j++)
            {
                sum += A[i][j] * x[j];
            }
            newB[i] = sum;
        }

        //calculate actual value
        double sumOfNewB = 0.0;
        for(int i = 0; i < n; i++)
        {
            sumOfNewB += newB[i] * newB[i];
        }
        actualValue = sqrt(sumOfNewB);

        //calculate error
        error = fabs(checkValue - actualValue);

        iteration++;

        //-------Print out-------//
        std::cout << "checkValue = " << checkValue << std::endl;
        std::cout << "actual = " << actualValue << std::endl;
        std::cout << "error = " << error << std::endl;
        std::cout << "----------------- iteration = " << iteration << "-----------------" << std::endl;
        std::cout << std::endl;

        //fail safe
        if(iteration > 99)
            goto getout;
    }
    while(error > tolerableError);

    getout:

    int rootN = sqrt(n);
    fout.open("q13_solution.txt");
    fout << x[0] << " ";
    for(int i = 1; i < n; i++)
    {
        if(i%rootN == 0)
            fout << std::endl;
        fout << x[i] << " ";
    }

    return 0;
}
