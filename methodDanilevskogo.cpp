#include <iostream>
#include <vector>
#include <iomanip>
#include <stdio.h>

using namespace std;

const int N = 5;
const vector<vector<double>> E =
{
	{1, 0, 0, 0, 0},
	{0, 1, 0, 0, 0},
	{0, 0, 1, 0, 0},
	{0, 0, 0, 1, 0},
	{0, 0, 0, 0, 1}
};

void findSource(vector<vector<double>>& matrix);
void printMatrix(const vector<vector<double>> matrix);
void solve(vector<vector<double>>& matrix);
vector<vector<double>> findInverse(vector<vector<double>> matrix);
void multiplyMatrix(vector<vector<double>> matrix,
	vector<vector<double>> otherMatrix,vector<vector<double>>& result);


int main()
{
	vector<vector<double>> matrix(N);
	findSource(matrix);
	vector<vector<double>> sourceMatrix = matrix;
	printMatrix(matrix);
	cout << endl << endl;
	solve(matrix);
	
	return 0;
}


void multiplyMatrix(vector<vector<double>> matrix,
	vector<vector<double>> otherMatrix, vector<vector<double>>& result)
{
	for (int i = 0; i < N; ++i)
	{
		result[i].clear();
	}
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			result[i].push_back(0.0);
			for (int jj = 0; jj < N; ++jj)
			{
				result[i][j] += matrix[i][jj] * otherMatrix[jj][j];
			}
		}
	}
}


void findSource(vector<vector<double>>& matrix)
{
	vector<vector<double>> A =
	{
		{0.6897, -0.0908, 0.0182, 0.0363, 0.1271},
		{0.0944, 1.0799, 0.0000, -0.0726, 0.0726},
		{0.0545, 0.0000, 0.8776, -0.2541, 0.1452,},
		{-0.1089, 0.2287, 0.0000, 0.8531, -0.0363},
		{0.4538, 0.0000, 0.1634, 0.0182, 1.0164}
	};

	vector<vector<double>> AT(N);
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			AT[i].push_back(A[j][i]);
		}
	}
	multiplyMatrix(AT, A, matrix);
}


void printMatrix(const vector<vector<double>> matrix)
{
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			cout << matrix[i][j] << setw(10) << setprecision(5);
		}
		cout << endl;
	}
}


vector<vector<double>> findInverse(vector<vector<double>> matrix)
{
	vector<vector<double>> inverse(N);
	vector<vector<double>> tmp = matrix;
	for (int ii = 0; ii < N; ++ii)
	{
		matrix = tmp;
		for (int i = 0; i < N; ++i)
		{
			matrix[i].push_back(0);
		}
		matrix[ii][N] = 1;

		////////////////////////////////////////nulify column
		for (int i = 0; i < N; ++i) 
		{
			for (int j = i + 1; j < N; ++j) 
			{
				double scalingFactors = matrix[j][i] / matrix[i][i];
				for (int jj = i; jj < N + 1; ++jj) 
				{
					matrix[j][jj] -= matrix[i][jj] * scalingFactors;
				}
			}
		}
		/////////////////////////////////////////

		/////////////////////////////////////////findX
		vector<double> col(N);
		for (int i = N - 1; i >= 0; --i) 
		{
			double x = matrix[i][N] / matrix[i][i];
			col[i] = x;
			for (int j = i - 1; j >= 0; --j) 
			{
				matrix[j][N] -= x * matrix[j][i];
			}
		}
		//////////////////////////////////////////

		for (int i = 0; i < N; ++i) 
		{
			inverse[i].push_back(col[i]);
		}
	}
	return inverse;
}


void solve(vector<vector<double>>& matrix)
{
	vector<vector<double>> tmp(N);
	for (int ii = 0; ii < N - 1; ++ii)
	{
		vector<vector<double>> M = E;
		int num = N - ii - 1;
		for (int i = 0; i < N; ++i)
		{
			M[num-1][i] = -matrix[num][i] / matrix[num][num-1];
		}
		M[num-1][num-1] = 1 / matrix[num][num-1];
		multiplyMatrix(matrix, M, tmp);
		matrix = tmp;

		vector<vector<double>> inverse = findInverse(M);
		multiplyMatrix(inverse, matrix, tmp);
		matrix = tmp;
	}
	printMatrix(matrix);
	cout << "Polynomial = -h^5";
	for (int i = N - 1, j = 0; i >= 0; --i, ++j)
	{
		if (matrix[0][j] < 0)
		{
			cout << " + ";
		}
		cout << " " << -1 * matrix[0][j] << "h^" << i << " ";
	}
}