#include <iostream>
#include <random>
#include <array>
#include <vector>
#include <thread>
#include <ppl.h>

#include "HighResStopwatch.h"

using namespace std;
using namespace concurrency;

const size_t matSize = 256; //Careful of stack size :D
typedef array<array<double, matSize>, matSize> Matrix;

// Computes the product of two square matrices.
template<class T, size_t arrSize>
array<T, arrSize> matrix_multiply(const array<T, arrSize> &m1, const array<T, arrSize> &m2)
{
	auto result = m1;

	for (size_t i = 0; i < arrSize; i++)
	{
		for (size_t j = 0; j < arrSize; j++)
		{
			double temp = 0;

			for (size_t k = 0; k < arrSize; k++)
			{
				temp += m1[i][k] * m2[k][j];
			}
			result[i][j] = temp;
		}
	}

	return result;
}

// Computes the product of two square matrices in parallel.
template<class T, size_t arrSize>
array<T, arrSize> parallel_matrix_multiply(const array<T, arrSize> &m1, const array<T, arrSize> &m2)
{
	auto result = m1;

	parallel_for(size_t(0), arrSize, [&](size_t i)
	{
		for (size_t j = 0; j < arrSize; j++)
		{
			double temp = 0;

			for (size_t k = 0; k < arrSize; k++)
			{
				temp += m1[i][k] * m2[k][j];
			}
			result[i][j] = temp;
		}
	});

	return result;
}

// Computes the product of two square matrices.
template<class T>
void matrix_multiply(vector<T> &vMList, vector<T> &vResult)
{
	for (auto vm : vMList)
	{
		auto tempMat = vm;
		auto result = vm = {};

		for (size_t i = 0; i < 256; i++)
		{
			for (size_t j = 0; j < 256; j++)
			{
				double temp = 0;

				for (size_t k = 0; k < 256; k++)
				{
					temp += tempMat[i][k] * tempMat[k][j];
				}
				result[i][j] = temp;
			}
		}

		vResult.push_back(result);
	}
}

// Computes the product of two square matrices in parallel.
template<class T, size_t arrSize>
array<T, arrSize> pdarallel_matrix_multiply(const array<T, arrSize> &m1, const array<T, arrSize> &m2)
{
	auto result = m1;

	parallel_for(size_t(0), arrSize, [&](size_t i)
	{
		for (size_t j = 0; j < arrSize; j++)
		{
			double temp = 0;

			for (size_t k = 0; k < arrSize; k++)
			{
				temp += m1[i][k] * m2[k][j];
			}
			result[i][j] = temp;
		}
	});

	return result;
}

void main()
{
	/* Steves Timer */
	HighResStopwatch NormTimer, ParaTimer, listTimer, paraListTimer; 

	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<> dis(0, 1);
	
	Matrix Mat1, Mat2;

	/* Generate random doubles 0..10 */
	for (auto row = Mat1.begin(); row != Mat1.end(); row++)
		std::generate(row->begin(), row->end(), [&]() {	return dis(gen); }); 
	
	for (auto row = Mat2.begin(); row != Mat2.end(); row++)
		std::generate(row->begin(), row->end(), [&]() {	return dis(gen); });
	

	/************************************************************************/
	/* Part 2                                                               */
	/************************************************************************/

	vector<Matrix> MatrixList;
	vector<Matrix> Results;

	for (size_t i = 0; i < 2; i++)
	{
		array<array<double, matSize>, matSize> newMat = {};

		/* Generate random doubles 0..10 */
		for (auto row = newMat.begin(); row != newMat.end(); row++)
			std::generate(row->begin(), row->end(), [&]() {	return dis(gen); });

		MatrixList.push_back(newMat);
	}
	
	cout << "Single core calculation for Matrix list..." << endl;
	this_thread::sleep_for(chrono::milliseconds(200));

	listTimer.StartTimer();

	matrix_multiply(MatrixList, Results);

	listTimer.StopTimer();
	listTimer.Report();

	/************************************************************************/
	/* Matrix Multiplication with single & multi core. 
	   Using timers and threads accuracy and calculation.                   */
	/************************************************************************/

	/**************************************************/

	cout << "Single core calculation..." << endl;
	this_thread::sleep_for(chrono::milliseconds(200));

	NormTimer.StartTimer();

	auto resultMat = matrix_multiply(Mat1, Mat2);
	
	NormTimer.StopTimer();
	NormTimer.Report();

	/**************************************************/

	cout << "Multi core calculation..." << endl;
	this_thread::sleep_for(chrono::milliseconds(200));

	ParaTimer.StartTimer();

	auto resultPMat = parallel_matrix_multiply(Mat1, Mat2);

	ParaTimer.StopTimer();
	ParaTimer.Report();

	/**************************************************/

	cout << "End simulation" << endl;
}