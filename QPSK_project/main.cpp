#include <iostream>
#include <cmath>
#include <random>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const int N = 20;		// Vector dimention
const int M = 1000000;		// Bits to Transmit
const double Amp = 10;	// Eb = 100;
const double PI = 3.1415926535897932384626433832795028;

int main(int argc, char* argv[])
{
	random_device rd;
	mt19937_64 generator(rd());
	uniform_int_distribution<int> idis(0, 1);
	normal_distribution<double> gaussian(0, 100000.);

	/*
		Build Biases
	*/
	VectorXd c0t(N), c1t(N);
	for (int i = 0; i < N; i++)
	{
		double f = 1.;
		double t = (double)i / N;
		c0t(i) = sqrt(2. / (double)N) * cos(2. * PI * f * t);
		c1t(i) = sqrt(2. / (double)N) * sin(2. * PI * f * t);
	}
	cout << "c0t * c0t = " << c0t.dot(c0t) << endl;
	cout << "c0t * c1t = " << c0t.dot(c1t) << endl;
	cout << "c1t * c0t = " << c1t.dot(c0t) << endl;
	cout << "c1t * c1t = " << c1t.dot(c1t) << endl << endl;

	/*
		Sender
	*/

	VectorXd data(M);
	VectorXd* Send, * Received;
	Send = new VectorXd[M / 2];
	for (int i = 0; i < M / 2; i++)
	{
		data(i * 2)		= idis(generator) * 2 - 1;
		data(i * 2 + 1) = idis(generator) * 2 - 1;
		Send[i] = Amp * (c0t * data(i * 2) + c1t * data(i * 2 + 1));
	}

	/*
		Channel
	*/

	Received = new VectorXd[M / 2];
	for (int i = 0; i < M / 2; i++)
	{
		Received[i] = Send[i];
		for (auto& A : Received[i])
		{
			A += gaussian(generator);
		}
	}

	/*
		Receiver
	*/
	double EBR = 0.0;

	VectorXd Decode(M);
	for (int i = 0; i < M / 2; i++)
	{
		Decode(i * 2)	  = Received[i].dot(c0t) / Amp > 0 ? 1 : -1;
		Decode(i * 2 + 1) = Received[i].dot(c1t) / Amp > 0 ? 1 : -1;

		if (Decode(i * 2)	  != data(i * 2))	  EBR++;
		if (Decode(i * 2 + 1) != data(i * 2 + 1)) EBR++;
	}
	EBR = EBR / M;
	cout << "EBR = " << EBR << endl << endl;

	delete[] Send, Received;
	return 0;
}