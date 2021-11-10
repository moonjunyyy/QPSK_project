#pragma once

#include <iostream>
#include <random>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

const double PI = 3.1415926535897932384626433832795028;

class QPSK
{
	int Dim = 0;
	int Bit = 0;
	int Amp = 0;
	double N0 = 0.;
	double ebr = 0.;

	VectorXd c0t, c1t, data, decode;
	vector<VectorXd> send, received;

public :
	QPSK() {}
	QPSK(int Dim, int Bits, int Amp, double N0)
		: Dim(Dim), Bit(Bits), Amp(Amp), N0(N0)
	{
		c0t.resize(Dim), c1t.resize(Dim);
		data.resize(Bit), decode.resize(Bit);
		ebr = 0.;
	}
	void init(int Dim, int Bits, int Amp, double N0);
	void build_Bias();
	void build_Data(mt19937_64& rnd);
	void send_Data();
	void recieve_Channel(double sigma, mt19937_64& rnd);
	void decode_Data();

	void execute_Simulation(mt19937_64& rnd);
	void test_Orthonomality(ostream& os);
	void print_Send(ostream& os);
	void print_Receive(ostream& os);

	double EBR();
	double Eb();
};

class QAM16
{
	int Dim = 0;
	int Bit = 0;
	int Amp = 0;
	double N0 = 0.;
	double ebr = 0.;

	VectorXd c0t, c1t, data, decode;
	vector<VectorXd> send, received;

public:
	QAM16() {}
	QAM16(int Dim, int Bits, int Amp, double N0)
		: Dim(Dim), Bit(Bits), Amp(Amp), N0(N0)
	{
		c0t.resize(Dim), c1t.resize(Dim);
		data.resize(Bit), decode.resize(Bit);
		ebr = 0.;
	}
	void init(int Dim, int Bits, int Amp, double N0);
	void build_Bias();
	void build_Data(mt19937_64& rnd);
	char graycode(char x1, char x2);
	void send_Data();
	void recieve_Channel(double sigma, mt19937_64& rnd);
	void decode_Data();

	void execute_Simulation(mt19937_64& rnd);
	void test_Orthonomality(ostream& os);
	void print_Send(ostream& os);
	void print_Receive(ostream& os);

	double EBR();
	double Eb();
};