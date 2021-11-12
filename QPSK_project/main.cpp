#include <iostream>
#include <fstream>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "Modulation.h"

using namespace std;
using namespace Eigen;

const int N = 20;		// Vector dimention
const int M = 1E6;		// Bits to Transmit
const double Amp = 10;	// Eb = 100;

int main(int argc, char* argv[])
{
	random_device rd;
	mt19937_64 generator(rd());
	
	fstream fio("QPSK20bits.csv", ios::out);
	if (fio.fail()) { cout << "Cannot Open QPSK File!" << endl; return 1; }

	QPSK qpskSim(N, 20, Amp, 10);
	qpskSim.execute_Simulation(generator);
	qpskSim.print_Send(fio);
	fio << endl;
	qpskSim.print_Receive(fio);
	cout << "EBR =\t" << qpskSim.EBR() << endl
		<< "Es / N0 =\t" << 10 * log10(qpskSim.Eb() / 10) << endl << endl;
	fio.close();
	
	qpskSim.init(N, 2000, Amp, 10);
	qpskSim.execute_Simulation(generator);
	fio.open("QPSK_Constellation.csv", ios::out);
	qpskSim.decode_Data();
	qpskSim.print_Decode(fio);
	qpskSim.evaluate_Data();
	fio.close();
	
	fio.open("QPSK_Eb_over_N0.csv", ios::out);
	if (fio.fail()) { cout << "Cannot Open SNR File!" << endl; return 1; }
	for (double dB = 0; dB < 35; dB += 1./3.)
	{
		double N0 = qpskSim.Eb() / pow(10., dB / 10.);

		cout << "for  N0 = " << N0 << "..." << endl;
		qpskSim.init(N, M, Amp, N0);
		qpskSim.execute_Simulation(generator);
		fio << dB << ", " << qpskSim.EBR() << endl;
	}
	fio.close();
	
	fio.open("16QAM20bits.csv", ios::out);
	if (fio.fail()) { cout << "Cannot Open 16QAM File!" << endl; return 1; }

	QAM16 qam16Sim(N, 40, Amp, 0.1);
	qam16Sim.execute_Simulation(generator);
	qam16Sim.print_Send(fio);
	fio << endl;
	qam16Sim.print_Receive(fio);
	cout << "EBR =\t" << qam16Sim.EBR() << endl
		<< "Es / N0 =\t" << 10 * log10(qam16Sim.Eb() / 0.5) << endl << endl;
	fio.close();

	qam16Sim.init(N, 4000, Amp / 4, 1.0);
	qam16Sim.execute_Simulation(generator);
	fio.open("16QAM_Constellation.csv", ios::out);
	qam16Sim.decode_Data();
	qam16Sim.print_Decode(fio);
	qam16Sim.evaluate_Data();
	fio.close();
	
	fio.open("16QAM_Eb_over_N0.csv", ios::out);
	if (fio.fail()) { cout << "Cannot Open SNR File!" << endl; return 1; }
	for (double dB = 0; dB < 35; dB += 1. / 3.)
	{
		double N0 = qam16Sim.Eb() / pow(10., dB / 10.);

		cout << "for  N0 = " << N0 << "..." << endl;
		qam16Sim.init(N, M, Amp, N0);
		qam16Sim.execute_Simulation(generator);
		fio << dB << ", " << qam16Sim.EBR() << endl;
	}
	fio.close();
	

	fio.open("16PSK20bits.csv", ios::out);
	if (fio.fail()) { cout << "Cannot Open 16QAM File!" << endl; return 1; }

	PSK16 psk16Sim(N, 40, Amp, 0.10);
	psk16Sim.execute_Simulation(generator);
	psk16Sim.print_Send(fio);
	fio << endl;
	psk16Sim.print_Receive(fio);
	cout << "EBR =\t" << psk16Sim.EBR() << endl
		<< "Es / N0 =\t" << 10 * log10(psk16Sim.Eb() / 0.1) << endl << endl;
	fio.close();
	 
	psk16Sim.init(N, 4000, Amp / 4, 0.01);
	psk16Sim.execute_Simulation(generator);
	fio.open("16PSK_Constellation.csv", ios::out);
	psk16Sim.decode_Data();
	psk16Sim.print_Decode(fio);
	psk16Sim.evaluate_Data();
	fio.close();
	
	fio.open("16PSK_Eb_over_N0.csv", ios::out);
	if (fio.fail()) { cout << "Cannot Open SNR File!" << endl; return 1; }
	for (double dB = 0; dB < 35; dB += 1./3.)
	{
		double N0 = psk16Sim.Eb() / pow(10., dB / 10.);

		cout << "for  N0 = " << N0 << "...\t";
		psk16Sim.init(N, M, Amp, N0);
		psk16Sim.execute_Simulation(generator);
		cout << psk16Sim.EBR() << endl;
		fio << dB << ", " << psk16Sim.EBR() << endl;
	}
	fio.close();
	
	return 0;
}