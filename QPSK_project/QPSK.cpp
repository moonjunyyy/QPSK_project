#include "QPSK.h"

void QPSK::init(int Dim, int Bits, int Amp, double N0)
{
	this->Dim = Dim, this->Bit = Bits,
		this->Amp = Amp, this->N0 = N0;
	ebr = 0.;
	c0t.resize(Dim), c1t.resize(Dim);
	data.resize(Bit), decode.resize(Bit);
	send.clear(), received.clear();
}

void QPSK::build_Bias()
{
	for (int i = 0; i < Dim; i++)
	{
		double f = 1.;
		double t = (double)i / Dim;
		c0t(i) = sqrt(2. / (double)Dim) * cos(2. * PI * f * t);
		c1t(i) = sqrt(2. / (double)Dim) * sin(2. * PI * f * t);
	}
}

void QPSK::build_Data(mt19937_64& rnd)
{
	uniform_int_distribution<int> idis(0, 1);
	for (auto& A : data) A = idis(rnd) * 2 - 1;
}

void QPSK::send_Data()
{
	for (int i = 0; i < Bit / 2; i++)
	{
		auto s_dat = Amp * (c0t * data(i * 2) + c1t * data(i * 2 + 1));
		send.push_back(s_dat);
	}
}

void QPSK::recieve_Channel(double sigma, mt19937_64& rnd)
{
	normal_distribution<double> gaussian(0, sigma);
	for (auto& V : send)
	{
		auto R = V;
		for (auto& B : R)
		{
			B += gaussian(rnd);
		}
		received.push_back(R);
	}
}



void QPSK::decode_Data()
{
	for (int i = 0; i < Bit / 2; i++)
	{
		decode(i * 2) = received[i].dot(c0t) / Amp;
		decode(i * 2 + 1) = received[i].dot(c1t) / Amp;
	}
}

void QPSK::evaluate_Data()
{
	for (int i = 0; i < Bit; i++)
	{
		decode(i) = decode(i) > 0 ? 1 : -1;
		if (decode(i) != data(i)) ebr++;
	}
	ebr = ebr / (double)Bit;
}

void QPSK::execute_Simulation(mt19937_64& rnd)
{
	build_Bias();
	build_Data(rnd);
	send_Data();
	recieve_Channel(sqrt(N0 / 2), rnd);
	decode_Data();
	evaluate_Data();
}

void QPSK::test_Orthonomality(ostream& os)
{
	os
		<< "C0t * C0t = " << c0t.dot(c0t) << endl
		<< "C0t * C1t = " << c0t.dot(c1t) << endl
		<< "C1t * C0t = " << c1t.dot(c0t) << endl
		<< "C1t * C1t = " << c1t.dot(c1t) << endl;

}

void QPSK::print_Send(ostream& os)
{
	int i = 0;
	for (auto& V : send)
	{
		for (auto& D : V)
		{
			os << i++ << ", " << D << endl;
		}
	}
}

void QPSK::print_Receive(ostream& os)
{
	int i = 0;
	for (auto& V : received)
	{
		for (auto& D : V)
		{
			os << i++ << ", " << D << endl;
		}
	}
}

void QPSK::print_Decode(ostream& os)
{
	for (int i = 0; i < Bit / 2; i++)
	{
		os << decode(i * 2) << ", " << decode(i * 2 + 1) << endl;
	}
}

double QPSK::EBR() { return ebr; }

double QPSK::Eb()
{
	/*
	* 	2 / Bit per Signal = 1
	*/
	return Amp * Amp; 
}

void QAM16::init(int Dim, int Bits, int Amp, double N0)
{
	this->Dim = Dim, this->Bit = Bits,
		this->Amp = Amp, this->N0 = N0;
	ebr = 0.;
	c0t.resize(Dim), c1t.resize(Dim);
	data.resize(Bit), decode.resize(Bit);
	send.clear(), received.clear();
}

void QAM16::build_Bias()
{
	for (int i = 0; i < Dim; i++)
	{
		double f = 1.;
		double t = (double)i / Dim;
		c0t(i) = sqrt(2. / (double)Dim) * cos(2. * PI * f * t);
		c1t(i) = sqrt(2. / (double)Dim) * sin(2. * PI * f * t);
	}
}

void QAM16::build_Data(mt19937_64& rnd)
{
	uniform_int_distribution<int> idis(0, 1);
	for (auto& A : data) A = idis(rnd);
}

char QAM16::graycode(char x1, char x2)
{
	x1 = x1 << 1;
	x1 += x2;
	x2 = (x1 >> 1);
	x1 = x2 ^ (x1 % 2);
	x2 = x2 << 1;
	x2 += x1;

	return x2;
}

void QAM16::send_Data()
{
	for (int i = 0; i < Bit / 4; i++)
	{
		data(i * 4) > 0, data(i * 4 + 1) > 0;
		double d1 = graycode(data(i * 4) > 0, data(i * 4 + 1) > 0) * 2. - 3.,
			d2 = graycode(data(i * 4 + 2) > 0, data(i * 4 + 3) > 0) * 2. - 3.;
		auto s_dat = Amp * (d1 * c0t + d2 * c1t);
		send.push_back(s_dat);
	}
}

void QAM16::recieve_Channel(double sigma, mt19937_64& rnd)
{
	normal_distribution<double> gaussian(0, sigma);
	for (auto& V : send)
	{
		auto R = V;
		for (auto& B : R)
		{
			B += gaussian(rnd);
		}
		received.push_back(R);
	}
}

void QAM16::decode_Data()
{
	for (int i = 0; i < Bit / 4; i++)
	{
		double r1 = (received[i].dot(c0t) / Amp);
		double r2 = (received[i].dot(c1t) / Amp);
		
		decode(i * 4) = r1;
		decode(i * 4 + 1) = r1;
		decode(i * 4 + 2) = r2;
		decode(i * 4 + 3) = r2;
	}
}

void QAM16::evaluate_Data()
{
	for (int i = 0; i < Bit / 4; i++)
	{
		double r1 = decode(i * 4);
		double r2 = decode(i * 4 + 2);

		if (r1 > 2.)
		{
			decode(i * 4)	  = 1;
			decode(i * 4 + 1) = 0;
		}
		else if (r1 > 0.)
		{
			decode(i * 4)	  = 1;
			decode(i * 4 + 1) = 1;
		}
		else if (r1 > -2.)
		{
			decode(i * 4)	  = 0;
			decode(i * 4 + 1) = 1;
		}
		else
		{
			decode(i * 4)	  = 0;
			decode(i * 4 + 1) = 0;
		}
		if (r2 > 2.)
		{
			decode(i * 4 + 2) = 1;
			decode(i * 4 + 3) = 0;
		}
		else if (r2 > 0.)
		{
			decode(i * 4 + 2) = 1;
			decode(i * 4 + 3) = 1;
		}
		else if (r2 > -2.)
		{
			decode(i * 4 + 2) = 0;
			decode(i * 4 + 3) = 1;
		}
		else
		{
			decode(i * 4 + 2) = 0;
			decode(i * 4 + 3) = 0;
		}
	}
	for (int i = 0; i < Bit; i++)
	{
		if (decode(i) != data(i)) ebr++;
	}
	ebr /= Bit;
}

void QAM16::execute_Simulation(mt19937_64& rnd)
{
	build_Bias();
	build_Data(rnd);
	send_Data();
	recieve_Channel(sqrt(N0 / 2), rnd);
	decode_Data();
}

void QAM16::test_Orthonomality(ostream& os)
{
	os
		<< "C0t * C0t = " << c0t.dot(c0t) << endl
		<< "C0t * C1t = " << c0t.dot(c1t) << endl
		<< "C1t * C0t = " << c1t.dot(c0t) << endl
		<< "C1t * C1t = " << c1t.dot(c1t) << endl;

}

void QAM16::print_Send(ostream& os)
{
	int i = 0;
	for (auto& V : send)
	{
		for (auto& D : V)
		{
			os << i++ << ", " << D << endl;
		}
	}
}

void QAM16::print_Receive(ostream& os)
{
	int i = 0;
	for (auto& V : received)
	{
		for (auto& D : V)
		{
			os << i++ << ", " << D << endl;
		}
	}
}

void QAM16::print_Decode(ostream& os)
{
	for (int i = 0; i < Bit / 4; i++)
	{
		os << decode(i * 4) << ", " << decode(i * 4 + 2) << endl;
	}
}

double QAM16::EBR() { return ebr; }

double QAM16::Eb()
{
	/*
	*
	*/
	return 240 * Amp * Amp;
}