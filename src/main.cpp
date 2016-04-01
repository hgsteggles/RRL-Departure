#include <iostream>
#include <fstream>
#include <iomanip>

#include "Departure/DepartureCoeffs.hpp"


int main (int argc, char** argv) {
	DepartureCoeffs coeffs(0, 0);
	std::array<int, 75> MVAL = std::array<int, 75>{5,6,7,8,9,10,11,12,13,14,
		15,17,19,21,23,26,29,32,35,38,41,44,48,52,56,60,64,68,72,77,82,87,
		92,97,102,107,112,118,124,130,136,142,148,154,160,167,174,181,188,
		195,202,210,218,226,234,242,250,259,268,277,286,295,305,315,325,
		335,345,355,366,377,388,399,410,421,432};
	double tem = 4.0e3;
	double den = 1.0e4;
	double nmin = 2;
	int icyc = 1;
	int IC = 75;
	int IR = 2;
	int NFIT = 4;
	int nplo = 1;
	int nphi = 0;

	//coeffs.run(tem, den, nmin, MVAL, IC, IR, NFIT, nplo, nphi);

	// Creating filename.
	std::ofstream file("bsubn.txt", std::ios_base::out);
	if (!file)
		throw std::runtime_error("main: unable to open bsubn.txt");

	// Writing data to file.
	file << std::setprecision(10) << std::scientific;

	double logT_Min = 0;
	double logT_Max = 9;
	double logD_Min = -1;
	double logD_Max = 8.5;
	int nT = 100;
	int nD = 100;

	for (int i = 0; i < nT; ++i) {
		double iT = std::pow(10.0, logT_Min + i*(logT_Max - logT_Min)/double(nT - 1));
		for (int j = 0; j < nD; ++j) {
			double iD = std::pow(10.0, logD_Min + j*(logD_Max - logD_Min)/double(nD - 1));
			std::array<double, 2> b = coeffs.bsubnnp1(66, iT, iD, nmin, MVAL, IC, IR, NFIT);
			file << iT << '\t' << iD << '\t' << b[0] << '\t' << b[1] << '\n';
		}
	}

	file.close();

	return 0;
}
