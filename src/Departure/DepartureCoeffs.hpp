/*
 * DepartureCoeffs.hpp
 *
 *  Created on: 8 Oct 2015
 *      Author: harry
 */

#ifndef DEPARTURECOEFFS_HPP_
#define DEPARTURECOEFFS_HPP_

#include <array>
#include <vector>
#include <cmath>
#include <iostream>

enum class RecombinationCase : unsigned int {CASE_A, CASE_B, UNKNOWN};

class DepartureCoeffs {
public:
	DepartureCoeffs(double bgTemp, double bgEmiss);

	void run(double tem, double den, int NMIN, const std::array<int, 75>& MVAL, int IC, int IR, int nfit, int NPLO, int NPHI);
	std::array<double, 2> bsubnnp1(int nlevel, double tem, double den, int NMIN, const std::array<int, 75>& MVAL, int IC, int IR, int nfit);

	void initColIon(int IONZ, double CTE);
	void initialiseGaunt();

	/**
	 *  Calls appropriate routines for calculation of atomic data for array
	 *  SK.
	 *
	 *  N = initial level, NDASH = final level, IS = subscript which identifies
	 *  value of N in condensed matrix.
	 */
	double BK(int N, int NDASH, int IS);
	double BK2(int N, int NDASH, int IS);

	/**
	 *  Computes free-free absorption coefficient given electron temperature
	 *  T, TL = log(T), T32 = T^1.5, and frequency F in GHz.
	 */
	double CAPPA(double T, double TL, double T32, double F);

	/**
	 *  Calculates collision rates from level N to higher level NDASH at
	 *  electron temperature T. TE12 is sqrt(T). Uses Gauss-Laguerre
	 *  integration of cross-sections (function CROSS) over Maxwell Distribution.
	 *  Function COLGL is used for values outside the region of validity of
	 *  function COLRAT.
	 */
	double COLGL(int N, int NDASH, double T, double TE12);

	/**
	 *  Calculates rate of collisions from level N to higher level NP at
	 *  electron temperature T. TE12 is sqrt(T). Sets rate=0 for (NP-N)>40
	 *  but this is easily modified.
	 *
	 *  This function must be initialised by being called with N=0 before
	 *  any collision rates are computed.
	 *  Theory: Gee, Percival, Lodge and Richards, MNRAS 175, 209-215 (1976).
	 *
	 *  Range of validity of GPLR rates is 10^6/N^2 < T << 3e9. Outside this
	 *  range numerical integration of the GPLR cross-sections is resorted to.
	 *  These cross-sections are valid down to energies of 4/N^2 Rydbergs; the
	 *  cross-section formula can be used at lower energies for BN calculations,
	 *  the innacuracy in the cross-sections having little effect on the BN's.
	 */
	double COLRAT(int N, int NP, double T, double TE12);

	/**
	 *  Correction to rate of population of level N due to radiation field.
	 *  See Radio Recombination Lines (Gordon & Sorochenko). This function may be
	 *  replaced by the user if other than a dilute blackbody radiation field is
	 *  required.
	 *
	 *  This gets called once with ISW = 0 before starting the calculations. Any
	 *  initialisation (including the reading of data if required) should be
	 *  carried out during this first call.
	 *
	 *  Computes terms of
	 *  D(NU)RHO(NU)(N(N+1)B(N+1,N)+N(N-1)B(N-1,N)-N(N)(B(N,N-1)+B(N,N+1)))
	 *  with removed factor
	 *  (C/4PI)(H*H/2PI M KT)**3/2 N*N NE NI EXP(CHI1/N*N KT)
	 */
	double COR(int N, int ISW);

	double C2(double x, double y);

	/**
	 *  Computes cross-section for transition from level N to higher level NP
	 *  due to collision with electron of energy E.
	 *
	 *  The formula is valid for energies in the range 4/N^2 < E << 137^2.
	 *  This function does not check that E is within this range.
	 *  Theory: Gee, Percival, Lodge and Richards, MNRAS 175, 209-215 (1976).
	 */
	double CROSS(int N, int NP, double E);

	double DMJ1(int K, int KM);

	double DMJ2(int K);

	/**
	 *  Lagangian Interpolation.
	 */
	double DPHI(const std::array<int, 8>& IQ, int LG, int ITAU, int N);

	/**
	 *  Lagangian Interpolation.
	 */
	double PHI(const std::array<int, 8>& IQ, int LG, int ITAU, int N);

	double POL (int K, int KM);

	void BCPCH(int KB, int KC, double T, double DENS, double NMIN, int NPLO, int NPHI, int NDIM);

	void COLION(int N, int IONZ, double T, double& QI);

	void HELPME(double& AID, int KM);

	void INTERP(const std::array<int, 75>& M, const std::array<double, 75>& CO, std::array<double, 507>& VAL, std::array<double, 507>& DVAL, int IC, int IR);

	void JMD(std::array<std::array<double, 75>, 75>& SK, std::array<double, 75>& CO, const std::array<int, 75>& MVAL, int IC);

	/**
	 *  Solves simultaneous equations if L=1.
	 *  Inverts matrix if L=0.
	 */
	template <unsigned int NDA>
	void MATINV(std::array<std::array<double, NDA>, NDA>& A, int N, std::array<double, NDA>& B, int L, double& D, int& IRROR, std::array<int, NDA>& IPIV, std::array<std::array<double, 2>, NDA>& IND) {
		int M = std::abs(L);
		D = 1.0;

		for (int i = 1; i <= N; ++i)
			IPIV[i-1] = 0;

		for (int i = 1; i <= N; ++i) {
			double AMAX = 0.0;
			int ICOL = 0, IROW = 0;

			for (int j = 1; j <= N; ++j) {
				if (IPIV[j-1] < 0) {
					IRROR = 1;
					return;
				}
				else if (IPIV[j-1] > 0) {
					continue;
				}
				else {
					for (int k = 1; k <= N; ++k) {
						if (std::abs(A[j-1][k-1]) < 1.0e-38)
							A[j-1][k-1] = 0.0;
						if (IPIV[k-1] - 1 > 0) {
							IRROR = 1;
							return;
						}
						else if (IPIV[k-1] - 1 == 0)
							continue;
						else {
							if (std::abs(A[j-1][k-1]) - AMAX > 0.0) {
								IROW = j;
								ICOL = k;
								AMAX = std::abs(A[j-1][k-1]);
							}
						}
					}
				}
			}

			IPIV[ICOL-1] += 1;

			if (AMAX - 1.0e-38 <= 0) {
				IRROR = 1;
				return;
			}
			if (IROW - ICOL != 0) {
				D = -D;

				for (int j = 1; j <= N; ++j) {
					AMAX = A[IROW-1][j-1];
					A[IROW-1][j-1] = A[ICOL-1][j-1];
					A[ICOL-1][j-1] = AMAX;
				}

				if (M > 0) {
					AMAX = B[IROW-1];
					B[IROW-1] = B[ICOL-1];
					B[ICOL-1] = AMAX;
				}
			}

			IND[i-1][0] = IROW;
			IND[i-1][1] = ICOL;
			AMAX = A[ICOL-1][ICOL-1];
			A[ICOL-1][ICOL-1] = 1.0;

			for (int j = 1; j <= N; ++j)
				A[ICOL-1][j-1] /= AMAX;

			if (M > 0) {
				B[ICOL-1] /= AMAX;
			}

			for (int j = 1; j <= N; ++j) {
				if (j - ICOL == 0)
					continue;
				AMAX = A[j-1][ICOL-1];
				A[j-1][ICOL-1] = 0.0;

				for (int k = 1; k <= N; ++k)
					A[j-1][k-1] -= A[ICOL-1][k-1]*AMAX;

				if (M > 0)
					B[j-1] -= B[ICOL-1]*AMAX;
			}
		}

		if (L > 0) {
			IRROR = 1;
			return;
		}

		for (int i = 1; i <= N; ++i) {
			int j = N + 1 - i;

			if (IND[j-1][0] - IND[j-1][1] == 0)
				continue;

			int IROW = IND[j-1][0];
			int ICOL = IND[j-1][1];

			for (int k = 1; k <= N; ++k) {
				double AMAX = A[k-1][IROW-1];
				A[k-1][IROW-1] = A[k-1][ICOL-1];
				A[k-1][ICOL-1] = AMAX;
			}
		}
	}

	/**
	 *  Computes Einstein radiative rates A and Gaunt factors G for transitions
	 *  from NDASH to N.
	 */
	void RAD(double& A, int N, int NDASH, double& G);


	void RADCOL(double T, const std::array<int, 75>& MVAL, int IC, int NMIN);

	/**
	 *  Computes recombination coefficient, ALPHA, ont level N1 for ions of
	 *  effective charge Z at electron temperature ETEMP.
	 */
	void RECOMB(double Z, double ETEMP, int N1, double& ALPHA);

	/**
	 *  Given a set of integers
	 *  	M[i-1], i = 1,IC  s.t.
	 *  	1) M[i] = M[i-1] + 1 for i < IA
	 *  	where IA >= 1 and
	 *  	2) (M[i] - M[i-1])) > 1 for i >= IA,
	 *  and given a function BK which calculates the elements of a large
	 *  M[IC x IC] matrix, this function uses Lagrange interpolation
	 *  of order 2*(IR+1) to calculate a smaller IC x IC matrix SK.
	 *
	 *  Requires function PHI and IR must be <= (IA - 1).
	 */
	void REDUCE(const std::array<int, 75>& M, int IC, int IR, std::array<std::array<double, 75>, 75>& SK);

	/**
	 *  Computes the right hand side of equations (2.7) of
	 *  Brocklehurst, MNRAS 148, 417 (1970)
	 */
	void RHS(std::array<double, 75>& CO, const std::array<int, 75>& MVAL, int IC);

private:
	// PARMS
	double DENS;
	double T, T1;
	int ITM;

	std::array<double, 506> AL18S4;
	std::array<double, 506> S23TRM;
	double bgT;
	double bgE;
	bool bgRadField = false;

	// TDEP
	double TE32, TE12, CTE;

	// FITDAT
	std::array<std::array<double, 4>, 4> AFIT;
	std::array<int, 4> IVAL;
	int NFIT;

	// HIGHER
	std::array<std::array<double, 4>, 224> STORE1;
	std::array<double, 224> STORE2;
	std::array<double, 24> JMDVAL;
	double B;
	std::array<std::array<double, 5>, 224> STORE3;
	std::array<double, 24> RTVAL;
	int LIMIT;

	// GAUSS
	std::array<double, 12> VALUE;

	// GAUNTS
	std::array<std::array<double, 22>, 50> GAUNT;
	std::array<int, 12> IXV;

	// RCMB
	std::array<double, 99> SV0, SV1, SV2;

	// EXPDAT
	std::array<double, 707> CXP;
	int MAXN;

	// RCRATS
	std::array<double, 75> RADTOT;
	std::array<double, 75> COLTOT;

	// COR
	std::array<double, 508> DX;
	std::array<double, 508> DILT;

	// COLION
	double CONS;
	std::array<double, 507> EXPX;

	bool flag = false;
};




#endif /* DEPARTURECOEFFS_HPP_ */
