/*
 * DepartureCoeffs.cpp
 *
 *  Created on: 8 Oct 2015
 *      Author: harry
 */
#include "DepartureCoeffs.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>

DepartureCoeffs::DepartureCoeffs(double bgTemp, double bgEmiss) {
	T1 = 0.0;

	initialiseGaunt();
	for (int i = 1; i <= 506; ++i) {
		int Z = i;
		AL18S4[i-1] = std::log(18.0*Z)/(4.0*Z);
		S23TRM[i-1] = 0.184 - 0.04*std::pow(Z, -2.0/3.0);
	}

	bgT = bgTemp;
	bgE = bgEmiss;

	if (bgE == 0 || bgT == 0) {
		bgRadField = false;
		std::cout << "No background radiation field." << std::endl << std::endl;
	}

	if (bgRadField) {
		for (int i = 1; i <= 508; ++i)
			DILT[i-1] = 0.5;

		double C15 = 15.778e4/bgT;
		int MAXP = MAXN+1;

		for (int i = 1; i <= MAXP; ++i) {
			int AI = i*(i+1);
			double ARG = C15*(2*i + 1)/double(AI*AI);
			DX[i-1] = std::expm1(ARG);
		}

		if (bgE < 0.9999e10) {
			double bgT32 = bgT*std::sqrt(bgT);
			double TL = std::log10(bgT);

			for (int i = 1; i <= MAXP; ++i) {
				double FREQG = std::pow(6.58e6/i, 3);
				double TOW = bgE*CAPPA(bgT, TL, bgT32, FREQG);

				if (TOW <= 20)
					DILT[i-1] = -0.5*std::expm1(-TOW);

			}
		}
	}
}

void DepartureCoeffs::run(double tem, double den, int NMIN, const std::array<int, 75>& MVAL, int IC, int IR, int nfit, int NPLO, int NPHI) {
	T = tem;
	DENS = den;
	NFIT = nfit;

	IVAL = std::array<int, 4>{75, 72, 69, 66};

	MAXN = MVAL[IC-1];

	int ILIM = 100;

	// Loop starts here by reading values...

	if (T <= 0 || DENS <= 0)
		throw std::runtime_error("DepartureCoeffs::run: Temperature or Density <= 0");

	if (NMIN <= 0)
		NMIN = 2;

	RecombinationCase icase = RecombinationCase::UNKNOWN;
	if (NMIN == 1)
		icase = RecombinationCase::CASE_A;
	else if (NMIN == 2)
		icase = RecombinationCase::CASE_B;
	NPLO = std::max(NPLO, MVAL[0]);
	NPHI = std::min(NPHI, MVAL[IC-1]);
	int ND = NPHI - NPLO + 1;

	// Reduce calculations depending only on T.
	if (T != T1) {
		ITM = 1;

		if (T >= 1000.0)
			ITM = 3;
		TE12 = std::sqrt(T);
		TE32 = T*TE12;
		CTE = 15.778e4/T;

		for (int i = 1; i <= 707; ++i) {
			double CX = 0.0;
			double ARG = CTE/double(i*i);

			if (ARG <= 165.0)
				CX = std::exp(-ARG);

			CXP[i-1] = CX;
		}

		initColIon(1, T);

		RADCOL(T, MVAL, IC, NMIN);
	}

	std::array<std::array<double, 75>, 75> SK;
	REDUCE(MVAL, IC, IR, SK);

	std::array<double, 75> CO;
	RHS(CO, MVAL, IC);

	JMD(SK, CO, MVAL, IC);

	int IRROR = 0;
	double D;
	std::array<int, 75> IPIV;
	std::array<std::array<double, 2>, 75> IND;
	MATINV<75>(SK, IC, CO, 1, D, IRROR, IPIV, IND);

	std::array<double, 507> VAL, DVAL;
	INTERP(MVAL, CO, VAL, DVAL, IC, 2);

	int j = MVAL[0];
	int k = MVAL[IC-1];

	if (ILIM != 0)
		k = ILIM;

	int IPUN = 0;

	std::array<double, 507> KBOUT, KCOUT;

	for (int i = j; i <= k; ++i) {
		double RATIO = DVAL[i-1]/VAL[i-1];
		double H = T*std::pow(i, 3)*RATIO/3.158e5;
		double HH = 1.0 - H;
		double XXXI = std::exp(15.778e4/(i*i*T))*VAL[i-1];
		double X = -HH*XXXI*100.0/T;

		//std::cout << i << " " << VAL[i-1] << " " << HH << " " << DVAL[i-1] << " " << RATIO << " " << H << " " << X << std::endl;

		if (NPHI != 0 && i >= NPLO && i <= NPHI) {
			IPUN += 1;
			KBOUT[IPUN-1] = VAL[i-1]*1.0e4;
			KCOUT[IPUN-1] = std::log10(RATIO)*1.0e4;
		}
	}

	T1 = T;

	if (NPHI != 0) {
		// PRINT.
	}

	// GOTO 10.

}

std::array<double, 2> DepartureCoeffs::bsubnnp1(int nlevel, double tem, double den, int NMIN, const std::array<int, 75>& MVAL, int IC, int IR, int nfit) {
	T = tem;
	DENS = den;
	NFIT = nfit;

	IVAL = std::array<int, 4>{75, 72, 69, 66};

	double T1 = 0.0;
	MAXN = MVAL[IC-1];

	int ILIM = 100;

	// Loop starts here by reading values...

	if (T <= 0 || DENS <= 0)
		throw std::runtime_error("DepartureCoeffs::run: Temperature or Density <= 0");

	if (NMIN <= 0)
		NMIN = 2;

	RecombinationCase icase = RecombinationCase::UNKNOWN;
	if (NMIN == 1)
		icase = RecombinationCase::CASE_A;
	else if (NMIN == 2)
		icase = RecombinationCase::CASE_B;

	// Reduce calculations depending only on T.
	if (T != T1) {
		ITM = 1;

		if (T >= 1000.0)
			ITM = 3;
		TE12 = std::sqrt(T);
		TE32 = T*TE12;
		CTE = 15.778e4/T;

		for (int i = 1; i <= 707; ++i) {
			double CX = 0.0;
			double ARG = CTE/double(i*i);

			if (ARG <= 165.0)
				CX = std::exp(-ARG);

			CXP[i-1] = CX;
		}

		initColIon(1, T);

		RADCOL(T, MVAL, IC, NMIN);
	}

	std::array<std::array<double, 75>, 75> SK;
	REDUCE(MVAL, IC, IR, SK);

	std::array<double, 75> CO;
	RHS(CO, MVAL, IC);

	JMD(SK, CO, MVAL, IC);

	int IRROR = 0;
	double D;
	std::array<int, 75> IPIV;
	std::array<std::array<double, 2>, 75> IND;
	MATINV<75>(SK, IC, CO, 1, D, IRROR, IPIV, IND);

	std::array<double, 507> VAL, DVAL;
	INTERP(MVAL, CO, VAL, DVAL, IC, 2);

	T1 = T;

	return std::array<double, 2>{VAL[nlevel-1], VAL[nlevel]};
}

void DepartureCoeffs::initColIon(int IONZ, double CTE) {
	// Malcolm Walmsley has modified this program to consider collisional
	// ionisation from level 10.
	// The original program considered levels only down to n = 20.
	CONS = double(IONZ*IONZ)*CTE;
	int nlow = 9;
	int nlow1 = nlow + 1;
	for (int i = nlow1; i <= 507; ++i) {
		EXPX[i-1] = std::exp(-CONS/double(i*i));
	}
}

void DepartureCoeffs::initialiseGaunt() {
	std::array<double, 50> A1 = std::array<double, 50>{.7166,.7652,.7799,.7864,.7898,.7918,.7931,.7940,.7946,.7951,
	.7954,.7957,.7959,.7961,.7963,.7964,.7965,.7966,.7966,.7967,.7968,
	.7968,.7968,.7969,.7969,.7969,.7970,.7970,.7970,.7970,.7970,.7971,
	.7971,.7971,.7971,.7971,.7971,.7971,.7971,.7971,.7972,.7972,.7972,
	.7972,.7972,.7972,.7972,.7972,.7972,.7972};

	std::array<double, 50> A2 = std::array<double, 50>{.7566,.8217,.8441,.8549,.8609,.8647,.8672,.8690,.8702,.8712,
	.8720,.8726,.8730,.8734,.8737,.8740,.8742,.8744,.8746,.8747,.8749,
	.8750,.8751,.8751,.8752,.8753,.8753,.8754,.8755,.8755,.8755,.8756,
	.8756,.8756,.8757,.8757,.8757,.8757,.8758,.8758,.8758,.8758,.8758,
	.8759,.8759,.8759,.8759,.8759,.8759,.8759};

	std::array<double, 50> A3 = std::array<double, 50>{.7674,.8391,.8653,.8784,.8861,.8910,.8944,.8968,.8986,.9000,
	.9011,.9019,.9026,.9032,.9037,.9041,.9044,.9047,.9049,.9052,.9054,
	.9055,.9057,.9058,.9059,.9060,.9061,.9062,.9063,.9064,.9064,.9065,
	.9066,.9066,.9067,.9067,.9067,.9068,.9068,.9068,.9069,.9069,.9069,
	.9070,.9070,.9070,.9070,.9070,.9071,.9071};

	std::array<double, 50> A4 = std::array<double, 50>{.7718,.8467,.8750,.8896,.8984,.9041,.9081,.9110,.9132,.9149,
	.9163,.9173,.9182,.9190,.9196,.9201,.9205,.9209,.9213,.9215,.9218,
	.9220,.9222,.9224,.9226,.9227,.9228,.9230,.9231,.9232,.9233,.9233,
	.9234,.9235,.9235,.9236,.9237,.9237,.9238,.9238,.9238,.9239,.9239,
	.9240,.9240,.9240,.9240,.9241,.9241,.9241};

	std::array<double, 50> A5 = std::array<double, 50>{.7741,.8507,.8804,.8960,.9055,.9118,.9162,.9195,.9220,.9240,
	.9255,.9268,.9278,.9287,.9294,.9300,.9306,.9310,.9314,.9318,.9321,
	.9324,.9326,.9329,.9331,.9332,.9334,.9335,.9337,.9338,.9339,.9340,
	.9341,.9342,.9343,.9344,.9344,.9345,.9345,.9346,.9347,.9347,.9348,
	.9348,.9348,.9349,.9349,.9349,.9350,.9350};

	std::array<double, 50> A6 = std::array<double, 50>{.7753,.8531,.8837,.8999,.9099,.9167,.9215,.9251,.9278,.9300,
	.9317,.9331,.9343,.9352,.9361,.9368,.9374,.9379,.9384,.9388,.9392,
	.9395,.9398,.9400,.9403,.9405,.9407,.9408,.9410,.9412,.9413,.9414,
	.9415,.9416,.9417,.9418,.9419,.9420,.9420,.9421,.9422,.9422,.9423,
	.9423,.9424,.9424,.9425,.9425,.9426,.9426};

	std::array<double, 50> A7 = std::array<double, 50>{.7761,.8547,.8858,.9025,.9130,.9200,.9251,.9289,.9318,.9342,
	.9360,.9376,.9389,.9399,.9408,.9416,.9423,.9429,.9434,.9439,.9443,
	.9447,.9450,.9453,.9455,.9458,.9460,.9462,.9464,.9466,.9467,.9468,
	.9470,.9471,.9472,.9473,.9474,.9475,.9476,.9477,.9477,.9478,.9479,
	.9479,.9480,.9480,.9481,.9481,.9482,.9482};

	std::array<double, 50> A8 = std::array<double, 50>{.7767,.8558,.8873,.9044,.9151,.9224,.9277,.9317,.9348,.9373,
	.9393,.9409,.9423,.9434,.9444,.9453,.9460,.9467,.9472,.9477,.9482,
	.9486,.9489,.9493,.9496,.9498,.9501,.9503,.9505,.9507,.9509,.9510,
	.9512,.9513,.9514,.9515,.9517,.9518,.9519,.9519,.9520,.9521,.9522,
	.9522,.9523,.9524,.9524,.9525,.9525,.9526};

	std::array<double, 50> A9 = std::array<double, 50>{.7771,.8565,.8884,.9058,.9167,.9242,.9297,.9338,.9370,.9396,
	.9417,.9434,.9449,.9461,.9472,.9481,.9489,.9496,.9502,.9507,.9512,
	.9517,.9520,.9524,.9527,.9530,.9533,.9535,.9537,.9539,.9541,.9543,
	.9545,.9546,.9548,.9549,.9550,.9551,.9552,.9553,.9554,.9555,.9556,
	.9557,.9557,.9558,.9559,.9559,.9560,.9561};

	std::array<double, 50> A10 = std::array<double, 50>{.7773,.8571,.8892,.9068,.9179,.9256,.9312,.9354,.9388,.9414,
	.9436,.9454,.9470,.9482,.9494,.9503,.9512,.9519,.9526,.9531,.9537,
	.9541,.9545,.9549,.9553,.9556,.9559,.9561,.9564,.9566,.9568,.9570,
	.9572,.9573,.9575,.9576,.9577,.9579,.9580,.9581,.9582,.9583,.9584,
	.9585,.9585,.9586,.9587,.9588,.9588,.9589};

	std::array<double, 50> A11 = std::array<double, 50>{.7775,.8575,.8898,.9076,.9188,.9267,.9324,.9367,.9402,.9429,
	.9452,.9470,.9486,.9500,.9511,.9521,.9530,.9538,.9545,.9551,.9556,
	.9561,.9566,.9570,.9573,.9577,.9580,.9583,.9585,.9588,.9590,.9592,
	.9594,.9595,.9597,.9599,.9600,.9601,.9603,.9604,.9605,.9606,.9607,
	.9608,.9609,.9610,.9610,.9611,.9612,.9612};

	std::array<double, 50> A14 = std::array<double, 50>{.7779,.8583,.8910,.9091,.9207,.9288,.9347,.9393,.9429,.9459,
	.9483,.9503,.9520,.9535,.9548,.9559,.9569,.9578,.9585,.9592,.9598,
	.9604,.9609,.9614,.9618,.9622,.9625,.9629,.9632,.9634,.9637,.9639,
	.9642,.9644,.9646,.9647,.9649,.9651,.9652,.9654,.9655,.9656,.9657,
	.9659,.9660,.9661,.9662,.9662,.9663,.9664};

	std::array<double, 50> A17 = std::array<double, 50>{.7781,.8587,.8916,.9099,.9217,.9300,.9361,.9408,.9446,.9476,
	.9502,.9523,.9541,.9557,.9571,.9582,.9593,.9602,.9611,.9618,.9625,
	.9631,.9637,.9642,.9646,.9651,.9655,.9658,.9662,.9665,.9668,.9670,
	.9673,.9675,.9677,.9679,.9681,.9683,.9685,.9686,.9688,.9689,.9691,
	.9692,.9693,.9694,.9695,.9696,.9697,.9698};

	std::array<double, 50> A20 = std::array<double, 50>{.7782,.8590,.8920,.9105,.9224,.9307,.9370,.9418,.9456,.9488,
	.9514,.9536,.9555,.9571,.9586,.9598,.9609,.9619,.9628,.9636,.9643,
	.9649,.9655,.9661,.9666,.9670,.9675,.9679,.9682,.9686,.9689,.9692,
	.9694,.9697,.9699,.9702,.9704,.9706,.9708,.9709,.9711,.9713,.9714,
	.9716,.9717,.9718,.9719,.9721,.9722,.9723};

	std::array<double, 50> A25 = std::array<double, 50>{.7784,.8592,.8924,.9110,.9230,.9315,.9378,.9428,.9467,.9500,
	.9527,.9550,.9569,.9586,.9601,.9614,.9626,.9637,.9646,.9655,.9662,
	.9669,.9676,.9682,.9687,.9692,.9697,.9701,.9705,.9709,.9712,.9715,
	.9718,.9721,.9724,.9727,.9729,.9731,.9733,.9735,.9737,.9739,.9741,
	.9742,.9744,.9745,.9747,.9748,.9749,.9751};

	std::array<double, 50> A30 = std::array<double, 50>{.7784,.8594,.8926,.9113,.9234,.9319,.9383,.9433,.9474,.9507,
	.9534,.9558,.9578,.9595,.9611,.9625,.9637,.9648,.9657,.9666,.9674,
	.9682,.9689,.9695,.9701,.9706,.9711,.9715,.9720,.9724,.9727,.9731,
	.9734,.9737,.9740,.9743,.9745,.9748,.9750,.9752,.9754,.9756,.9758,
	.9760,.9762,.9763,.9765,.9766,.9768,.9769};

	std::array<double, 50> A40 = std::array<double, 50>{.7785,.8595,.8928,.9116,.9237,.9324,.9389,.9440,.9480,.9514,
	.9542,.9567,.9587,.9606,.9622,.9636,.9649,.9660,.9670,.9680,.9688,
	.9696,.9703,.9710,.9716,.9722,.9727,.9732,.9737,.9741,.9745,.9749,
	.9752,.9756,.9759,.9762,.9765,.9768,.9770,.9773,.9775,.9777,.9779,
	.9781,.9783,.9785,.9787,.9788,.9790,.9791};

	std::array<double, 50> A50 = std::array<double, 50>{.7785,.8596,.8929,.9117,.9239,.9326,.9391,.9443,.9484,.9518,
	.9547,.9571,.9592,.9611,.9627,.9642,.9655,.9666,.9677,.9687,.9696,
	.9704,.9711,.9718,.9724,.9730,.9736,.9741,.9746,.9751,.9755,.9759,
	.9763,.9766,.9770,.9773,.9776,.9779,.9781,.9784,.9786,.9789,.9791,
	.9793,.9795,.9797,.9799,.9801,.9803,.9804};

	std::array<double, 50> A100 = std::array<double, 50>{.7785,.8597,.8931,.9119,.9242,.9329,.9395,.9447,.9489,.9523,
	.9553,.9578,.9599,.9619,.9635,.9651,.9664,.9676,.9687,.9698,.9707,
	.9716,.9724,.9731,.9738,.9744,.9750,.9756,.9761,.9766,.9771,.9775,
	.9779,.9783,.9787,.9791,.9794,.9797,.9801,.9803,.9806,.9809,.9812,
	.9814,.9817,.9819,.9821,.9823,.9825,.9827};

	std::array<double, 50> A150 = std::array<double, 50>{.7786,.8597,.8931,.9119,.9242,.9330,.9396,.9448,.9490,.9525,
	.9554,.9579,.9601,.9620,.9637,.9652,.9666,.9678,.9690,.9700,.9710,
	.9718,.9726,.9734,.9741,.9747,.9754,.9759,.9765,.9770,.9775,.9779,
	.9783,.9787,.9791,.9795,.9799,.9802,.9805,.9808,.9811,.9814,.9817,
	.9819,.9822,.9824,.9827,.9829,.9831,.9833};

	std::array<double, 50> A225 = std::array<double, 50>{.7786,.8597,.8931,.9120,.9243,.9330,.9396,.9448,.9490,.9525,
	.9554,.9580,.9602,.9621,.9638,.9653,.9667,.9680,.9691,.9701,.9711,
	.9720,.9728,.9735,.9742,.9749,.9755,.9761,.9766,.9772,.9776,.9781,
	.9785,.9790,.9793,.9797,.9801,.9804,.9808,.9811,.9814,.9817,.9819,
	.9822,.9825,.9827,.9829,.9832,.9834,.9836};

	std::array<double, 50> A500 = std::array<double, 50>{.7786,.8597,.8931,.9120,.9243,.9330,.9396,.9448,.9490,.9525,
	.9555,.9580,.9602,.9621,.9639,.9654,.9668,.9680,.9692,.9702,.9712,
	.9720,.9729,.9736,.9743,.9750,.9756,.9762,.9768,.9773,.9778,.9782,
	.9787,.9791,.9795,.9799,.9802,.9806,.9809,.9812,.9815,.9818,.9821,
	.9824,.9827,.9829,.9831,.9834,.9836,.9838};

	for (int i = 0; i < 50; ++i) {
		GAUNT[i][0] = A1[i];
		GAUNT[i][1] = A2[i];
		GAUNT[i][2] = A3[i];
		GAUNT[i][3] = A4[i];
		GAUNT[i][4] = A5[i];
		GAUNT[i][5] = A6[i];
		GAUNT[i][6] = A7[i];
		GAUNT[i][7] = A8[i];
		GAUNT[i][8] = A9[i];
		GAUNT[i][9] = A10[i];
		GAUNT[i][10] = A11[i];
		GAUNT[i][11] = A14[i];
		GAUNT[i][12] = A17[i];
		GAUNT[i][13] = A20[i];
		GAUNT[i][14] = A25[i];
		GAUNT[i][15] = A30[i];
		GAUNT[i][16] = A40[i];
		GAUNT[i][17] = A50[i];
		GAUNT[i][18] = A100[i];
		GAUNT[i][19] = A150[i];
		GAUNT[i][20] = A225[i];
		GAUNT[i][21] = A500[i];
	}

	IXV = std::array<int, 12>{11,14,17,20,25,30,40,50,100,150,225,507};

	VALUE = std::array<double, 12>{.4975936099985107,.4873642779856548,.4691372760013664,
	.4432077635022005,.4100009929869515,.3700620957892772,
	.3240468259684878,.2727107356944198,.2168967538130226,
	.1575213398480817,.0955594337368082,.0320284464313028};

	std::array<double, 33> SV0A = std::array<double, 33>{.06845,.07335,
		.07808,.08268,.08714,.09148,.09570,.09982,.10385,.10778,.11163,
		.1209,.1297,.1382,.1462,.1540,.1615,.1687,.1757,.1824,.1889,
		.1953,.2015,.2133,.2245,.2352,.2454,.2552,.2646,.2736,.2823,
		.2906,.2987};

	std::array<double, 33> SV0B = std::array<double, 33>{.3140,.3284,.3419,
		.3547,.3668,.3783,.3892,.3996,.4096,.4191,.4413,.4615,.4799,.4968,
		.5124,.5269,.5404,.5530,.5648,.5759,.5864,.5963,.6146,.6311,.6461,
		.6598,.6724,.6840,.6947,.7047,.7140,.7226,.7384};

	std::array<double, 33> SV0C = std::array<double, 33>{.7524,.7649,.7761,
		.7862,.7955,.8039,.8117,.8188,.8254,.8399,.8521,.8626,.8716,.8795,
		.8865,.8927,.8982,.9032,.9077,.9118,.9156,.9223,.9279,.9328,.9370,
		.9408,.9441,.9471,.9498,.9522,.9544,.9544,.9544};

	std::array<double, 33> SV1A = std::array<double, 33>{.00417,.00444,
		.00469,.00493,.00516,.00538,.00558,.00578,.00597,.00614,.00631,
		.0067,.0070,.0073,.0076,.0078,.0080,.0082,.0083,.0085,.0086,
		.0087,.0087,.0088,.0088,.0088,.0087,.0086,.0084,.0082,.0080,
		.0077,.0074};

	std::array<double, 33> SV1B = std::array<double, 33>{.0068,.0060,.0053,
		.0044,.0035,.0025,.0016,+.0005,-.0005,-.0016,-.0044,-.0072,-.0101,
		-.0130,-.0160,-.0190,-.0220,-.0250,-.0279,-.0308,-.0337,-.0366,
		-.0422,-.0478,-.0532,-.0586,-.0638,-.0689,-.0739,-.0788,-.0835,
		-.0882,-.0972};

	std::array<double, 33> SV1C = std::array<double, 33>{-.1058,-.1141,
		-.1220,-.1295,-.1368,-.1438,-.1506,-.1571,-.1634,-.1784,-.1923,
		-.2053,-.2174,-.2289,-.2397,-.2499,-.2596,-.2689,-.2777,-.2862,
		-.2944,-.3099,-.3242,-.3376,-.3502,-.3622,-.3736,-.3844,-.3948,
		-.4047,-.4147,-.4147,-.4147};

	std::array<double, 33> SV2A = std::array<double, 33>{-.00120,-.00131,
		-.00142,-.00153,-.00164,-.00175,-.00186,-.00197,-.00207,-.00218,
		-.00228,-.0025,-.0028,-.0030,-.0033,-.0035,-.0038,-.0040,-.0043,
		-.0045,-.0047,-.0050,-.0052,-.0056,-.0061,-.0065,-.0070,-.0074,
		-.0078,-.0082,-.0086,-.0091,.0095};

	std::array<double, 33> SV2B = std::array<double, 33>{-.0103,-.0110,
		-.0118,-.0126,-.0133,-.0141,-.0148,-.0155,-.0162,-.0170,-.0187,
		-.0204,-.0221,-.0237,-.0253,-.0269,-.0284,-.0299,-.0314,-.0329,
		-.0344,-.0359,-.0388,-.0416,-.0444,-.0471,-.0497,-.0523,-.0549,
		-.0575,-.0600,-.0625,-.0674};

	std::array<double, 33> SV2C = std::array<double, 33>{-.0722,-.0768,
		-.0814,-.0859,-.0904,-.0947,-.0989,-.1031,-.1072,-.1174,-.1272,
		-.1367,-.146,-.155,-.1638,-.1723,-.1807,-.189,-.197,-.2049,
		-.2127,-.228,-.243,-.257,-.272,-.285,-.299,-.312,-.325,-.337,
		-.35,-.35,-.35};

	for (int i = 0; i < 33; ++i) {
		SV0[i] = SV0A[i];
		SV0[i+33] = SV0B[i];
		SV0[i+66] = SV0C[i];
		SV1[i] = SV1A[i];
		SV1[i+33] = SV1B[i];
		SV1[i+66] = SV1C[i];
		SV2[i] = SV2A[i];
		SV2[i+33] = SV2B[i];
		SV2[i+66] = SV2C[i];
	}
}

// Depends on: COLRAT, COLION, COR
double DepartureCoeffs::BK (int N, int NDASH, int IS) {
	if (N - NDASH < 0) {
		double A, G;
		RAD(A, N, NDASH, G);
		double C = COLRAT(N, NDASH, T, TE12);
		int AN = N;
		int ANDASH = NDASH;

		double TEMP;
		if (NDASH > 707 || CXP[NDASH-1] < 1.0e-30 || CXP[N-1] < 1.0e-30) {
			TEMP = std::exp(-CTE*(1.0/double(AN*AN) - 1.0/double(ANDASH*ANDASH)));
		}
		else
			TEMP = CXP[N-1]/CXP[NDASH-1];
		TEMP *= std::pow(ANDASH/double(AN), 2);

		double bk = (A + DENS*C)*TEMP;

		if (N > 20 && NDASH == N + 1)
			bk += COR(N, 1);

		return bk;
	}
	else if (N - NDASH == 0) {
		double RT;
		COLION(N, 1, T, RT);

		double bk = -RADTOT[IS-1] - (COLTOT[IS-1] + RT)*DENS;

		if (N > 20)
			bk += COR(N, 3);

		return bk;
	}
	else {
		double C = COLRAT(NDASH, N, T, TE12);
		double bk = C*DENS;

		if (NDASH == N - 1 && N > 20)
			bk += COR(N, 2);

		return bk;
	}
}

double DepartureCoeffs::BK2 (int N, int NDASH, int IS) {
	if (N - NDASH < 0) {
		double A, G;
		RAD(A, N, NDASH, G);
		double C = COLRAT(N, NDASH, T, TE12);
		int AN = N;
		int ANDASH = NDASH;

		double TEMP;
		if (NDASH > 707 || CXP[NDASH-1] < 1.0e-30 || CXP[N-1] < 1.0e-30) {
			TEMP = std::exp(-CTE*(1.0/double(AN*AN) - 1.0/double(ANDASH*ANDASH)));
		}
		else
			TEMP = CXP[N-1]/CXP[NDASH-1];
		TEMP *= std::pow(ANDASH/double(AN), 2);

		double bk = (A + DENS*C)*TEMP;

		if (N > 20 && NDASH == N + 1)
			bk += COR(N, 1);

		return bk;
	}
	else if (N - NDASH == 0) {
		double RT;
		COLION(N, 1, T, RT);

		double bk = -RADTOT[IS-1] - (COLTOT[IS-1] + RT)*DENS;

		if (N > 20)
			bk += COR(N, 3);

		return bk;
	}
	else {
		double C = COLRAT(NDASH, N, T, TE12);
		double bk = C*DENS;

		if (NDASH == N - 1 && N > 20)
			bk += COR(N, 2);

		return bk;
	}
}

double DepartureCoeffs::CAPPA(double T, double TL, double T32, double F) {
	double V = 0.6529 + 0.6666667*std::log10(F) - TL;
	double ALIV;

	if (V > -2.6) {
		if (V > -0.25)
			ALIV = -1.0842*V + 0.1359;
		else
			ALIV = -1.2326*V + 0.0987;
	}
	else {
		ALIV = -1.1249*V + 0.3788;
	}

	return 4.646*std::expm1(.047993*F/T)*std::exp(2.302585*ALIV)/(std::pow(F, 2.33333)*T32);
}

double DepartureCoeffs::COLGL(int N, int NDASH, double T, double TE12) {
	std::array<double, 10> XGL = {0.4157746, 2.2942800, 6.289945000, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	std::array<double, 10> WGL = {0.7110930, 0.2785177, 1.038926e-2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	int NGL = 3;
	double BETA = 1.58e5/T;
	int EN2 = N*N;
	int END2 = NDASH*NDASH;
	double DE = 1.0/EN2 - 1.0/END2;

	double colGL = 0;

	for (int i = 1; i <= NGL; ++i) {
		double E = XGL[i-1]/BETA + DE;
		colGL += WGL[i-1]*CROSS(N, NDASH, E)*E*BETA;
	}

	colGL *= 6.21241e5*TE12*EN2/END2;

	return colGL;
}

double DepartureCoeffs::COLRAT(int N, int NP, double T, double TE12) {
	double IS = NP - N;
	if (IS > 40)
		return 0;
	int S = IS;
	double EN2 = N*N;
	if (T < 1e6/EN2) {
		// GPLR collision rate formulas invalid at low temperatures.
		// Integrate cross-sections.
		return COLGL(N, NP, T, TE12);
	}
	int EN = N;
	int ENP = NP;
	int IPOW = 1+IS+IS;
	int POW = IPOW;
	int ENNP = N*NP;
	double BETA = 1.58e5/T;
	double BETA1 = 1.4*std::sqrt(ENNP);
	double BETRT = BETA1/BETA;
	double BETSUM = BETA1 + BETA;
	double F1 = 0.2*S/ENNP;
	F1 = F1 > 0.02 ? std::pow(1.0 - F1, IPOW) : 1.0 - POW*F1;

	double A = (2.666667/S)*std::pow(ENP/double(S*EN), 3)*S23TRM[IS-1]*F1;
	double L = 0.85/BETA;
	L = std::log((1.0 + 0.53*L*L*ENNP)/(1.0 + 0.4*L));
	double J1 = 1.333333*A*L*BETRT/BETSUM;
	double DRT = std::sqrt(2.0 - double(N*N)/double(NP*NP));
	F1 = 0.3*S/ENNP;
	F1 = F1 > 0.02 ? std::pow(1.0-F1, IPOW) : 1.0 - POW*F1;

	double J2 = 0;
	double ARG = BETA/BETA1;
	if (ARG <= 150) {
		J2 = 1.777778*F1*std::pow(ENP*(DRT+1.0)/double((EN+ENP)*S), 3)*std::exp(-ARG)/(BETA/(1.0 - AL18S4[IS-1]));
	}

	double XI = 2.0/(EN2*(DRT - 1.0));
	double Z = 0.75*XI*BETSUM;

	double EXPZ = 0;
	if (Z <= 150)
		EXPZ = std::exp(-Z);

	double J4 = 2.0/(Z*(2.0 + Z*(1.0 + EXPZ)));
	double J3 = 0.25*std::pow(EN2*XI/ENP, 3)*J4*std::log(1.+0.5*BETA*XI)/BETSUM;

	return EN2*EN2*(J1+J2+J3)/std::pow(TE12, 3);
}

double DepartureCoeffs::COR(int N, int ISW) {
	// TODO: Some setting up is done here.

	if (bgRadField) {
		if (N > MAXN) {
			throw std::runtime_error("COR called with N too large. N = " + std::to_string(N));
		}
		else if (ISW == 3) {
			double AM, AP, G, EX;
			RAD(AM, N-1, N, G);
			RAD(AP, N, N+1, EX);
			EX = DX[N-2];
			G = DX[N-1];

			return DILT[N-1]*(AM/EX+((N+1)*(N+1)/double(N*N))*AP/G);
		}
		else if (ISW == 2) {
			double A, G;
			RAD(A, N-1, N, G);
			double EX = DX[N-2];

			return DILT[N-1]*A/(EX/(EX+1.0));
		}
		else {
			double A, G;
			RAD(A, N, N+1, G);
			double EX = DX[N-1];

			return ((N+1)*(N+1)/double(N*N))*DILT[N-1]*A/(EX*(EX+1.0));
		}
	}
	else {
		return 0;
	}
}

double DepartureCoeffs::C2(double x, double y) {
	return x*x*std::log(1.0 + 2.0*x/3.0)/(2*y + 1.5*x);
}

double DepartureCoeffs::CROSS(int N, int NP, double E) {
	int EN = N;
	int ENP = NP;
	int IS = NP-N;
	int IPOW = 1 + 2*IS;
	int POW = IPOW;
	int S = IS;
	int ENNP = N*NP;
	double EENNP = E*E*ENNP;
	int EN2 = N*N;
	double ENN = E*EN2;
	double D = 0.2*S/ENNP;
	D = D > 0.02 ? std::pow(1.0 - D, IPOW) : 1.0 - POW*D;

	double A = (2.666667/S)*std::pow(ENP/(S*EN), 3)*S23TRM[IS-1]*D;
	D = 0;

	double ARG = 1.0/EENNP;
	if (ARG < 150)
		D = std::exp(-ARG);

	double L = std::log((1.0 + 0.53*EENNP)/(1.0 + 0.4*E));
	double F = std::pow(1.0 - 0.3*S*D/ENNP, IPOW);
	double G = 0.5*std::pow(ENN/ENP, 3);
	double Y = 1./(1.0 - D*AL18S4[IS-1]);
	double DRT = std::sqrt(2.0 - N*N/double(NP*NP));
	double XP = 2.0/(ENN*(DRT + 1.0));
	double XM = 2.0/(ENN*(DRT - 1.0));
	double H = C2(XM, Y) - C2(XP,Y);

	return 8.797016e-17*(EN2*EN2/E)*(A*D*L + F*G*H);
}

double DepartureCoeffs::DMJ1(int K, int KM) {
	double value = 0;
	for (int i = 1; i <= NFIT; ++i)
		value += AFIT[i-1][KM-1]*STORE1[K-1][i-1];
	return value;
}

double DepartureCoeffs::DMJ2(int K) {
	double tot = 0;
	for (int i = 1; i <= NFIT; ++i) {
		double sum = 0;
		for (int j = 1; j <= NFIT; ++j) {
			sum += AFIT[i-1][j-1];
		}
		tot += sum*STORE1[K-1][i-1];
	}

	return 1.0 - tot;
}

double DepartureCoeffs::DPHI(const std::array<int, 8>& IQ, int LG, int ITAU, int N) {
	double value = 0;

	for (int i = 1; i <= LG; ++i) {
		if (i == ITAU)
			continue;
		double A = 1.0;
		for (int j = 1; j <= LG; ++j) {
			if (j == ITAU || j == i)
				continue;
			A *= N - IQ[j-1];
		}
		value += A;
	}

	return value;
}

/**
 *  Lagangian Interpolation.
 */
double DepartureCoeffs::PHI(const std::array<int, 8>& IQ, int LG, int ITAU, int N) {
	double value = 1.0;

	for (int i = 1; i <= LG; ++i) {
		if (i == ITAU)
			continue;
		value *= N - IQ[i-1];
	}

	return value;
}

double DepartureCoeffs::POL(int K, int KM) {
	return RTVAL[K-1]*STORE3[K+LIMIT-1][KM-1]*STORE2[K+LIMIT-1];
}

void DepartureCoeffs::BCPCH(int KB, int KC, double T, double DENS, double NMIN, int NPLO, int NPHI, int NDIM) {
	// TODO: printing...
}

void DepartureCoeffs::COLION (int N, int IONZ, double T, double& QI) {
	QI = 0;

	if (N == 0)
		throw std::runtime_error("COLION: N cannot equal 0!");

	if (N <= 507)
		return;

	double X = CONS/(N*N);
	double DXP = EXPX[N-1];

	if (N > 507)
		DXP = std::exp(-X);

	double E1;
	if (X <= 1)
		E1 = -0.57721566 + X*(0.9999193 + X*(-0.24991055 + X*(0.5519968e-1 + X*(-0.976004e-2 + X*0.107857e-2)))) - std::log(X);
	else
		E1 = ((0.250621 + X*(2.334733 + X))/(1.681534 + X*(3.330657 + X)))*DXP/X;

	if (E1 < 1.0e-100)
		return;

	double EI = DXP*(1.666667 - (2.0/3.0)*X)/X + E1*((2.0/3.0)*X - 1.0) - 0.5*E1*E1/DXP;
	QI = 5.444089*EI/double(TE32);
}

void DepartureCoeffs::HELPME(double& AID, int KM) {
	double SUM = 0;
	double C = 0;
	for (int i = 1; i <= LIMIT; ++i) {
		C = STORE2[i-1]*STORE3[i-1][KM-1];
		SUM += C;
	}
	SUM -= 0.5*C;

	double Y;
	Y =  0.6170614899993600e-2*(POL(1,  KM) + POL(2,  KM));
	Y += 0.1426569431446683e-1*(POL(3,  KM) + POL(4,  KM));
	Y += 0.2213871940870990e-1*(POL(5,  KM) + POL(6,  KM));
	Y += 0.2964929245771839e-1*(POL(7,  KM) + POL(8,  KM));
	Y += 0.3667324070554015e-1*(POL(9,  KM) + POL(10, KM));
	Y += 0.4309508076597664e-1*(POL(11, KM) + POL(12, KM));
	Y += 0.4880932605205694e-1*(POL(13, KM) + POL(14, KM));
	Y += 0.5372213505798282e-1*(POL(15, KM) + POL(16, KM));
	Y += 0.5775283402686280e-1*(POL(17, KM) + POL(18, KM));
	Y += 0.6083523646390170e-1*(POL(19, KM) + POL(20, KM));
	Y += 0.6291872817341415e-1*(POL(21, KM) + POL(22, KM));
	Y += 0.6396909767337608e-1*(POL(23, KM) + POL(24, KM));

	AID = SUM + B*Y;
}

void DepartureCoeffs::INTERP (const std::array<int, 75>& M, const std::array<double, 75>& CO, std::array<double, 507>& VAL, std::array<double, 507>& DVAL, int IC, int IR) {
	int LG = 2*(IR+1);
	int IB = IC - IR;
	int IBB = IB - 1;
	int ICC = IC - 1;
	int K = M[IC-1];

	std::array<int, 8> IQ;

	for (int i = 1; i <= K; ++i) {
		VAL[i-1] = 0;
		DVAL[i-1] = 0;
	}

	int IA = 0;
	for (int  i = 1; i <= IC; ++i) {
		IA = i;
		int MIT = M[i-1];
		VAL[MIT-1] = CO[i-1];

		if (i < IC && M[i] - MIT > 1)
			break;
	}

	if (IA == IC)
		return;

	if (IA != IB) {
		for (int i = IA; i <= IBB; ++i) {
			int N1 = M[i-1] + 1;
			int N2 = M[i];
			int ITR1 = i - IR - 1;

			for (int itau = 1; itau <= LG; ++itau) {
				int IND = ITR1 + itau;
				IQ[itau-1] = M[IND-1];
			}

			for (int itau = 1; itau <= LG; ++itau) {
				double PHITAU = 1.0/PHI(IQ, LG, itau, IQ[itau-1]);
				for (int N = N1; N <= N2; ++N) {
					double FL = PHI(IQ, LG, itau, N)*PHITAU;
					double DFL = DPHI(IQ, LG, itau, N)*PHITAU;

					int IND = ITR1 + itau;
					double COT = CO[IND-1];
					VAL[N-1] += FL*COT;
					DVAL[N-1] += DFL*COT;
				}
			}
		}
	}

	if (IR == 0)
		return;

	int ICLG = IC - LG;
	for (int itau = 1; itau <= LG; ++itau) {
		int IND = ICLG + itau;
		IQ[itau-1] = M[IND-1];
	}

	for (int itau = 1; itau <= LG; ++itau) {
		double PHITAU = 1.0 / PHI(IQ, LG, itau, IQ[itau-1]);

		for (int i = IB; i <= ICC; ++i) {
			int N1 = M[i-1] + 1;
			int N2 = M[i];

			for (int N = N1; N <= N2; ++N) {
				double FL = PHI(IQ, LG, itau, N)*PHITAU;
				double DFL = DPHI(IQ, LG, itau, N)*PHITAU;
				int IND = ICLG + itau;
				double COT = CO[IND-1];
				VAL[N-1] += FL*COT;
				DVAL[N-1] += DFL*COT;
			}
		}
	}
}

double SOS(int i, double a, double itm) {
	return std::pow(std::sqrt(-a), (2*i + itm))/std::log(-a);
}

/**/
void DepartureCoeffs::JMD (std::array<std::array<double, 75>, 75>& SK, std::array<double, 75>& CO, const std::array<int, 75>& MVAL, int IC) {
	LIMIT = 200;
	int NG = 24;

	for (int j = 1; j <= NFIT; ++j) {
		int k = IVAL[j-1];
		double AJ = -1.0/double(MVAL[k-1]*MVAL[k-1]);

		for (int i = 1; i <= NFIT; ++i)
			AFIT[j-1][i-1] = SOS(i, AJ, ITM);
	}

	std::array<double, 4> AZ;
	std::array<int, 4> IPIV;
	std::array<std::array<double, 2>, 4> IND;
	int IRROR = 0;
	double D;
	MATINV<4>(AFIT, NFIT, AZ, 0, D, IRROR, IPIV, IND);

	B = 1.0/std::pow(double(MVAL[IC-1] + LIMIT), 2.0);
	double A = -0.5*B;
	int NH = NG/2;

	for (int k = 1; k <= NH; ++k) {
		JMDVAL[2*k - 2] = A + VALUE[k-1]*B;
		JMDVAL[2*k - 1] = A - VALUE[k-1]*B;
	}

	for (int k = 1; k <= LIMIT; ++k) {
		A = MVAL[IC-1] + k;
		double AC = -1.0/(A*A);

		for (int j = 1; j <= NFIT; ++j) {
			STORE1[k-1][j-1] = SOS(j, AC, ITM);
		}
	}

	for (int k = 1; k <= NG; ++k) {
		for (int j = 1; j <= NFIT; ++j) {
			int INP = k + LIMIT;
			STORE1[INP-1][j-1] = SOS(j, JMDVAL[k-1], ITM);
		}
	}

	int KK = LIMIT + NG;
	for (int k = 1; k <= KK; ++k) {
		for (int j = 1; j <= NFIT; ++j) {
			STORE3[k-1][j-1] = DMJ1(k, j);
		}
		STORE3[k-1][NFIT] = DMJ2(k);
	}

	for (int j = 1; j <= IC; ++j) {
		int i = MVAL[j-1];

		for (int k = 1; k <= LIMIT; ++k) {
			KK = MVAL[IC-1] + k;
			int AKK = KK;
			//if (k <= 40)
				//STORE2[k-1] = BK2(i, KK, 0);
			//else
				STORE2[k-1] = BK(i, KK, 0);
		}

		for (int k = 1; k <= NG; ++k) {
			double AK = std::sqrt(-1.0/JMDVAL[k-1]);
			RTVAL[k-1] = std::pow(AK, 3)*0.5;
			KK = AK;
			int INP = k + LIMIT;
			//if (INP <= 40)
				//STORE2[INP-1] = BK2(i, KK, 0);
			//else
				STORE2[INP-1] = BK(i, KK, 0);
		}

		double AID;
		HELPME(AID, NFIT+1);
		CO[j-1] -= AID;

		for (int km = 1; km <= NFIT; ++km) {
			HELPME(AID, km);
			int L = IVAL[km-1];
			SK[j-1][L-1] += AID;
		}
	}

}


void DepartureCoeffs::RAD(double& A, int N, int NDASH, double& G) {
	int EN = N;
	int EN2 = N*N;
	int END2 = NDASH*NDASH;
	int NDMN = NDASH - N;

	if (NDMN > 50) {
		double ALF2 = EN2/double(END2);
		double TERM1 = 1.0 - ALF2;
		TERM1 = std::pow(TERM1*EN, 2.0/3.0);
		double TERM2 = 0.1728*(1.0 + ALF2)/TERM1;
		double TERM3 = 0.04960*(1.0 - 1.3333333333333*ALF2 + ALF2*ALF2)/(TERM1*TERM1);
		G = 1.0 - TERM2 - TERM3;
	}
	else {
		if (N > 10) {
			bool smaller = false;

			for (int j = 2; j <= 12; ++j) {
				double IVJ = IXV[j-1];
				double IVJ1 = IXV[j-2];
				smaller = N <= IVJ;

				if (smaller) {
					double P = 1.0/(IVJ - IVJ1);
					double Q = (IVJ - N)*P;
					P *= (N - IVJ1);
					G = Q*GAUNT[NDMN-1][j+9-1] + P*GAUNT[NDMN-1][j+10-1];
					break;
				}
			}
			if (!smaller)
				G = GAUNT[NDMN-1][N-1];
		}
		else {
			G = GAUNT[NDMN-1][N-1];
		}
	}

	A = G*15.7457e9/(double(EN)*(END2-EN2)*NDASH*END2);
}


void DepartureCoeffs::RADCOL(double T, const std::array<int, 75>& MVAL, int IC, int NMIN) {
	 // Radiation cascade coefficients.
	 // Sum from N to all levels down to NMIN (=1 or 2).
	for (int m = 1; m <= IC; ++m) {
		int j = MVAL[m-1];
		double TOT = 0.0;

		if (j > NMIN) {
			int k = j - 1;

			for (int i = NMIN; i <= k; ++i) {
				double A, G;
				RAD(A, i, j, G);
				TOT += A;
			}
		}
		RADTOT[m-1] = TOT;
	}

	// Collisional rate totals on to level N.
	// Summed from NMIN to INFINITY.
	//
	// Here also, Malcom Walmsley has changed the lower principal quantum
	// number from 20 to 1 (CASE A) or 2 (CASE B) for the rates.
	int NLO = NMIN;
	for (int m = 1; m <= IC; ++m) {
		int N = MVAL[m-1];
		int AN = N;
		double TOT = 0.0;
		int L = N - 1;

		for (int kk = NLO; kk <= L; ++kk) {
			double C = COLRAT(kk, N, T, TE12);
			TOT += C;
		}

		L = N + 1;
		int MAX = N + 40;

		for (int kk = L; kk <= MAX; ++kk) {
			int AKK = kk;
			double C = COLRAT(N, kk, T, TE12);
			double CX = CXP[kk-1];

			double Q;
			if (CX < 1e-30 || CXP[N-1] < 1.0e-30)
				Q = std::exp(15.778e4*(1.0/double(AKK*AKK) - 1.0/double(AN*AN))/T);
			else
				Q = CXP[N-1]/CX;

			C *= (AKK/double(AN))*(AKK/double(AN))*Q;
			TOT += C;
		}

		COLTOT[m-1] = TOT;
	}
}

/**
 *  Computes recombination coefficient, ALPHA, ont level N1 for ions of
 *  effective charge Z at electron temperature ETEMP.
 */
void DepartureCoeffs::RECOMB(double Z, double ETEMP, int N1, double& ALPHA) {
	std::array<double, 99> XV;

	int Z2 = Z*Z;
	double TE = ETEMP*1.0e-4;
	double CONST = 15.778/TE;
	double FL = 1.0/std::pow(CONST*Z2, 1.0/3.0);
	CONST = 5.197e-14*Z2*std::sqrt(CONST);
	XV[0] = 0.02;

	for (int N = 2; N <= 11; ++N)
		XV[N-1] = XV[N-2] + 0.002;
	for (int N = 12; N <= 23; ++N)
		XV[N-1] = XV[N-2] + 0.005;
	for (int N = 24; N <= 33; ++N)
		XV[N-1] = XV[N-2] + 0.01;
	for (int N = 34; N <= 65; ++N)
		XV[N-1] = XV[N-33]*10.0;
	for (int N = 66; N <= 97; ++N)
		XV[N-1] = XV[N-33]*10.0;

	int F = N1;
	double X = 15.778*Z2/(F*F*TE);

	double S0, S1, S2;
	if (X < 0.02) {
		double XXX = std::pow(X, 1.0/3.0);
		S0 = X*std::exp(X)*(-std::log(X) - 0.5772 + X);
		S1 = 0.46290*X*(1.0 + 4.0*X) - 1.03680*std::pow(XXX, 4)*(1.0 + 1.875*X);
		S2 = -0.06720*X*(1.0 + 3.0*X) + 0.14880*std::pow(XXX, 5)*(1.0 + 1.80*X);
	}
	else if (X > 20.0) {
		double U = 1.0/X;
		double V = U/3.0;
		double XXX = std::pow(X, 1.0/3.0);
		S0 = 1.0 - U*(1.0 - 2.0*U*(1.0 - 3.0*U*(1.0 - 4.0*U)));
		S1 = -0.17280*XXX*(1.0 - V*(8.0 - V*(70.0 - V*(800.0 - V*11440.0))));
		S2 = -0.04960*XXX*XXX*(1.0 - V*(3.0 - V*(32.0 - V*448.0)));
	}
	else {
		for (int i = 2; i <= 99; ++i) {
			double XVI = XV[i-1];
			if (X > XVI)
				continue;
			int IM1 = i - 1;
			double XVI1 = XV[IM1-1];
			double P = 1.0/(XVI - XVI1);
			double Q = (XVI - X)*P;
			P *= (X - XVI1);
			S0 = P*SV0[i-1] + Q*SV0[IM1-1];
			S1 = P*SV1[i-1] + Q*SV1[IM1-1];
			S2 = P*SV2[i-1] + Q*SV2[IM1-1];
			break;
		}
	}

	double Y = (S0 + FL*(S1 + FL*S2))/F;
	ALPHA = CONST*Y;
}

void DepartureCoeffs::REDUCE(const std::array<int, 75>& M, int IC, int IR, std::array<std::array<double, 75>, 75>& SK) {
	std::array<int, 8> IQ;
	std::array<int, 8> STORE_A, STORE_B;

	int LG = 2*(IR + 1);
	int IB = IC - IR;
	int IBB = IB - 1;
	int ICC = IC - 1;

	for (int i = 1; i <= IC; ++i) {
		for (int j = 1; j <= IC; ++j)
			SK[i-1][j-1] = BK(M[i-1], M[j-1], i);
	}

	int IA = 0;
	for (int i = 1; i <= IC; ++i) {
		IA = i;
		if (M[i] - M[i-1] > 1)
			break;
	}

	if (IA == IC)
		return;
	if (IA != IB) {
		for (int i = IA; i <= IBB; ++i) {
			int N1 = M[i-1] + 1;
			int N2 = M[i] - 1;

			for (int itau = 1; itau <= LG; ++itau)
				IQ[itau-1] = M[i-IR-1+itau-1];
			for (int itau = 1; itau <= LG; ++itau) {
				STORE_A[itau-1] = PHI(IQ, LG, itau, IQ[itau-1]);
			}

			for (int N = N1; N <= N2; ++N) {
				for (int itau = 1; itau <= LG; ++itau)
					STORE_B[itau-1] = PHI(IQ, LG, itau, N);
				for (int j = 1; j <= IC; ++j) {
					double DUCKIT = BK(M[j-1], N, j);
					for (int itau = 1; itau <= LG; ++itau) {
						double FL = STORE_B[itau-1]/double(STORE_A[itau-1]);
						int IND = i - IR - 1 + itau;
						SK[j-1][IND-1] += DUCKIT*FL;
					}
				}
			}
		}
	}

	if (IR == 0)
		return;

	for (int itau = 1; itau <= LG; ++itau) {
		int IND = IC - LG + itau;
		IQ[itau-1] = M[IND-1];
	}

	for (int itau = 1; itau <= LG; ++itau) {
		double PHITAU = 1.0/PHI(IQ, LG, itau, IQ[itau-1]);

		for (int i = IB; i <= ICC; ++i) {
			int N1 = M[i-1] + 1;
			int N2 = M[i] - 1;

			for (int N = N1; N <= N2; ++N) {
				double FL = PHI(IQ, LG, itau, N)*PHITAU;

				for (int j = 1; j <= IC; ++j) {
					int IND = IC - LG + itau;
					SK[j-1][IND-1] += BK(M[j-1], N, j-1)*FL;
				}
			}
		}
	}
}

void DepartureCoeffs::RHS(std::array<double, 75>& CO, const std::array<int, 75>& MVAL, int IC) {
	for (int i = 1; i <= IC; ++i) {
		int j = MVAL[i-1];
		double RT, ALFA;
		COLION(j, 1, T, RT);
		RECOMB(1.0, T, j, ALFA);
		CO[i-1] = -ALFA*CXP[j-1]*TE32*0.24146879e16/double(j*j) - RT*DENS;
	}
}



