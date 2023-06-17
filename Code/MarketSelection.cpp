#include "MarketSelection.h"
#include <iostream>
#include <algorithm>
#include <sstream>


void MarketSelection::SetDimensions(int m, int t)
{
	T = t;
	SetupCost.resize(T);
	ProductionCost.resize(T);
	HoldingCost.resize(T);

	Markets.resize(M);
	for (int i = 0; i < m; ++i)
	{
		Markets[i].Index = i;
		Markets[i].Demand.resize(T);
	}
}

void MarketSelection::Resize(int NewM, int NewT)
{
	auto ResizeFunction = [](auto& Vect, int NewSize)
	{
		int CurrentSize = Vect.size();
		if (NewSize < CurrentSize)
			Vect.resize(NewSize);
		else if (NewSize > CurrentSize)
		{
			Vect.reserve(NewSize);
			for (int i = 0; i < NewSize - CurrentSize; ++i)
				Vect.push_back(Vect[i % CurrentSize]);
		}
	};

	ResizeFunction(Markets, NewM);
	ResizeFunction(SetupCost, NewT);
	ResizeFunction(ProductionCost, NewT);
	ResizeFunction(HoldingCost, NewT);
	for (auto& m : Markets)
		ResizeFunction(m.Demand, NewT);

	if (NewT > T)
	{
		for (auto& m : Markets)
			m.Revenue *= double(NewT) / T;
	}

	M = NewM;
	T = NewT;
}

void MarketSelection::ScaleRevenueCost(double lambda)
{
	for (auto& m : Markets)
		m.Revenue *= lambda;

	for (int t = 0; t < T; ++t)
	{
		SetupCost[t] *= (1 - lambda);
		ProductionCost[t] *= (1 - lambda);
		HoldingCost[t] *= (1 - lambda);
	}

	CalculateStatistics();
}

void MarketSelection::ReadData(ifstream& file)
{
	const int MAX_LINE_LENGTH = 500;
	char DummyLine[MAX_LINE_LENGTH];
	string DummyToken;

	file.getline(DummyLine, MAX_LINE_LENGTH);
	file >> DummyToken >> T;
	file >> DummyToken >> M;

	SetDimensions(M, T);

	file.getline(DummyLine, MAX_LINE_LENGTH); file.getline(DummyLine, MAX_LINE_LENGTH);

	for (int t = 0; t < T; ++t)
		file >> SetupCost[t] >> ProductionCost[t] >> HoldingCost[t];

	file.getline(DummyLine, MAX_LINE_LENGTH); file.getline(DummyLine, MAX_LINE_LENGTH);

	for (int m = 0; m < M; ++m)
		file >> Markets[m].Revenue;

	file.getline(DummyLine, MAX_LINE_LENGTH); file.getline(DummyLine, MAX_LINE_LENGTH);
	for (int m = 0; m < M; ++m)
		for (int t = 0; t < T; ++t)
			file >> Markets[m].Demand[t];
}

void MarketSelection::CalculateStatistics()
{
	MaxDemand = -DBL_MAX;
	TotalDemand = 0;
	TotalRevenue = 0;
	AggregateDemand.resize(T, 0);
	for (int m = 0; m < M; ++m)
	{
		auto& Market = Markets[m];
		Market.TotalDemand = 0;
		for (int t = 0; t < T; ++t)
		{
			AggregateDemand[t] += Market.Demand[t];
			Market.TotalDemand += Market.Demand[t];
			if (Market.Demand[t] > MaxDemand)
				MaxDemand = Market.Demand[t];
		}
		TotalDemand += Market.TotalDemand;
		TotalRevenue += Market.Revenue;
	}

	{
		double ConstantH;
		IntegerVector Ones(T);
		DoubleVector Duals(T, 0);
		int counter;
		WWgeneralPD(T, AggregateDemand.data(), SetupCost.data(), ProductionCost.data(), HoldingCost.data(), TotalCost, Ones.data(), Duals.data(), ConstantH, counter);
	}

	for (int m = 0; m < M; ++m)
	{
		auto& Market = Markets[m];
		Market.DemandRatio = Market.TotalDemand / TotalDemand;

		IntegerVector Ones(T);
		DoubleVector Duals(T, 0);
		int counter;
		DoubleVector ZeroSetup(T, 0);
		WWgeneralPD(T, Market.Demand.data(), ZeroSetup.data(), ProductionCost.data(), HoldingCost.data(), Market.VariableCost, Ones.data(), Duals.data(), Market.ConstantH, counter);

		WWgeneralPD(T, Market.Demand.data(), SetupCost.data(), ProductionCost.data(), HoldingCost.data(), Market.TotalCost, Ones.data(), Duals.data(), Market.ConstantH, counter);

		DoubleVector OtherDemand = AggregateDemand;
		for (int t = 0; t < T; ++t)
			OtherDemand[t] -= Market.Demand[t];
		double OtherCost, ConstantH;
		WWgeneralPD(T, OtherDemand.data(), SetupCost.data(), ProductionCost.data(), HoldingCost.data(), OtherCost, Ones.data(), Duals.data(), ConstantH, counter);
		Market.IncrementalCost = TotalCost - OtherCost;
	}
	{
		MinSetupCost = *min_element(begin(SetupCost), end(SetupCost));
	}
}

double MarketSelection::GetCumulativeCost(int i, int j)
{
	double Cost = ProductionCost[i];
	for (int t = i; t < j; ++t)
		Cost += HoldingCost[t];
	return Cost;
}

double MarketSelection::CalculateRevenue(const IntegerVector& Selected)
{
	double Revenue = 0;
	for (int m = 0; m < M; ++m)
		if (Selected[m])
			Revenue += Markets[m].Revenue;
	return Revenue;
}

double MarketSelection::CalculateCost(const IntegerVector& Selected)
{
	DoubleVector Demand(T, 0);
	for (int m = 0; m < M; ++m)
	{
		if (m < Selected.size() && Selected[m])
		{
			auto& Market = Markets[m];
			for (int t = 0; t < T; ++t)
				Demand[t] += Market.Demand[t];
		}
	}

	double ConstantH;
	double TotalCost;

	IntegerVector Ones(T);
	DoubleVector Duals(T, 0);
	int counter;
	WWgeneralPD(T, Demand.data(), SetupCost.data(), ProductionCost.data(), HoldingCost.data(), TotalCost, Ones.data(), Duals.data(), ConstantH, counter);

	return TotalCost;
}

// Wagner Whitin algorithm which computes both primal and dual 
void MarketSelection::WWgeneralPD(int T, double* d, double* K1, double* p1, double* h1, double& opt, int* Ones, double* mu, double& constH, int& counter)
{
	double* f, * c, * K, * p, * h, eps = 0.0001, *v, objD = 0.0;
	int* s, fstNZ = 0;				// first non-zero demand
	bool printOutput = false;

	f = new double[T + 1]; f[0] = 0;	// f[t]: optimal cost for period 1,...,t
	c = new double[T];					// c[t]: se below
	counter = 0;						// counts the number of iterations needed
	K = new double[T];
	p = new double[T];
	h = new double[T];
	s = new int[T + 1];
	s[0] = 0;
	for (int i = 0; i < T; i++)
	{
		p[i] = p1[i];
		h[i] = h1[i];
		K[i] = K1[i];
	}

	// compute first non-zero demand period
	while ((fstNZ < T) && (d[fstNZ] < eps))
		fstNZ++;
	s[fstNZ] = fstNZ;
	f[fstNZ] = 0;

	/////////////////////////////////
	// Primal problem computations //
	/////////////////////////////////

	// Redefine cost parameters -> transformed problem without holding cost
	p[T - 1] = p[T - 1] + h[T - 1];
	opt = -h[T - 1] * d[T - 1];
	for (int i = T - 2; i >= 0; i--)
	{
		h[i] = h[i] + h[i + 1];
		opt -= h[i] * d[i];
		p[i] += h[i];
	}
	constH = opt;

	// start the forward recursion
	int period = 1;
	double tmp, min;
	for (int i = fstNZ; i < T; i++)
	{
		min = pow(10.0, 10);

		for (int j = fstNZ; j <= i; j++)
		{
			counter++;
			// Calculate cost parameters
			// N.B. c[j-1]:=cost of periods j till i
			if (j == i)
				c[j] = K[j] + p[j] * d[j];
			else
				c[j] += p[j] * d[i];

			// Find the optimum cost for period i
			tmp = f[j] + c[j];
			if (tmp < min)
			{
				min = tmp;
				period = j;
			}
		}

		f[i + 1] = min;
		s[i + 1] = period + 1;	// s[i+1]: optimal previous production period for period i+1
		//cout << i+1 << ' ' << s[i+1] << endl;
	}

	opt += f[T];

	// set the setup periods
	int per = T;
	do
	{
		Ones[s[per] - 1] = 1;
		for (int i = s[per]; i <= per - 1; i++)
		{
			Ones[i] = 0;
		}
		per = s[per] - 1;
	} while (per > fstNZ);

	for (int i = 0; i < fstNZ; i++)
		Ones[i] = 0;

	if (printOutput)
	{
		cout << "Setup periods: ";
		for (int i = 0; i < T; i++)
		{
			cout << Ones[i] << " ";
		}
		cout << endl;
	}

	///////////////////////////////
	// Dual problem computations //
	///////////////////////////////

	// Compute duals of transformed model
	v = new double[T]; // dual variables

	// initialize duals 
	for (int i = 0; i < T; i++)
	{
		v[i] = 0.0;
	}

	if (fstNZ < T) // else all duals are zero
	{
		// Compute duals and objective
		for (int i = fstNZ; i < T; i++)
		{
			if (d[i] > eps)
			{
				v[i] = (f[i + 1] - f[i]) / d[i];
			}
			else
			{
				v[i] = 0;
			}
			objD += v[i] * d[i];
		}

		// Compute objective to check correctness
		if (isnan(objD) || isinf(objD) || abs(objD - f[T]) > eps)
		{
			cout << "Primal not equal to (transformed) dual objectve!" << endl;
		}

		// Check dual feasibility to check correctness
		double lhs = 0.0;
		for (int i = 0; i < T; i++)
		{
			lhs = 0.0;
			for (int t = i; t < T; t++)
			{
				lhs += d[t] * fmax(v[t] - p[i], 0.0);
			}
			if (lhs > K[i] + eps)
			{
				cout << "Dual (transformed) is infeasible!" << endl;
			}
		}
	}

	if (printOutput)
	{
		// Print duals
		cout << "Dual transformed: ";
		for (int i = 0; i < T; i++)
		{
			cout << v[i] << " ";
		}
		cout << endl;
	}

	// Dual problem computations: second approach (non-transformed model)

	// precompute the FL cost
	double** cst, slack, tmpMu, objD2 = 0.0;
	cst = new double* [T];
	for (int i = 0; i < T; i++)
		cst[i] = new double[T];

	// compute unit cost [i,t]
	for (int i = 0; i < T; i++)
	{
		cst[i][i] = p1[i];
		for (int t = i + 1; t < T; t++)
		{
			cst[i][t] = cst[i][t - 1] + h1[t - 1];
		}
	}

	// compute duals mu
	for (int j = 0; j < T; j++)
	{
		if (d[j] > eps)
		{
			mu[j] = cst[j][j] + K1[j] / d[j]; // case l = j
			for (int l = 0; l <= j - 1; l++)
			{
				slack = K1[l];
				for (int t = l; t < j; t++)
				{
					slack -= d[t] * fmax(mu[t] - cst[l][t], 0);
				}
				tmpMu = cst[l][j] + slack / d[j];


				if (tmpMu < mu[j])
					mu[j] = tmpMu;
			}
		}
		else
		{
			mu[j] = 0;
		}

	}

	// compute objective dual to check correctness
	for (int i = 0; i < T; i++)
	{
		objD2 += mu[i] * d[i];
	}
	if (isnan(objD2) || isinf(objD2) || abs(objD2 - opt) > eps)
	{
		cout << "Primal not equal to dual (simple plant) objective!" << endl;
	}

	if (printOutput)
	{
		// Print duals
		cout << "Dual simple plant: ";
		for (int i = 0; i < T; i++)
		{
			cout << mu[i] << " ";
		}
		cout << endl << endl;
	}


	delete[] f;
	delete[] c;
	delete[] K;
	delete[] p;
	delete[] h;
	delete[] s;
	delete[] v;
	for (int i = 0; i < T; i++)
		delete[] cst[i];
	delete[] cst;
	//delete[] mu;
}

void MarketSelection::GetMu(int noT, int noM, const DoubleMatrix& cCost, const DoubleVector& fixedCost, const DoubleMatrix& demand, const DoubleVector& zBar, double& Objective, DoubleMatrix& muValue)
{
	for (int t = 0; t < noT; t++) {
		for (int m = 0; m < noM; m++) {
			if (zBar[m] == 0 || demand[m][t] == 0)
				muValue[m][t] = 0;
			else {
				double minVal = DBL_MAX;
				for (int l = 0; l <= t; l++) {
					int parantIns = cCost[l][t];
					double pay = fixedCost[l];
					for (int lp = l; lp < t; lp++) {
						for (int mp = 0; mp < noM; mp++) {
							pay -= demand[mp][lp] * max<double>(0, muValue[mp][lp] - cCost[l][lp]);
						}
					}
					for (int mp = 0; mp < m; mp++) {
						pay -= demand[mp][t] * max<double>(0, muValue[mp][t] - cCost[l][t]);
					}
					pay = pay / demand[m][t];

					if (parantIns + pay < minVal)
						minVal = parantIns + pay;
				}
				muValue[m][t] = minVal;
			}
		}
	}
	int sum = 0;
	for (int m = 0; m < noM; m++) {
		for (int t = 0; t < noT; t++) {
			sum += muValue[m][t] * demand[m][t];
		}
	}
	Objective = sum;
}

string to_string(const IntegerVector& vect)
{
	stringstream str;
	str << "[";
	for (auto i : vect)
		str << i << " ";
	str << "]";
	return str.str();
}
