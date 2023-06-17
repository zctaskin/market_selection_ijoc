#include "IterativeHeuristic.h"
#include "MarketSelection.h"
#include <chrono>

IterativeHeuristic::IterativeHeuristic(MarketSelection& MSIn, ParameterMap& PM, double lambdaIn) : MS(MSIn), Parameters(PM), lambda(lambdaIn)
{
	BestM.resize(MS.M);
	BestT.resize(MS.T);
	LB = -1;
	CPU = 0;
}

bool IterativeHeuristic::Solve()
{
	// Inputs
	DoubleVector Revenue;
	vector<double*> d;
	for (auto& m : MS.Markets)
	{
		Revenue.push_back(m.Revenue * lambda);
		d.push_back(m.Demand.data());
	}

	DoubleVector SetupCost = MS.SetupCost;
	DoubleVector ProductionCost = MS.ProductionCost;
	DoubleVector HoldingCost = MS.HoldingCost;
	for (int t = 0; t < MS.T; ++t)
	{
		SetupCost[t] *= (1 - lambda);
		ProductionCost[t] *= (1 - lambda);
		HoldingCost[t] *= (1 - lambda);
	}

	// Outputs
	double cntIA = 0;
	int mcntIA = 0;

	auto startTime = chrono::high_resolution_clock::now();
	IAT2(MS.T, MS.M, Revenue.data(), d.data(), SetupCost.data(), ProductionCost.data(), HoldingCost.data(), BestM.data(), BestT.data(), LB, cntIA, mcntIA);
	CPU = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

	return LB >= 0;
}

double IterativeHeuristic::GetRevenue()
{
	double Revenue = 0;
	for (int m = 0; m < MS.M; ++m)
		Revenue += BestM[m] * MS.Markets[m].Revenue;
	return Revenue;
}

double IterativeHeuristic::GetCost()
{
	double ScaledRevenue = 0;
	for (int m = 0; m < MS.M; ++m)
		ScaledRevenue += BestM[m] * MS.Markets[m].Revenue * lambda;
	return (ScaledRevenue - LB)/(1-lambda);
}

size_t IterativeHeuristic::GetSelectedMarketCount()
{
	int nSelected = 0;
	for (int m = 0; m < MS.M; ++m)
		if (BestM[m])
			++nSelected;
	return nSelected;
}

void IterativeHeuristic::IAT2(int T, int nM, double* R, double** d, double* K, double* p, double* h, int* bestM, int* bestT, double& best, double& cntIA, int& mcntIA)
{
	int* OnesM, * OnesT, ind, cnt, nI = 0;
	bool* candT;
	double Prof, eps = 0.0001, * cD;
	OnesM = new int[nM];
	OnesT = new int[T];
	candT = new bool[T];
	cD = new double[T];
	cntIA = 0;
	mcntIA = 0;

	// Trivial solution
	best = 0;
	for (int i = 0; i < T; i++)
	{
		bestT[i] = 0;
		candT[i] = false;
	}
	for (int m = 0; m < nM; m++)
	{
		bestM[m] = 0;
	}

	//Find candidate start setup periods
	for (int m = 0; m < nM; m++)
	{
		ind = 0;
		while (d[m][ind] < eps && ind < T - 1)
			ind++;
		candT[ind] = true;
	}

	// Apply the IA starting with T^2 solutions
	for (int n = 0; n < T; n++)
	{
		if (candT[n])
		{
			for (int i = 0; i < (T - n); i++)
			{
				nI++;
				for (int j = 0; j < T; j++)
					OnesT[j] = 0;
				for (int j = 0; j <= i; j++)
				{
					ind = n + int(j * (T - n) / (i + 1));
					OnesT[ind] = 1;
				}

				IterateT(T, nM, R, d, K, p, h, OnesM, OnesT, Prof, cnt);
				cntIA += cnt;
				if (cnt > mcntIA)
					mcntIA = cnt;

				if (Prof > best)
				{
					best = Prof;
					for (int j = 0; j < T; j++)
					{
						bestT[j] = OnesT[j];
					}
					for (int m = 0; m < nM; m++)
					{
						bestM[m] = OnesM[m];
					}
				}
			}
		}
	}
	cntIA = cntIA / nI;

	// Print the final demand vector
	for (int j = 0; j < T; j++)
		cD[j] = 0;
	for (int m = 0; m < nM; m++)
	{
		if (bestM[m] == 1)
		{
			for (int j = 0; j < T; j++)
				cD[j] += d[m][j];
		}
	}
	cout << "IA demand vector:" << endl;
	for (int i = 0; i < T; i++)
		cout << cD[i] << ' ';
	cout << endl;

	delete[] OnesM;
	delete[] OnesT;
	delete[] candT;
	delete[] cD;
}

void IterativeHeuristic::IterateT(int T, int nM, double* R, double** d, double* K, double* p, double* h, int* OnesM, int* OnesT, double& Prof, int& cnt)
{
	double pr_old, cost;
	bool flag = true, firstIt = true, type1 = false, type2 = false; // type 1 (2): change from 0 to 1 (1 to 0)
	int* oldT, * nChanges;
	oldT = new int[T];
	nChanges = new int[nM]; // counts the number of changes of a partcilar market
	cnt = 0;

	for (int i = 0; i < T; i++)
	{
		oldT[i] = OnesT[i];
	}

	pr_old = -pow(10.0, 10);
	do
	{	// find the best market selection given a prod plan
		GivenT(oldT, nM, R, T, d, K, p, h, OnesM, Prof);

		if (Prof < pr_old - 0.001)
		{
			cout << "Error" << endl;
		}
		pr_old = Prof;

		// find the best prod plan given a market selection
		Prof = 0;
		for (int m = 0; m < nM; m++)
		{
			if (OnesM[m] == 1)
			{
				Prof += R[m];
			}
		}
		GivenM(OnesM, nM, T, d, K, p, h, OnesT, cost);
		Prof -= cost;

		AllM[IntegerVector{OnesM, OnesM + nM}] = Prof;

		if (Prof < pr_old - 0.001)
		{
			cout << "Error" << endl;
		}
		pr_old = Prof;

		// check whether a market selection changes like 0-1-0 or 1-0-1
		if (firstIt)
		{
			for (int m = 0; m < nM; m++)
				nChanges[m] = OnesM[m];
			firstIt = false;
		}
		else
		{
			for (int m = 0; m < nM; m++)
			{
				if ((nChanges[m] % 2 == 0) && (OnesM[m] == 1))
				{
					nChanges[m] += 3;
					type1 = true;
				}
				else if ((nChanges[m] % 2 == 1) && (OnesM[m] == 0))
				{
					nChanges[m] += 1;
					type2 = true;
				}
			}
		}

		// check whether two consecutive production plans are identical
		flag = true;
		for (int i = 0; i < T; i++)
		{
			if (oldT[i] != OnesT[i])
			{
				flag = false;
				oldT[i] = OnesT[i];
			}
		}
		cnt++;

	} while (flag == false);
	//cout << endl;

	delete[] oldT;
	delete[] nChanges;
}

// Given the production plan, find the optimal subset of markets
void IterativeHeuristic::GivenT(int* OnesT, int nM, double* R, int T, double** d, double* K, double* p, double* h, int* OnesM, double& Prof)
{
	double* Rev, totSetup, c, eps = 0.0001;
	int start;
	bool* valid;
	valid = new bool[nM];
	Rev = new double[nM];

	for (int m = 0; m < nM; m++)
	{
		Rev[m] = R[m];
		valid[m] = true;
	}
	totSetup = 0;
	Prof = 0;

	start = 0; // first setup period
	while (OnesT[start] == 0)
	{
		for (int m = 0; m < nM; m++)
		{
			if (d[m][start] > eps)
				valid[m] = false;
		}
		start++;
	}

	for (int i = start; i < T; i++)
	{
		if (OnesT[i] == 1)
		{
			c = p[i];
			totSetup += K[i];
		}
		else
		{
			c += h[i - 1];
		}

		for (int m = 0; m < nM; m++)
		{
			Rev[m] -= c * d[m][i];
		}
	}

	for (int m = 0; m < nM; m++)
	{
		if ((Rev[m] > 0) && valid[m])
		{
			OnesM[m] = 1;
			Prof += Rev[m];
		}
		else
		{
			OnesM[m] = 0;
		}
	}

	Prof -= totSetup;

	delete[] Rev;
	delete[] valid;
}

// Given the markets, find the optimal production plan
void IterativeHeuristic::GivenM(int* OnesM, int nM, int T, double** d, double* K, double* p, double* h, int* OnesT, double& cost)
{
	double* cD;
	int cnt;
	cD = new double[T];

	//Calculate demands
	for (int j = 0; j < T; j++)
	{
		cD[j] = 0;
	}
	for (int m = 0; m < nM; m++)
	{
		if (OnesM[m] == 1)
		{
			for (int j = 0; j < T; j++)
			{
				cD[j] += d[m][j];
			}
		}
	}

	//Calculate costs
	WWgeneralBW(T, cD, K, p, h, cost, OnesT, cnt);

	delete[] cD;
}

void IterativeHeuristic::WWgeneralBW(int T, double* d, double* K1, double* p1, double* h1, double& opt, int* Ones, int& counter)
{
	double* f, c, * K, * p, * h, eps = 0.0001;
	int* s, fstNZ = 0; // first non-zero demand
	counter = 0; // counts the number of iterations needed

	f = new double[T + 1]; f[T] = 0;	// f[t]: optimal cost for period t,...,T	
	K = new double[T];
	p = new double[T];
	h = new double[T];
	s = new int[T + 1];
	s[T] = T;
	for (int i = 0; i < T; i++)
	{
		p[i] = p1[i];
		h[i] = h1[i];
		K[i] = K1[i];
	}

	while ((fstNZ < T) && (d[fstNZ] < eps))
		fstNZ++;
	s[fstNZ] = fstNZ;
	f[fstNZ] = 0;

	// Redefine cost parameters -> problem without holding cost
	p[T - 1] = p[T - 1] + h[T - 1];
	opt = -h[T - 1] * d[T - 1];
	for (int i = T - 2; i >= 0; i--)
	{
		h[i] = h[i] + h[i + 1];
		opt -= h[i] * d[i];
		p[i] += h[i];
	}

	// start the recursion
	int period = 1;
	double tmp, min;
	for (int j = T - 1; j >= 0; j--)
	{
		min = pow(10.0, 10);

		for (int i = j; i < T; i++)
		{
			counter++;

			// Calculate cost parameters
			// c:=cost from period i,...,j
			if (i == j)
				c = K[i] + p[i] * d[i];
			else
				c += p[j] * d[i];

			// Find the optimum cost for period j
			tmp = c + f[i + 1];
			if (tmp < min)
			{
				min = tmp;
				period = i + 1;
			}
		}

		f[j] = min;
		s[j] = period;	// s[j]: optimal next production period after period j		
	}

	// find the optimal first production period
	int per = fstNZ;
	double tmpOpt = f[fstNZ];
	for (int i = fstNZ - 1; i >= 0; i--)
	{
		if (f[i] < tmpOpt)
		{
			tmpOpt = f[i];
			per = i;
		}
	}
	opt += tmpOpt;

	// assign the production periods
	for (int i = 0; i < per; i++)
		Ones[i] = 0;

	while (per < T)
	{
		Ones[per] = 1;
		for (int i = per + 1; i <= s[per] - 1; i++)
		{
			Ones[i] = 0;
			if (i >= T)
				cout << "Error" << endl;
		}
		per = s[per];
	}

	delete[] f;
	delete[] K;
	delete[] p;
	delete[] h;
	delete[] s;
}