#ifndef MARKET_SELECTION
#define MARKET_SELECTION

#include <vector>
#include <fstream>
#include <map>
#include <set>
#include <unordered_map>

using namespace std;

#define VALID_INEQUALITY_TOTAL_COST 1
//#define VALID_INEQUALITY_NON_TRIVIAL 2
#define VALID_INEQUALITY_FULL_SELECTION 4
#define VALID_INEQUALITY_RECOVER_LOSS 8

typedef vector<int> IntegerVector;
typedef vector<double> DoubleVector;
typedef set<double> DoubleSet;
typedef vector<DoubleVector> DoubleMatrix;
typedef map<IntegerVector, double> SolutionMap;

typedef map<double, IntegerVector> CostSolutionMap;
typedef map<double, CostSolutionMap> RevenueCostSolutionDetailMap;

typedef map<double, pair<double, IntegerVector>> RevenueCostSolutionMap;
typedef map<double, pair<double, IntegerVector>> CostRevenueSolutionMap;

struct Market
{
	int Index;
	double Revenue;
	DoubleVector Demand;

	double TotalDemand;
	double DemandRatio;

	double ConstantH;
	double TotalCost;		//Minimum setup + production + inventory cost of this market
	double VariableCost;	//Minimum production + inventory cost of this market
	double IncrementalCost; //Cost of satisfying all markets - (all except this one)
};

struct MarketSelection
{
	int T;
	int M;

	double MaxDemand;
	double TotalDemand;
	double TotalRevenue;
	double TotalCost;
	double MinSetupCost;

	vector<Market> Markets;
	DoubleVector AggregateDemand;
	DoubleVector SetupCost;
	DoubleVector ProductionCost;
	DoubleVector HoldingCost;

	void SetDimensions(int m, int t);
	void ReadData(ifstream& file);
	void Resize(int NewM, int NewT);
	void ScaleRevenueCost(double lambda);

	void CalculateStatistics();
	double GetCumulativeCost(int i, int j);		//Cumulative Production + Holding cost for production at i for demand at j

	double CalculateCost(const IntegerVector& Selected);
	double CalculateRevenue(const IntegerVector& Selected);


	void WWgeneralPD(int T, double* d, double* K1, double* p1, double* h1, double& opt, int* Ones, double* mu, double& constH, int& counter);
	void GetMu(int noT, int noM, const DoubleMatrix& cCost, const DoubleVector& fixedCost, const DoubleMatrix& demand, const DoubleVector& zBar, double& Objective, DoubleMatrix& muValue);
};

string to_string(const IntegerVector& vect);
#endif // !MARKET_SELECTION

