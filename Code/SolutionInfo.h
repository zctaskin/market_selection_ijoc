#pragma once

#include <vector>
using namespace std;

typedef vector<int> IntegerVector;

struct SolutionInfo
{
	SolutionInfo(double lambda, double objVal, double revenue, double cost, const IntegerVector& selected);

	double Lambda;
	double ObjVal;				// weighted objective value
	double Revenue;				// revenue objective value
	double Cost;				// cost objective value

	double CPUTime_Heuristic;
	double CPUTime_MIP;
	double CPUTime_AutoBenders;
	double CPUTime_Decomposition;
	double CPUTime_DecompositionReuse;

	IntegerVector SelectedMarkets;

	double GetIntercept();
	double GetSlope();
	void Print();
};
