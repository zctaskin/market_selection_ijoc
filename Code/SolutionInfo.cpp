#include "SolutionInfo.h"
#include <iostream>

SolutionInfo::SolutionInfo(double lambda, double objVal, double revenue, double cost, const IntegerVector& selected)
{
	Lambda = lambda;
	ObjVal = objVal;
	Revenue = revenue;
	Cost = cost;

	CPUTime_Heuristic = 0;
	CPUTime_MIP = 0;
	CPUTime_AutoBenders = 0;
	CPUTime_Decomposition = 0;
	CPUTime_DecompositionReuse = 0;

	SelectedMarkets = selected;
}

double SolutionInfo::GetIntercept()
{
	return -Cost;
}

double SolutionInfo::GetSlope()
{
	return Revenue + Cost;
}

void SolutionInfo::Print()
{
	cout << "Lambda: " << Lambda << " Revenue: " << Revenue << " Cost: " << Cost << " ObjVal: " << ObjVal << endl;
}