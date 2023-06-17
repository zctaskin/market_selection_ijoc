#ifndef MIP_H
#define MIP_H

#include "Common.h"
#include "MarketSelection.h"

class MIP
{
	MarketSelection& MS;
	ParameterMap& Parameters;

	const double Epsilon = 0.0001;
	double CPU = 0;

	IloEnv env;
	IloModel model;
	IloObjective obj;
	IloCplex cplex;
	bool Solved = false;

	int nFixedMarkets = 0;

	IloNumVarArray z;
	IloNumVarArray y;
	NumVarArray3 x;
	IloNumVar Cost;
	IloNumVar Revenue;
	//IloNumVar nontrivial;

public:
	MIP(MarketSelection& MSIn, ParameterMap& PM);
	~MIP();
	void SetupModel();
	void SetObjective(double lambda);
	void SetRevenueBounds(double LB, double UB);
	void Preprocess(double lambda);
	void AddValidInequalities(double LowerBound, double lambda, bool MultiObj);
	bool Solve(double timeLimit, bool autoBenders = false);
	double GetCPUTime() { return CPU; }
	double GetLB(); 
	double GetUB();
	IntegerVector GetSelectedMarkets(); 
	int GetSelectedMarketCount();
	double GetCost();
	double GetRevenue();
	int GetBendersCutCount();
};
#endif