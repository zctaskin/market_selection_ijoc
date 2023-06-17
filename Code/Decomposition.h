#ifndef Decomposition_H
#define Decomposition_H

#include "Common.h"
#include "MarketSelection.h"
#include <unordered_map>
#include <list>
#include <mutex>

struct Decomposition;

/** The DecompositionCallback thread-local class. */

struct Worker
{
	Decomposition* pDecomposition;
	size_t CallCount;
	size_t CutCount;
	double CPU;

	Worker(Decomposition* pDecomposition);

	virtual bool separate(MarketSelection& MS, const IloNum thetaVal, const IloNumArray& zSol, double& OptimalCost, IloExpr& cutLhs, bool applyLifting) = 0;
};

struct WorkerWW : public Worker
{
	WorkerWW(Decomposition* pDecomposition);

	bool separate(MarketSelection& MS, const IloNum thetaVal, const IloNumArray& zSol, double& OptimalCost, IloExpr& cutLhs, bool applyLifting);
};

struct DecompositionCallback : public IloCplex::Callback::Function
{
	Decomposition* pDecomposition;
	vector<Worker*> workers;

	DecompositionCallback(Decomposition* pDecomposition);
	~DecompositionCallback();

	void invoke(const IloCplex::Callback::Context& context);
};

typedef pair<IntegerVector, double> Solution;

struct Decomposition
{
	mutex Mutex;
	unique_ptr<DecompositionCallback> pCallback;
	double BestBoundThreshold;

	const double Epsilon = 0.001;
	double CPU = 0;

	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloObjective obj;
	bool Solved = false;

	int nFixedMarkets = 0;

	IloBoolVarArray z;
	IloNumVar Cost;
	IloNumVar Revenue;

	MarketSelection& MS;
	ParameterMap& Parameters;

	IntegerVector InitialSolution;

	int IntegerSolutionLimit = INT_MAX;
	bool ReuseCuts = false;
	list <pair<IntegerVector, IloConstraintArray>> CutPool;

	list<IloConstraint> GeneratedCuts;

	Decomposition(MarketSelection& MSIn, ParameterMap& PM);
	~Decomposition();

	void SetupModel();
	void SetObjective(double lambda);
	void SetRevenueBounds(double LB, double UB);
	void SetCostThreshold(double threshold);
	void Preprocess(double lambda);
	void AddValidInequalities(double LowerBound, double lambda, bool MultiObj);
	bool AddInitialSolution(const IntegerVector& SelectedMarkets);
	double GenerateCut(const IntegerVector& SelectedMarkets);
	void AddGeneratedCuts();
	void SetReuseCuts(bool b) { ReuseCuts = b; }
	void SetIntegerSolutionLimit(int i) { IntegerSolutionLimit = i; }
	bool Solve(double timeLimit);
	double GetCPUTime() { return CPU; }
	double GetLB(); 
	double GetUB();
	double GetCost();
	double GetRevenue();
	double GetCallbackCPU();
	size_t GetCallCount();
	size_t GetCutCount();
	size_t GetFixedMarketCount() { return nFixedMarkets; }
	IntegerVector GetSelectedMarkets();
	int GetSelectedMarketCount();
};
#endif