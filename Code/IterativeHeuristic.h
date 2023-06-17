#ifndef ITERATIVE_HEURISTIC
#define ITERATIVE_HEURISTIC

#include "Common.h"
#include "MarketSelection.h"

class IterativeHeuristic
{
	MarketSelection& MS;
	ParameterMap& Parameters;

	SolutionMap AllM;
	IntegerVector BestM;
	IntegerVector BestT;
	double LB;
	double CPU;
	double lambda;

	void IAT2(int T, int nM, double* R, double** d, double* K, double* p, double* h, int* bestM, int* bestT, double& best, double& cntIA, int& mcntIA);
	void IterateT(int T, int nM, double* R, double** d, double* K, double* p, double* h, int* OnesM, int* OnesT, double& Prof, int& cnt);
	void GivenT(int* OnesT, int nM, double* R, int T, double** d, double* K, double* p, double* h, int* OnesM, double& Prof);
	void GivenM(int* OnesM, int nM, int T, double** d, double* K, double* p, double* h, int* OnesT, double& cost);
	void WWgeneralBW(int T, double* d, double* K1, double* p1, double* h1, double& opt, int* Ones, int& counter);

public:
	IterativeHeuristic(MarketSelection& MSIn, ParameterMap& PM, double lambdaIn);
	bool Solve();
	double GetCPUTime() { return CPU; }
	double GetLB() { return LB; }
	IntegerVector GetSelectedMarkets() { return BestM; }
	size_t GetSelectedMarketCount();
	double GetRevenue();
	double GetCost();
	const SolutionMap& GetAllSolutions() { return AllM; }
};

#endif // !ITERATIVE_HEURISTIC

