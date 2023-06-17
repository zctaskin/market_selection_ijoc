#pragma once

#include "Common.h"
#include "MarketSelection.h"
#include <chrono>

struct Decomposition;

class FrontierCalculator
{
	MarketSelection& MS;
	ParameterMap& Parameters;

	double CPU = 0;

public:

	FrontierCalculator(MarketSelection& MSIn, ParameterMap& PM);

	RevenueCostSolutionMap CalculateRevenueCostSolutionMap(Decomposition* pDecomposition = nullptr);
	size_t EliminateDominatedSolutions(RevenueCostSolutionMap& All);

	double GetCPUTime() { return CPU; }

};
