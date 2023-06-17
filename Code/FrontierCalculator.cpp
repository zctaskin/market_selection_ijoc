#include "FrontierCalculator.h"
#include "Decomposition.h"
#include <iostream>

using namespace std;

FrontierCalculator::FrontierCalculator(MarketSelection& MSIn, ParameterMap& PM) : MS(MSIn), Parameters(PM)
{
}

// Revenues are exact, costs are upper bounds
RevenueCostSolutionMap FrontierCalculator::CalculateRevenueCostSolutionMap(Decomposition* pDecomposition)
{
	auto startTime = chrono::high_resolution_clock::now();

	cout << "Calculating possible revenues" << endl;
	RevenueCostSolutionMap All; 	// [Revenue -> (Cost, Selection)]
	All[0] = make_pair(0, IntegerVector(MS.M));

	for (int m = 0; m < MS.M; ++m)
	{
		double Revenue = MS.Markets[m].Revenue;
		auto Current = All;
		for (auto& [d, p] : All)
		{
			IntegerVector Selected = p.second;
			Selected[m] = 1;

			double Cost = pDecomposition ? pDecomposition->GenerateCut(Selected) : MS.CalculateCost(Selected);

			auto it = Current.find(d + Revenue);
			if (it == Current.end())
				Current[d + Revenue] = make_pair(Cost, Selected);
			else if (Cost < it->second.first)
				it->second = make_pair(Cost, Selected);
		}
		All.swap(Current);
		cout << "m " << m << " Number of possible revenues : " << All.size() << endl;
	}
	cout << "Number of possible revenues: " << All.size() << endl;

	CPU += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

	return All;
}

size_t FrontierCalculator::EliminateDominatedSolutions(RevenueCostSolutionMap& All)
{
	auto startTime = chrono::high_resolution_clock::now();

	size_t OriginalSize = All.size();

	cout << "Eliminating dominated solutions. OriginalSize: " << OriginalSize << endl;

	if (!All.empty())
	{
		auto itNext = prev(All.end());
		while (true)
		{
			if (itNext == All.begin())
				break;
			auto itCurrent = prev(itNext);
			if (itCurrent->second.first < itNext->second.first)
				itNext = itCurrent;
			else
				All.erase(itCurrent);
		}
	}

	cout << "Finished eliminating dominated solutions. Number of Pareto efficient solutions: " << All.size() << endl;

	CPU += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

	return OriginalSize - All.size();
}
