#include "Common.h"
#include "MarketSelection.h"
#include "IterativeHeuristic.h"
#include "MIP.h"
#include "Decomposition.h"
#include "FrontierCalculator.h"
#include "SolutionInfo.h"

int main(int nParams, char* params[])
{
	string inputFileName = (nParams > 1 ? params[1] : "../../Data/Results/SetPosProd_M40_T40/M40_T40_alpha2_9.txt");

	string inputFileNameOnly;
	auto fileNameBegin = inputFileName.find_last_of("/\\");
	if (fileNameBegin == string::npos)
		inputFileNameOnly = inputFileName;
	else
		inputFileNameOnly = inputFileName.substr(fileNameBegin + 1);

	auto fileNameEnd = inputFileNameOnly.find_last_of(".");
	if (fileNameEnd != string::npos)
		inputFileNameOnly = inputFileNameOnly.substr(0, fileNameEnd);

	string outputFileName = (nParams > 2 ? params[2] : inputFileNameOnly + "_result.csv");
	string summaryFileName = (nParams > 3 ? params[3] : "summary.csv");

	string parameterFileName = (nParams > 4 ? params[4] : "../Run/Parameters/MultiObj_Hybrid_MIP.txt");

	int mOverride = stoi((nParams > 5 ? params[5] : "0"));
	int tOverride = stoi((nParams > 6 ? params[6] : "0"));

	ParameterMap Parameters;
	ReadParameterMapFromFile(Parameters, parameterFileName);

	ofstream summaryFile;
	if (nParams > 3)
	{
		summaryFile.open(summaryFileName.c_str(), ios::app);
		if (!summaryFile)
			cerr << "Unable to open summary file in append mode " << summaryFileName << endl;
	}
	else
	{
		summaryFile.open(summaryFileName.c_str(), ios::trunc);
		if (!summaryFile)
			cerr << "Unable to create summary file " << summaryFileName << endl;
		else if (GetParameterValue(Parameters, "HYBRID_MULTIOBJ"))
		{
			summaryFile << "Name,Parameters,M,T,nPareto,nDiscovered,nEliminated,nImproved,"
				"FrontierCalculator_CPU,FrontierCalculator_Initial,FrontierCalculator_Pareto,"
				"MIP_CPU,MIP_Iterations,MIP_Improved,AutoBenders_CPU,AutoBenders_Iterations,AutoBenders_Improved,AutoBenders_CutCount,"
				"Decomposition_CPU,Decomposition_Iterations,Decomposition_Improved,Decomposition_CallbackCPU,Decomposition_CallCount,Decomposition_CutCount,"
				"DecompositionReuse_CPU,DecompositionReuse_Iterations,DecompositionReuse_Improved,DecompositionReuse_CallbackCPU,DecompositionReuse_CallCount,DecompositionReuse_CutCount"
				<< endl;
		}
		else
		{
			summaryFile << "Name,Parameters,M,T,Lambda,IsPareto,Heuristic_CPU,Heuristic_LB,Heuristic_Revenue,Heuristic_Cost,Heuristic_nSelected,"
				"MIP_CPU,MIP_LB,MIP_UB,MIP_Revenue,MIP_Cost,MIP_nSelected,"
				"AutoBenders_CPU,AutoBenders_LB,AutoBenders_UB,AutoBenders_Revenue,AutoBenders_Cost,AutoBenders_nSelected,AutoBenders_CutCount,"
				"Decomposition_CPU,Decomposition_LB,Decomposition_UB,Decomposition_Revenue,Decomposition_Cost,Decomposition_nFixed,Decomposition_nSelected,Decomposition_CallbackCPU,Decomposition_CallCount,Decomposition_CutCount"
				<< endl;
		}
	}

	MarketSelection MS;
	ifstream file(inputFileName);
	if (file)
	{
		cout << endl << endl << "Reading " << inputFileName << endl;
		MS.ReadData(file);
		if (mOverride > 0 && tOverride > 0)
			MS.Resize(mOverride, tOverride);
		MS.CalculateStatistics();

		cout << "Parameter " << parameterFileName << endl;

		int timeLimit = GetParameterValue(Parameters, "TIME_LIMIT");

		string paramFileName = parameterFileName;
		if (summaryFile)
		{
			auto pos = parameterFileName.find_last_of("/\\");
			if (pos != string::npos)
				paramFileName = parameterFileName.substr(pos + 1);
		}

		if (GetParameterValue(Parameters, "HYBRID_MULTIOBJ"))
		{
			double lambda = 0;
			const double Epsilon = 0.001;

			size_t nDiscovered = 0;
			size_t nEliminated = 0;
			size_t nImproved = 0;
			size_t nImproved_MIP = 0;
			size_t nImproved_AutoBenders = 0;
			size_t nImproved_Decomposition = 0;
			size_t nImproved_DecompositionReuse = 0;


			double totalTime_MIP = 0;
			double totalTime_AutoBenders = 0;
			double totalTime_Decomposition = 0;
			double totalTime_DecompositionReuse = 0;

			size_t iterations_MIP = 0;
			size_t iterations_AutoBenders = 0;
			size_t iterations_Decomposition = 0;
			size_t iterations_DecompositionReuse = 0;

			size_t bendersCutCount_AutoBenders = 0;

			double callbackTime_Decomposition = 0;
			double callbackTime_DecompositionReuse = 0;

			size_t callbackCallCount_Decomposition = 0;
			size_t callbackCallCount_DecompositionReuse = 0;

			size_t callbackCutCount_Decomposition = 0;
			size_t callbackCutCount_DecompositionReuse = 0;

			FrontierCalculator FC(MS, Parameters);
			RevenueCostSolutionMap ParetoMap = FC.CalculateRevenueCostSolutionMap();
			size_t nFC_Initial = ParetoMap.size();
			FC.EliminateDominatedSolutions(ParetoMap);
			size_t nFC_Pareto = ParetoMap.size();
			cout << "FrontierCalculator generated " << ParetoMap.size() << " Pareto efficient solutions in " << FC.GetCPUTime() << " seconds " << endl;

			Decomposition decompositionReuse(MS, Parameters);
			if (GetParameterValue(Parameters, "DECOMPOSITION_REUSE"))
			{
				decompositionReuse.SetupModel();
				decompositionReuse.AddValidInequalities(0, 0, true);
				decompositionReuse.SetObjective(lambda);

				for (auto& solution : ParetoMap)
					decompositionReuse.GenerateCut(solution.second.second);

				decompositionReuse.AddGeneratedCuts();
			}

			auto it = prev(ParetoMap.end());
			while (it != ParetoMap.begin())
			{
				double RevenueUB = it->first;
				double Cost = it->second.first;

				 //Check if the current point is (still) efficient
				if (auto itNext = next(it); itNext != ParetoMap.end())
				{
					if (Cost >= itNext->second.first)
					{
						++nEliminated;
						it = ParetoMap.erase(it);
						continue;
					}
				}

				auto itPrev = prev(it);
				double RevenueLB = itPrev->first + 1;

				// Only possible if there is a numerical issue. Skip and avoid infeasibility. 
				if (RevenueLB > RevenueUB)
				{
					--it;
					continue;
				}

				bool timeLimitReached = true;
				SolutionInfo solution(0, 0, 0, 0, IntegerVector(MS.M));

				MIP mip(MS, Parameters);
				MIP autoBenders(MS, Parameters);
				Decomposition decomposition(MS, Parameters);

				if (GetParameterValue(Parameters, "MIP"))
				{
					if (totalTime_MIP < timeLimit)
					{
						timeLimitReached = false;
						cout << "Started solving MIP" << endl;
						mip.SetupModel();
						mip.SetObjective(lambda);
						mip.SetRevenueBounds(RevenueLB, RevenueUB);
						mip.AddValidInequalities(0, 0, true);
						mip.Solve(timeLimit - totalTime_MIP);
						cout << "Finished solving MIP. UB: " << mip.GetUB() << " LB: " << mip.GetLB() << endl;

						solution.CPUTime_MIP = mip.GetCPUTime();
						if (mip.GetLB() < solution.ObjVal)
						{
							solution.ObjVal = mip.GetLB();
							solution.Cost = mip.GetCost();
							solution.Revenue = mip.GetRevenue();
							solution.SelectedMarkets = mip.GetSelectedMarkets();
						}
						if (mip.GetCost() < Cost - Epsilon)
							++nImproved_MIP;
						totalTime_MIP += mip.GetCPUTime();
						++iterations_MIP;
					}
				}

				if (GetParameterValue(Parameters, "AUTO_BENDERS"))
				{
					if (totalTime_AutoBenders < timeLimit)
					{
						timeLimitReached = false;
						cout << "Started solving AUTO_BENDERS" << endl;
						autoBenders.SetupModel();
						autoBenders.SetObjective(lambda);
						autoBenders.SetRevenueBounds(RevenueLB, RevenueUB);
						autoBenders.AddValidInequalities(0, 0, true);
						autoBenders.Solve(timeLimit - totalTime_AutoBenders, true);
						cout << "Finished solving AUTO_BENDERS. UB: " << autoBenders.GetUB() << " LB: " << autoBenders.GetLB() << endl;

						solution.CPUTime_AutoBenders = autoBenders.GetCPUTime();
						if (autoBenders.GetLB() < solution.ObjVal)
						{
							solution.ObjVal = autoBenders.GetLB();
							solution.Cost = autoBenders.GetCost();
							solution.Revenue = autoBenders.GetRevenue();
							solution.SelectedMarkets = autoBenders.GetSelectedMarkets();
						}
						if (autoBenders.GetCost() < Cost - Epsilon)
							++nImproved_AutoBenders;
						totalTime_AutoBenders += autoBenders.GetCPUTime();
						bendersCutCount_AutoBenders += autoBenders.GetBendersCutCount();
						++iterations_AutoBenders;
					}
				}

				if (GetParameterValue(Parameters, "DECOMPOSITION"))
				{
					if (totalTime_Decomposition < timeLimit)
					{
						timeLimitReached = false;
						cout << "Started solving Decomposition" << endl;
						decomposition.SetupModel();
						decomposition.SetObjective(lambda);
						decomposition.SetRevenueBounds(RevenueLB, RevenueUB);
						decomposition.AddValidInequalities(0, 0, true);
						decomposition.Solve(timeLimit - totalTime_Decomposition);
						cout << "Finished solving Decomposition. UB: " << decomposition.GetUB() << " LB: " << decomposition.GetLB() << endl;

						solution.CPUTime_Decomposition = decomposition.GetCPUTime();
						if (decomposition.GetLB() < solution.ObjVal)
						{
							solution.ObjVal = decomposition.GetLB();
							solution.Cost = decomposition.GetCost();
							solution.Revenue = decomposition.GetRevenue();
							solution.SelectedMarkets = decomposition.GetSelectedMarkets();
						}
						if (decomposition.GetCost() < Cost - Epsilon)
							++nImproved_Decomposition;

						totalTime_Decomposition += decomposition.GetCPUTime();
						callbackTime_Decomposition += decomposition.GetCallbackCPU();
						callbackCallCount_Decomposition += decomposition.GetCallCount();
						callbackCutCount_Decomposition += decomposition.GetCutCount();
						++iterations_Decomposition;
					}
				}
				
				if (GetParameterValue(Parameters, "DECOMPOSITION_REUSE"))
				{
					if (totalTime_DecompositionReuse < timeLimit)
					{
						timeLimitReached = false;
						decompositionReuse.SetRevenueBounds(RevenueLB, RevenueUB);
						decompositionReuse.AddInitialSolution(it->second.second);
						decompositionReuse.AddGeneratedCuts();
						decompositionReuse.Solve(timeLimit - totalTime_DecompositionReuse);

						solution.CPUTime_DecompositionReuse = decompositionReuse.GetCPUTime();
						if (decompositionReuse.GetLB() < solution.ObjVal)
						{
							solution.ObjVal = decompositionReuse.GetLB();
							solution.Cost = decompositionReuse.GetCost();
							solution.Revenue = decompositionReuse.GetRevenue();
							solution.SelectedMarkets = decompositionReuse.GetSelectedMarkets();
						}
						if (decompositionReuse.GetCost() < Cost - Epsilon)
							++nImproved_DecompositionReuse;

						totalTime_DecompositionReuse += decompositionReuse.GetCPUTime();
						callbackTime_DecompositionReuse += decompositionReuse.GetCallbackCPU();
						callbackCallCount_DecompositionReuse += decompositionReuse.GetCallCount();
						callbackCutCount_DecompositionReuse += decompositionReuse.GetCutCount();
						++iterations_DecompositionReuse;
					}
				}

				if (timeLimitReached)
				{
					cout << "Time limit reached for all algorithms. Skipping remaining iterations." << endl;
					break;
				}
				else if (solution.Cost < Cost - Epsilon)
				{
					if (solution.Revenue < RevenueUB)
					{
						ParetoMap[solution.Revenue] = make_pair(solution.Cost, solution.SelectedMarkets);
						++nDiscovered;
					}
					else
					{
						it->second.first = solution.Cost;
						it->second.second = solution.SelectedMarkets;
						++nImproved;

						--it;
					}
				}
				else
					--it;
			}

			cout << "Finished hybrid multiobj algorithm " << endl;

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< paramFileName << ","
					<< MS.M << ","
					<< MS.T << ","
					<< ParetoMap.size() << ","
					<< nDiscovered << ","
					<< nImproved << ","
					<< nEliminated << ","
					<< FC.GetCPUTime() << ","
					<< nFC_Initial << ","
					<< nFC_Pareto << ","
					<< totalTime_MIP << ","
					<< iterations_MIP << ","
					<< nImproved_MIP << ","
					<< totalTime_AutoBenders << ","
					<< iterations_AutoBenders << ","
					<< nImproved_AutoBenders << ","
					<< bendersCutCount_AutoBenders << ","
					<< totalTime_Decomposition << ","
					<< iterations_Decomposition << ","
					<< nImproved_Decomposition << ","
					<< callbackTime_Decomposition << ","
					<< callbackCallCount_Decomposition << ","
					<< callbackCutCount_Decomposition << ","
					<< totalTime_DecompositionReuse << ","
					<< iterations_DecompositionReuse << ","
					<< nImproved_DecompositionReuse << ","
					<< callbackTime_DecompositionReuse << ","
					<< callbackCallCount_DecompositionReuse << ","
					<< callbackCutCount_DecompositionReuse << ","
					<< endl;
			}

			ofstream outputFile(outputFileName.c_str());
			if (outputFile)
			{
				outputFile << "Revenue,Cost,SelectedMarkets" << endl;
				for (auto& Solution : ParetoMap)
				{
					outputFile << Solution.first << ","
						<< Solution.second.first<< ","
						<< to_string(Solution.second.second) 
						<< endl;
				}
			}
		}
		else
		{
			IterativeHeuristic IH(MS, Parameters, 0.5);
			IH.Solve();
			cout << "Heuristic LB: " << IH.GetLB() << endl;

			MIP mip(MS, Parameters);
			MIP autoBenders(MS, Parameters);
			Decomposition decomposition(MS, Parameters);

			if (GetParameterValue(Parameters, "MIP"))
			{
				cout << "Started solving MIP" << endl;
				mip.SetupModel();
				mip.SetObjective(0.5);
				if (GetParameterValue(Parameters, "USE_PREPROCESSING"))
					mip.Preprocess(0.5);
				mip.AddValidInequalities(IH.GetLB(), 0.5, false);
				mip.Solve(timeLimit);
				cout << "Finished solving MIP. UB: " << mip.GetUB() << " LB: " << mip.GetLB() << endl;
			}

			if (GetParameterValue(Parameters, "AUTO_BENDERS"))
			{
				cout << "Started solving AUTO_BENDERS" << endl;
				autoBenders.SetupModel();
				autoBenders.SetObjective(0.5);
				if (GetParameterValue(Parameters, "USE_PREPROCESSING"))
					autoBenders.Preprocess(0.5);
				autoBenders.AddValidInequalities(IH.GetLB(), 0.5, false);
				autoBenders.Solve(timeLimit, true);
				cout << "Finished solving AUTO_BENDERS. UB: " << autoBenders.GetUB() << " LB: " << autoBenders.GetLB() << endl;
			}

			if (GetParameterValue(Parameters, "DECOMPOSITION"))
			{
				cout << "Started solving Decomposition" << endl;
				decomposition.SetupModel();
				decomposition.SetObjective(0.5);
				if (GetParameterValue(Parameters, "USE_PREPROCESSING"))
					decomposition.Preprocess(0.5);
				decomposition.AddValidInequalities(IH.GetLB(), 0.5, false);
				decomposition.AddInitialSolution(IH.GetSelectedMarkets());
				decomposition.Solve(timeLimit);
				cout << "Finished solving Decomposition. UB: " << decomposition.GetUB() << " LB: " << decomposition.GetLB() << endl;
			}

			if (summaryFile)
			{
				summaryFile << outputFileName << ","
					<< paramFileName << ","
					<< MS.M << ","
					<< MS.T << ","
					<< 0.5 << ","
					<< 0 << ","
					<< IH.GetCPUTime() << ","
					<< IH.GetLB() << ","
					<< IH.GetRevenue() << ","
					<< IH.GetCost() << ","
					<< IH.GetSelectedMarketCount() << ","
					<< mip.GetCPUTime() << ","
					<< mip.GetLB() << ","
					<< mip.GetUB() << ","
					<< mip.GetRevenue() << ","
					<< mip.GetCost() << ","
					<< mip.GetSelectedMarketCount() << ","
					<< autoBenders.GetCPUTime() << ","
					<< autoBenders.GetLB() << ","
					<< autoBenders.GetUB() << ","
					<< autoBenders.GetRevenue() << ","
					<< autoBenders.GetCost() << ","
					<< autoBenders.GetSelectedMarketCount() << ","
					<< autoBenders.GetBendersCutCount() << ","
					<< decomposition.GetCPUTime() << ","
					<< decomposition.GetLB() << ","
					<< decomposition.GetUB() << ","
					<< decomposition.GetRevenue() << ","
					<< decomposition.GetCost() << ","
					<< decomposition.GetFixedMarketCount() << ","
					<< decomposition.GetSelectedMarketCount() << ","
					<< decomposition.GetCallbackCPU() << ","
					<< decomposition.GetCallCount() << ","
					<< decomposition.GetCutCount() << ","
					<< endl;
			}
		}

	}
	else
	{
		cout << "Input file could not be opened!" << endl;
		return -1;
	}

	return 0;
}
