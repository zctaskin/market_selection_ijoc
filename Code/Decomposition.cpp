#include "Decomposition.h"
#include "MarketSelection.h"
#include <chrono>
#include <numeric>

Decomposition::Decomposition(MarketSelection& MSIn, ParameterMap& PM) : MS(MSIn), Parameters(PM), model(env), cplex(env)
{
	BestBoundThreshold = DBL_MAX;
}

Decomposition::~Decomposition()
{
	env.end();
}

void Decomposition::SetCostThreshold(double threshold)
{
	BestBoundThreshold = threshold;
}

void Decomposition::SetRevenueBounds(double LB, double UB)
{
	Revenue.setLB(LB);
	Revenue.setUB(UB);
}

void Decomposition::Preprocess(double lambda)
{
	IntegerVector MarketStates(MS.M);

	nFixedMarkets = 0;
	// Preprocessing
	{
		for (int m = 0; m < MS.M; ++m)
		{
			auto& Market = MS.Markets[m];
			if (Market.Revenue * lambda >= Market.TotalCost * (1 - lambda))
			{
				cout << "Trivial market found " << m << " Revenue: " << Market.Revenue * lambda << " TotalCost: " << Market.TotalCost * (1 - lambda) << endl;
				z[m].setBounds(1, 1);
				++nFixedMarkets;

				MarketStates[m] = 1;
			}
			else if (Market.Revenue * lambda < Market.VariableCost * (1 - lambda))
			{
				cout << "Trivial market found " << m << " Revenue: " << Market.Revenue * lambda << " VariableCost: " << Market.VariableCost * (1 - lambda) << endl;
				z[m].setBounds(0, 0);
				++nFixedMarkets;

				MarketStates[m] = -1;
			}
			else
			{
				z[m].setBounds(0, 1);

				MarketStates[m] = 0;
			}
		}
		cout << "Preprocessing fixed " << nFixedMarkets << " / " << MS.M << " markets" << endl;
	}

	if (ReuseCuts)
	{
		cplex.clearUserCuts();

		cout << "Processing cut pool: " << CutPool.size() << endl;
		for (auto& cp : CutPool)
		{
			bool isCompatible = true;
			IntegerVector& PrevMarketStates = cp.first;
			for (int m = 0; m < MS.M; ++m)
			{
				if (MarketStates[m] == 1 && PrevMarketStates[m] == -1 || MarketStates[m] == -1 && PrevMarketStates[m] == 1)
				{
					isCompatible = false;
					break;
				}
			}
			if (isCompatible)
			{
				cout << "Adding compatible cuts from cut pool: " << cp.second.getSize() << endl;
				cplex.addUserCuts(cp.second);
			}
		}

		CutPool.emplace_back(MarketStates, IloConstraintArray(env));
	}
}

void Decomposition::AddValidInequalities(double LowerBound, double lambda, bool MultiObj)
{
	int SelectedVIs = GetParameterValue(Parameters, "USE_VALID_INEQUALITIES");
	// Valid inequality: Cost >= min total cost of selected markets
	if (SelectedVIs & VALID_INEQUALITY_TOTAL_COST)
	{
		for (int m = 0; m < MS.M; ++m)
		{
			auto& Market = MS.Markets[m];
			IloExpr Expr(env);
			for (int m2 = 0; m2 < MS.M; m2++)
			{
				auto& Market2 = MS.Markets[m2];
				if (m == m2)
					Expr += MS.Markets[m].TotalCost * z[m2];
				else
					Expr += MS.Markets[m].VariableCost * z[m2];
			}
			model.add(Cost >= Expr);
			Expr.end();
		}
	}

	// theta >= cost of satisfying all demand - \sum market cost * (1-z) 
	if (SelectedVIs & VALID_INEQUALITY_FULL_SELECTION)
	{
		IloExpr Expr(env);
		Expr += MS.TotalCost;

		for (int m = 0; m < MS.M; ++m)
		{
			auto& Market = MS.Markets[m];
			Expr -= Market.TotalCost * (1 - z[m]);
		}
		model.add(Cost >= Expr);
		Expr.end();
	}

	// Valid inequality: loss of market >= sum(revenue - variable cost from others) 
	if (SelectedVIs & VALID_INEQUALITY_RECOVER_LOSS)
	{
		if (!MultiObj && !ReuseCuts)
		{
			for (int m = 0; m < MS.M; ++m)
			{
				auto& Market = MS.Markets[m];
				IloExpr Expr(env);
				for (int m2 = 0; m2 < MS.M; m2++)
				{
					auto& Market2 = MS.Markets[m2];
					if (m == m2)
						Expr += z[m2] * (Market2.Revenue * lambda - Market2.TotalCost * (1 - lambda));
					else
						Expr += z[m2] * (Market2.Revenue * lambda - Market2.VariableCost * (1 - lambda));
				}
				model.add(Expr >= LowerBound);
				Expr.end();
			}
		}
	}
}

void Decomposition::AddGeneratedCuts()
{
	if (!GeneratedCuts.empty())
	{
		cout << "Adding " << GeneratedCuts.size() << " generated cuts" << endl;
		IloConstraintArray cuts(env, GeneratedCuts.size());

		int i = 0;
		for (auto& cut : GeneratedCuts)
			cuts[i++] = cut;

		cplex.clearUserCuts();
		cplex.addUserCuts(cuts);
		GeneratedCuts.clear();
	}
}

double Decomposition::GenerateCut(const IntegerVector& SelectedMarkets)
{
	WorkerWW worker(this);
	double OptimalCost;
	IloNumArray zVal(env, MS.M);
	IloExpr cutLhs(env);

	for (int m = 0; m < MS.M; ++m)
		zVal[m] = SelectedMarkets[m];
	IloBool sepStat = worker.separate(MS, -DBL_MAX, zVal, OptimalCost, cutLhs, false);
	if (sepStat)
	{
		IloRange r(env, 0, cutLhs, IloInfinity);
		GeneratedCuts.push_back(r);
	}

	zVal.end();
	cutLhs.end();

	return OptimalCost;
}

void Decomposition::SetupModel()
{
	int nM = MS.M;
	int T = MS.T;
	double maxD = MS.MaxDemand;

	z = CreateBoolVarArray(env, nM, "z");
	Cost = IloNumVar(env, 0, IloInfinity, "cost");
	Revenue = IloNumVar(env, 0, IloInfinity, "revenue");

	Cost.setStringProperty("Type", "cost");
	Cost.setIntProperty("Market", -1);
	Revenue.setStringProperty("Type", "revenue");
	Revenue.setIntProperty("Market", -1);
	for (int m = 0; m < nM; ++m)
	{
		z[m].setStringProperty("Type", "z");
		z[m].setIntProperty("Market", m);
	}

	model.add(z);
	model.add(Cost);
	model.add(Revenue);

	IloExpr revenueObj(env);
	// Revenue
	for (int m = 0; m < nM; m++)
		revenueObj += MS.Markets[m].Revenue * z[m];

	model.add(Revenue == revenueObj);
	revenueObj.end();

	obj = IloMaximize(env);
	model.add(obj);

	cplex = IloCplex(model);
}

void Decomposition::SetObjective(double lambda)
{
	obj.setLinearCoef(Revenue, lambda);
	obj.setLinearCoef(Cost, lambda - 1);
}

bool Decomposition::AddInitialSolution(const IntegerVector& SelectedMarkets)
{
	InitialSolution = SelectedMarkets;
	return true;
}

bool Decomposition::Solve(double timeLimit)
{
	cplex.setParam(IloCplex::ClockType, 2);
	auto startTime = chrono::high_resolution_clock::now();

	if (timeLimit > 0)
		cplex.setParam(IloCplex::TiLim, timeLimit);

	cplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));

	try
	{
		pCallback = make_unique<DecompositionCallback>(this);
		CPXLONG contextmask = IloCplex::Callback::Context::Id::Candidate
			| IloCplex::Callback::Context::Id::GlobalProgress;
		cplex.use(pCallback.get(), contextmask);

		// Initial solution
		if (GetParameterValue(Parameters, "USE_INITIAL_SOLUTION") && !InitialSolution.empty())
		{
			IloNumVarArray vars(env);
			IloNumArray values(env);
			for (int m = 0; m < MS.M; ++m)
			{
				vars.add(z[m]);
				values.add(InitialSolution[m]);
			}

			if (cplex.getNMIPStarts())
				cplex.deleteMIPStarts(0);
			cplex.addMIPStart(vars, values, IloCplex::MIPStartEffort::MIPStartSolveFixed);

			vars.end();
			values.end();
		}

		cplex.setParam(IloCplex::IntSolLim, IntegerSolutionLimit);
		Solved = cplex.solve();
	}
	catch (IloException & ex)
	{
		cout << ex.getMessage() << endl;
	}

	CPU = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

	return Solved;
}

double Decomposition::GetLB()
{
	return Solved ? cplex.getObjValue() : 0;
}

double Decomposition::GetUB()
{
	return Solved ? cplex.getBestObjValue() : DBL_MAX;
}

double Decomposition::GetCost()
{
	if (!Solved)
		return 0;

	IntegerVector Selected = GetSelectedMarkets();
	return MS.CalculateCost(Selected);
}

double Decomposition::GetRevenue()
{
	if (!Solved)
		return 0;

	IntegerVector Selected = GetSelectedMarkets();
	return MS.CalculateRevenue(Selected);
}

IntegerVector Decomposition::GetSelectedMarkets()
{
	IntegerVector Result(MS.M, 0);
	if (Solved)
	{
		IloNumArray zVal(env);
		cplex.getValues(z, zVal);
		for (int m = 0; m < MS.M; ++m)
			Result[m] = int(zVal[m] + Epsilon);
		zVal.end();
	}
	return Result;
}

int Decomposition::GetSelectedMarketCount()
{
	if (!Solved)
		return -1;

	int nSelected = 0;
	for (int m = 0; m < MS.M; ++m)
	{
		double Selected = cplex.getValue(z[m]);
		if (Selected > Epsilon)
			++nSelected;
	}
	return nSelected;
}

double Decomposition::GetCallbackCPU()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CPU;
	return Total;
}

size_t Decomposition::GetCallCount()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CallCount;
	return Total;
}

size_t Decomposition::GetCutCount()
{
	if (!Solved)
		return 0;

	double Total = 0;
	for (auto& Worker : pCallback->workers)
		Total += Worker->CutCount;
	return Total;
}

DecompositionCallback::DecompositionCallback(Decomposition* pB) : pDecomposition(pB)
{
	int nThreads = GetParameterValue(pDecomposition->Parameters, "THREAD_COUNT");
	workers.resize(nThreads);

	for (int i = 0; i < nThreads; ++i)
		workers[i] = new WorkerWW(pDecomposition);
}

DecompositionCallback::~DecompositionCallback()
{}

void DecompositionCallback::invoke(const IloCplex::Callback::Context& context)
{
	auto startTime = chrono::high_resolution_clock::now();

	int const threadNo = context.getIntInfo(IloCplex::Callback::Context::Info::ThreadId);

	// Get the right worker
	auto& worker = workers[threadNo];

	IloEnv env = context.getEnv();
	IloNumArray zVal(env, pDecomposition->MS.M);
	IloNum thetaVal;

	// Get the current z solution
	switch (context.getId()) {
	case IloCplex::Callback::Context::Id::Candidate:
		if (!context.isCandidatePoint()) // The model is always bounded
			throw IloCplex::Exception(-1, "Unbounded solution");
		context.getCandidatePoint(pDecomposition->z, zVal);
		thetaVal = context.getCandidatePoint(pDecomposition->Cost);
		break;
	case IloCplex::Callback::Context::Id::Relaxation:
		context.getRelaxationPoint(pDecomposition->z, zVal);
		thetaVal = context.getRelaxationPoint(pDecomposition->Cost);
		break;
	case IloCplex::Callback::Context::Id::GlobalProgress:
	{
		double BestBound = context.getDoubleInfo(IloCplex::Callback::Context::Info::BestBound);
		if (abs(BestBound) >= pDecomposition->BestBoundThreshold)
		{
			cout << "Aborting callback. BestBound: " << BestBound << " Threshold: " << pDecomposition->BestBoundThreshold << endl;
			context.abort();
		}
		return;
	}
	default:
		// Free memory
		zVal.end();
		throw IloCplex::Exception(-1, "Unexpected contextID");
	}

	IloExpr cutLhs(env);
	double OptimalCost;
	bool applyLifting = context.isCandidatePoint();
	IloBool sepStat = worker->separate(pDecomposition->MS, thetaVal, zVal, OptimalCost, cutLhs, applyLifting);

	if (context.getId() == IloCplex::Callback::Context::Id::Candidate)
	{
		IloNum CandidateObjective = context.getCandidateObjective();
		IloNum IncumbentObjective = context.getIncumbentObjective();
		IloNum ActualCandidateObjective = CandidateObjective + thetaVal - OptimalCost;
	}

	zVal.end();

	if (sepStat) {
		// Add the cut
		IloRange r(env, 0, cutLhs, IloInfinity);

		switch (context.getId()) {
		case IloCplex::Callback::Context::Id::Candidate:
			context.rejectCandidate(r);
			break;
		case IloCplex::Callback::Context::Id::Relaxation:
			if (OptimalCost - thetaVal > 1)
			{
				cout << r << endl;
				context.addUserCut(r,
					IloCplex::UseCutPurge,
					IloFalse);
			}
			break;
		default:
			r.end();
			throw IloCplex::Exception(-1, "Unexpected contextID");
		}

		{
			unique_lock lock(pDecomposition->Mutex);
			pDecomposition->GeneratedCuts.push_back(r);
		}
	}
	++worker->CallCount;
	if (sepStat)
		++worker->CutCount;
	worker->CPU += chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;
}

Worker::Worker(Decomposition* pB) : pDecomposition(pB)
{
	CallCount = 0;
	CutCount = 0;
	CPU = 0;
}

WorkerWW::WorkerWW(Decomposition* pB) : Worker(pB)
{
}

bool WorkerWW::separate(MarketSelection& MS, const IloNum thetaVal, const IloNumArray& zVal, double& OptimalCost, IloExpr& cutLhs, bool applyLifting)
{
	//int nSelected = 0;
	DoubleVector Demand(MS.T, 0);
	double TotalIndividualCost = 0;
	double TotalConstantH = 0;
	vector<bool> IsSelected(MS.M);
	for (int m = 0; m < MS.M; ++m)
	{
		auto& Market = MS.Markets[m];

		for (int t = 0; t < MS.T; ++t)
			Demand[t] += Market.Demand[t] * zVal[m];
	}

	IntegerVector Ones(MS.T);
	DoubleVector Duals(MS.T, 0);
	double ConstantH;

	int counter;

	MS.WWgeneralPD(MS.T, Demand.data(), MS.SetupCost.data(), MS.ProductionCost.data(), MS.HoldingCost.data(), OptimalCost, Ones.data(), Duals.data(), ConstantH, counter);

	for (auto d : Duals)
	{
		if (isnan(d) || isnan(ConstantH))
		{
			cout << "Numerical trouble detected" << endl;
			MS.WWgeneralPD(MS.T, Demand.data(), MS.SetupCost.data(), MS.ProductionCost.data(), MS.HoldingCost.data(), OptimalCost, Ones.data(), Duals.data(), ConstantH, counter);
			return false;
		}
	}

	double check = 0;
	if (OptimalCost > thetaVal + pDecomposition->Epsilon)
	{
		cutLhs.clear();
		cutLhs += pDecomposition->Cost;

		for (int m = 0; m < MS.M; ++m)
		{
			auto& Market = MS.Markets[m];
			if (zVal[m] > pDecomposition->Epsilon)
			{
				double Coefficient = 0;
				for (int t = 0; t < MS.T; ++t)
					Coefficient += Duals[t] * Market.Demand[t] * zVal[m];

				cutLhs -= pDecomposition->z[m] * Coefficient;
				check += (zVal[m] * Coefficient);
			}
			else if (applyLifting)
			{
				// Apply lifting
				cutLhs -= pDecomposition->z[m] * (Market.VariableCost);
			}
		}

		return true;
	}
	return false;
}
