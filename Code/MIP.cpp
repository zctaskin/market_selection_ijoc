#include "MIP.h"
#include "MarketSelection.h"
#include <chrono>

MIP::MIP(MarketSelection& MSIn, ParameterMap& PM) : MS(MSIn), Parameters(PM), model(env), cplex(env)
{
}

MIP::~MIP()
{
	env.end();
}


void MIP::SetupModel()
{
	int nM = MS.M;
	int T = MS.T;
	double maxD = MS.MaxDemand;

	x = CreateNumVarArray3(env, nM, T, T, "x", 0, maxD);
	y = CreateNumVarArray(env, T, "y", 0, 1);
	z = CreateNumVarArray(env, nM, "z", 0, 1);
	Revenue = IloNumVar(env, 0, IloInfinity, "revenue");
	Cost = IloNumVar(env, 0, IloInfinity, "cost");

	model.add(IloConversion(env, z, ILOBOOL));

	for (int m = 0; m < nM; m++)
	{
		for (int j = 0; j < T; j++)
		{
			double Demand = MS.Markets[m].Demand[j];
			if (Demand > Epsilon)
			{
				IloExpr v(env);
				for (int i = 0; i < j + 1; ++i)
					v += x[m][i][j];
				model.add(v == Demand * z[m]);
				v.end();
			}
		}
	}

	for (int m = 0; m < nM; m++)
	{
		for (int j = 0; j < T; j++)
		{
			double Demand = MS.Markets[m].Demand[j];
			if (Demand > Epsilon)
			{
				for (int i = 0; i < j + 1; ++i)
					model.add(x[m][i][j] <= Demand * y[i]);
			}
		}
	}

	// Objective function terms
	IloExpr costObj(env);

	// Setup cost		
	for (int i = 0; i < T; i++)
		costObj += MS.SetupCost[i] * y[i];

	// Production + holding cost
	for (int m = 0; m < nM; m++)
	{
		for (int j = 0; j < T; j++)
		{
			double Demand = MS.Markets[m].Demand[j];
			if (Demand > Epsilon)
			{
				for (int i = 0; i < j + 1; i++)
				{
					double CumulativeCost = MS.GetCumulativeCost(i, j);
					costObj += CumulativeCost * x[m][i][j];
				}
			}
		}
	}
	model.add(Cost == costObj);
	costObj.end();

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

void MIP::SetObjective(double lambda)
{
	obj.setLinearCoef(Revenue, lambda);
	obj.setLinearCoef(Cost, lambda - 1);
}

void MIP::SetRevenueBounds(double LB, double UB)
{
	Revenue.setLB(LB);
	Revenue.setUB(UB);
}

void MIP::Preprocess(double lambda)
{
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
			}
			else if (Market.Revenue * lambda < Market.VariableCost * (1 - lambda))
			{
				cout << "Trivial market found " << m << " Revenue: " << Market.Revenue * lambda << " VariableCost: " << Market.VariableCost * (1 - lambda) << endl;
				z[m].setBounds(0, 0);
				++nFixedMarkets;
			}
			else
			{
				z[m].setBounds(0, 1);
			}
		}
		cout << "Preprocessing fixed " << nFixedMarkets << " / " << MS.M << " markets" << endl;
	}
}

void MIP::AddValidInequalities(double LowerBound, double lambda, bool MultiObj)
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
		if (!MultiObj)
		{
			IloNumArray Revenue(env);
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

bool MIP::Solve(double timeLimit, bool autoBenders)
{
	cplex.setParam(IloCplex::ClockType, 2);
	auto startTime = chrono::high_resolution_clock::now();

	if (timeLimit > 0)
		cplex.setParam(IloCplex::TiLim, timeLimit);

	cplex.setParam(IloCplex::Threads, GetParameterValue(Parameters, "THREAD_COUNT"));
	if (autoBenders)
		cplex.setParam(IloCplex::Param::Benders::Strategy, 3);

	try
	{
		cplex.exportModel("MIP.lp");
		Solved = cplex.solve();
	}
	catch (IloException & ex)
	{
		cout << ex.getMessage() << endl;
	}

	CPU = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - startTime).count() / 1000.0;

	return Solved;
}

double MIP::GetLB()
{
	return Solved ? cplex.getObjValue() : DBL_MAX;
}

double MIP::GetUB()
{
	return Solved ? cplex.getBestObjValue() : DBL_MAX;
}

IntegerVector MIP::GetSelectedMarkets()
{
	IntegerVector Result(MS.M, 0);
	if (Solved)
	{
		IloNumArray zVal(env);
		cplex.getValues(z, zVal);
		for (int m = 0; m < MS.M; ++m)
			Result[m] = zVal[m];
		zVal.end();
	}
	return Result;
}

int MIP::GetSelectedMarketCount()
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

double MIP::GetCost()
{
	return Solved ? cplex.getValue(Cost) : DBL_MAX;
}

double MIP::GetRevenue()
{
	return Solved ? cplex.getValue(Revenue) : DBL_MAX;
}

int MIP::GetBendersCutCount()
{
	return cplex.getNcuts(IloCplex::CutType::CutBenders);
}