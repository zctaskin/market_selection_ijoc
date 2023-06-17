#ifndef COMMON_H
#define COMMON_H

#include <climits>
#include <cfloat>
#include <cstring>
#include <ilcplex/ilocplex.h>

#include <vector>
#include <algorithm>
#include <fstream>
#include <map>

using namespace std;

typedef map<string, IloTimer> TimerMap;
typedef map<string, double> TimePointMap;

class TimeManager;
ostream& operator << (ostream& os, const TimeManager& tm);

class TimeManager
{
	TimerMap timerMap;
	TimePointMap timePointMap;
	IloEnv env;

	TimerMap::iterator GetTimer(const string& timerCode)
	{
		TimerMap::iterator timerIt = timerMap.find(timerCode);
		if (timerIt == timerMap.end())
		{
			IloTimer timer = IloTimer(env);
			pair<TimerMap::iterator, bool> insertReturn = timerMap.insert(pair<string, IloTimer>(timerCode, timer));
			timerIt = insertReturn.first;
		}
		return timerIt;
	}

public:
	
	TimeManager(const IloEnv& envIn): env(envIn) {}

	double MarkCurrentTime(const string& timeCode)
	{
		double currentTime = env.getTime();
		if (timeCode != "")
			timePointMap[timeCode] = currentTime;
		return currentTime;
	}

	double GetMarkedTimePoint(const string& timeCode)
	{
		TimePointMap::iterator timeIt = timePointMap.find(timeCode);
		if (timeIt == timePointMap.end())
			return -1;
		else
			return timeIt->second;
	}

	void StartTimer(const string& timerCode)
	{
		TimerMap::iterator timerIt = GetTimer(timerCode);
		timerIt->second.start();
	}

	void StopTimer(const string& timerCode)
	{
		TimerMap::iterator timerIt = GetTimer(timerCode);
		timerIt->second.stop();
	}

	void ResetTimer(const string& timerCode)
	{
		TimerMap::iterator timerIt = GetTimer(timerCode);
		timerIt->second.reset();
	}

	double GetElapsedTime(const string& timerCode)
	{
		TimerMap::iterator timerIt = GetTimer(timerCode);
		return timerIt->second.getTime();
	}

	friend ostream& operator << (ostream& os, const TimeManager& tm);

};


typedef map<string, int> ParameterMap;
void ReadParameterMapFromFile(ParameterMap& map, string& fileName);
int GetParameterValue(ParameterMap& map, const char* paramName);



typedef IloArray<IloNumArray> NumArray2;
typedef IloArray<IloBoolArray> BoolArray2;
typedef IloArray<IloIntArray> IntArray2;

typedef IloArray<IntArray2> IntArray3;
typedef IloArray<NumArray2> NumArray3;
typedef IloArray<BoolArray2> BoolArray3;

typedef IloArray<IloNumVarArray> NumVarArray2;
typedef IloArray<IloBoolVarArray> BoolVarArray2;
typedef IloArray<IloIntVarArray> IntVarArray2;

typedef IloArray<NumVarArray2> NumVarArray3;
typedef IloArray<BoolVarArray2> BoolVarArray3;
typedef IloArray<IntVarArray2> IntVarArray3;

typedef IloArray<IntVarArray3> IntVarArray4;
typedef IloArray<BoolVarArray3> BoolVarArray4;

typedef IloArray<IloRangeArray> RangeArray2;

IloRangeArray CreateRangeArray(IloEnv env, IloInt m, const char* prefix, IloNum LB = -IloInfinity, IloNum UB = IloInfinity);
RangeArray2 CreateRangeArray2(IloEnv env, IloInt m, IloInt n, const char* prefix, IloNum LB = -IloInfinity, IloNum UB = IloInfinity);

NumArray2 CreateNumArray2(IloEnv env, IloInt m, IloInt n);
IntArray2 CreateIntArray2(IloEnv env, IloInt m, IloInt n);
BoolArray2 CreateBoolArray2(IloEnv env, IloInt m, IloInt n);

IloNumVarArray CreateNumVarArray(IloEnv env, IloInt m, const char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
IloIntVarArray CreateIntVarArray(IloEnv env, IloInt m, const char* prefix, IloInt LB = 0, IloInt UB = IloInfinity);
IloBoolVarArray CreateBoolVarArray(IloEnv env, IloInt m, const char* prefix);

NumVarArray2 CreateNumVarArray2(IloEnv env, IloInt m, IloInt n, const char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
IntVarArray2 CreateIntVarArray2(IloEnv env, IloInt m, IloInt n, const char* prefix, IloInt LB = 0, IloInt UB = IloInfinity);
BoolVarArray2 CreateBoolVarArray2(IloEnv env, IloInt m, IloInt n, const char* prefix);

NumVarArray3 CreateNumVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, const char* prefix, IloNum LB = 0, IloNum UB = IloInfinity);
IntVarArray3 CreateIntVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, const char* prefix, IloInt LB = 0, IloInt UB = IloInfinity);
BoolVarArray3 CreateBoolVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, const char* prefix);

IntArray3 CreateIntArray3(IloEnv env, IloInt m, IloInt n, IloInt k);
NumArray3 CreateNumArray3(IloEnv env, IloInt m, IloInt n, IloInt k);
BoolArray3 CreateBoolArray3(IloEnv env, IloInt m, IloInt n, IloInt k);

IntVarArray4 CreateIntVarArray4(IloEnv env, IloInt m, IloInt n, IloInt k, IloInt l, const char* prefix, IloInt LB = 0, IloInt UB = IloInfinity);
BoolVarArray4 CreateBoolVarArray4(IloEnv env, IloInt m, IloInt n, IloInt k, IloInt l, const char* prefix);

void SetName(IloExtractable obj, const char* prefix, int i);
void SetName2(IloExtractable obj, const char* prefix, int i, int j);
void SetName3(IloExtractable obj, const char* prefix, int i, int j, int k);
void SetName4(IloExtractable obj, const char* prefix, int i, int j, int k, int l);

IloNumArray ReadNumArrayFromFile(IloEnv env, IloInt n, ifstream& file);

IntArray2 ReadIntArray2FromFile(IloEnv env, IloInt m, IloInt n, ifstream& file);
NumArray2 ReadNumArray2FromFile(IloEnv env, IloInt m, IloInt n, ifstream& file);

#endif

