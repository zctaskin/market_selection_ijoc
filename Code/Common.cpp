#include "Common.h"
#include <list>

ILOSTLBEGIN

void SetName(IloExtractable obj, const char* prefix, int i)
{
	string name = prefix;
	name += '_' + to_string(i);
	obj.setName(name.c_str());
}

void SetName2(IloExtractable obj, const char* prefix, int i, int j)
{
	string name = prefix;
	name += '_' + to_string(i);
	name += '_' + to_string(j);
	obj.setName(name.c_str());
}

void SetName3(IloExtractable obj, const char* prefix, int i, int j, int k)
{
	string name = prefix;
	name += '_' + to_string(i);
	name += '_' + to_string(j);
	name += '_' + to_string(k);
	obj.setName(name.c_str());

}

void SetName4(IloExtractable obj, const char* prefix, int i, int j, int k, int l)
{
	string name = prefix;
	name += '_' + to_string(i);
	name += '_' + to_string(j);
	name += '_' + to_string(k);
	name += '_' + to_string(l);
	obj.setName(name.c_str());
}

IloBoolVarArray CreateBoolVarArray(IloEnv env, IloInt m, const char* prefix)
{
	IloBoolVarArray array = IloBoolVarArray(env, m);
	for (int i = 0; i < m; i++)
		SetName(array[i], prefix, i);
	return array;
}

IloNumVarArray CreateNumVarArray(IloEnv env, IloInt m, const char* prefix, IloNum LB, IloNum UB)
{
	IloNumVarArray array = IloNumVarArray(env, m, LB, UB);
	for (int i = 0; i < m; i++)
		SetName(array[i], prefix, i);
	return array;
}

IloIntVarArray CreateIntVarArray(IloEnv env, IloInt m, const char* prefix, IloInt LB, IloInt UB)
{
	IloIntVarArray array = IloIntVarArray(env, m, LB, UB);
	for (int i = 0; i < m; i++)
		SetName(array[i], prefix, i);
	return array;
}

NumArray2 CreateNumArray2(IloEnv env, IloInt m, IloInt n)
{
	NumArray2 array =  NumArray2(env, m);
	for (int i = 0; i < m; i++)
		array[i] = IloNumArray(env, n);
	return array;
}

IntArray2 CreateIntArray2(IloEnv env, IloInt m, IloInt n)
{
	IntArray2 array =  IntArray2(env, m);
	for (int i = 0; i < m; i++)
		array[i] = IloIntArray(env, n);
	return array;
}

BoolArray2 CreateBoolArray2(IloEnv env, IloInt m, IloInt n)
{
	BoolArray2 array =  BoolArray2(env, m);
	for (int i = 0; i < m; i++)
		array[i] = IloBoolArray(env, n);
	return array;
}

IntArray2 ReadIntArray2FromFile(IloEnv env, IloInt m, IloInt n, ifstream& file)
{
	IntArray2 array = CreateIntArray2(env, m, n);
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			file >> array[i][j];
	return array;
}

NumArray2 ReadNumArray2FromFile(IloEnv env, IloInt m, IloInt n, ifstream& file)
{
  NumArray2 array = CreateNumArray2(env, m, n);
  for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++)
      file >> array[i][j];
  return array;
}

IloNumArray ReadNumArrayFromFile(IloEnv env, IloInt n, ifstream& file)
{
  IloNumArray array = IloNumArray(env, n);
    for (int j = 0; j < n; j++)
      file >> array[j];
  return array;
}

NumVarArray2 CreateNumVarArray2(IloEnv env, IloInt m, IloInt n, const char* prefix, IloNum LB, IloNum UB)
{
	NumVarArray2 array =  NumVarArray2(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = IloNumVarArray(env, n, LB, UB);
		for (int j = 0; j < n; j++)
			SetName2(array[i][j], prefix, i, j);
	}
	return array;
}

IntVarArray2 CreateIntVarArray2(IloEnv env, IloInt m, IloInt n, const char* prefix, IloInt LB, IloInt UB)
{
	IntVarArray2 array =  IntVarArray2(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = IloIntVarArray(env, n, LB, UB);
		for (int j = 0; j < n; j++)
			SetName2(array[i][j], prefix, i, j);
	}
	return array;
}

IloRangeArray CreateRangeArray(IloEnv env, IloInt m, const char* prefix, IloNum LB, IloNum UB)
{
	IloRangeArray array = IloRangeArray(env, m, LB, UB);
	for (int i = 0; i < m; i++)
		SetName(array[i], prefix, i);
	return array;
}
		
RangeArray2 CreateRangeArray2(IloEnv env, IloInt m, IloInt n, const char* prefix, IloNum LB, IloNum UB)
{
	RangeArray2 array = RangeArray2(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = IloRangeArray(env, n, LB, UB);
		for (int j = 0; j < n; j++)
			SetName2(array[i][j], prefix, i, j);
	}
	return array;
}

BoolVarArray2 CreateBoolVarArray2(IloEnv env, IloInt m, IloInt n, const char* prefix)
{
	BoolVarArray2 array =  BoolVarArray2(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = IloBoolVarArray(env, n);
		for (int j = 0; j < n; j++)
			SetName2(array[i][j], prefix, i, j);
	}
	return array;
}

NumVarArray3 CreateNumVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, const char* prefix, IloNum LB, IloNum UB)
{
	NumVarArray3 array = NumVarArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateNumVarArray2(env, n, k, "", LB, UB);
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				SetName3(array[i][j][c], prefix, i, j, c);
	}
	return array;
}

NumArray3 CreateNumArray3(IloEnv env, IloInt m, IloInt n, IloInt k)
{
	NumArray3 array = NumArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateNumArray2(env, n, k);
	}
	return array;
}

BoolArray3 CreateBoolArray3(IloEnv env, IloInt m, IloInt n, IloInt k)
{
	BoolArray3 array = BoolArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateBoolArray2(env, n, k);
	}
	return array;
}

IntArray3 CreateIntArray3(IloEnv env, IloInt m, IloInt n, IloInt k)
{
	IntArray3 array = IntArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateIntArray2(env, n, k);
	}
	return array;
}

BoolVarArray3 CreateBoolVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, const char* prefix)
{
	BoolVarArray3 array = BoolVarArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateBoolVarArray2(env, n, k, "");
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				SetName3(array[i][j][c], prefix, i, j, c);
	}
	return array;
}

BoolVarArray4 CreateBoolVarArray4(IloEnv env, IloBool m, IloBool n, IloBool k, IloBool l, const char* prefix)
{
	BoolVarArray4 array = BoolVarArray4(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateBoolVarArray3(env, n, k, l, "");
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				for (int d = 0; d < l; d++)
					SetName4(array[i][j][c][d], prefix, i, j, c, d);
	}
	return array;
}

IntVarArray3 CreateIntVarArray3(IloEnv env, IloInt m, IloInt n, IloInt k, const char* prefix, IloInt LB, IloInt UB)
{
	IntVarArray3 array = IntVarArray3(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateIntVarArray2(env, n, k, "", LB, UB);
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				SetName3(array[i][j][c], prefix, i, j, c);
	}
	return array;
}

IntVarArray4 CreateIntVarArray4(IloEnv env, IloInt m, IloInt n, IloInt k, IloInt l, const char* prefix, IloInt LB, IloInt UB)
{
	IntVarArray4 array = IntVarArray4(env, m);
	for (int i = 0; i < m; i++)
	{
		array[i] = CreateIntVarArray3(env, n, k, l, "", LB, UB);
		for (int j = 0; j < n; j++)
			for (int c = 0; c < k; c++)
				for (int d = 0; d < l; d++)
					SetName4(array[i][j][c][d], prefix, i, j, c, d);
	}
	return array;
}

void ReadParameterMapFromFile(ParameterMap& map, string& fileName)
{
	ifstream file(fileName.c_str());
	if(file)
	{
		string name;
		int value;
		char dummy[1000];
		while(!file.eof())
		{
			char ch = file.peek();
			if (ch == '/' || ch == '\n' || ch == '\r')
			{
				// This is a comment line, just disregard it
				file.getline(dummy, 1000);
				continue;
			}
			else
			{
				file >> name >> value;
				map[name] = value;
				// Read the rest of the line as comments
				file.getline(dummy, 1000);
			}
		}
	}
	else
	{
		cout << "Input parameter file could not be opened! " << fileName << endl;
		exit(-1);
	}
}

int GetParameterValue(ParameterMap& map, const char* paramName)
{
	string paramNameStr(paramName);
	auto it = map.find(paramNameStr);
	if (it != map.end())
		return it->second;
	else
	{
		cout << "Unknown parameter " << paramNameStr << endl;
		throw paramNameStr;
	}

	return -1;	// Can never come here	
}

bool CompareMarkedTimePoints(const pair<string, double>& tp1, const pair<string, double>& tp2)
{
	return tp1.second < tp2.second;
}

typedef list <pair<string, double> > TimePointList;

ostream& operator << (ostream& os, const TimeManager& tm)
{
	if (tm.timePointMap.size())
	{
		TimePointList timePointList;
		for (TimePointMap::const_iterator timeIt = tm.timePointMap.begin(); timeIt != tm.timePointMap.end(); ++timeIt)
			timePointList.push_back(pair<string, double>(timeIt->first, timeIt->second));
		timePointList.sort(CompareMarkedTimePoints);
		
		os << "Marked time points: " << endl;
		for (TimePointList::const_iterator timeIt = timePointList.begin(); timeIt != timePointList.end(); ++timeIt)
			os << timeIt->first << " : " << timeIt->second << endl;
	}
	else
		os << "No marked time points available." << endl; 

	if (tm.timerMap.size())
	{
		os << endl << "Elapsed time for timers " << endl;
		for (TimerMap::const_iterator timerIt = tm.timerMap.begin(); timerIt != tm.timerMap.end(); ++timerIt)
			os << timerIt->first << " : " << timerIt->second.getTime() << endl;
	}
	else
		os << "No timers available." << endl; 

	os << "Total elapsed time : " << tm.env.getTime();
	return os;
}
