#pragma once

#include "GlobalDefinitions.h"
#include <vector>
using namespace std;

class AdjustedRandIndex
{
  vector<double> &partitionVariables;
  int nClasses;
  typedef vector<int> contingencyTableElement;
  vector<contingencyTableElement> contingencyTable;
  double mChosenk(int m,int k=2);
  long factorial(int m);
  void makeContingencyTable(int clusters,vector<int> ** &algorithmPartition);
  long getContingencyTableElement(int i,int j);
  double getA();  
public:
	AdjustedRandIndex(vector<double> &partitionVariables,int nClasses);	
	double getAdjustedRandIndex(int clusters,vector<int> **&algorithmPartition);
	~AdjustedRandIndex(void);
};
