#ifndef STATISTICAL_SUMMARY
  #define STATISTICAL_SUMMARY
#include <vector>
#include "GlobalDefinitions.h"
#include <math.h>
#include <string>
using namespace std;
class StatisticalSummary{

	vector<double> values,bestValAssociationsVal,worstValAssociationsVal;
	vector<string> bestValAssociationsDesc,worstValAssociationsDesc;
	string description;
	double bestValue,worstValue,mean;
    bool bestIsMax,bestHasChanged,worstHasChanged;
public:
  
  StatisticalSummary(string desc, bool bestIsMax);
  
  ~StatisticalSummary();
  
  double getBestValue();

  double getLastValueAdded();

  double getPreivouslyLastValueAdded();

  double getWorstValue();

  double getMean();

  double getStdDeviation();
  
  double getVariance();

  string getValueDescription();

  void addValue(double value);

  bool getBestCriteria();

  vector<string> getBestValueAssociationsDescription();

  vector<double> getBestValueAssociationsValues();

  double getBestValueAssociationsValue(int idx);

  //index start in zero
  void setBestValueAssociationDescription(int idx,string desc);

  //index start in zero
  void setBestValueAssociationsValue(int idx,double val);
  
  vector<string> getWorstValueAssociationsDescription();
  
  vector<double> getWorstValueAssociationsValues();

  double getWorstValueAssociationsValue(int idx);
  
  void setWorstValueAssociationDescription(int idx,string desc);
  
  void setWorstValueAssociationsValue(int idx,double val);
};

#endif