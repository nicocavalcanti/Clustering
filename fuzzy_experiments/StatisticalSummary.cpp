#include "StatisticalSummary.h"

StatisticalSummary::StatisticalSummary(string desc, bool bestIsMax):description(desc),bestIsMax(bestIsMax){
	  this->mean=0;
	  double best = (bestIsMax?MIN_DOUBLE:MAX_DOUBLE);
      this->bestValue = best;
      double worst = (bestIsMax?MAX_DOUBLE:MIN_DOUBLE);
	  this->worstValue = worst;
	  this->bestHasChanged = false;		  
	  this->worstHasChanged = false;
  };
  
  StatisticalSummary::~StatisticalSummary(){
  };
  
  double StatisticalSummary::getBestValue(){
    return this->bestValue;
  };
  double StatisticalSummary::getWorstValue(){
    return this->worstValue;
  };

  double StatisticalSummary::getMean(){
    return this->mean/this->values.size();
  };

  bool StatisticalSummary::getBestCriteria(){
	  return this->bestIsMax;
  };

  double StatisticalSummary::getStdDeviation(){
	 double variance = this->getVariance();
	 double standardDeviation = pow(variance,0.5);
	 return standardDeviation;
  };
  
  double StatisticalSummary::getVariance(){
    double variance=0.0;
	double mean = this->getMean();
    for (int i=0;i<this->values.size();i++)
  	  variance += pow(fabs(this->values[i]-mean),2);
    variance /= this->values.size()-1;  
    return variance;
  }

  string StatisticalSummary::getValueDescription(){
    return this->description;
  };

  void StatisticalSummary::addValue(double value){
    this->values.push_back(value);
	this->bestHasChanged = false;		  
	this->worstHasChanged = false;
	if (bestIsMax){
		if (value>this->bestValue){
			this->bestValue = value;
			this->bestHasChanged = true;		  
		}
		if (value<this->worstValue){
			this->worstValue = value;
			this->worstHasChanged = true;		  
		}
	}
	else {
		if (value<this->bestValue){
			this->bestValue = value;
			this->bestHasChanged = true;		  
		}
		if (value>this->worstValue){
			this->worstValue = value;
			this->worstHasChanged = true;		  
		}
	}
	this->mean += value;
  };
  //pre-conditions: at least one value has to be previously added
  double StatisticalSummary::getLastValueAdded(){
	  return this->values[this->values.size()-1]; 
  };

  //pre-conditions: at least two values have to be previously added
  double StatisticalSummary::getPreivouslyLastValueAdded(){
	  if ((long)this->values.size()-2>=0)
	    return this->values[this->values.size()-2];
	  return this->worstValue;
  };

  vector<string> StatisticalSummary::getBestValueAssociationsDescription(){
    return this->bestValAssociationsDesc;		
  };

  vector<double> StatisticalSummary::getBestValueAssociationsValues(){
	  return this->bestValAssociationsVal;
  };

  double StatisticalSummary::getBestValueAssociationsValue(int idx){
	  if (idx>((long)this->bestValAssociationsVal.size()-1) || (idx<0))
		  return INDF;
	  return this->bestValAssociationsVal[idx];
  };

  //index start in zero
  void StatisticalSummary::setBestValueAssociationDescription(int idx,string desc){
      if (this->bestValAssociationsDesc.size()<=idx+1)
		  this->bestValAssociationsDesc.resize(idx+1);
	  this->bestValAssociationsDesc[idx]=desc;
  };

  //index start in zero
  void StatisticalSummary::setBestValueAssociationsValue(int idx,double val){
	  if (!this->bestHasChanged)
		  return;
	  if (this->bestValAssociationsVal.size()<=idx+1)
		  this->bestValAssociationsVal.resize(idx+1);
	  this->bestValAssociationsVal[idx]=val;
  };
  
  vector<string> StatisticalSummary::getWorstValueAssociationsDescription(){
    return this->worstValAssociationsDesc;
  };
  
  vector<double> StatisticalSummary::getWorstValueAssociationsValues(){
    return this->worstValAssociationsVal;
  };

  double StatisticalSummary::getWorstValueAssociationsValue(int idx){
    if (idx>((long)this->worstValAssociationsVal.size()-1) || (idx<0))
		  return INDF;
	return this->worstValAssociationsVal[idx];
  }
  
  void StatisticalSummary::setWorstValueAssociationDescription(int idx,string desc){
      if (this->worstValAssociationsDesc.size()<=idx+1)
		  this->worstValAssociationsDesc.resize(idx+1);
	  this->worstValAssociationsDesc[idx]=desc;
  };
  
  void StatisticalSummary::setWorstValueAssociationsValue(int idx,double val){
	  if (!this->worstHasChanged)
		  return;
	  if (this->worstValAssociationsVal.size()<=idx+1)
		  this->worstValAssociationsVal.resize(idx+1);
	  this->worstValAssociationsVal[idx]=val;
  };