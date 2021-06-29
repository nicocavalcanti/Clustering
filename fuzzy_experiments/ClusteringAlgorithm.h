// ClusteringAlgorithm.h: interface for the ClusteringAlgorithm class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ClusteringAlgorithm_H__62839795_566F_44CD_9EC4_9C45C8D09889__INCLUDED_)
#define AFX_ClusteringAlgorithm_H__62839795_566F_44CD_9EC4_9C45C8D09889__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <iostream>
#include <math.h>
#include "GlobalDefinitions.h"
#include "AdjustedRandIndex.h"
#include "StatisticalSummary.h"
using namespace std;

#define EPSLON 0.0000001


#ifdef DEBUGTIME
  #include <sys/timeb.h>
  #include <time.h>
#endif


class Partitioner;
//  In those methods which ask code the convencion is the following
//   1 = FuzzyCMeans
//   2 = FuzzyCMeansAdaptive
//This class provides a bridge to the classes which make clustering according to some particular algorithm

class ClusteringAlgorithm  {
public:
	ClusteringAlgorithm();
	static Partitioner * getAlgorithm(short algorithm);
	virtual ~ClusteringAlgorithm();

};

class Partitioner{
protected:

	string algDescription;

	vector<tVecDouble> centerCluster,centerClusterForBestJ,centerClusterTmp,centerClusterInitialForBestJ;
	//Indices
	//indices for each cluster and each feature
	tVecDouble3 T_CV, B_CV;
	tVecDouble2 J_CV;
	//do i have to work out statiscals for this indices?
	//Cluster heterogeneity with respect to single variables
	tVecDouble3 COR_CV,CTR_CV, CE_CV;	
    //indices for variables contribution
    tVecDouble2 COR_V,CTR_V;
    //indices for relative cluster heterogeneity
    tVecDouble2 T_C,B_C;
	tVecDouble J_C;
	//to measure time
    #ifdef DEBUGTIME
      #ifndef _WIN32 
         #define TIME timeb
      #else
         #define TIME _timeb
      #endif
      struct TIME timebufferBeforeRun,timebufferAfterRun;
	  struct TIME timebuffer;
    #endif
    
    //global center for new indices
    tVecDouble globalCenterZ;
	//
	bool bDoCleanUp;
	vector<tVecDouble > *tabData;
    StatisticalSummary **indexs;
    tVecDouble labelFeatureData;
	tVecDouble2 globalCenter;
    int nIndividuals,nFeatures,nClusters,labelFeature,nClasses;
	string fileOutName;
	FILE *out,*outInte,*outAux;
	FILE *outClusteringResult;
	FILE *outCustomReport;
	FILE *outDataStatistics;
    AdjustedRandIndex *aRand;
	int nIterations,nRuns;
	double avgIte;
	//index and statistic informations
	double *idxAndStatistic;

		vector<tVecInt> hardPartition;
		//methods
	void initializeCentroids();
	virtual void startUp()=0;
	virtual void cleanUp()=0;
    double getAdjustedRandIndex(vector<tVecInt> &hardPartition);
	virtual double objectiveFunction()=0;
    virtual void workOutCVIndices();
   virtual void workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime);
	virtual void workOutBTRPrime(double &Bp,double &Tp, double &Rp);
    virtual void workOutIndexAndStatistics(bool bCalculateJustJ);	
	virtual void makeInterpretationReport();
	virtual void makeOutPutReportFooter();
	virtual void makeOutPutCustomReport();
    virtual void makeOutPutDataStatistics();
	virtual void printIndividualDistributionMatrix(vector<tVecInt> &hardPartiion,string Header, FILE *out);
	virtual void printCentroids(vector<tVecDouble> centerCluster,string Header, FILE *out);
	virtual void makeOutPutReportAboutRun(int run);
    virtual void makeOutPutReportAboutIteration(int iteration);
	virtual void makeOutPutReportHeader(string fileInName,int runs, int iterations);

public:
	void doCleanUp();

    virtual void execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs=true,bool bWillUsed=false)=0; 

	double getBestJ();
    double getWorstJ();
	double getRCForBestJ();
	double getBestRC();
	double getJForBestRC();
	double getWorstRC();
	double getAvgIte();
	StatisticalSummary ** getStatisticalIndexs();
	vector<double> getBestJAssociatedValues();
    vector<string> getBestJAssociatedValuesDescription();
	virtual vector<string> getConfiguration()=0;
};
//Generalizes all classes which implements some kind of partition algorithm
class FuzzyPartitioner:  public Partitioner{

	void saveInitialCenterCluster();
	
protected:

   
	double ** membershipRaisedToM;
	double *** difIndFromCenter;
	vector<tVecDouble> membership,tmpInitialMembershipJ,bestInitialMembershipJ,bestMembershipJ;

	
	double parameterM;



    //methods
	void initializeMembership();
	void createHardPartition(vector<tVecDouble> &membership,vector<tVecInt> &hardPartiion);
	
	virtual double objectiveFunction();
	virtual void workOutCVIndices();
	inline virtual double getAdaptiveValue(int &c,int &p);
    virtual void workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime);
	virtual void workOutBTRPrime(double &Bp,double &Tp, double &Rp);
    virtual void workOutIndexAndStatistics(bool bCalculateJustJ);
	void saveMembershipForBestJ();
	virtual void makeOutPutReportFooter();
	virtual void makeOutPutCustomReport();
    void printMembershipMatrix(vector<tVecDouble> membership,string Header, FILE *out);
	void printHardPartition(vector<tVecDouble> membership,string Header, FILE *out);
	
	virtual void makeOutPutReportHeader(string fileInName,int runs, int iterations);
	inline double getMembershipRaisedToM(int c,int i);
	virtual void upDateCenterClusterForAdaptiveFuzzy();
	virtual void upDateCenterClustersAndRelated()=0;
    virtual void upDateMembershipAndRelated()=0;
	virtual void upDateLambda();
    virtual void saveLambdaForBestJ();

	virtual void startUp()=0;
	virtual void cleanUp()=0;
public:
	FuzzyPartitioner();

	virtual ~FuzzyPartitioner();
	virtual void execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs=true,bool bWillUsed=false); 
};

//CODIGO=1
//It partitione a dataset using FuzzyCMeans as proposed by Bezdek
class FuzzyCMeans: public FuzzyPartitioner{
  //methods
  FuzzyCMeans();
  void startUp();
  void cleanUp();
  //void initializeMembership();
  void upDateCenterClustersAndRelated();
  void upDateMembershipAndRelated();
  
public:
//  FuzzyCMeans();

  virtual ~FuzzyCMeans();
  //void execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs=true,bool bWillUsed=false);
  virtual vector<string> getConfiguration();
  friend class ClusteringAlgorithm;
};

//CODIGO=2
//It implements our approach to fuzzy partitioner. This approach consists of including a Adaptive value to each feature and here 
//we have a first of other models
class FuzzyCMeansAdaptive: public FuzzyPartitioner{
  tVecDouble lambda,bestJLambda;
  //methods
  FuzzyCMeansAdaptive();
  void startUp();
  void cleanUp();
  void upDateCenterClustersAndRelated();
  void upDateMembershipAndRelated();
  void printLambda(string header,FILE *out);
  void makeOutPutCustomReport();
  void makeOutPutReportFooter();
  void upDateLambda();
  void saveLambdaForBestJ();
protected:
  inline virtual double getAdaptiveValue(int &c,int &p);
  virtual double objectiveFunction();
  virtual void workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime);
  virtual void workOutBTRPrime(double &Bp,double &Tp, double &Rp);
public:
  virtual ~FuzzyCMeansAdaptive();
  //void execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs=true,bool bWillUsed=false);
  virtual vector<string> getConfiguration();
  friend class ClusteringAlgorithm;
};

//CODIGO=3
//It implements our approach to fuzzy partitioner. This approach consists of including a Adaptive value to each feature and here 
//we have a first of other models
class FuzzyCMeansAdaptive2: public FuzzyPartitioner{
  tVecDouble2 lambda,bestJLambda;

  //methods
  FuzzyCMeansAdaptive2();
  void startUp();
  void cleanUp();
  void upDateCenterClustersAndRelated();
  void upDateLambda();
  void upDateMembershipAndRelated();
  void printLambda(string header,FILE *out);
  void makeOutPutCustomReport();
  void makeOutPutReportFooter();
  void saveLambdaForBestJ();
protected:
  inline virtual double getAdaptiveValue(int &c,int &p);
  virtual double objectiveFunction();
  virtual void workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime);
  virtual void workOutBTRPrime(double &Bp,double &Tp, double &Rp);
public:
  virtual ~FuzzyCMeansAdaptive2();
  //void execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs=true,bool bWillUsed=false);
  virtual vector<string> getConfiguration();
  friend class ClusteringAlgorithm;
};


class KMeans: public Partitioner{

	vector<tVecInt> hardPartitionTmp,hardPartitionInitial,hardPartitionFinal;
		//methods
	KMeans();
  void saveCentroids(vector<tVecDouble> &centerCluster,vector<tVecDouble> &centerClusterO);
  void saveHardPartition(vector<tVecInt> &hardPartition, vector<tVecInt> &hardPartitionO);

  void startUp();
  void cleanUp();
  void upDateCenterClustersAndRelated();
  void allocate();
  void printHardPartition(vector<tVecInt> hardPartition,string Header, FILE *out);
  void makeOutPutCustomReport();
  void makeOutPutReportFooter();
protected:
  virtual double objectiveFunction();
public:
  virtual ~KMeans();
  virtual vector<string> getConfiguration();
	    virtual void execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs=true,bool bWillUsed=false); 
  friend class ClusteringAlgorithm;
};

#endif // !defined(AFX_ClusteringAlgorithm_H__62839795_566F_44CD_9EC4_9C45C8D09889__INCLUDED_)
