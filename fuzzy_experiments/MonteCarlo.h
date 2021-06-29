// MonteCarlo.h: interface for the MonteCarlo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MONTECARLO_H__42ECD131_0666_44D6_8675_FE260C872EA4__INCLUDED_)
#define AFX_MONTECARLO_H__42ECD131_0666_44D6_8675_FE260C872EA4__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "UserInterface.h"
#include "GlobalDefinitions.h"
#include "ClusteringAlgorithm.h"
#include "StatisticalSummary.h"
#include <iostream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <math.h>

#ifdef DEBUGTIME
  #include <sys/timeb.h>
  #include <time.h>
#endif

using namespace std;
class MonteCarlo;

class MonteCarloUserInterface: public UserInterface{
	void readOutPutFileName();
	void readClassParameters();
	void readVariableNumber();
	void readIntervals();
	void readRepetitions();
    MonteCarlo *mc;
public:
	MonteCarloUserInterface(MonteCarlo *mc);
	void readParameters();
};


class MonteCarlo  
{
  //tbl &tabData;
  vector<tVecDouble> tabData;
  vector<double> individualLabels;
  vector<double *>partitionVariables;
  vector<char *>partitionFeature;
  MonteCarloUserInterface *usrInt;
  FILE *out;
  FILE *mc;
  FILE *mcAUX;
  FILE *table;
  FILE *symbolicTable;
  FILE *sdsFile;
  string fileOutPreffixName;
  static int mcRepetitions;
  static int numBase;
  string fileSDSName;
  typedef struct{
     double lowerBound,upperBound;
  }interval;

      #ifdef DEBUGTIME
      #ifndef _WIN32 
         #define TIME timeb
      #else
         #define TIME _timeb
      #endif
      struct TIME timebufferBeforeRun,timebufferAfterRun;
	  struct TIME timebuffer;
    #endif

  typedef vector<double> vecDoubleType;
  typedef vector<vecDoubleType> vecDoubleType2;
  vector<vecDoubleType2> tempValue;
  vector<interval> intervals;
  vector<int> individuals;
  vector<double> correlation;
  vector<vecDoubleType> average;
  vector<vecDoubleType> standardDeviation;
  StatisticalSummary **indexs;
  int nVariables;
  int nClasses;
  int nIndividuals;

  void memoryAllocation();
  void createNormals(bool bMakeReport);
  void createSymbolics(double lower, double upper,bool bMakeReport);
  void createSODASFile();
  void freePartitionVector();
  void createPartitionVector();
  void printSummaryTitle();  
  void printOutStatistics(double lower, double upper, double bestRC,double worstRC, double ofForBestRC, double averageRC, double stdDeviationRC, double bestOF, double bestOFRC, double averageOF, double stdDeviationOF,double &NumIte,vector<double>& bestJAssociatedValues,vector<string>& bestJAssociatedValueDescriptions);
  void printOutPartitionerAlgorithmCfg(Partitioner &algorithm);
public:
	MonteCarlo(int &nClasses);
	MonteCarloUserInterface * getUserInterface();
	void initialize();
	void finalize();
	void execute(Partitioner &algorithm,string fileOutName,int runs,int iterations, double parameterM);
	void setSDSFileName(string sdsFileName);
	void setClassesNumber(int n);
	void setVariableNumber(int n);	
	void setIntervalNumber(int n);
	void setInterval(int interval,double lower,double upper);
	void setIndividuals(int cluster,int number);
	void setCorrelation(int cluster,double correlation);
	void setAverage(int cluster,int variable,double average);
	void setStandardDeviation(int cluster,int variable, double standardDeviation);
	double getNumClasses();
	double getNumVariables();
	static void setMCRepetitions(int repetitions);
	//tbl &getTabData();
	//void setTabData(tbl &tabData);
	vector<double *> getPartitionVariable();
    vector<char *> getPartitionFeature();
	string getBaseName();
	void printOutInformation(string str);
	virtual ~MonteCarlo();
    friend class MonteCarloUserInterface;
};

#endif // !defined(AFX_MONTECARLO_H__42ECD131_0666_44D6_8675_FE260C872EA4__INCLUDED_)
