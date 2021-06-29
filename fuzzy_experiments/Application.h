// Application.h: interface for the Application class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_APPLICATION_H__94A0A29D_498D_43AF_BFE7_F554520BA5B5__INCLUDED_)
#define AFX_APPLICATION_H__94A0A29D_498D_43AF_BFE7_F554520BA5B5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "GlobalDefinitions.h"
#include "MonteCarlo.h"
#include "UserInterface.h"
#include <limits>
using namespace std;
class Application  {
	typedef void (Application::*tFuncProcessToken)(char *); 
	MonteCarlo *mc;
	bool bNeedParamM;
	UserInterface *usrInterface;
	string fileInName;
	string fileOutName;
    vector<bool> individuals;
	vector<bool> features;
	vector<tVecDouble> tabData;
	vector<double> individualLabels;
	short algorithm;
	string algorithmDesc;
	double parameterM;
	int runs;
	int iterations;
	int labelFeature;
	int nClasses;
	int nClusters;
    int appKind;
    //methods
	void readApplicationKind();
	void readParameters();
	void readIndividuals();
	void readFeatures();
	void readDataFileName();
	void readAlgorithm();
	void readMParameter();
	void readRuns();
	void readIterations();
	void readFileOutName();
	void readLabelFeature();
	void readNumClasses();
	void readClusterNumber();
    void loadData();
	void readTokens( string tokens,tFuncProcessToken func);
	void processIndividualTokens(char *token);
	void processFeatureTokens(char *token);
public:
	Application();
	virtual ~Application();
	void execute();
};

#endif // !defined(AFX_APPLICATION_H__94A0A29D_498D_43AF_BFE7_F554520BA5B5__INCLUDED_)
