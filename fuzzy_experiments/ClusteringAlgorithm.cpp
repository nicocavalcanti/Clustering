// ClusteringAlgorithm.cpp: implementation of the ClusteringAlgorithm class.
//
//////////////////////////////////////////////////////////////////////
#include <limits>
#include "ClusteringAlgorithm.h"


//Defines for debuging
//#define DEBUG_MEMBERSHIP
//#define DEBUG_LAMBDA
//#define DEBUG_CLUSTER_CENTER
#define DEBUGTIME
//#define PRINT_STATISTICS




//////////////////////////////////////////////////////////////////////
// ClusteringAlgorithm methods
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// ClusteringAlgorithm Construction/Destruction
//////////////////////////////////////////////////////////////////////

ClusteringAlgorithm::ClusteringAlgorithm(){

}

ClusteringAlgorithm::~ClusteringAlgorithm(){

}

//////////////////////////////////////////////////////////////////////
// Other ClusteringAlgorithm methods
//////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1) The object returned by this method must be destroyed by the caller of the method
//   (n)
// Coments:
//    XXX
// Parameters:
//   -experimentCode:
//       A code which is a unique identifier of the classes
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
Partitioner * ClusteringAlgorithm::getAlgorithm(short algorithm){
	switch (algorithm){
	  case 1:
		return new FuzzyCMeans();
	  case 2:
		return new FuzzyCMeansAdaptive();
	  case 3:
		return new FuzzyCMeansAdaptive2(); 
	  case 4:
		return new KMeans(); 
	};
	return NULL;
};



//////////////////////////////////////////////////////////////////////
// FuzzyCMeans methods
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// FuzzyCMeans Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyCMeans::FuzzyCMeans():FuzzyPartitioner(){
	//FuzzyPartitioner::FuzzyPartitioner();
	this->algDescription = "Fuzzy c-means";
};//FuzzyCMeans::FuzzyCMeans

FuzzyCMeans::~FuzzyCMeans(){
	if (this->bDoCleanUp)
		this->cleanUp();
};//FuzzyCMeans::~FuzzyCMeans


//////////////////////////////////////////////////////////////////////
// Other FuzzyCMeans methods
//////////////////////////////////////////////////////////////////////

//depend on having the individuals number known
void FuzzyCMeans::startUp(){
	FuzzyPartitioner::startUp();
};//FuzzyCMeans::startUp

void FuzzyCMeans::cleanUp(){
	FuzzyPartitioner::cleanUp();

};//method FuzzyCMeans::cleanUp



//post-condition: The membership matrix generated has to be all positive

//it updates the center cluster and the global center
void FuzzyCMeans::upDateCenterClustersAndRelated(){
	for (int c=0;c<this->nClusters;c++){
		for (int f=0;f<this->nFeatures;f++){
		  double deno = 0;
		  double num = 0;
                  int i;
		  for (i=0;i<this->nIndividuals;i++){
			double memb = this->membership[c][i];
			double membershipM = this->membershipRaisedToM[c][i];
			deno += membershipM;
			num += membershipM*(*this->tabData)[i][f];

#ifndef NO_INDF_CHECK
			if ( (deno<=INDF)  || (deno>=-INDF) || (num<=INDF)  || (num>=-INDF))
             cerr<<"\nERROR1 - Indefinition -> upDateCenterClusters()";
#endif
		  }//for i

#ifndef NOT_CHECK_DIVISION_BY_ZERO
		  if (deno==0)
			  cerr<<"\nERROR2 - Division by Zero -> upDateCenterClusters()";
#endif
		  this->centerCluster[c][f] = num/deno;
		  //Updates the difference from individuals to centercluster[c][f]
		  for (i=0;i<this->nIndividuals;i++){
			  this->difIndFromCenter[i][c][f] = pow((*this->tabData)[i][f]-this->centerCluster[c][f],2.0);
		  }//for i
		}//for f
	}//for c

#ifndef NO_INTERPRETATION_REPORT
	int f;
	for ( f=0;f<this->nFeatures;f++){
		double val=0;
		for (int c = 0;c<this->nClusters;c++){
			val += this->centerCluster[c][f];
		}//for c
		this->globalCenter[0][f] = val/this->nClusters;
	}//for f

	//Update globacenter prime
    for (f=0;f<this->nFeatures;f++){
		double tmpGCPNum=0,tmpGCPDeno=0;
		for (int i=0;i<this->nClusters;i++){
			for (int k=0;k<this->nIndividuals;k++){
              double tmp = this->membershipRaisedToM[i][k];
			  tmpGCPNum += tmp*(*this->tabData)[k][f];
			  tmpGCPDeno += tmp;
			}//for k
		}//for i
		this->globalCenter[1][f] = tmpGCPNum/tmpGCPDeno;
	}//for f
#endif

};//void FuzzyCMeans::upDateCenterClusters(){

void FuzzyCMeans::upDateMembershipAndRelated(){
	
	vector<double> dist;
	vector<bool> distZero;

	for (int i=0;i<this->nIndividuals;i++){

		int distZeros=0,c;
		dist.clear();
		distZero.clear();
		dist.resize(this->nClusters,0.0);		
		distZero.resize(this->nClusters,false);
		for ( c=0;c<this->nClusters;c++){
			dist[c] = 0.0;
			for (int f=0;f<this->nFeatures;f++){
			  dist[c] += this->difIndFromCenter[i][c][f];
			}//for f
			//dist[c] = pow(dist[c],0.5);
			if (dist[c]==0.0){
				distZero[c] = true;
				distZeros++;
			}
		}//for c
		double sum=0;
		bool bHadDistZero = distZeros>0;
		for (c=0;c<this->nClusters;c++){			
			if (bHadDistZero){
				if (distZero[c]){
					this->membership[c][i] = 1/distZeros;
					this->membershipRaisedToM[c][i] = pow(this->membership[c][i],this->parameterM);
				}
				else{
                    this->membership[c][i] = 0;
					this->membershipRaisedToM[c][i] = 0;
				}
			}
			else{
				double tmpDist=0;
				for (int j=0;j<this->nClusters;j++){

#ifndef NOT_CHECK_DIVISION_BY_ZERO
					if (dist[j]==0)
			          cerr<<"\nERROR - Division by Zero ->1  upDateMembership()";
#endif
					double div = dist[c]/dist[j];
					tmpDist += pow(div,1.0/(this->parameterM-1)); 
				}//for j
				this->membership[c][i]=1/tmpDist;
				this->membershipRaisedToM[c][i] = pow(this->membership[c][i],this->parameterM);

#ifndef NOT_CHECK_DIVISION_BY_ZERO
				if (this->membership[c][i]<0)
                   cerr<<"\nERROR - Membership cant be negative ->2  upDateMembership()";
				if (tmpDist==0)
			      cerr<<"\nERROR - Division by Zero ->3  upDateMembership()";
#endif
				sum += this->membership[c][i];
			}
			
		}//for c
		/*double one=1.0;
		if (sum>one){
			cout<<endl<<"ALERT - For cluster "<<(c+1)<<" membership sum differs from 1";
		}*/
		if (!bHadDistZero){
		  if (sum==0)
			cerr<<"\nERROR - Division by Zero ->4  upDateMembership()";
		  for (c=0;c<this->nClusters;c++){
		    this->membership[c][i] /= sum;
		  }//for c
		}
    
	}//for i
};//void FuzzyCMeans::upDateMembership(){

/*
void FuzzyCMeans::execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs,bool bWillUsed){
//initialize fields
	this->tabData = &tabData;
	this->nClusters = clusters;
	this->nIndividuals = tabData.size();
	this->nFeatures = tabData.begin()->size();
	//this->labelFeature = ;
	this->parameterM = parameterM;
	this->fileOutName = fileOutName;
	this->labelFeatureData = labelFeatureData;
	this->nClasses = nClasses;
	this->labelFeature = labelFeature;
	this->nIterations = iterations;
	this->nRuns = runs;
//First of all
	this->startUp();
    int run;
	int ite;
	long NumIte=0;
	double oldJ=MAX_DOUBLE,bestJ=MAX_DOUBLE;
	bool bConvergence=false;

	//It allows measure the time used to do certain operation in the algorithm
	#ifdef DEBUGTIME
	   struct _timeb timebufferBeforeRun,timebufferAfterRun;
	   struct _timeb timebuffer;
	#endif

	if (this->fileOutName.size()>0)
		this->makeOutPutReportHeader(fileInName,"Fuzzy c-means",runs,iterations);
	
    #ifdef DEBUGTIME
      _ftime( &timebufferBeforeRun );           
    #endif

	for (run=0;run<runs;run++){
		cout<<endl<<"*Initializing run "<<(run+1)<<endl;
		this->initializeMembership();
		bConvergence = false;
		oldJ=MAX_DOUBLE;
		if (this->fileOutName.size()>0)
		  this->makeOutPutReportAboutRun(run+1);
		for (ite=0;ite<iterations && !bConvergence;ite++){
          
		  //core elements of the clustering process
		    upDateCenterClusters();

          this->workOutIndexAndStatistics(!bAllwaysCalculateIndexs);
		  //Since J value is allways decreasing in each iterations there is no necessity to save the membership in each iterations so that you can save process time and save it just at the end of the iterations
		  //if (bAllwaysCalculateIndexs)
		  // this->saveMembershipForBestJ();
		

          if (fabs(oldJ-this->idxAndStatistic[0])<EPSLON){
			  bConvergence=true;
			  cout<<endl<<"Algorithm is stable";
		  }
	      // J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
	      // 0,1,2, 3, 4,      5,            6,        7, 8, 9
		  
		  if (oldJ<this->idxAndStatistic[0])
			  cout<<endl<<"ALERT - J increses in this iteration in: "<<(this->idxAndStatistic[0]-oldJ);
          oldJ = this->idxAndStatistic[0]; 
		  if (bAllwaysCalculateIndexs){  
			  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
				  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
			  }
			  
			  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
				  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
			  }
			  #ifdef PRINT_STATISTICS
			     printf("\n****Done iteration %d \tJ=%5.6f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%f R=%f  B/(J+B)=%f NF=%f D'=%f",ite+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
			 	   this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
              #endif
		  }//if

		  if (this->fileOutName.size()>0)
		    this->makeOutPutReportAboutIteration(ite+1);


          //The membership must be updated after the indixes and things like that had been worked out
		  if (!bConvergence && ite<iterations-1)
		    upDateMembership();
		  NumIte++;
          //cout<<endl<<"****Done iteration "<<(ite+1)<<"\tJ="<<J<<" R="<<R<<" T="<<T<<" B="<<B<<endl;
		}//for ite

        //As the algorithm makes J decrease over the iterations, J now has the smallest value
       // if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){
		
        if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){
		  this->saveMembershipForBestJ();
		  this->workOutCVIndices();
          if (run==0)
			bestJ = this->idxAndStatistic[0];//bestJ has not been initialized yet
		  else 
		    bestJ = this->indexs[0]->getBestValue();
        }
				

		if (!bAllwaysCalculateIndexs){
          this->workOutIndexAndStatistics(false);
		  
		  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
			  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
		  }
			  
		  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
			  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
		  }
			  
		  printf("\n****Done Run %d \tJ=%5.10f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%5.10f R=%f  B/(J+B)=%f NF=%f D'=%f",run+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
				  this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
		}
		printf("\n****Done Run %d \tJ=%5.10f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%5.10f R=%f  B/(J+B)=%f NF=%f D'=%f",run+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
				  this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
	}//for run

	this->avgIte = (double)NumIte/this->nRuns;
    #ifdef DEBUGTIME           
	   //It should be measured just after the end of the runs in order to be more accurated
	   _ftime( &timebufferAfterRun );    
    #endif

//creates the output file
    if (this->fileOutName.size()>0){
      this->makeOutPutReportFooter();
	  this->makeInterpretationReport();
	}
	#ifdef DEBUGTIME           
	   // _ftime( &timebufferAfterRun );           
	   timebuffer.millitm = timebufferAfterRun.millitm - timebufferBeforeRun.millitm;
	   (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?timebuffer.millitm += 1000 : 0; 
	   printf("\n    Time used to execute the runs: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);	
	   if (this->fileOutName.size()>0){
		   fprintf(this->out,"\n\n\n[Algorithm execution information]");
		   fprintf(this->out,"\n-Time used to execute the runs: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);
		   fprintf(this->out,"\n-Avg Iterations = %3.3f",(double)NumIte/this->nRuns);
       }
    #endif
//After all
	if (!bWillUsed)
	  this->cleanUp();

};//method FuzzyCMeans::execute
*/

vector<string> FuzzyCMeans::getConfiguration(){
    vector<string> strTmp;
	string strAux;
	strAux.assign("\nPartitioner Algorithm = Fuzzy C-Means");
	strAux.append("\nDistance Measure = Euclidian");
	strTmp.push_back(strAux);        
	strAux.assign("\nParameter M = ");
    char buffer[30];
	gcvt(this->parameterM,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	strAux.assign("\nRuns = ");
    gcvt(this->nRuns,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	strAux.assign("\nIterations per run = ");
    gcvt(this->nIterations,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	return strTmp;
};


//////////////////////////////////////////////////////////////////////
// FuzzyCMeansAdaptive methods
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// FuzzyCMeansAdaptive Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyCMeansAdaptive::FuzzyCMeansAdaptive():FuzzyPartitioner(){
	//FuzzyPartitioner::FuzzyPartitioner();
	this->lambda.clear();
	this->bestJLambda.clear();
this->algDescription = "Fuzzy c-means adaptive 1";
};//FuzzyCMeansAdaptive::FuzzyCMeansAdaptive


FuzzyCMeansAdaptive::~FuzzyCMeansAdaptive(){
	if (this->bDoCleanUp)
		this->cleanUp();
};//FuzzyCMeansAdaptive::~FuzzyCMeansAdaptive

//////////////////////////////////////////////////////////////////////
// Other FuzzyCMeansAdaptiveAdaptive methods
//////////////////////////////////////////////////////////////////////

void FuzzyCMeansAdaptive::startUp(){
	FuzzyPartitioner::startUp();
	this->lambda.resize(this->nFeatures,0.0);
	this->bestJLambda.resize(this->nFeatures,0.0);

};//FuzzyCMeansAdaptive::startUp

void FuzzyCMeansAdaptive::cleanUp(){
	FuzzyPartitioner::cleanUp();
	this->lambda.clear();
	this->bestJLambda.clear();
};//method FuzzyCMeansAdaptive::cleanUp


//it updates the center cluster and the global center
//it's exactally a copy of void FuzzyCMeans::upDateCenterClusters()
void FuzzyCMeansAdaptive::upDateCenterClustersAndRelated(){
	FuzzyPartitioner::upDateCenterClusterForAdaptiveFuzzy();
};//void FuzzyCMeans::upDateCenterClusters(){


void FuzzyCMeansAdaptive::upDateLambda(){
    
	int j,i,k,f;
	vector<double> sumDistIndToCenter;

    sumDistIndToCenter.clear();
    sumDistIndToCenter.resize(this->lambda.size(),0.0);

	for (j=0;j<this->lambda.size();j++){
	  double dist=0.0;
      for (k=0;k<this->nIndividuals;k++){		
		for (i=0;i<this->nClusters;i++){
			dist += this->membershipRaisedToM[i][k]*this->difIndFromCenter[k][i][j];
		}//for k
	  }//for i
      sumDistIndToCenter[j] = dist;
	}//for j

    double prod = 1;
	for (f=0;f<this->lambda.size();f++){
		prod *= sumDistIndToCenter[f];
	}//for f
	prod = pow(prod,1.0/this->lambda.size());

	for (j=0;j<this->lambda.size();j++){		
        //This calcutation above doesnt change when j change so that it can be worked out just one time and used here
		/*
		double prod = 1;
		for (f=0;f<this->lambda.size();f++){
			prod *= pow(sumDistIndToCenter[f],1.0/this->lambda.size());
		}//for f
		*/
		this->lambda[j] = prod/sumDistIndToCenter[j];

#ifndef NO_INDF_CHECK
        double dis = sumDistIndToCenter[j];
		if ( (this->lambda[j]<=INDF)  || (this->lambda[j]>=-INDF) ){
			cout<<endl<<"ALERT - Lambda got a indefined value";
		}
#endif
	}//for j

	//normalization, the product of each lambdaj has to be equals to 1
	//it's no longer necessary due to the lamda's formula fullfils the constrain above
/*   double prod = 1;
   for (j=0;j<this->lambda.size();j++){
	   double lamb = this->lambda[j];
	   if (j==this->lambda.size()-1){
		 this->lambda[j] *= (1/this->lambda[j])*prod;
		 break;
	   }
	   prod *= (2/(this->lambda[j]*this->lambda[j]));
	   this->lambda[j] *= this->lambda[j]/2;
   }//for j
*/
   #ifdef DEBUG_LAMBDA
     cout<<endl<<"Lambdas"<<endl;
	 double prodC=1;
     for (j=0;j<this->lambda.size();j++){
	   prodC *= this->lambda[j];
	   cout<<"L"<<(j+1)<<"="<<this->lambda[j]<<"    ";
	 }//for j     
	 cout<<"|Prod="<<prodC<<endl;
   #endif

};

void FuzzyCMeansAdaptive::upDateMembershipAndRelated(){
	
	vector<double> dist;
	vector<bool> distZero;

	for (int i=0;i<this->nIndividuals;i++){
		
		int distZeros=0,c;
		dist.clear();
		distZero.clear();
		dist.resize(this->nClusters,0.0);		
		distZero.resize(this->nClusters,false);
		for ( c=0;c<this->nClusters;c++){
			dist[c] = 0.0;
			for (int f=0;f<this->nFeatures;f++){
			  dist[c] += this->lambda[f]*this->difIndFromCenter[i][c][f];
			}//for f
			//dist[c] = pow(dist[c],0.5);
			if (dist[c]==0.0){
				distZero[c] = true;
				distZeros++;
			}
		}//for c
		double sum=0;
		bool bHadDistZero = distZeros>0;
		for (c=0;c<this->nClusters;c++){			
			if (bHadDistZero){
				if (distZero[c]){
					this->membership[c][i] = 1/distZeros;
					this->membershipRaisedToM[c][i] = pow(this->membership[c][i],this->parameterM);
				}
				else{
                    this->membership[c][i] = 0;
					this->membershipRaisedToM[c][i] = 0;
				}
			}
			else{
				double tmpDist=0;
				for (int j=0;j<this->nClusters;j++){
					if (dist[j]==0)
			          cerr<<"\nERROR - Division by Zero ->1  upDateMembership()";
					double div = dist[c]/dist[j];
					tmpDist += pow(div,1.0/(this->parameterM-1)); 
				}//for j
				this->membership[c][i]=1/tmpDist;
				this->membershipRaisedToM[c][i] = pow(this->membership[c][i],this->parameterM);

#ifndef NOT_CHECK_DIVISION_BY_ZERO
				if (this->membership[c][i]<0)
                   cerr<<"\nERROR - Membership cant be negative ->2  upDateMembership()";
				if (tmpDist==0)
			      cerr<<"\nERROR - Division by Zero ->3  upDateMembership()";
#endif
				//sum += this->membership[c][i];
			}//else
		}//for c
		if (!bHadDistZero){
		 /*f (sum==0)
			cerr<<"\nERROR - Division by Zero ->4  upDateMembership()";
		  for (c=0;c<this->nClusters;c++){
		    this->membership[c][i] /= sum;
		  }//for c
		  */
		} 
	}//for i
};

void FuzzyCMeansAdaptive::workOutBTRPrime(double &Bp,double &Tp, double &Rp){
 
	double retRP=0,retBP=0,retTP=0;
	int i,k,p;

//It calculates NonFuz, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += this->lambda[p]*pow((*this->tabData)[k][p]-this->globalCenter[1][p],2.0);
		  }//for p
		  retTP += this->membershipRaisedToM[i][k]*dist;		  
	  }//for k
  }//for i

    //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += this->lambda[p]*pow(this->centerCluster[i][p]-this->globalCenter[1][p],2.0);
	  }//for p
	  double u=0;
      for (k=0;k<this->nIndividuals;k++){
		  u += this->membershipRaisedToM[i][k];
	  }//for k
      retBP += dist*u;
  }//for i

  
  //the results are put into the parameters

  retRP = retBP/retTP;
  Rp=retRP;
  Bp=retBP;
  Tp=retTP;
  

};

void FuzzyCMeansAdaptive::workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime){

  double retR=0,retB=0,retT=0,retNonFuz=0,retDPrime=0;
  int i,k,p;
  
  //It calculates NonFuz, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += this->lambda[p]*pow((*this->tabData)[k][p]-this->globalCenter[0][p],2.0);
		  }//for p
		  retT += this->membershipRaisedToM[i][k]*dist;
		  //This is the Dunn's partition coefficient (1976), which is defined as the sum
		  //of squares of all the membership coefficients, divided by the number of objects, i.e
		  //  Fk(U)= S(i=1 to n)S(v=1 to k)pow(uik,2)/n
		  // Due to optimization questions the division by n will be done out of the sums we can do that since n neither depends on i nor k
          retNonFuz += pow(this->membership[i][k],2.0);
		  
	  }//for k
  }//for i

  retNonFuz /= (this->nIndividuals+0.0);//the nIndividual is added to 0 in order to assure that a division which results a double value
  //The coefficient above can be normalized to vary from 1 (hard cluster) to 0 (entirely fuzzy), independently of the number of clusters,
  //by the following transformation:
  //  Fk '(U) = (Fk(U)-(1/k))/( 1 - (1/k)) = (kFk(U)-1)/(k-1)
  //This normalized coefficient has sometimes been called "nonfuzziness index" (Roubens,1982).
  retNonFuz = (this->nClusters*retNonFuz - 1.0 ) / (this->nClusters -1.0);

  //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += this->lambda[p]*pow(this->centerCluster[i][p]-this->globalCenter[0][p],2.0);
	  }//for p
	  double u=0;
      for (k=0;k<this->nIndividuals;k++){
		  u += this->membershipRaisedToM[i][k];
	  }//for k
      retB += dist*u;
  }//for i

  //it calculates Dprime
  for (k=0;k<this->nIndividuals;k++){
	  int hardClus=-1;
	  double maxU=MIN_DOUBLE;
	  for (i=0;i<this->nClusters;i++){
		 if (this->membership[i][k]>maxU){
			 maxU = this->membership[i][k];
			 hardClus = i;
		 }
	  }//for i

	  for (i=0;i<this->nClusters;i++){
		  double uHard = 0.0;
		  if (i==hardClus)
            uHard = 1.0;
		  //It's a function of the membership, defined as Dk(U)=S(i=1 to N)S(v=1 to K)pow(Wiv-Uiv,2)/N
          //In this definition, W is the closest hard representation of the fuzzy U. Hence Dk(U) represents the average squared error of 
		  //a fuzzy clustering with respect to the closes hard clustering. It can be shown to vary between 0 (hard clustering) and 1-1/k (completely fuzzy)
		   // Due to optimization questions the division by n will be done out of the sums we can do that since n neither depends on i nor k
 		  retDPrime += pow(uHard-this->membership[i][k],2.0);
	  }//for i

  }//for k  
  retDPrime /= (this->nIndividuals+0.0);//the nIndividual is added to 0 in order to assure that a division which results a double value
  //Normalizing this function (in the same wa as Dunn's partition coefficient), we obtain: Dk'(U) = Dk(U)/(1-(1/K))=KDk(U)/(K-1)
  //To details see:Finding Groups in Data - an introdutory to cluster analysis. Kaufman, L and Rousseeuw, P. J.
  retDPrime = this->nClusters*retDPrime/(this->nClusters-1.0);

  //the results are put into the parameters

  retR = retB/retT;
  R=retR;
  B=retB;
  T=retT;
  nonFuz = retNonFuz;
  dPrime = retDPrime;
};


vector<string> FuzzyCMeansAdaptive::getConfiguration(){
vector<string> strTmp;
	string strAux;
	strAux.assign("\nPartitioner Algorithm = Fuzzy C-Means Adaptive 1");
	strAux.append("\nDistance Measure = Euclidian");
	strTmp.push_back(strAux);        
	strAux.assign("\nParameter M = ");
    char buffer[30];
	gcvt(this->parameterM,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	strAux.assign("\nRuns = ");
    gcvt(this->nRuns,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	strAux.assign("\nIterations per run = ");
    gcvt(this->nIterations,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	return strTmp;
};

double FuzzyCMeansAdaptive::getAdaptiveValue(int &c,int &p){
  return this->lambda[p];
};

double  FuzzyCMeansAdaptive::objectiveFunction(){
	double ret=0;
	for (int k=0;k<this->nIndividuals;k++){	
		for (int i=0;i<this->nClusters;i++){
			double dist=0;
			for (int f=0;f<this->nFeatures;f++){
			  	dist += this->lambda[f]*this->difIndFromCenter[k][i][f];
                //dist += pow(fabs((*this->tabData)[k][f]-this->centerCluster[i][f]),2.0);
			}//for f
			//dist = pow(dist,0.5);
			//ret += pow(this->membership[i][k],this->parameterM)*pow(dist,2.0);
			ret += this->membershipRaisedToM[i][k]*dist;
		}//for k
	}//for i
	return ret;
};

void FuzzyCMeansAdaptive::makeOutPutCustomReport(){

	FuzzyPartitioner::makeOutPutCustomReport();

	this->printLambda("Lambda values for the minimal J obtained",this->outCustomReport);
};

void FuzzyCMeansAdaptive::printLambda(string header,FILE *out){

	fprintf(out,"\n\n------------------- %s ------------------- ",header.c_str());
	for (int i=0;i<this->bestJLambda.size();i++){
		fprintf(out,"\n%25s%-3d = %.8f","Lambda for feature",i+1,this->bestJLambda[i]);
	}//for i
	fflush(out);
};

void FuzzyCMeansAdaptive::makeOutPutReportFooter(){
	
    
	this->printLambda("Lambda values for the minimal J obtained",this->out);
	FuzzyPartitioner::makeOutPutReportFooter();
};

void FuzzyCMeansAdaptive::saveLambdaForBestJ(){
	//The test is made in execute
	//if (this->idxAndStatistic[0]<this->indexs[0]->getPreivouslyLastValueAdded()){
	  for (int c=0;c<this->bestJLambda.size();c++){
		this->bestJLambda[c] = this->lambda[c];		
	  }//for c
	//}//if
};

/*
void FuzzyCMeansAdaptive::execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs,bool bWillUsed){
//initialize fields
	this->tabData = &tabData;
	this->nClusters = clusters;
	this->nIndividuals = tabData.size();
	this->nFeatures = tabData.begin()->size();
	//this->labelFeature = ;
	this->parameterM = parameterM;
	this->fileOutName = fileOutName;
	this->labelFeatureData = labelFeatureData;
	this->nClasses = nClasses;
	this->labelFeature = labelFeature;
	this->nIterations = iterations;
	this->nRuns = runs;


	long NumIte=0;
//First of all
	this->startUp();
    int run;
	int ite;
	double oldJ=MAX_DOUBLE,bestJ=MAX_DOUBLE;
	bool bConvergence=false;



	  //It allows measure the time used to do certain operation in the algorithm
	  #ifdef DEBUGTIME
		struct _timeb timebufferBeforeRun,timebufferAfterRun;
		struct _timeb timebuffer;
	  #endif

	if (this->fileOutName.size()>0)
		this->makeOutPutReportHeader(fileInName,"Adaptive Fuzzy c-means 1",runs,iterations);

	#ifdef DEBUGTIME
      _ftime( &timebufferBeforeRun );           
    #endif

	for (run=0;run<runs;run++){
		cout<<endl<<"*Initializing run "<<(run+1)<<endl;
		
		this->initializeMembership();
		//initialize  the center prototypes with the randomically generated U matrix
		this->upDateCenterClusters();
		bConvergence = false;
		oldJ=MAX_DOUBLE;

		if (this->fileOutName.size()>0)
		  this->makeOutPutReportAboutRun(run+1);

		for (ite=0;ite<iterations && !bConvergence;ite++){
          //double J,R,T,B;      

		  //core elements of the clustering process
            this->upDateLambda();

		  
          this->workOutIndexAndStatistics(!bAllwaysCalculateIndexs);
		  
		  //Since J value is allways decreasing in each iterations there is no necessity to save the membership in each iterations so that you can save process time and save it just at the end of the iterations
		  //if (bAllwaysCalculateIndexs){
		    //this->saveMembershipForBestJ();
		  //}
		  
          if (fabs(oldJ-this->idxAndStatistic[0])<EPSLON){
			  bConvergence=true;
			  cout<<endl<<"Algorithm is stable";
		  }
	      // J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
	      // 0,1,2, 3, 4,      5,            6,        7, 8, 9
		  if (oldJ<this->idxAndStatistic[0])
			  cout<<endl<<"ALERT - J increses in this iteration in: "<<(this->idxAndStatistic[0]-oldJ);
		  if (bAllwaysCalculateIndexs){  
			  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
				  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
			  }
			  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
				  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
			  }
              #ifdef PRINT_STATISTICS
			     printf("\n****Done iteration %d \tJ=%5.6f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%f R=%f  B/(J+B)=%f NF=%f D'=%f",ite+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
			     this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
              #endif
		  }//if
		  oldJ = this->idxAndStatistic[0];

		  if (this->fileOutName.size()>0)
		    this->makeOutPutReportAboutIteration(ite+1);
          
          //The membership must be updated after the indixes and things like that had been worked out
		  if (!bConvergence && ite<iterations-1){
		    this->upDateMembership();
		    this->upDateCenterClusters();
		  }
		  NumIte++;
          //cout<<endl<<"****Done iteration "<<(ite+1)<<"\tJ="<<J<<" R="<<R<<" T="<<T<<" B="<<B<<endl;
		}//for ite

		this->avgIte = (double)NumIte/this->nRuns;

		//As the algorithm makes J decrease over the iterations, J now has the smallest value
        if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){
		  this->saveMembershipForBestJ();
		  this->workOutCVIndices();
          this->saveLambdaForBestJ();
		  if (run==0)
			bestJ = this->idxAndStatistic[0];//bestJ has not been initialized yet
		  else 
		    bestJ = this->indexs[0]->getBestValue();
        }

		if (!bAllwaysCalculateIndexs){
          this->workOutIndexAndStatistics(false);		  
		  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
			  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
		  }
			  
		  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
			  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
		  }
			  
		  printf("\n****Done Run %d \tJ=%5.6f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%f R=%f  B/(J+B)=%f NF=%f D'=%f",run+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
				  this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
		}
	}//for run
    if (this->fileOutName.size()>0){
      this->makeOutPutReportFooter();
	  this->makeInterpretationReport();
	}
	#ifdef DEBUGTIME           
	   _ftime( &timebufferAfterRun );           
	   timebuffer.millitm = timebufferAfterRun.millitm - timebufferBeforeRun.millitm;
	   (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?timebuffer.millitm += 1000 : 0; 
	   printf("\n    Time used to execute the runs: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);	
	   if (this->fileOutName.size()>0){
		   fprintf(this->out,"\n\n\n[Algorithm execution information]");
		   fprintf(this->out,"\n-Time used to execute the runs: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);
		   fprintf(this->out,"\n-Avg Iterations = %3.3f",(double)NumIte/this->nRuns);
       }
    #endif
//After all
	if  (!bWillUsed) 
	  this->cleanUp();

};//method FuzzyCMeansAdaptive::execute
*/

//////////////////////////////////////////////////////////////////////
// FuzzyCMeansAdaptive2 methods
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// FuzzyCMeansAdaptive2 Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyCMeansAdaptive2::FuzzyCMeansAdaptive2():FuzzyPartitioner(){
	//FuzzyPartitioner:;
	this->lambda.clear();
	this->bestJLambda.clear();
	this->algDescription = "Fuzzy c-means adaptive 2";
};//FuzzyCMeansAdaptive2::FuzzyCMeansAdaptive2


FuzzyCMeansAdaptive2::~FuzzyCMeansAdaptive2(){
	if (this->bDoCleanUp)
		this->cleanUp();
};//FuzzyCMeansAdaptive2::~FuzzyCMeansAdaptive2

//////////////////////////////////////////////////////////////////////
// Other FuzzyCMeansAdaptive2Adaptive methods
//////////////////////////////////////////////////////////////////////

void FuzzyCMeansAdaptive2::startUp(){
	FuzzyPartitioner::startUp();
	this->lambda.resize(this->nClusters);
	this->bestJLambda.resize(this->nClusters);
	for (int i=0;i<this->nClusters;i++){
		this->lambda[i].resize(this->nFeatures,0.0);
		this->bestJLambda[i].resize(this->nFeatures,0.0);
	}//for i
  
};//FuzzyCMeansAdaptive2::startUp

void FuzzyCMeansAdaptive2::cleanUp(){
	FuzzyPartitioner::cleanUp();
	for (int i=0;i<this->nClusters;i++){
		this->lambda[i].clear();
		this->bestJLambda[i].clear();
	}//for i
	this->lambda.clear();
	this->bestJLambda.clear();

};//method FuzzyCMeansAdaptive2::cleanUp


//it updates the center cluster and the global center
//it's exactally a copy of void FuzzyCMeans::upDateCenterClusters()
void FuzzyCMeansAdaptive2::upDateCenterClustersAndRelated(){
	FuzzyPartitioner::upDateCenterClusterForAdaptiveFuzzy();
//	int h,j,i,k;
//	for (i=0;i<this->nClusters;i++){
//		for (j=0;j<this->nFeatures;j++){
//           double num=0,deno=0;
//		   for (k=0;k<this->nIndividuals;k++){
//			   double memberM = this->membershipRaisedToM[i][k];//pow (this->membership[i][k],this->parameterM);
//               num += memberM*(*this->tabData)[k][j];
//			   deno += memberM;
//		   }//for k
//		   this->centerCluster[i][j] = num/deno;
//		   if (deno==0)
//				cout<<endl<<"ALERT - In the prototype calculation appear a denominator equals to zero";
//		}//for k
//	}//for i
	//Update global center
//    int f;
//	for (f=0;f<this->nFeatures;f++){
//		double val=0;
//		for (int c = 0;c<this->nClusters;c++){
//			val += this->centerCluster[c][f];
//		}//for c
//		this->globalCenter[0][f] = val/this->nClusters;
//	}//for f
	//Update global center prime
	//Update globacenter prime
//    for (f=0;f<this->nFeatures;f++){
//		double tmpGCPNum=0,tmpGCPDeno=0;
//		for (int i=0;i<this->nClusters;i++){
//			for (int k=0;k<this->nIndividuals;k++){
//              double tmp = this->membershipRaisedToM[i][k];//pow(this->membership[i][k],this->parameterM);
//			  tmpGCPNum += tmp*(*this->tabData)[k][f];
//			  tmpGCPDeno += tmp;
//			}//for k
//		}//for i
//		this->globalCenter[1][f] = tmpGCPNum/tmpGCPDeno;
//	}//for f
#ifdef DEBUG_CLUSTER_CENTER
	cout<<endl<<"Cluster Centers";
	for (int c=0;c<this->nClusters;c++){
		cout<<endl<<"Centro "<<(c+1)<<" -> ";
		for (int f=0;f<this->nFeatures;f++){
			cout<<"C"<<(c+1)<<(f+1)<<"="<<this->centerCluster[c][f]<<",";
		}
	}
#endif
};//void FuzzyCMeans::upDateCenterClusters(){


void FuzzyCMeansAdaptive2::upDateLambda(){
    

	int j,i,k,f;
	tVecDouble2 sumDistIndToCenter;

    sumDistIndToCenter.clear();
    sumDistIndToCenter.resize(this->lambda.size());
	for (i=0;i<this->lambda.size();i++){
      sumDistIndToCenter[i].resize(this->lambda[i].size(),0.0);
	}//for i

	for (i=0;i<this->lambda.size();i++){
	  
	  for (j=0;j<this->lambda[i].size();j++){
        double dist=0.0;
        
		for (k=0;k<this->nIndividuals;k++){		
		  dist += this->membershipRaisedToM[i][k]*this->difIndFromCenter[k][i][j];
		}//for k
		
		sumDistIndToCenter[i][j] = dist;
	  }//for j
       
	}//for i

	for (i=0;i<this->lambda.size();i++){
      double prod = 1;
	  for (f=0;f<this->lambda[i].size();f++){
	 	prod *= sumDistIndToCenter[i][f];
	  }//for f
	  prod = pow(prod,1.0/this->lambda[i].size());

	  for (j=0;j<this->lambda[i].size();j++){		
        //This calcutation above doesnt change when j change so that it can be worked out just one time and used here
		//
		//double prod = 1;
		//for (f=0;f<this->lambda.size();f++){
		//	prod *= pow(sumDistIndToCenter[f],1.0/this->lambda.size());
		//}//for f
		//
		this->lambda[i][j] = prod/sumDistIndToCenter[i][j];

#ifndef NO_INDF_CHECK
        double dis = sumDistIndToCenter[i][j];
		if ( (this->lambda[i][j]<=INDF)  || (this->lambda[i][j]>=-INDF) ){
			cout<<endl<<"ALERT - Lambda got a indefined value";
		}
#endif
	  }//for j
	}//for i

   #ifdef DEBUG_LAMBDA
     cout<<endl<<"Lambdas"<<endl;
	 double prodC=1;
	 for (i=0;i<this->lambda.size();i++){
	   prodC = 1;
	   cout<<endl<<"Cluster "<<(i+1)<<": "<<endl;
       for (j=0;j<this->lambda[i].size();j++){
	     prodC *= this->lambda[i][j];
	     cout<<"L"<<(j+1)<<"="<<this->lambda[i][j]<<"    ";
	   }//for j     
	   cout<<" | Prod="<<prodC<<endl;
	   if (fabs(1-prodC)>EPSLON)
		   cout<<endl<<"ALERT - Produt is different from one";
	}//for i
	 
   #endif

};

void FuzzyCMeansAdaptive2::upDateMembershipAndRelated(){

	vector<double> dist;
	vector<bool> distZero;

	for (int k=0;k<this->nIndividuals;k++){
		
		int distZeros=0,i;
		dist.clear();
		distZero.clear();
		dist.resize(this->nClusters,0.0);		
		distZero.resize(this->nClusters,false);
		for ( i=0;i<this->nClusters;i++){
			dist[i] = 0.0;
			for (int f=0;f<this->nFeatures;f++){
			  dist[i] += this->lambda[i][f]*this->difIndFromCenter[k][i][f];
			}//for f
			//dist[c] = pow(dist[c],0.5);
			if (dist[i]==0.0){
				distZero[i] = true;
				distZeros++;
			}
		}//for c
		double sum=0;
		bool bHadDistZero = distZeros>0;
		for (i=0;i<this->nClusters;i++){			
			if (bHadDistZero){
				if (distZero[i]){
					this->membership[i][k] = 1/distZeros;
					this->membershipRaisedToM[i][k] = pow(this->membership[i][k],this->parameterM);
				}
				else{
                    this->membership[i][k] = 0;
					this->membershipRaisedToM[i][k] = 0;
				}
			}
			else{
				double tmpDist=0;
//				for (int j=0;j<this->nClusters;j++){
//					if (dist[j]==0)
//			          cerr<<"\nERROR - Division by Zero ->1  upDateMembership()";
//					double div = dist[c]/dist[j];
//					tmpDist += pow(div,1.0/(this->parameterM-1)); 
//				}//for j
				int p=0;
				double Num=dist[i];
			
				for (int h=0;h<this->nClusters;h++){
  					double deno= dist[h];
                    tmpDist += pow(Num/deno,(1.0/(this->parameterM-1)));
				}//for h
				this->membership[i][k]=1.0/tmpDist;
				this->membershipRaisedToM[i][k] = pow(this->membership[i][k],this->parameterM);

#ifndef NOT_CHECK_DIVISION_BY_ZERO
				if (this->membership[i][k]<0)
                   cerr<<"\nERROR - Membership cant be negative ->2  upDateMembership()";
				if (tmpDist==0)
			      cerr<<"\nERROR - Division by Zero ->3  upDateMembership()";
#endif
				//sum += this->membership[i][k];
			}
		}//for i
	}//for k

#ifdef DEBUG_MEMBERSHIP
    cout<<endl<<endl<<"Membership Matrix"<<endl;
	cout<<"Individual -> Memberships to C1 | Cn : Sum"<<endl;
	for (k=0;k<this->nIndividuals;k++){
		cout<<(k+1)<<" -> ";
		double sum=0;
		for (int i=0;i<this->nClusters;i++){
			cout<<this->membership[i][k]<<" | ";
			sum += this->membership[i][k];
		}//for c
		cout<<" : "<<sum<<endl;
		if (fabs(1-sum)>EPSLON)
			cout<<endl<<"ALERT - Membership sum differs from one"<<endl;
	}//for k
	
  #endif
};

void FuzzyCMeansAdaptive2::workOutBTRPrime(double &Bp,double &Tp, double &Rp){
 
	double retRP=0,retBP=0,retTP=0;
	int i,k,p;

//It calculates NonFuz, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += this->lambda[i][p]*pow((*this->tabData)[k][p]-this->globalCenter[1][p],2.0);
		  }//for p
		  retTP += this->membershipRaisedToM[i][k]*dist;		  ;//pow(this->membership[i][k],this->parameterM)*dist;		  
	  }//for k
  }//for i

    //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += this->lambda[i][p]*pow(this->centerCluster[i][p]-this->globalCenter[1][p],2.0);
	  }//for p
	  double u=0;
      for (k=0;k<this->nIndividuals;k++){
		  u += this->membershipRaisedToM[i][k];//pow(this->membership[i][k],this->parameterM);
	  }//for k
      retBP += dist*u;
  }//for i

  
  //the results are put into the parameters

  retRP = retBP/retTP;
  Rp=retRP;
  Bp=retBP;
  Tp=retTP;
  

};

void FuzzyCMeansAdaptive2::workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime){

  double retR=0,retB=0,retT=0,retNonFuz=0,retDPrime=0;
  int i,k,p;
  
  //It calculates NonFuz, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += this->lambda[i][p]*pow((*this->tabData)[k][p]-this->globalCenter[0][p],2.0);
		  }//for p
		  retT += this->membershipRaisedToM[i][k]*dist;//pow(this->membership[i][k],this->parameterM)*dist;
		  //This is the Dunn's partition coefficient (1976), which is defined as the sum
		  //of squares of all the membership coefficients, divided by the number of objects, i.e
		  //  Fk(U)= S(i=1 to n)S(v=1 to k)pow(uik,2)/n
		  // Due to optimization questions the division by n will be done out of the sums we can do that since n neither depends on i nor k
          retNonFuz += pow(this->membership[i][k],2.0);
		  
	  }//for k
  }//for i

  retNonFuz /= (this->nIndividuals+0.0);//the nIndividual is added to 0 in order to assure that a division which results a double value
  //The coefficient above can be normalized to vary from 1 (hard cluster) to 0 (entirely fuzzy), independently of the number of clusters,
  //by the following transformation:
  //  Fk '(U) = (Fk(U)-(1/k))/( 1 - (1/k)) = (kFk(U)-1)/(k-1)
  //This normalized coefficient has sometimes been called "nonfuzziness index" (Roubens,1982).
  retNonFuz = (this->nClusters*retNonFuz - 1.0 ) / (this->nClusters -1.0);

  //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += this->lambda[i][p]*pow(this->centerCluster[i][p]-this->globalCenter[0][p],2.0);
	  }//for p
	  double u=0;
      for (k=0;k<this->nIndividuals;k++){
		  u += this->membershipRaisedToM[i][k];//pow(this->membership[i][k],this->parameterM);
	  }//for k
      retB += dist*u;
  }//for i

  //it calculates Dprime
  for (k=0;k<this->nIndividuals;k++){
	  int hardClus=-1;
	  double maxU=MIN_DOUBLE;
	  for (i=0;i<this->nClusters;i++){
		 if (this->membership[i][k]>maxU){
			 maxU = this->membership[i][k];
			 hardClus = i;
		 }
	  }//for i

	  for (i=0;i<this->nClusters;i++){
		  double uHard = 0.0;
		  if (i==hardClus)
            uHard = 1.0;
		  //It's a function of the membership, defined as Dk(U)=S(i=1 to N)S(v=1 to K)pow(Wiv-Uiv,2)/N
          //In this definition, W is the closest hard representation of the fuzzy U. Hence Dk(U) represents the average squared error of 
		  //a fuzzy clustering with respect to the closes hard clustering. It can be shown to vary between 0 (hard clustering) and 1-1/k (completely fuzzy)
		   // Due to optimization questions the division by n will be done out of the sums we can do that since n neither depends on i nor k
 		  retDPrime += pow(uHard-this->membership[i][k],2.0);
	  }//for i

  }//for k  
  retDPrime /= (this->nIndividuals+0.0);//the nIndividual is added to 0 in order to assure that a division which results a double value
  //Normalizing this function (in the same wa as Dunn's partition coefficient), we obtain: Dk'(U) = Dk(U)/(1-(1/K))=KDk(U)/(K-1)
  //To details see:Finding Groups in Data - an introdutory to cluster analysis. Kaufman, L and Rousseeuw, P. J.
  retDPrime = this->nClusters*retDPrime/(this->nClusters-1.0);

  //the results are put into the parameters

  retR = retB/retT;
  R=retR;
  B=retB;
  T=retT;
  nonFuz = retNonFuz;
  dPrime = retDPrime;

  
  //new indices calculation

};


vector<string> FuzzyCMeansAdaptive2::getConfiguration(){
vector<string> strTmp;
	string strAux;
	strAux.assign("\nPartitioner Algorithm = Fuzzy C-Means Adaptive 2");
	strAux.append("\nDistance Measure = Euclidian");
	strTmp.push_back(strAux);        
	strAux.assign("\nParameter M = ");
    char buffer[30];
	gcvt(this->parameterM,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	strAux.assign("\nRuns = ");
    gcvt(this->nRuns,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	strAux.assign("\nIterations per run = ");
    gcvt(this->nIterations,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	return strTmp;
};

double FuzzyCMeansAdaptive2::getAdaptiveValue(int &c,int &p){
  return this->lambda[c][p];
};

double  FuzzyCMeansAdaptive2::objectiveFunction(){
	double ret=0;
	for (int k=0;k<this->nIndividuals;k++){	
		for (int i=0;i<this->nClusters;i++){
			double dist=0;
			for (int f=0;f<this->nFeatures;f++){
			  	dist += this->lambda[i][f]*pow((*this->tabData)[k][f]-this->centerCluster[i][f],2.0);
                //dist += pow(fabs((*this->tabData)[k][f]-this->centerCluster[i][f]),2.0);
			}//for f
			//dist = pow(dist,0.5);
			//ret += pow(this->membership[i][k],this->parameterM)*pow(dist,2.0);
			ret += this->membershipRaisedToM[i][k]*dist;//pow(this->membership[i][k],this->parameterM)*dist;
		}//for k
	}//for i
	return ret;
};


void FuzzyCMeansAdaptive2::makeOutPutCustomReport(){

	FuzzyPartitioner::makeOutPutCustomReport();

	this->printLambda("Lambda values for the minimal J obtained",this->outCustomReport);
};

void FuzzyCMeansAdaptive2::printLambda(string header,FILE *out){

	fprintf(out,"\n\n------------------- %s ------------------- \n",header.c_str());
	fprintf(out,"  Cluster  ");
	int dash=11;
        int k;
	for (k=0;k<this->nFeatures;k++){
      fprintf(out,"  %7s%-3d  ","Feature",k+1);
	  dash += 14; 
	}//for k
	fprintf(out,"\n");
	for (k=0;k<dash;k++)
		fprintf(out,"-");

	for (int i=0;i<this->bestJLambda.size();i++){
		fprintf(out,"\n");
        fprintf(out,"  %7d  ",i+1);
		for (int f=0;f<this->bestJLambda[i].size();f++){
		  fprintf(out,"  %10.8f  ",this->bestJLambda[i][f]);
		}//for f
		
	}//for i
	fflush(out);

};

void FuzzyCMeansAdaptive2::makeOutPutReportFooter(){
    
	this->printLambda("Lambda values for the minimal J obtained",this->out);

	FuzzyPartitioner::makeOutPutReportFooter();
};

void FuzzyCMeansAdaptive2::saveLambdaForBestJ(){
	//The test is made in execute
	//if (this->idxAndStatistic[0]<this->indexs[0]->getPreivouslyLastValueAdded()){
	  for (int c=0;c<this->bestJLambda.size();c++){
		this->bestJLambda[c] = this->lambda[c];		
	  }//for c
	//}//if
};

/*
void FuzzyCMeansAdaptive2::execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs,bool bWillUsed){
//initialize fields
	this->tabData = &tabData;
	this->nClusters = clusters;
	this->nIndividuals = tabData.size();
	this->nFeatures = tabData.begin()->size();
	//this->labelFeature = ;
	this->parameterM = parameterM;
	this->fileOutName = fileOutName;
	this->labelFeatureData = labelFeatureData;
	this->nClasses = nClasses;
	this->labelFeature = labelFeature;
	this->nIterations = iterations;
	this->nRuns = runs;


	long NumIte=0;
//First of all
	this->startUp();
    int run;
	int ite;
	double oldJ=MAX_DOUBLE,bestJ=MAX_DOUBLE;
	bool bConvergence=false;



	  //It allows measure the time used to do certain operation in the algorithm
	  #ifdef DEBUGTIME
		struct _timeb timebufferBeforeRun,timebufferAfterRun;
		struct _timeb timebuffer;
	  #endif

	if (this->fileOutName.size()>0)
		this->makeOutPutReportHeader(fileInName,"Adaptive Fuzzy c-means 2",runs,iterations);

	#ifdef DEBUGTIME
      _ftime( &timebufferBeforeRun );           
    #endif

	for (run=0;run<runs;run++){
		cout<<endl<<"*Initializing run "<<(run+1)<<endl;
		
		this->initializeMembership();
		//initialize  the center prototypes with the randomically generated U matrix
		this->upDateCenterClusters();
		bConvergence = false;
		oldJ=MAX_DOUBLE;

		if (this->fileOutName.size()>0)
		  this->makeOutPutReportAboutRun(run+1);

		for (ite=0;ite<iterations && !bConvergence;ite++){
          //double J,R,T,B;        

		  cout<<endl<<"Iteration: "<<(ite+1)<<endl;

		  //core elements of the clustering process
            this->upDateLambda();

		  
          this->workOutIndexAndStatistics(!bAllwaysCalculateIndexs);
		  
		  //Since J value is allways decreasing in each iterations there is no necessity to save the membership in each iterations so that you can save process time and save it just at the end of the iterations
		  //if (bAllwaysCalculateIndexs){
		    //this->saveMembershipForBestJ();
		  //}
		  
          if (fabs(oldJ-this->idxAndStatistic[0])<EPSLON){
			  bConvergence=true;
			  cout<<endl<<"Algorithm is stable";
		  }
	      // J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
	      // 0,1,2, 3, 4,      5,            6,        7, 8, 9
		  if (oldJ<this->idxAndStatistic[0]){
			  int k,j,i;
			  cout<<endl<<"ALERT - J increses in this iteration in "<<(this->idxAndStatistic[0]-oldJ);
			  cout<<endl<<"Prototypes:"<<endl;
	          for (int c=0;c<this->nClusters;c++){
		        cout<<endl<<"Centro "<<(c+1)<<" -> ";
		        for (int f=0;f<this->nFeatures;f++){
			      cout<<"C"<<(c+1)<<(f+1)<<"="<<this->centerCluster[c][f]<<",";
				}
			  }
			  //
			  cout<<endl<<"Lambdas:"<<endl;
              double prodC=1;
	          for (i=0;i<this->lambda.size();i++){
	            prodC = 1;
	            cout<<endl<<"Cluster "<<(i+1)<<": "<<endl;
                for (j=0;j<this->lambda[i].size();j++){
	              prodC *= this->lambda[i][j];
	              cout<<"L"<<(j+1)<<"="<<this->lambda[i][j]<<"    ";
				}//for j     
	            cout<<" | Prod="<<prodC<<endl;
	            if (fabs(1-prodC)>EPSLON)
		          cout<<endl<<"ALERT - Produt is different from one";
			  }//for i
			  //
			  cout<<endl<<endl<<"Membership Matrix"<<endl;
	          cout<<"Individual -> Memberships to C1 | Cn : Sum"<<endl;
	          for (k=0;k<this->nIndividuals;k++){
		          cout<<(k+1)<<" -> ";
		          double sum=0;
		          for (int i=0;i<this->nClusters;i++){
			        cout<<this->membership[i][k]<<" | ";
			        sum += this->membership[i][k];
				  }//for c
		          cout<<" : "<<sum<<endl;
		          if (fabs(1-sum)>EPSLON)
			        cout<<endl<<"ALERT - Membership sum differs from one"<<endl;
			  }//for k
	          //
			  cout<<endl<<"Individuals and their features"<<endl;
	          for (i=0;i<this->tabData->size();i++){
		        tVecDouble vFeatures = tabData.at(i);
		        cout<<(i+1)<<"-> ";
		        for (int f=0;f<vFeatures.size();f++){
			      cout<<vFeatures[f]<<" | ";
				}//for f
		        cout<<endl;
			  }//for i
	          cout<<endl<<endl<<"Label Feature"<<endl;
	          cout<<"Individual | Class"<<endl;
	          for (int l=0;l<labelFeatureData.size();l++){
                cout<<(l+1)<<" | "<<this->labelFeatureData[l]<<endl;
			  }//for l
		  }
		  if (bAllwaysCalculateIndexs){  
			  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
				  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
			  }
			  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
				  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
			  }
              #ifdef PRINT_STATISTICS
			     printf("\n****Done iteration %d \tJ=%5.6f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%f R=%f  B/(J+B)=%f NF=%f D'=%f",ite+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
			     this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
              #endif
		  }//if
		  oldJ = this->idxAndStatistic[0];

		  if (this->fileOutName.size()>0)
		    this->makeOutPutReportAboutIteration(ite+1);
          
          //The membership must be updated after the indixes and things like that had been worked out
		  if (!bConvergence && ite<iterations-1){
		    this->upDateMembership();
		    this->upDateCenterClusters();
		  }
		  NumIte++;
          //cout<<endl<<"****Done iteration "<<(ite+1)<<"\tJ="<<J<<" R="<<R<<" T="<<T<<" B="<<B<<endl;
		}//for ite

		this->avgIte = (double)NumIte/this->nRuns;
		

		//As the algorithm makes J decrease over the iterations, J now has the smallest value
        if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){
		  this->saveMembershipForBestJ();
		  this->workOutCVIndices();
          this->saveLambdaForBestJ();
		  if (run==0)
			bestJ = this->idxAndStatistic[0];//bestJ has not been initialized yet
		  else 
		    bestJ = this->indexs[0]->getBestValue();
        }
		if (!bAllwaysCalculateIndexs){
          this->workOutIndexAndStatistics(false);		  
		  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
			  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
		  }
			  
		  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
			  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
		  }
			  
		  printf("\n****Done Run %d \tJ=%5.6f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%f R=%f  B/(J+B)=%f NF=%f D'=%f",run+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
				  this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
		}
	}//for run
    if (this->fileOutName.size()>0){
      this->makeOutPutReportFooter();
	  this->makeInterpretationReport();
	}
	#ifdef DEBUGTIME           
	   _ftime( &timebufferAfterRun );           
	   timebuffer.millitm = timebufferAfterRun.millitm - timebufferBeforeRun.millitm;
	   (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?timebuffer.millitm += 1000 : 0; 
	   printf("\n    Time used to execute the runs: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);	
	   if (this->fileOutName.size()>0){
		   fprintf(this->out,"\n\n\n[Algorithm execution information]");
		   fprintf(this->out,"\n-Time used to execute the runs: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);
		   fprintf(this->out,"\n-Avg Iterations = %3.3f",(double)NumIte/this->nRuns);
       }
    #endif
//After all
	if  (!bWillUsed) 
	  this->cleanUp();

};//method FuzzyCMeansAdaptive2::execute
*/


//K-Means

//
//
//
KMeans::KMeans():Partitioner(){
	this->algDescription = "K-means";
};//KMeans

//
//
//
void KMeans::startUp(){
  this->hardPartitionFinal.clear();
  this->hardPartitionInitial.clear();
  this->hardPartitionTmp.clear();
  this->hardPartitionFinal.resize(this->nClusters);
  this->hardPartitionInitial.resize(this->nClusters);
  this->hardPartitionTmp.resize(this->nClusters);
  Partitioner::startUp();
};//startUp
//
//
//
void KMeans::cleanUp(){
   Partitioner::cleanUp();
};//cleanUp
//
//
//
void KMeans::upDateCenterClustersAndRelated(){
	//compute centroids
    for (int c=0;c<this->nClusters;c++){		
			for (int f=0;f<this->nFeatures;f++){
				this->centerCluster[c][f] = 0.0;
				for (int i=0;i<this->hardPartition[c].size();i++){
				  this->centerCluster[c][f] += (*this->tabData)[this->hardPartition[c][i]][f]; 
				}//for i
                this->centerCluster[c][f] /= this->hardPartition[c].size();

				
			}//for f
		//this->hardPartition[c].clear();
	}//for c


#ifndef NO_INTERPRETATION_REPORT	
	//Update global center
    int f;
	for (f=0;f<this->nFeatures;f++){
		double val=0;
		for (int c = 0;c<this->nClusters;c++){
			val += this->centerCluster[c][f];
		}//for c
		this->globalCenter[0][f] = val/this->nClusters;
	}//for f
	//Update global center prime
	//Update globacenter prime
    for (f=0;f<this->nFeatures;f++){
		double tmpGCPNum=0,tmpGCPDeno=0;
		for (int i=0;i<this->nClusters;i++){
			for (int k=0;k<this->nIndividuals;k++){
              //double tmp = this->membershipRaisedToM[i][k];
			  tmpGCPNum += (*this->tabData)[k][f];//tmp*(*this->tabData)[k][f];
			  //tmpGCPDeno += tmp;
			}//for k
		}//for i
		this->globalCenter[1][f] = tmpGCPNum/this->nIndividuals;//tmpGCPDeno;
	}//for f
#endif

};//upDateCenterClustersAndRelated
//
//
//pre-condition: 
void KMeans::allocate(){
	
	for (int c=0;c<this->nClusters;c++){
      this->hardPartition[c].clear();
	}//for c

	for (int i=0;i<this->nIndividuals;i++){
		double minDist = MAX_DOUBLE;
		int cluster=-1;
        for (int c=0;c<this->nClusters;c++){
			double dist=0.0;
			for (int f=0;f<this->nFeatures;f++){
				dist += pow ((*this->tabData)[i][f]-this->centerCluster[c][f],2.0);
			}//for f
			if (dist<=minDist){
				minDist = dist;
				cluster = c;
			}
		}//for c
		this->hardPartition[cluster].push_back(i);
	}//for i
};//upDateMembershipAndRelated


void KMeans::printHardPartition(vector<tVecInt> hardPartition,string Header, FILE *out){

	int c;

	fprintf(out,"\n\n%s",Header.c_str());
	
    for ( c=0;c<hardPartition.size();c++){
        fprintf(out,"\nCluster %d -> ",c+1);
		int countIndLine=0;
		for (int j=0;j<hardPartition[c].size();j++,countIndLine++){
			fprintf(out,"%d,",hardPartition[c][j]+1);
			if (countIndLine==50){
			  fprintf(out,"\n");
			  countIndLine = 0;
			}
		}
		fprintf(out,"\nCluster size-> %d",hardPartition[c].size());
	}//for c

};


void KMeans::makeOutPutCustomReport(){


    this->printHardPartition(this->hardPartitionInitial,"Hard partition generated from Initial Centroids",this->outCustomReport);

	this->printCentroids(this->centerClusterInitialForBestJ,"Initial Center Cluster",this->outCustomReport);

	this->printIndividualDistributionMatrix(this->hardPartitionInitial,"Individuals Distribution into clusters according to Initial Centroids",this->outCustomReport);

    //Information related to Best J
	this->printHardPartition(this->hardPartitionFinal,"Hard partition for the minimal J obtained",this->outCustomReport);

	this->printCentroids(this->centerClusterForBestJ,"Center Cluster for the minimal J obtained",this->outCustomReport);

	this->printIndividualDistributionMatrix(this->hardPartitionFinal,"Individuals Distribution into clusters according to Final Centroids",this->outCustomReport);
};
//
//
//
void KMeans::makeOutPutReportFooter(){

	this->makeOutPutCustomReport();

	Partitioner::makeOutPutReportFooter();
};//makeOutPutReportFooter

//
//
//
double KMeans::objectiveFunction(){
	double J=0.0;
	for (int j=0;j<this->nClusters;j++){
		for (int i=0;i<this->hardPartition[j].size();i++){
			for (int f=0;f<this->nFeatures;f++){
				J += pow ((*this->tabData)[this->hardPartition[j][i]][f]-this->centerCluster[j][f],2.0);
			}//for f   
		}//for i
	}//for j

	return J;
};//objectiveFunction
//
//
//
KMeans::~KMeans(){
	if (this->bDoCleanUp)
		this->cleanUp();
};//~KMeans
//
//
//
vector<string> KMeans::getConfiguration(){
vector<string> strTmp;
	string strAux;
	strAux.assign("\nPartitioner Algorithm = K-means");
	strAux.append("\nDistance Measure = Euclidian");
	strTmp.push_back(strAux);        

	char buffer[30];

	strAux.assign("\nRuns = ");
    gcvt(this->nRuns,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	strAux.assign("\nIterations per run = ");
    gcvt(this->nIterations,10,buffer);
	strAux.append(buffer);
    strTmp.push_back(strAux);

	return strTmp;
};//getConfiguration

void KMeans::saveCentroids(vector<tVecDouble> &centerCluster,vector<tVecDouble> &centerClusterO){

	int c,f;
	for (c=0;c<this->nClusters;c++){
		for (f=0;f<this->nFeatures;f++){
			centerCluster[c][f] = centerClusterO[c][f];
		}//for f
	}//for c
};

void KMeans::saveHardPartition(vector<tVecInt> &hardPartition,vector<tVecInt> &hardPartitionO){

	int c,i;

	for (c=0;c<this->nClusters;c++){
		hardPartition[c].clear();
		for (i=0;i<hardPartitionO[c].size();i++){
			hardPartition[c].push_back(hardPartitionO[c][i]);
		}//for f
	}//for c

};

//
//
//
void KMeans::execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs,bool bWillUsed){
//initialize fields
	this->tabData = &tabData;
	this->nClusters = clusters;
	this->nIndividuals = tabData.size();
	this->nFeatures = tabData.begin()->size();
	//this->labelFeature = ;
	this->fileOutName = fileOutName;
	this->labelFeatureData = labelFeatureData;
	this->nClasses = nClasses;
	this->labelFeature = labelFeature;
	this->nIterations = iterations;
	this->nRuns = runs;
//First of all
	this->startUp();
    int run;
	int ite;
	long NumIte=0;
	double oldJ=MAX_DOUBLE,bestJ=MAX_DOUBLE;
	bool bConvergence=false;

	//It allows measure the time used to do certain operation in the algorithm
	#ifdef DEBUGTIME
	   //struct _timeb timebufferBeforeRun,timebufferAfterRun;
	   //struct _timeb timebuffer;
	#endif

	if (this->fileOutName.size()>0){
		this->makeOutPutReportHeader(fileInName,runs,iterations);
		this->makeOutPutDataStatistics();
	}
	
    #ifdef DEBUGTIME
      #ifndef _WIN32 
	    ftime( &timebufferBeforeRun );           
      #else
        _ftime( &timebufferBeforeRun );           

      #endif
    #endif
    cout<<endl<<this->algDescription<<endl;
	for (run=0;run<runs;run++){
		cout<<endl<<"*Initializing run "<<(run+1)<<endl;
		this->initializeCentroids();
		bConvergence = false;
		oldJ=MAX_DOUBLE;
        #ifndef NO_RUN_REPORT
		  if (this->fileOutName.size()>0)
		    this->makeOutPutReportAboutRun(run+1);
        #endif
		for (ite=0;ite<iterations && !bConvergence;ite++){
          
		  this->upDateCenterClustersAndRelated();

		  if (ite==0){
			  this->saveCentroids(this->centerClusterTmp,this->centerCluster);
			  this->saveHardPartition(this->hardPartitionTmp,this->hardPartition);
		  }//if

          this->allocate();
		  
		  

          this->workOutIndexAndStatistics(!bAllwaysCalculateIndexs);
		  //Since J value is allways decreasing in each iterations there is no necessity to save the membership in each iterations so that you can save process time and save it just at the end of the iterations
		  /*if (bAllwaysCalculateIndexs)
		    this->saveMembershipForBestJ();*/

          if (fabs(oldJ-this->idxAndStatistic[0])<EPSLON){
			  bConvergence=true;
			  cout<<endl<<"Algorithm is stable";
		  }
	      // J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
	      // 0,1,2, 3, 4,      5,            6,        7, 8, 9

#ifndef NO_OSCILATION_CHECK		  
		  if (oldJ<this->idxAndStatistic[0])
			  cout<<endl<<"ALERT - J increses in this iteration in: "<<(this->idxAndStatistic[0]-oldJ);
#endif
          oldJ = this->idxAndStatistic[0]; 
#ifndef NO_OSCILATION_CHECK	
		  if (bAllwaysCalculateIndexs){  
			  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
				  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
			  }
			  
			  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
				  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
			  }
			  #ifdef PRINT_STATISTICS
			     printf("\n****Done iteration %d \tJ=%5.6f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%f R=%f  B/(J+B)=%f NF=%f D'=%f",ite+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
			 	   this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
              #endif
		  }//if
#endif


		  #ifndef  NO_RUN_REPORT
		    if (this->fileOutName.size()>0)
		      this->makeOutPutReportAboutIteration(ite+1);
          #endif


          //The membership must be updated after the indixes and things like that had been worked out
		  //if (!bConvergence && ite<iterations-1){
		  NumIte++;
          //cout<<endl<<"****Done iteration "<<(ite+1)<<"\tJ="<<J<<" R="<<R<<" T="<<T<<" B="<<B<<endl;
		}//for ite

        //As the algorithm makes J decrease over the iterations, J now has the smallest value
       // if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){
		
        if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){

#ifndef NO_INTERPRETATION_INDICE
		  this->workOutCVIndices();
#endif
          if (run==0)
			bestJ = this->idxAndStatistic[0];//bestJ has not been initialized yet
		  else 
		    bestJ = this->indexs[0]->getBestValue();


          this->saveCentroids(this->centerClusterInitialForBestJ,this->centerClusterTmp);
		  this->saveHardPartition(this->hardPartitionInitial,this->hardPartitionTmp);

		  this->saveCentroids(this->centerClusterForBestJ,this->centerCluster);
		  this->saveHardPartition(this->hardPartitionFinal,this->hardPartition);


        }
				

		if (!bAllwaysCalculateIndexs){
          this->workOutIndexAndStatistics(false);
#ifndef NO_INTERPRETATION_INDICE		  
		  if (fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]))>EPSLON){
			  cout<<endl<<"ALERT - J+B differs from T. The difference is "<<fabs(this->idxAndStatistic[1]-(this->idxAndStatistic[0]+this->idxAndStatistic[4]));
		  }
			  
		  if (fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]))>EPSLON){
			  cout<<endl<<"ALERT - J+B' differs from T'. The difference is "<<fabs(this->idxAndStatistic[7]-(this->idxAndStatistic[0]+this->idxAndStatistic[8]));
		  }
#endif
          #ifndef DEBUGTIME  
		    printf("\n****Done Run %d \tJ=%5.10f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%5.10f R=%f  B/(J+B)=%f NF=%f D'=%f",run+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
				  this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
          #endif
		}
		#ifndef DEBUGTIME  
		  printf("\n****Done Run %d \tJ=%5.10f B=%f T=%f J+B=%f B'=%f T'=%f J+B'=%f R'=%f RC=%5.10f R=%f  B/(J+B)=%f NF=%f D'=%f",run+1,this->idxAndStatistic[0],this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[0]+this->idxAndStatistic[4],
				  this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[0]+this->idxAndStatistic[8],this->idxAndStatistic[9],this->idxAndStatistic[3],this->idxAndStatistic[2],this->idxAndStatistic[4]/(this->idxAndStatistic[0]+this->idxAndStatistic[4]),this->idxAndStatistic[5],this->idxAndStatistic[6]);
        #endif
	}//for run

	this->avgIte = (double)NumIte/this->nRuns;
    #ifdef DEBUGTIME           
	   //It should be measured just after the end of the runs in order to be more accurated
      #ifndef _WIN32 
	    ftime( &timebufferAfterRun );           
      #else
        _ftime( &timebufferAfterRun );           

      #endif   
    #endif



//creates the output file
    if (this->fileOutName.size()>0){
      this->makeOutPutReportFooter();

#ifndef NO_INTERPRETATION_INDICE
	  this->makeInterpretationReport();
#endif
	}

	#ifdef DEBUGTIME           
	   // _ftime( &timebufferAfterRun );   
	   //when generating output this is done
	   if (this->fileOutName.size()>0){
	     timebuffer.millitm = timebufferAfterRun.millitm - timebufferBeforeRun.millitm;
	     (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?timebuffer.millitm += 1000 : 0; 
	   }
	   printf("\n    Time used to execute the runs: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);	
	   
    #endif
//After all
	if (!bWillUsed)
	  this->cleanUp();
};//execute 