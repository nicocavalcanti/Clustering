#include "ClusteringAlgorithm.h"

//#define DEBUG_CREATE_HARD_PARTITION
//#define DEBUG_MEMBERSHIP
//#define DEBUG_CLUSTER_CENTER

#define NO_RUN_REPORT
//#define DEBUG_MEMBERSHIP_INTIALIZATION

#define NUM_OF_GLOBAL_PROTOTYPE 2
#define NUMBER_INDEX_OBSERVED 10

//pre-condition
//there must be more elements than clusters
void Partitioner::initializeCentroids(){
	
	int c,i;

	vector<bool> bAllocated;

	bAllocated.resize(this->nIndividuals,false);

	for (i=0;i<this->nIndividuals;i++){
		bAllocated[i] = false;
	}//for i


    
	randomize(); 
    
	//First step
	//To guarantee that each cluster has at least one element

	for (c=0;c<this->nClusters;c++){

        int ind;
        
		this->hardPartition[c].clear();

		do{
          ind = rand()%this->nIndividuals;
		  if (!bAllocated[ind]){
            bAllocated[ind] = true;
			break;
		  }
		}while(true);

		this->hardPartition[c].push_back(ind);

	}//for c

	//Finish allocating all the element into any of the clusters
    int lastCluster = -1;
    for (i=0;i<this->nIndividuals;i++){
		if (bAllocated[i])
			continue;


         int cluster;

		 do{
			 cluster = rand()%this->nClusters;
			 if (cluster!=lastCluster){
				 break;
			 }
		 }while(true);
         lastCluster = cluster;
		 this->hardPartition[cluster].push_back(i);

	}//for i

	//compute centroids

//    for (c=0;c<this->nClusters;c++){		
//			for (int f=0;f<this->nFeatures;f++){
//				this->centerCluster[c][f] = 0.0;
//				for (i=0;i<this->hardPartition[c].size();i++){
//				  this->centerCluster[c][f] += (*this->tabData)[this->hardPartition[c][i]][f]; 
//				}//for i
//                this->centerCluster[c][f] /= this->hardPartition[c].size();
//
//				this->centerClusterTmp[c][f] = this->centerCluster[c][f];
//			}//for f
//
//			this->hardPartition[c].clear();
//	}//for c

};

void Partitioner::workOutIndexAndStatistics(bool bCalculateJustJ){

		  this->idxAndStatistic[0] = this->objectiveFunction();
		  
		  if (bCalculateJustJ){              
			  //this->indexs[0]->addValue(this->idxAndStatistic[0]);
			  return;
		  }

		  


		  this->workOutBTR(this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[2],this->idxAndStatistic[5],this->idxAndStatistic[6]);

#ifndef NO_INTERPRETATION_INDICE
		  this->workOutBTRPrime(this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[9]);
#else
    this->idxAndStatistic[8]=this->idxAndStatistic[7]=this->idxAndStatistic[9]=-1;
#endif
		  
		  this->idxAndStatistic[3] = this->getAdjustedRandIndex(this->hardPartition);

          for (int i=0;i<NUMBER_INDEX_OBSERVED;i++){
            this->indexs[i]->addValue(this->idxAndStatistic[i]);
		  }//for i
          this->indexs[3]->setBestValueAssociationsValue(0,this->idxAndStatistic[0]);
          // T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
          //indexs in best J associated values
		  // 1,2,0 ,3,   4          ,  5         , 6 ,  7, 8
		  //RC for best J
		  //this->indexs[0]->setBestValueAssociationDescription(0,this->indexs[3]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(0,this->idxAndStatistic[3]);
          //Non Fuzziness for best J
		  //this->indexs[0]->setBestValueAssociationDescription(4,this->indexs[5]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(4,this->idxAndStatistic[5]);
		  //Dk(U) Prime for best J
		  //this->indexs[0]->setBestValueAssociationDescription(5,this->indexs[6]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(5,this->idxAndStatistic[6]);

#ifndef NO_INTERPRETATION_INDICE
		  //T for best J
		  //this->indexs[0]->setBestValueAssociationDescription(1,this->indexs[1]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(1,this->idxAndStatistic[1]);
		  //R for best J
		  //this->indexs[0]->setBestValueAssociationDescription(2,this->indexs[2]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(2,this->idxAndStatistic[2]);
		  //B for best J
		  //this->indexs[0]->setBestValueAssociationDescription(3,this->indexs[4]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(3,this->idxAndStatistic[4]);
		  
		  //Tp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(6,this->indexs[7]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(6,this->idxAndStatistic[7]);
		  //Bp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(7,this->indexs[8]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(7,this->idxAndStatistic[8]);
		  //Rp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(8,this->indexs[9]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(8,this->idxAndStatistic[9]);
#else
		  //in order to make the report look nice
		  //T for best J
		  //this->indexs[0]->setBestValueAssociationDescription(1,this->indexs[1]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(1,-1);
		  //R for best J
		  //this->indexs[0]->setBestValueAssociationDescription(2,this->indexs[2]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(2,-1);
		  //B for best J
		  //this->indexs[0]->setBestValueAssociationDescription(3,this->indexs[4]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(3,-1);
		  
		  //Tp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(6,this->indexs[7]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(6,-1);
		  //Bp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(7,this->indexs[8]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(7,-1);
		  //Rp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(8,this->indexs[9]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(8,-1);
#endif
		  
};//void FuzzyPartitioner::workOutIndexAndStatistics(){

//Pre-conditions:
// The best membership must be already defined
void Partitioner::workOutCVIndices(){
	
	int i,j,m;

	double *B;
	B = (double *)malloc(sizeof(double)*NUM_OF_GLOBAL_PROTOTYPE);

	//T_C, J_C,B_C
    for (i=0;i<this->nClusters;i++){
       for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
		   this->B_C[m][i] = 0;
		   this->T_C[m][i] = 0;
		   for (j=0;j<this->nFeatures;j++){
			   this->T_CV[m][i][j] = 0;
			   this->B_CV[m][i][j] = 0;			   
		   }
	   }//for m
	   for (j=0;j<this->nFeatures;j++)
		   this->J_CV[i][j] = 0;
	   this->J_C[i] = 0;
	}//for i

	for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++)
      B[m] = 0;  

	//T_CV,J_CV,B_CV
	for (i=0;i<this->nClusters;i++){
	
		for (j=0;j<this->nFeatures;j++){
          for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
            for (int k=0;k<this->nIndividuals;k++){
			   this->T_CV[m][i][j] += pow( ((*this->tabData)[k][j]-this->globalCenter[m][j]),2.0);
			}//for k 
            this->T_C[m][i] += this->T_CV[m][i][j];
		  }//for m
		  for (int k=0;k<this->nIndividuals;k++)
		    this->J_CV[i][j] += pow( ((*this->tabData)[k][j]-this->centerCluster[i][j]),2.0);
		  //J_C
		  this->J_C[i] += this->J_CV[i][j];		  
		  for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
		    this->B_CV[m][i][j] = pow( (this->centerCluster[i][j]-this->globalCenter[m][j]),2.0);			
			//B_C
			this->B_C[m][i] += this->B_CV[m][i][j];
		  }//for m
		}//for j
        
		for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++)
		  B[m] += this->B_C[m][i];
	}//for i

	for (j=0;j<this->nFeatures;j++){
	    for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
          double Tj=0,Bj=0;
	      for (i=0;i<this->nClusters;i++){
			Tj += this->T_CV[m][i][j];
			Bj += this->B_CV[m][i][j];
			
			//CTR(j,i)
			this->CTR_CV[m][i][j] = this->B_CV[m][i][j]/this->B_C[m][i];
			//CE(j,i) 
			this->CE_CV[m][i][j] = this->B_CV[m][i][j]/B[m];
		  }//for i 
		  for (i=0;i<this->nClusters;i++)		
		    this->COR_CV[m][i][j] = this->B_CV[m][i][j]/Tj;
          
		  this->COR_V[m][j] = Bj/Tj;
		  this->CTR_V[m][j] = Bj/B[m];
		}//for m				
		
	}//for j
	
	free(B);
	

};//void FuzzyPartitioner::workOutCVIndices(){

void Partitioner::workOutBTRPrime(double &Bp,double &Tp, double &Rp){


	double retRP=0,retBP=0,retTP=0;
	int i,k,p;

//It calculates, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += pow((*this->tabData)[k][p]-this->globalCenter[1][p],2.0);
		  }//for p
		  retTP += dist;		  
	  }//for k
  }//for i

    //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += pow(this->centerCluster[i][p]-this->globalCenter[1][p],2.0);
	  }//for p
      retBP += dist;
  }//for i

  
  //the results are put into the parameters

  retRP = retBP/retTP;
  Rp=retRP;
  Bp=retBP;
  Tp=retTP;

};

void Partitioner::workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime){
  double retR=0,retB=0,retT=0,retNonFuz=0,retDPrime=0;

  int i,k,p;
  
  //It calculates NonFuz, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += pow((*this->tabData)[k][p]-this->globalCenter[0][p],2.0);
		  }//for p
		  retT += dist;
		  //This is the Dunn's partition coefficient (1976), which is defined as the sum
		  //of squares of all the membership coefficients, divided by the number of objects, i.e
		  //  Fk(U)= S(i=1 to n)S(v=1 to k)pow(uik,2)/n
		  // Due to optimization questions the division by n will be done out of the sums we can do that since n neither depends on i nor k
		  
	  }//for k
  }//for i

  
#ifndef NO_INTERPRETATION_INDICE
  //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += pow(this->centerCluster[i][p]-this->globalCenter[0][p],2.0);
	  }//for p
      retB += dist;
  }//for i
#endif

 
  //the results are put into the parameters

  retR = retB/retT;
  R=retR;
  B=retB;
  T=retT;
  nonFuz = retNonFuz;
  dPrime = retDPrime;
};

void Partitioner::makeOutPutReportHeader(string fileInName,int runs, int iterations){
	fprintf(this->out,"Results:\n");
    fprintf(this->out,"\n%28s = %s","Algorithm",this->algDescription.c_str());
	if (fileInName.size()>0)
	  fprintf(this->out,"\n%28s = %s","Input data",fileInName.c_str());

	fprintf(this->out,"\n%28s = %d","Clusters",this->nClusters);	
	fprintf(this->out,"\n%28s = %d\n%28s = %d","Run Number",runs,"Maximun Iteration Number",iterations);
	fflush(this->out);
    //put this information in the interpretation output file
	fprintf(this->outInte,"Results:\n");
    fprintf(this->outInte,"\n%28s = %s","Algorithm",this->algDescription.c_str());
	if (fileInName.size()>0)
	  fprintf(this->outInte,"\n%28s = %s","Input data",fileInName.c_str());
	fprintf(this->outInte,"\n%28s = %d","Clusters",this->nClusters);	
	fprintf(this->outInte,"\n%28s = %d\n%28s = %d","Run Number",runs,"Maximun Iteration Number",iterations);
	fflush(this->outInte);

};

void Partitioner::makeOutPutCustomReport(){

	for (int i=0;i<this->nIndividuals;i++){
		fprintf(this->outCustomReport,"%d\t%d",(i+1),this->labelFeatureData[i]);
		for (int f=0;f<this->nFeatures;f++){
			fprintf(this->outCustomReport,"\t%f",(*this->tabData)[i][f]);
		}//for f
		fprintf(this->outCustomReport,"\n");

	}//for i
	fflush(this->outCustomReport);
};


void Partitioner::makeOutPutDataStatistics(){

	int c,i,f;

	vector<tVecDouble> mean,covariance;
	vector<int> classCounter;
	vector<int> classMapper;
	vector<tVecInt> classInds;

	classInds.resize(this->nClasses);

	//It might be usuful for cases such as when class attribute assumes discontinuous values like 0,2,6,7 in such case
	//we are given 4 classes but we can not assume that each individual into the dataset will have eather 1 or 2 or 3 or 4
	//as its class value
	//For the case cited one solution could be make
	//classMapper[0] = 1
	//classMapper[2] = 2
	//classMapper[6] = 3
	//classMapper[7] = 4
	classMapper.resize(this->nClasses+1);
	for (c=1;c<=this->nClasses;c++)
		classMapper[c] = c-1;

	covariance.resize(this->nFeatures);

	for (f=0;f<this->nFeatures;f++){
		covariance[f].resize(f+1,0.0);
	}//for f

	mean.resize(this->nClasses);
	for (c=0;c<this->nClasses;c++){
		mean[c].resize(this->nFeatures,0.0);
	}//for c

	classCounter.resize(this->nClasses,0);


	//for each class compute its mean vector	
		
	for (i=0;i<this->nIndividuals;i++){

		for (f=0;f<this->nFeatures;f++){
			mean[classMapper[this->labelFeatureData[i]]][f] += (*this->tabData)[i][f];
			classCounter[classMapper[this->labelFeatureData[i]]]++;
			classInds[classMapper[this->labelFeatureData[i]]].push_back(i);
		}//for f

	}//for i

    for (c=0;c<this->nClasses;c++){
		for (f=0;f<this->nFeatures;f++){
		   mean[c][f] /=  classCounter[c];
		}//for f

	}//for c

    for (c=0;c<this->nClasses;c++){

		//Compute this class covariance matrix

		for (f=0;f<this->nFeatures;f++){
			for (int p=0;p<covariance[f].size();p++){
				covariance[f][p]=0;
				for (i=0;i<classInds[c].size();i++){
					covariance[f][p] += ((*this->tabData)[classInds[c][i]][f]-mean[c][f])*((*this->tabData)[classInds[c][i]][p]-mean[c][p]);
				}//for i
				covariance[f][p] /= classInds[c].size();
			}//for p
		}//for f	 
	  

	  //Print out this class mean vector

		fprintf(this->outDataStatistics,"\n\n---------- %10s%02d ---------- ","Class",(c+1));
        
		fprintf(this->outDataStatistics,"\n\nMean Vector\n");

		for (f=0;f<this->nFeatures;f++){
			fprintf(this->outDataStatistics,"\n%10s%03d\t%15.8f","Feature",(f+1),mean[c][f]);
		}//for f


	  //Print out this class covariance matrix

		fprintf(this->outDataStatistics,"\n\nCovariance Matrix\n");


        fprintf(this->outDataStatistics,"\n%15s\t"," ");

		for (f=0;f<this->nFeatures;f++){
			fprintf(this->outDataStatistics,"%12s%03d\t","Feature",(f+1));
		}//for f

		for (f=0;f<this->nFeatures;f++){
			fprintf(this->outDataStatistics,"\n%12s%03d\t","Feature",(f+1));
			for (int p=0;p<this->nFeatures;p++){
				int r=f,c=p;
				if (f<p){
					r=p;
					c=f;
				}
				fprintf(this->outDataStatistics,"%15.8f\t",covariance[r][c]);
			}//for p

		}//for f




	   fflush(this->outDataStatistics);


	}//for c

	fclose(this->outDataStatistics);

	
};

void Partitioner::printCentroids(vector<tVecDouble> centerCluster,string Header, FILE *out){

	fprintf(out,"\n\n%s",Header.c_str());

	int c,f;
	
	fprintf(out,"\n\n%12s"," ");

	for (f=0;f<this->nFeatures;f++){
	  fprintf(out,"\t%15s%02d","F",(f+1));
	}//for f

	for (c=0;c<this->nClusters;c++){
		fprintf(out,"\n%10s%02d","C",(c+1));
		for (f=0;f<this->nFeatures;f++){
          fprintf(out,"\t%17.7f",centerCluster[c][f]);
		}//for f
	}//for c
};

void Partitioner::printIndividualDistributionMatrix(vector<tVecInt> &hardPartition,string Header, FILE *out){

	fprintf(out,"\n\n%s",Header.c_str());

	int p,c;

	fprintf(out,"\n\n%12s"," ");
    for (c=0;c<this->nClasses;c++){
      fprintf(out,"\t%10s%02d","C",(c+1));
	}//for c

	vector<int> totalColumn;
	totalColumn.resize(this->nClasses,0);

	for (p=0;p<this->nClusters;p++){
		fprintf(out,"\n%10s%02d","P",(p+1));
		int total = 0;
        for (c=0;c<this->nClasses;c++){
			int totalC =0;
			for (int i=0;i<hardPartition[p].size();i++){
				if (this->labelFeatureData[hardPartition[p][i]]==(c+1))
					totalC++;
			}//for i
            fprintf(out,"\t%12d",totalC);
			total += totalC;
			totalColumn[c] += totalC;
		}//for c
		fprintf(out,"\t%12d",total);
	}//for p

	int total =0;
    
	fprintf(out,"\n%12s"," ");
	
	for (c=0;c<this->nClasses;c++){			
			
            fprintf(out,"\t%12d",totalColumn[c]);
			total += totalColumn[c];
			
	}//for c

	fprintf(out,"\t%12d",total);



};

void Partitioner::makeOutPutReportFooter(){
//creates the output file

	
	int c,i;
	
  
  
   fprintf(this->out,"\n\n Avg Number of Iterations = %6.3f",this->avgIte);

    fprintf(this->out,"\n\nAbout J:");
	fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[0]->getBestValue(),"RC for the best J",this->indexs[0]->getBestValueAssociationsValue(0),"Worst",this->indexs[0]->getWorstValue(),
		     "Std",this->indexs[0]->getStdDeviation(),"Avg",this->indexs[0]->getMean());
    
	fprintf(this->out,"\n\t\t---- The values of the others indexs for the minimal J ----");
	vector<double> vecIdxJVal = this->indexs[0]->getBestValueAssociationsValues();
	vector<string> vecIdxJDesc = this->indexs[0]->getBestValueAssociationsDescription();
    for (c=0;c<vecIdxJVal.size();c++){
      fprintf(this->out,"\n\t\t\t->%-30s = %.8f",vecIdxJDesc[c].c_str(),vecIdxJVal[c]);
	}//for c

	fprintf(this->out,"\nAbout RC:");
	fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[3]->getBestValue(),"J for the best RC",this->indexs[3]->getBestValueAssociationsValue(0),"Worst",this->indexs[3]->getWorstValue(),
		     "Std",this->indexs[3]->getStdDeviation(),"Avg",this->indexs[3]->getMean());

	fprintf(this->out,"\nAbout T:");
	fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[1]->getBestValue(),"Worst",this->indexs[1]->getWorstValue(),
		     "Std",this->indexs[1]->getStdDeviation(),"Avg",this->indexs[1]->getMean());

	fprintf(this->out,"\nAbout R:");
	fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[2]->getBestValue(),"Worst",this->indexs[2]->getWorstValue(),
		     "Std",this->indexs[2]->getStdDeviation(),"Avg",this->indexs[2]->getMean());

	fprintf(this->out,"\nAbout B:");
	fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[4]->getBestValue(),"Worst",this->indexs[4]->getWorstValue(),
		     "Std",this->indexs[4]->getStdDeviation(),"Avg",this->indexs[4]->getMean());

	fprintf(this->out,"\nAbout %s:",this->indexs[5]->getValueDescription().c_str());
	fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[5]->getBestValue(),"Worst",this->indexs[5]->getWorstValue(),
		     "Std",this->indexs[5]->getStdDeviation(),"Avg",this->indexs[5]->getMean());

	fprintf(this->out,"\nAbout %s:",this->indexs[6]->getValueDescription().c_str());
	fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[6]->getBestValue(),"Worst",this->indexs[6]->getWorstValue(),
		     "Std",this->indexs[6]->getStdDeviation(),"Avg",this->indexs[6]->getMean());

	for (i=7;i<NUMBER_INDEX_OBSERVED;i++){
      fprintf(this->out,"\nAbout %s:",this->indexs[i]->getValueDescription().c_str());
	  fprintf(this->out,"\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f\n\t\t->%-20s = %f","Best",this->indexs[i]->getBestValue(),"Worst",this->indexs[i]->getWorstValue(),
		     "Std",this->indexs[i]->getStdDeviation(),"Avg",this->indexs[i]->getMean());
	}//for i

	fflush(this->out);

	//auxiliar output
    fprintf(this->outAux,"Number of runs = %d\nNumber of iterations per run = %d",this->nRuns,this->nIterations);
    vector<string> strTmp = this->getConfiguration();
    for (vector<string>::iterator itr=strTmp.begin();itr!=strTmp.end();itr++){
      fprintf(this->outAux,itr->c_str());
	}
	//RC
    fprintf(this->outAux,"\n%08.6f %08.6f %08.6f %08.6f %08.6f ",this->indexs[3]->getBestValue(),this->indexs[3]->getBestValueAssociationsValue(0),this->indexs[3]->getWorstValue(),this->indexs[3]->getMean(),this->indexs[3]->getStdDeviation());
	//J
	fprintf(this->outAux,"%016.7f %016.7f %016.7f %016.7f ",this->indexs[0]->getBestValue(),this->indexs[0]->getWorstValue(),this->indexs[0]->getMean(),this->indexs[0]->getStdDeviation());
	//other indexs in relation to the minimal J
// T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
//indexs in best J associated values
// 1,2,0 ,3,   4          ,  5         , 6 ,  7, 8
	fprintf(this->outAux,"%016.7f %016.7f %016.7f %016.7f ",this->indexs[0]->getBestValueAssociationsValue(0),this->indexs[0]->getBestValueAssociationsValue(1),this->indexs[0]->getBestValueAssociationsValue(2),this->indexs[0]->getBestValueAssociationsValue(3));
	fprintf(this->outAux,"%016.7f %016.7f %016.7f %016.7f ",this->indexs[0]->getBestValueAssociationsValue(4),this->indexs[0]->getBestValueAssociationsValue(5),this->indexs[0]->getBestValueAssociationsValue(6),this->indexs[0]->getBestValueAssociationsValue(7));
	fprintf(this->outAux,"%016.7f ",this->indexs[0]->getBestValueAssociationsValue(8));
    fprintf(this->outAux,"%016.7f ",this->avgIte);

	#ifdef DEBUGTIME           
	   // _ftime( &timebufferAfterRun );           
	   timebuffer.millitm = timebufferAfterRun.millitm - timebufferBeforeRun.millitm;
	   (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?timebuffer.millitm += 1000 : 0; 
	   fprintf(this->outAux,"%011fS%dMil ",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);
	   
	  fprintf(this->out,"\n\n\n[Algorithm execution information]");
	  fprintf(this->out,"\n-Time used to execute the runs: %.3fS%dMil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);
	  fprintf(this->out,"\n-Avg Iterations = %3.3f",(double)this->avgIte);
       
    #endif

    fflush(this->outAux);	

	//Create the file with the final clustering, each line corresponds to a pattern data plus its class
	long line=0;
	for (c=0;c<this->hardPartition.size();c++){
		for (i=0;i<this->hardPartition[c].size();i++){
			fprintf(this->outClusteringResult,"L%08d: %08d ",++line,c+1);
			for (int f=0;f<this->nFeatures;f++){
              fprintf(this->outClusteringResult,"%f ",(*this->tabData)[this->hardPartition[c][i]][f]);
			}//for f
			fprintf(this->outClusteringResult,"\n");
		}//for i
		
	}//for i
     fprintf(this->outClusteringResult,"Cluster  |  Size\n");
	for (c=0;c<this->hardPartition.size();c++){
		fprintf(this->outClusteringResult,"%08d %08d\n",c+1,this->hardPartition[c].size());
	}//for c
	fflush(this->outClusteringResult);

};//void FuzzyPartitioner::makeOutPutReport(){


void Partitioner::makeOutPutReportAboutRun(int run){
  fprintf(this->out,"\n==============================================================================================================================================================================================================");
  fprintf(this->out,"\nRun: %d",run);   
  fprintf(this->out,"\nIteration      |      J           |      T           |      R           |      B           |      RC          |  Non Fuzziness    |       Dk(U)'        |      Tp          |     Bp           |  Rp");
  fprintf(this->out,"\n---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
  fflush(this->out);
};

void Partitioner::makeOutPutReportAboutIteration(int iteration){
	             //0,1,2,3, 4, 5            , 6          , 7 ,  8, 9																
	 			// J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
  fprintf(this->out,"\n%14d   %16.8f   %16.8f   %16.8f   %16.8f   %16.8f    %16.8f      %16.8f   %16.8f   %16.8f   %16.8f",iteration,this->idxAndStatistic[0],this->idxAndStatistic[1],this->idxAndStatistic[2],this->idxAndStatistic[4],
	     this->idxAndStatistic[3],this->idxAndStatistic[5],this->idxAndStatistic[6],this->idxAndStatistic[7],this->idxAndStatistic[8],this->idxAndStatistic[9]);
  fflush(this->out);
};

//Pre-conditions:
// All the indices have to be already calculated

void Partitioner::makeInterpretationReport(){

	int i,j,m=1,d;
	double T=0,J=0,B=0;
  
	fprintf(this->outInte,"\n\n/-------------------------------------  Interpretation -------------------------------------/");
    fprintf(this->outInte,"\n\n/------------------------------  Global Inertia - T ------------------------------/");
    fprintf(this->outInte,"%-30s=%15.8f","\nData set global intertia",this->indexs[0]->getBestValueAssociationsValue(6) );

	//
   fprintf(this->outInte,"\n\nThe global inertia for each cluster" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   CLUSTER   |              VALUE      |"); 
   //                         X 012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   
   for (i=0;i<this->nClusters;i++){
	   fprintf(this->outInte,"\n| %9d      %17.8f      |",i+1,this->T_C[m][i]);
	   T += this->T_C[m][i];
   }//for i
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n  %-9s      %17.8f","Sum",T);

//
   fprintf(this->outInte,"\n\nThe global inertia for each variable" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   VARIABLE   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double Ti=0; T=0; 
   for (j=0;j<this->nFeatures;j++){
	   Ti=0;
	   for (i=0;i<this->nClusters;i++)
		   Ti += this->T_CV[m][i][j];
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",j+1,Ti);
	   T += Ti;
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",T);

//
   fprintf(this->outInte,"\n\nThe global inertia for each variable and cluster" );
   int dashes=18;/* the 21 spaces less the / and \ and final -*/

   for (i=0;i<this->nClusters;i++)
	   dashes += 19;

   fprintf(this->outInte,"\n\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

   fprintf(this->outInte,  "\n|   VAR\\CLUSTER    |"); //it takes 21 spaces to be shown
   //                             0123456789012345
   for (i=0;i<this->nClusters;i++){
     fprintf(this->outInte," CLUSTER %8d |",i+1);//it takes 19 spaces to be shown
   }//for i
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");

   fprintf(this->outInte,"\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

  
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|   VARIABLE %5d |",j+1);
	   for (i=0;i<this->nClusters;i++)
		   fprintf(this->outInte," %16.8f |",this->T_CV[m][i][j]);
	     
	    
   }//for j
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");
//should I put the sum of each element above?


//
    fprintf(this->outInte,"\n\n/------------------------------  Within-cluster Inertia - J ------------------------------/");
    fprintf(this->outInte,"%-30s=%15.8f","\nData set within-cluster intertia",this->indexs[0]->getBestValue() );

	//
   fprintf(this->outInte,"\n\nThe within-cluster inertia for each cluster" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   CLUSTER   |              VALUE      |"); 
   //                         X 012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   
   for (i=0;i<this->nClusters;i++){
	   fprintf(this->outInte,"\n| %9d      %17.8f      |",i+1,this->J_C[i]);
	   J += this->J_C[i];
   }//for i
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n  %-9s      %17.8f","Sum",J);

//
   fprintf(this->outInte,"\n\nThe within-cluster inertia for each variable" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   VARIABLE   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double Ji=0; J=0; 
   for (j=0;j<this->nFeatures;j++){
	   Ji=0;
	   for (i=0;i<this->nClusters;i++)
		   Ji += this->J_CV[i][j];
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",j+1,Ji);
	   J += Ji;
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",J);

//
   fprintf(this->outInte,"\n\nThe within-cluster inertia for each variable and cluster" );
   dashes=18;/* the 21 spaces less the / and \ */

   for (i=0;i<this->nClusters;i++)
	   dashes += 19;

   fprintf(this->outInte,"\n\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

   fprintf(this->outInte,  "\n|   VAR\\CLUSTER    |"); //it takes 21 spaces to be shown
   //                             0123456789012345
   for (i=0;i<this->nClusters;i++){
     fprintf(this->outInte," CLUSTER %8d |",i+1);//it takes 19 spaces to be shown
   }//for i
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");

   fprintf(this->outInte,"\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

  
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|   VARIABLE %5d |",j+1);
	   for (i=0;i<this->nClusters;i++)
		   fprintf(this->outInte," %16.8f |",this->J_CV[i][j]);
	     
	    
   }//for j
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");
//should I put the sum of each element above?


   //
    fprintf(this->outInte,"\n\n/------------------------------  Between-cluster Inertia - B ------------------------------/");
    fprintf(this->outInte,"%-30s=%15.8f","\nData set between-cluster intertia",this->indexs[0]->getBestValueAssociationsValue(7) );

	//
   fprintf(this->outInte,"\n\nThe between-cluster inertia for each cluster" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   CLUSTER   |              VALUE      |"); 
   //                         X 012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   
   for (i=0;i<this->nClusters;i++){
	   fprintf(this->outInte,"\n| %9d      %17.8f      |",i+1,this->B_C[m][i]);
	   B += this->B_C[m][i];
   }//for i
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n  %-9s      %17.8f","Sum",B);

//
   fprintf(this->outInte,"\n\nThe between-cluster inertia for each variable" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   VARIABLE   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double Bi=0; B=0; 
   for (j=0;j<this->nFeatures;j++){
	   Bi=0;
	   for (i=0;i<this->nClusters;i++)
		   Bi += this->B_CV[m][i][j];
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",j+1,Bi);
	   B += Bi;
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",B);

//
   fprintf(this->outInte,"\n\nThe between-cluster inertia for each variable and cluster" );
   dashes=18;/* the 21 spaces less the / and \ */

   for (i=0;i<this->nClusters;i++)
	   dashes += 19;

   fprintf(this->outInte,"\n\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

   fprintf(this->outInte,  "\n|   VAR\\CLUSTER    |"); //it takes 21 spaces to be shown
   //                             0123456789012345
   for (i=0;i<this->nClusters;i++){
     fprintf(this->outInte," CLUSTER %8d |",i+1);//it takes 19 spaces to be shown
   }//for i
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");

   fprintf(this->outInte,"\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

  
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|   VARIABLE %5d |",j+1);
	   for (i=0;i<this->nClusters;i++)
		   fprintf(this->outInte," %16.8f |",this->B_CV[m][i][j]);
	     
	    
   }//for j
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");
//should I put the sum of each element above?

  //
    fprintf(this->outInte,"\n\n/------------------------------  General Index - R ------------------------------/");
    fprintf(this->outInte,"%-60s=%15.8f","\n\nThe proportion of inertia explaines by the cluster",this->indexs[0]->getBestValueAssociationsValue(8) );
//
    fprintf(this->outInte,"\n\n/------------------------------  Variables contribution ------------------------------/");
//
   fprintf(this->outInte,"\n\nThe proportion of inertia of the variables taken into account by the cluster" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   VARIABLE   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double COR=0; 
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",j+1,this->COR_V[m][j]);
	   COR += this->COR_V[m][j];
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",COR);

//
   fprintf(this->outInte,"\n\nThe relative contribution of the variables to the betwwen-class inertia" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|   VARIABLE   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double CTR=0; 
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",j+1,this->CTR_V[m][j]);
	   CTR += this->CTR_V[m][j];
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",CTR);
//
    fprintf(this->outInte,"\n\n/------------------------------  Cluster description ------------------------------/");
//
   fprintf(this->outInte,"\n\nThe proportion of the global inertia explained by clusters" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|    CLUSTER   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double Taux=0; 
   for (i=0;i<this->nClusters;i++){
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",i+1,this->T_C[m][i]/T);
	   Taux += this->T_C[m][i]/T;
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",Taux);
//
   fprintf(this->outInte,"\n\nThe relative contribution of clusters to the between-cluster inertia" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|    CLUSTER   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double Baux=0; 
   for (i=0;i<this->nClusters;i++){
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",i+1,this->B_C[m][i]/B);
	   Baux += this->B_C[m][i]/B;
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",Baux);
//
   fprintf(this->outInte,"\n\nThe relative contribution of clusters to the within-cluster inertia" );
   fprintf(this->outInte,"\n\n/---------------------------------------\\");
   fprintf(this->outInte,  "\n|    CLUSTER   |              VALUE     |"); 
   //                         X  012345678      01234567890123456 
   fprintf(this->outInte,  "\n\\---------------------------------------/");
   fprintf(this->outInte,  "\n/---------------------------------------\\");
   double Jaux=0; 
   for (i=0;i<this->nClusters;i++){
	   fprintf(this->outInte,"\n|  %9d      %17.8f     |",i+1,this->J_C[i]/J);
	   Jaux += this->J_C[i]/J;
   }//for j
   fprintf(this->outInte, "\n\\---------------------------------------/");
   fprintf(this->outInte,"\n   %-9s      %17.8f","Sum",Jaux);
//
    fprintf(this->outInte,"\n\n/------------------------------  Cluster description by variables------------------------------/");
//   
   fprintf(this->outInte,"\n\nThe proportion of the discriminant power of variables taken into account by clusters" );
   dashes=18;/* the 21 spaces less the / and \ */

   for (i=0;i<this->nClusters;i++)
	   dashes += 19;

   fprintf(this->outInte,"\n\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

   fprintf(this->outInte,  "\n|   VAR\\CLUSTER    |"); //it takes 21 spaces to be shown
   //                             0123456789012345
   for (i=0;i<this->nClusters;i++){
     fprintf(this->outInte," CLUSTER %8d |",i+1);//it takes 19 spaces to be shown
   }//for i
   
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");

   fprintf(this->outInte,"\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

  
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|   VARIABLE %5d |",j+1);
	   for (i=0;i<this->nClusters;i++)
		   fprintf(this->outInte," %16.8f |",this->COR_CV[m][i][j]);
	     
	    
   }//for j
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");
//should I put the sum of each element above?

//   
   fprintf(this->outInte,"\n\nThe relative contribution of variables to the between-cluster inertia explained by clusters" );
   dashes=18;/* the 21 spaces less the / and \ */

   for (i=0;i<this->nClusters;i++)
	   dashes += 19;

   fprintf(this->outInte,"\n\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

   fprintf(this->outInte,  "\n|   VAR\\CLUSTER    |"); //it takes 21 spaces to be shown
   //                             0123456789012345
   for (i=0;i<this->nClusters;i++){
     fprintf(this->outInte," CLUSTER %8d |",i+1);//it takes 19 spaces to be shown
   }//for i
   
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");

   fprintf(this->outInte,"\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

  
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|   VARIABLE %5d |",j+1);
	   for (i=0;i<this->nClusters;i++)
		   fprintf(this->outInte," %16.8f |",this->CTR_CV[m][i][j]);
	     
	    
   }//for j
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");
//should I put the sum of each element above?
//   
   fprintf(this->outInte,"\n\nThe relative contribution of variables and clusters to the between-cluster inertia" );
   dashes=18;/* the 21 spaces less the / and \ */

   for (i=0;i<this->nClusters;i++)
	   dashes += 19;

   fprintf(this->outInte,"\n\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

   fprintf(this->outInte,  "\n|   VAR\\CLUSTER    |"); //it takes 21 spaces to be shown
   //                             0123456789012345
   for (i=0;i<this->nClusters;i++){
     fprintf(this->outInte," CLUSTER %8d |",i+1);//it takes 19 spaces to be shown
   }//for i
   
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");

   fprintf(this->outInte,"\n/");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"\\");

  
   for (j=0;j<this->nFeatures;j++){
	   fprintf(this->outInte,"\n|   VARIABLE %5d |",j+1);
	   for (i=0;i<this->nClusters;i++)
		   fprintf(this->outInte," %16.8f |",this->CE_CV[m][i][j]);
	     
	    
   }//for j
   fprintf(this->outInte,"\n\\");
   for (d=0;d<dashes;d++)
     fprintf(this->outInte,"-");
   fprintf(this->outInte,"/");
//should I put the sum of each element above?
};

StatisticalSummary ** Partitioner::getStatisticalIndexs(){
  return this->indexs;
};


// J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
double Partitioner::getBestJ(){
	return this->indexs[0]->getBestValue();
};

double Partitioner::getWorstJ(){
	return this->indexs[0]->getWorstValue();
};

double Partitioner::getRCForBestJ(){
	return this->indexs[0]->getBestValueAssociationsValue(0);
};

double Partitioner::getBestRC(){
	return this->indexs[3]->getBestValue();
};

double Partitioner::getJForBestRC(){
	return this->indexs[3]->getBestValueAssociationsValue(0);
};

double Partitioner::getWorstRC(){
	return this->indexs[3]->getWorstValue();
};

vector<double> Partitioner::getBestJAssociatedValues(){
  return this->indexs[0]->getBestValueAssociationsValues();
};

vector<string> Partitioner::getBestJAssociatedValuesDescription(){
  return this->indexs[0]->getBestValueAssociationsDescription();
};

double Partitioner::getAvgIte(){
  return this->avgIte;
};

//pre-conditions: - the hardPartition structure has to contain a valid hard partition, createHardPartition is able to create such a structure
double Partitioner::getAdjustedRandIndex(vector<tVecInt> &hardPartition){
  if ( (this->labelFeature==-1) && (this->labelFeatureData.size()==0) ){
		return -1;
  }
  vector<int>** part;
  part = (vector<int> **)malloc(sizeof(vector<int>* )*this->nClusters);  
  int i;
  for (i=0;i<this->nClusters;i++)
	  part[i] = new vector<int>;
  for (int cc=0;cc<hardPartition.size();cc++){//for 1
	  //int s = hardPartition[cc].size();
	  for(int i=0;i<hardPartition[cc].size();i++){//for 2
		  part[cc]->push_back(hardPartition[cc][i]);
	  }//for i
  }//for cc
  double ret = this->aRand->getAdjustedRandIndex(this->nClusters,part);
  for (i=0;i<this->nClusters;i++){
	  part[i]->clear();
	  delete part[i];
  }
  free(part);
  return ret;
};



void Partitioner::startUp(){
	int i,j,m;

	
	this->avgIte = 0;

	this->hardPartition.resize(this->nClusters);
 	this->centerCluster.resize(this->nClusters);
	this->centerClusterForBestJ.resize(this->nClusters);
	this->centerClusterInitialForBestJ.resize(this->nClusters);
	this->centerClusterTmp.resize(this->nClusters);
	for (int c=0;c<this->nClusters;c++){
		this->centerCluster[c].resize(this->nFeatures,0.0);
		this->centerClusterForBestJ[c].resize(this->nFeatures,0.0);
	    this->centerClusterInitialForBestJ[c].resize(this->nFeatures,0.0);
	    this->centerClusterTmp[c].resize(this->nFeatures,0.0);
	}//for c
	this->globalCenter.clear();
    this->globalCenter.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->globalCenter[0].resize(this->nFeatures,0.0);
	this->globalCenter[1].resize(this->nFeatures,0.0);

    //Indexs:
	// J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
	// 0,1,2, 3, 4,      5,            6,        7, 8, 9
	char *desc[]= {"Objective Function","T","R","Adjustaded Rand Index","B","Non Fuzziness","Dk(U) Prime","T'","B'","R'" };
	bool isMax[] = {false,true,true,true,true,true,false,true,true,true};
	//it should be implemented as a vector later in order to allow subclasses to add their own indices
	this->indexs = (StatisticalSummary **)malloc(sizeof(StatisticalSummary *)*NUMBER_INDEX_OBSERVED);
	//It remain the indexs and statistics more recently computed
	this->idxAndStatistic = (double *)malloc(sizeof(double)*NUMBER_INDEX_OBSERVED);
    for (i =0;i<NUMBER_INDEX_OBSERVED;i++){
		this->indexs[i] = new StatisticalSummary(desc[i],isMax[i]);
		this->idxAndStatistic[i] = -1;
	}//for i
	//RC for best J
		  this->indexs[0]->setBestValueAssociationDescription(0,this->indexs[3]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(0,this->idxAndStatistic[3]);
		  //T for best J
		  this->indexs[0]->setBestValueAssociationDescription(1,this->indexs[1]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(1,this->idxAndStatistic[1]);
		  //R for best J
		  this->indexs[0]->setBestValueAssociationDescription(2,this->indexs[2]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(2,this->idxAndStatistic[2]);
		  //B for best J
		  this->indexs[0]->setBestValueAssociationDescription(3,this->indexs[4]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(3,this->idxAndStatistic[4]);
		  //Non Fuzziness for best J
		  this->indexs[0]->setBestValueAssociationDescription(4,this->indexs[5]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(4,this->idxAndStatistic[5]);
		  //Dk(U) Prime for best J
		  this->indexs[0]->setBestValueAssociationDescription(5,this->indexs[6]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(5,this->idxAndStatistic[6]);
		  //Tp for best J
		  this->indexs[0]->setBestValueAssociationDescription(6,this->indexs[7]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(6,this->idxAndStatistic[7]);
		  //Bp for best J
		  this->indexs[0]->setBestValueAssociationDescription(7,this->indexs[8]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(7,this->idxAndStatistic[8]);
		  //Rp for best J
		  this->indexs[0]->setBestValueAssociationDescription(8,this->indexs[9]->getValueDescription());
		  //this->indexs[0]->setBestValueAssociationsValue(8,this->idxAndStatistic[9]);
	
	//treats the output
	if (this->fileOutName.size()>0){
	  this->out = fopen(this->fileOutName.c_str(),"w");
	  string fileOutName = this->fileOutName;
	  int pos;
	  string strSuf;
	  string tmpFile = fileOutName;
	  strSuf.assign("_AUX_DATA.wri");
	  this->outAux = fopen(fileOutName.replace((pos = fileOutName.find_last_of(".",fileOutName.length()))==string::npos?fileOutName.length():pos,strSuf.size(),strSuf).c_str(),"w");			    
	  fileOutName = tmpFile;
	  strSuf.assign("_ClusterResult.wri");
	  this->outClusteringResult= fopen(fileOutName.replace((pos = fileOutName.find_last_of(".",fileOutName.length()))==string::npos?fileOutName.length():pos,strSuf.size(),strSuf).c_str(),"w");			    
	  fileOutName = tmpFile;
	  strSuf.assign("_Interpretation.wri");
	  this->outInte = fopen(fileOutName.replace((pos = fileOutName.find_last_of(".",fileOutName.length()))==string::npos?fileOutName.length():pos,strSuf.size(),strSuf).c_str(),"w");			    
	  fileOutName = tmpFile;
	  strSuf.assign("_CustomReport.wri");
	  this->outCustomReport = fopen(fileOutName.replace((pos = fileOutName.find_last_of(".",fileOutName.length()))==string::npos?fileOutName.length():pos,strSuf.size(),strSuf).c_str(),"w");			    
	  fileOutName = tmpFile;
	  strSuf.assign("_DataStatistics.wri");
	  this->outDataStatistics = fopen(fileOutName.replace((pos = fileOutName.find_last_of(".",fileOutName.length()))==string::npos?fileOutName.length():pos,strSuf.size(),strSuf).c_str(),"w");			    

	}
	//Create the class which will calculate RC
    this->aRand = new AdjustedRandIndex(this->labelFeatureData,this->nClasses) ;
	//Indices
	this->B_CV.clear();
	this->T_CV.clear();
	this->J_CV.clear();
	this->CE_CV.clear();
	this->COR_CV.clear();
	this->CTR_CV.clear();
    // each indice will be calculated according to 2 different kinds of global prototype
	this->B_CV.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->T_CV.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->CE_CV.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->COR_CV.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->CTR_CV.resize(NUM_OF_GLOBAL_PROTOTYPE);

    for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){

		this->B_CV[m].resize(this->nClusters);
	    this->T_CV[m].resize(this->nClusters);
	    this->CE_CV[m].resize(this->nClusters);
	    this->COR_CV[m].resize(this->nClusters);
	    this->CTR_CV[m].resize(this->nClusters);

	}//for m
	this->J_CV.resize(this->nClusters);

	for (i=0;i<this->nClusters;i++){
	   for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
		    this->B_CV[m][i].resize(this->nFeatures,0.0);	        
			this->T_CV[m][i].resize(this->nFeatures,0.0);
	        this->CE_CV[m][i].resize(this->nFeatures,0.0);
	        this->COR_CV[m][i].resize(this->nFeatures,0.0);
	        this->CTR_CV[m][i].resize(this->nFeatures,0.0);
	   }//for m
    	this->J_CV[i].resize(this->nFeatures,0.0);		
	}// for i

	this->B_C.clear();
	this->T_C.clear();
	this->J_C.clear();
	this->COR_V.clear();
	this->CTR_V.clear();
    // each indice will be calculated according to 2 different kinds of global prototype
	this->B_C.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->T_C.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->COR_V.resize(NUM_OF_GLOBAL_PROTOTYPE);
	this->CTR_V.resize(NUM_OF_GLOBAL_PROTOTYPE);
    for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){

		this->B_C[m].resize(this->nClusters,0);
	    this->T_C[m].resize(this->nClusters,0);
	    this->COR_V[m].resize(this->nFeatures,0);
	    this->CTR_V[m].resize(this->nFeatures,0);

	}//for m
    this->J_C.resize(this->nClusters,0);


};


void Partitioner::cleanUp(){
	this->bDoCleanUp = false;
	this->centerCluster.clear();
		this->hardPartition.clear();
	this->globalCenter.clear();
	for (int i =0;i<NUMBER_INDEX_OBSERVED;i++){
		delete this->indexs[i];
	}//for i
    
	this->B_CV.clear();
	this->T_CV.clear();
	this->J_CV.clear();
	this->CE_CV.clear();
	this->COR_CV.clear();
	this->CTR_CV.clear();

	this->B_C.clear();
	this->T_C.clear();
	this->J_C.clear();
	this->COR_V.clear();
	this->CTR_V.clear();

/*
	for (int c=0;c<this->nClusters;c++){
      free(this->membershipRaisedToM[c]);
	}//for c

    free(this->membershipRaisedToM);


	*/

/*
	
	for (i=0;i<this->nIndividuals;i++){
		
		for (int c=0;c<this->nClusters;c++){
			free(this->difIndFromCenter[i][c]);
		}//for f
		free(this->difIndFromCenter[i]);
	}//for i
	free(this->difIndFromCenter);
	*/

	
	delete this->aRand;
	this->aRand = NULL;
	free(this->indexs);
	free(this->idxAndStatistic);
	this->indexs = NULL;
	if (this->fileOutName.size()>0){
	  fclose(this->out);
	  fclose(this->outAux);
	  fclose(this->outClusteringResult);
	  fclose(this->outInte);
	  fclose(this->outCustomReport);
	  fclose(this->outDataStatistics);
	}
};


void Partitioner::doCleanUp(){
  this->cleanUp();
};

//////////////////////////////////////////////////////////////////////
// FuzzyPartitioner methods
//////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////
// FuzzyPartitioner Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyPartitioner::FuzzyPartitioner():Partitioner(){
	this->bDoCleanUp = true;
	this->membership.clear();
	//this->membershipRaisedToM.clear();
	this->bestMembershipJ.clear();
	this->bestInitialMembershipJ.clear();
	this->tmpInitialMembershipJ.clear();
	this->globalCenter.clear();
	this->centerCluster.clear();
	this->hardPartition.clear();
	this->indexs = NULL;
	this->aRand = NULL;
    this->idxAndStatistic = NULL;
    
	this->membershipRaisedToM = NULL;
	this->difIndFromCenter = NULL;
	
};


FuzzyPartitioner::~FuzzyPartitioner(){
	int i;
    
	for (int c=0;c<this->nClusters;c++){
      free(this->membershipRaisedToM[c]);
	}//for c
	free(this->membershipRaisedToM);
	
	for (i=0;i<this->nIndividuals;i++){
		
		for (int c=0;c<this->nClusters;c++){
			free(this->difIndFromCenter[i][c]);
		}//for f
		free(this->difIndFromCenter[i]);
	}//for i
	free(this->difIndFromCenter);	
};


//////////////////////////////////////////////////////////////////////
// Other FuzzyPartitioner methods
//////////////////////////////////////////////////////////////////////

//pre-conditions:-the number of clusters has to be known, - the membership matrix has to be initialized 
void FuzzyPartitioner::createHardPartition(vector<tVecDouble> &membership,vector<tVecInt> &hardPartiion){
	for (int c=0;c<hardPartition.size();c++)
	  hardPartition[c].clear();
	for (int i=0;i<this->nIndividuals;i++){
		double maxMembership=-HUGE_VAL;
		int cluster = 0;
		for (int c=0;c<this->nClusters;c++){
			if (maxMembership<membership[c][i]){
				cluster = c;
				maxMembership = membership[c][i];
			}			
		}//for c
		hardPartition[cluster].push_back(i);
	}//for i
#ifdef DEBUG_CREATE_HARD_PARTITION
	for (c=0;c<this->nClusters;c++){
		
		if (this->hardPartition[c].size()<1)
			cout<<endl<<"Cluster "<<c+1<<" vazio!";
		else{
			cout<<endl<<"Cluster "<<c+1<<"-> ";
			for (int k=0;k<this->hardPartition[c].size();k++){
				cout<<this->hardPartition[c][k]<<", ";
			}//for k
			cout<<endl;
		}
	}//for c
#endif
};




void FuzzyPartitioner::startUp(){
	int i,j,m;

	bool bAllocate = false;

	Partitioner::startUp();

	this->membership.resize(this->nClusters);
	//this->membershipRaisedToM.resize(this->nClusters);
	if (this->membershipRaisedToM==NULL){
	  this->membershipRaisedToM = (double **)malloc(sizeof(double *)*this->nClusters);
	  bAllocate = true;
	}
	this->bestMembershipJ.resize(this->nClusters);
	this->bestInitialMembershipJ.resize(this->nClusters);
	this->tmpInitialMembershipJ.resize(this->nClusters);
	for (m=0;m<this->nClusters;m++){
		this->membership[m].resize(this->nIndividuals,0.0);
		this->bestMembershipJ[m].resize(this->nIndividuals,0.0);
        this->bestInitialMembershipJ[m].resize(this->nIndividuals,0.0);
		this->tmpInitialMembershipJ[m].resize(this->nIndividuals,0.0);
		if (bAllocate)
		   this->membershipRaisedToM[m] = (double *)malloc(sizeof(double)*this->nIndividuals);
		//this->membershipRaisedToM[m].resize(this->nIndividuals,0.0);
	}//for m

    if (this->difIndFromCenter==NULL){
		this->difIndFromCenter = (double ***)malloc(sizeof(double **)*this->nIndividuals);
		for (i=0;i<this->nIndividuals;i++){
			this->difIndFromCenter[i] = (double **)malloc(sizeof(double *)*this->nClusters);
			for (int c=0;c<this->nClusters;c++){
				this->difIndFromCenter[i][c] = (double *)malloc(sizeof(double)*this->nFeatures);
			}//for f
		}//for i
	}


   
};




double FuzzyPartitioner::getMembershipRaisedToM(int c,int i){
  return this->membershipRaisedToM[c][i];
};


//pre-conditions: parameterM must be setted
void  FuzzyPartitioner::initializeMembership(){
  
  
  
  vector<bool> choosen;
  vector<bool> clusterChoosen;
  choosen.resize(this->nIndividuals);
  clusterChoosen.resize(this->nClusters);
  double mults[6] = {1.5,1.4,1.3,1.2,1.1,1};
  bool doOrNot[2] = {true,false};
  int w;
  int somaPercentual = 0;
  double aux,piece,reservePiece=0;
  int i,c;
  int timesI=0,timesC=0;
  double cota, resto=0;

  randomize(); 

  for (w=0;w<choosen.size();w++)
	choosen[w] = false;
  
  for (w=0;w<clusterChoosen.size();w++)
	clusterChoosen[w] = false;
  
  
  cota = (double)this->nIndividuals/this->nClusters;
  
  
  c=0;
  
  for (int clu=0;clu<this->membership.size();clu++){
	double ind = cota +resto;	
	int idx=0;	
	int k = ind;
	timesC++;
	do{
	  c = rand()%this->nClusters;//random(this->nClusters);
	}while (clusterChoosen[c]);
	clusterChoosen[c]=true;
	if (timesC==this->membership.size())
	  k = this->nIndividuals-timesI;
	for (;k>0;k--){	  
	  timesI++;
	  do{
		idx = rand()%this->nIndividuals;//random(this->nIndividuals);
	  }while (choosen[idx]);
	  choosen[idx]=true;
	
		double minPiece =0;
		//it makes minPiece grows til minPiece be the major part of 100 but not equals to 100
		while ( (minPiece<=(100-minPiece)) || (minPiece>=100) ){
		  minPiece+=(100/this->nClusters)*mults[rand()%6];//random(6)];		
		  if (minPiece>=100)
			  minPiece -=(100/this->nClusters)*mults[rand()%6];//random(6)];		
		}//while
        int id1 = rand()%6;//random(6);
		double id2 = rand()%100;//random(100);
		if (doOrNot[rand()%2]){//random(2)]){
		  double inc = mults[rand()%6]/mults[rand()%6]*(id2/100)*(1.0/mults[id1])*((100-minPiece)/this->nClusters);
		  if (minPiece+inc<100)
		    minPiece += inc; 
		}
        //this membership is surely positive
		membership[c][idx] = minPiece/100;
		this->membershipRaisedToM[c][idx] = pow(membership[c][idx],this->parameterM); 
		//sorteia a pertinencia entre os outros clusters
		double sum=0;
		reservePiece=0;
		piece = (100-minPiece)/(this->nClusters-1);
		int i=0;		
		for (;i<membership.size();i++){   
		  if ( (c==membership.size()-1) && (i==membership.size()-2))
			break;
		  if (i==membership.size()-1)
			break;
		  if (i==c)
			continue;
		  //aux will be at least 0+1=1 and at most piece-1+1=piece
		  aux = (rand()%(long)(piece*100))/100+1;//random(piece*100)/100+1;
		  if (aux>piece)
			  aux=piece;
		  //aux = ( (aux+somaPercentual)>100?100-somaPercentual:aux);
		  //somaPercentual += aux;	  
		  membership[i][idx] = (aux+reservePiece)/100;
		  this->membershipRaisedToM[i][idx] = pow(membership[i][idx],this->parameterM);
		  reservePiece = (aux<piece)?piece-aux:0.0;	  
		  sum += membership[i][idx];		  
		  if (membership[i][idx]<0)
			  cerr<<endl<<"Error1 - Membership must not be negative - initializeMembership()";
		}//for i			
		membership[i][idx] = ((100-minPiece)/100)-sum;		
		this->membershipRaisedToM[i][idx] = pow(membership[i][idx],this->parameterM);
		if (membership[i][idx]<0)
			  cerr<<endl<<"Error2 - Membership must not be negative - initializeMembership()";
		
	}//for k
	resto = ind-(long)ind;	
  } 
  
  #ifdef DEBUG_MEMBERSHIP
    cout<<endl<<endl<<"Membership Matrix"<<endl;
	cout<<"Individual -> Memberships to C1 | Cn : Sum"<<endl;
	for (i=0;i<this->nIndividuals;i++){
		cout<<(i+1)<<" -> ";
		double sum=0;
		for (c=0;c<this->nClusters;c++){
			cout<<this->membership[c][i]<<" | ";
			sum += this->membership[c][i];
		}//for c
		cout<<" : "<<sum<<endl;
	}//for i
  #endif


//DEBUG
#ifdef DEBUG_MEMBERSHIP_INTIALIZATION
	if (this->fileOutName.size()>0){
	int i;
fprintf(this->out,"\n\n------------------- Membership matrix initial ------------------- ");
    
	fprintf(this->out,"\n Individuali");

	for (c=0;c<this->nClusters;c++){
      fprintf(this->out,"\t%10s%02d","Mi",c+1);
	}//for c
    fprintf(this->out,"\t%12s","Hard Cluster");
	if (this->labelFeature!=-1){
		fprintf(this->out,"\t%12s","A priori partition");
	}//if

	for (i=0;i<this->nIndividuals;i++){
		fprintf(this->out,"\n%12d",i+1);
		double maxM=MIN_DOUBLE;
		int clustHard=0;
		for (int c=0;c<this->nClusters;c++){
			fprintf(this->out,"\t%12.8f",this->membership[c][i]);
			if (maxM<this->membership[c][i]){
				maxM = this->membership[c][i];
				clustHard = c;
			}
		}//for c
		fprintf(this->out,"\t%12d",clustHard+1);
		if (this->labelFeature!=-1){
			fprintf(this->out,"\t%12.0f",this->labelFeatureData[i]);
		}//if
	}//for i

	}//if (this->fileOutName.size()>0)
#endif
  
   //copy the initial membership matrix into a tmp strucutre
  for (i=0;i<this->nIndividuals;i++){
	for (c=0;c<this->nClusters;c++){
	  this->tmpInitialMembershipJ[c][i] = this->membership[c][i];
	}//for c
  }//for i

};

double FuzzyPartitioner::objectiveFunction(){
	double ret=0;
	for (int i=0;i<this->nClusters;i++){	
		for (int k=0;k<this->nIndividuals;k++){
			double dist=0;
			for (int f=0;f<this->nFeatures;f++){
			  	dist += this->difIndFromCenter[k][i][f];
			}//for f
			//dist = pow(dist,0.5);
			//ret += pow(this->membership[i][k],this->parameterM)*pow(dist,2.0);
			ret += this->membershipRaisedToM[i][k]*dist;
		}//for k
	}//for i
	return ret;
};//double FuzzyCMeans::objectiveFunction();

void FuzzyPartitioner::workOutBTRPrime(double &Bp,double &Tp, double &Rp){


	double retRP=0,retBP=0,retTP=0;
	int i,k,p;

//It calculates, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += pow((*this->tabData)[k][p]-this->globalCenter[1][p],2.0);
		  }//for p
		  retTP += this->membershipRaisedToM[i][k]*dist;		  
	  }//for k
  }//for i

    //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += pow(this->centerCluster[i][p]-this->globalCenter[1][p],2.0);
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

//c = cluster
//p=variable
double FuzzyPartitioner::getAdaptiveValue(int &c,int &p){
	return 1.0;
};//double FuzzyPartitioner::getAdaptiveValue(int &i,int &c,int &p){
//Pre-conditions:
// The best membership must be already defined
void FuzzyPartitioner::workOutCVIndices(){
	
	int i,j,m;

	double *B;
	B = (double *)malloc(sizeof(double)*NUM_OF_GLOBAL_PROTOTYPE);

	//T_C, J_C,B_C
    for (i=0;i<this->nClusters;i++){
       for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
		   this->B_C[m][i] = 0;
		   this->T_C[m][i] = 0;
		   for (j=0;j<this->nFeatures;j++){
			   this->T_CV[m][i][j] = 0;
			   this->B_CV[m][i][j] = 0;			   
		   }
	   }//for m
	   for (j=0;j<this->nFeatures;j++)
		   this->J_CV[i][j] = 0;
	   this->J_C[i] = 0;
	}//for i

	for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++)
      B[m] = 0;  

	//T_CV,J_CV,B_CV
	for (i=0;i<this->nClusters;i++){
		double Ui=0;
		for (int k=0;k<this->nIndividuals;k++){
			Ui += this->membershipRaisedToM[i][k];
		}//for k
		for (j=0;j<this->nFeatures;j++){
          for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
            for (int k=0;k<this->nIndividuals;k++){
			   this->T_CV[m][i][j] += this->getAdaptiveValue(i,j)*this->membershipRaisedToM[i][k]*pow( ((*this->tabData)[k][j]-this->globalCenter[m][j]),2.0);
			}//for k 
            this->T_C[m][i] += this->T_CV[m][i][j];
		  }//for m
		  for (int k=0;k<this->nIndividuals;k++)
		    this->J_CV[i][j] += this->getAdaptiveValue(i,j)*this->membershipRaisedToM[i][k]*pow( ((*this->tabData)[k][j]-this->centerCluster[i][j]),2.0);
		  //J_C
		  this->J_C[i] += this->J_CV[i][j];		  
		  for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
		    this->B_CV[m][i][j] = this->getAdaptiveValue(i,j)*Ui*pow( (this->centerCluster[i][j]-this->globalCenter[m][j]),2.0);			
			//B_C
			this->B_C[m][i] += this->B_CV[m][i][j];
		  }//for m
		}//for j
        
		for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++)
		  B[m] += this->B_C[m][i];
	}//for i

	for (j=0;j<this->nFeatures;j++){
	    for (m=0;m<NUM_OF_GLOBAL_PROTOTYPE;m++){
          double Tj=0,Bj=0;
	      for (i=0;i<this->nClusters;i++){
			Tj += this->T_CV[m][i][j];
			Bj += this->B_CV[m][i][j];
			
			//CTR(j,i)
			this->CTR_CV[m][i][j] = this->B_CV[m][i][j]/this->B_C[m][i];
			//CE(j,i) 
			this->CE_CV[m][i][j] = this->B_CV[m][i][j]/B[m];
		  }//for i 
		  for (i=0;i<this->nClusters;i++)		
		    this->COR_CV[m][i][j] = this->B_CV[m][i][j]/Tj;
          
		  this->COR_V[m][j] = Bj/Tj;
		  this->CTR_V[m][j] = Bj/B[m];
		}//for m				
		
	}//for j
	
	free(B);
	

};//void FuzzyPartitioner::workOutCVIndices(){



void FuzzyPartitioner::upDateLambda(){
};

void FuzzyPartitioner::saveLambdaForBestJ(){
};

void FuzzyPartitioner::saveInitialCenterCluster(){
	
	int c,f;

	for (c=0;c<this->nClusters;c++){
		for (f=0;f<this->nFeatures;f++){
			this->centerClusterTmp[c][f] = this->centerCluster[c][f];
		}//for f
	}//for c
};

void FuzzyPartitioner::execute(vector<tVecDouble>&tabData,vector<double> &labelFeatureData,int nClasses,int labelFeature,string fileInName,string fileOutName,int runs,int iterations,int clusters, double parameterM,bool bAllwaysCalculateIndexs,bool bWillUsed){
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
		this->initializeMembership();
		bConvergence = false;
		oldJ=MAX_DOUBLE;
        #ifndef NO_RUN_REPORT
		  if (this->fileOutName.size()>0)
		    this->makeOutPutReportAboutRun(run+1);
        #endif
		for (ite=0;ite<iterations && !bConvergence;ite++){
          
		  //core elements of the clustering process
		  this->upDateCenterClustersAndRelated();

		  if (ite==0)
			  this->saveInitialCenterCluster();

		  //core of adaptive methods
		  this->upDateLambda();

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
		  if (!bConvergence){
		    this->upDateMembershipAndRelated();
		  }
		  NumIte++;
          //cout<<endl<<"****Done iteration "<<(ite+1)<<"\tJ="<<J<<" R="<<R<<" T="<<T<<" B="<<B<<endl;
		}//for ite

        //As the algorithm makes J decrease over the iterations, J now has the smallest value
       // if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){
		
        if ( (bestJ>this->indexs[0]->getBestValue()) ||  (run==0) ){
		  this->saveMembershipForBestJ();

#ifndef NO_INTERPRETATION_INDICE
		  this->workOutCVIndices();
#endif
		  this->saveLambdaForBestJ();
          if (run==0)
			bestJ = this->idxAndStatistic[0];//bestJ has not been initialized yet
		  else 
		    bestJ = this->indexs[0]->getBestValue();
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

};//method FuzzyPartitioner::execute

void FuzzyPartitioner::workOutIndexAndStatistics(bool bCalculateJustJ){

		  this->idxAndStatistic[0] = this->objectiveFunction();
		  
		  if (bCalculateJustJ){              
			  //this->indexs[0]->addValue(this->idxAndStatistic[0]);
			  return;
		  }

		  


		  this->workOutBTR(this->idxAndStatistic[4],this->idxAndStatistic[1],this->idxAndStatistic[2],this->idxAndStatistic[5],this->idxAndStatistic[6]);

#ifndef NO_INTERPRETATION_INDICE
		  this->workOutBTRPrime(this->idxAndStatistic[8],this->idxAndStatistic[7],this->idxAndStatistic[9]);
#else
    this->idxAndStatistic[8]=this->idxAndStatistic[7]=this->idxAndStatistic[9]=-1;
#endif
		  this->createHardPartition(this->membership,this->hardPartition);

		  this->idxAndStatistic[3] = this->getAdjustedRandIndex(this->hardPartition);

          for (int i=0;i<NUMBER_INDEX_OBSERVED;i++){
            this->indexs[i]->addValue(this->idxAndStatistic[i]);
		  }//for i
          this->indexs[3]->setBestValueAssociationsValue(0,this->idxAndStatistic[0]);
          // T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
          //indexs in best J associated values
		  // 1,2,0 ,3,   4          ,  5         , 6 ,  7, 8
		  //RC for best J
		  //this->indexs[0]->setBestValueAssociationDescription(0,this->indexs[3]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(0,this->idxAndStatistic[3]);
          //Non Fuzziness for best J
		  //this->indexs[0]->setBestValueAssociationDescription(4,this->indexs[5]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(4,this->idxAndStatistic[5]);
		  //Dk(U) Prime for best J
		  //this->indexs[0]->setBestValueAssociationDescription(5,this->indexs[6]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(5,this->idxAndStatistic[6]);

#ifndef NO_INTERPRETATION_INDICE
		  //T for best J
		  //this->indexs[0]->setBestValueAssociationDescription(1,this->indexs[1]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(1,this->idxAndStatistic[1]);
		  //R for best J
		  //this->indexs[0]->setBestValueAssociationDescription(2,this->indexs[2]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(2,this->idxAndStatistic[2]);
		  //B for best J
		  //this->indexs[0]->setBestValueAssociationDescription(3,this->indexs[4]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(3,this->idxAndStatistic[4]);
		  
		  //Tp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(6,this->indexs[7]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(6,this->idxAndStatistic[7]);
		  //Bp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(7,this->indexs[8]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(7,this->idxAndStatistic[8]);
		  //Rp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(8,this->indexs[9]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(8,this->idxAndStatistic[9]);
#else
		  //in order to make the report look nice
		  //T for best J
		  //this->indexs[0]->setBestValueAssociationDescription(1,this->indexs[1]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(1,-1);
		  //R for best J
		  //this->indexs[0]->setBestValueAssociationDescription(2,this->indexs[2]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(2,-1);
		  //B for best J
		  //this->indexs[0]->setBestValueAssociationDescription(3,this->indexs[4]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(3,-1);
		  
		  //Tp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(6,this->indexs[7]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(6,-1);
		  //Bp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(7,this->indexs[8]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(7,-1);
		  //Rp for best J
		  //this->indexs[0]->setBestValueAssociationDescription(8,this->indexs[9]->getValueDescription());
		  this->indexs[0]->setBestValueAssociationsValue(8,-1);
#endif
		  
};//void FuzzyPartitioner::workOutIndexAndStatistics(){

//pre-condition: must be called after workOutIndexAndStatistics
void FuzzyPartitioner::saveMembershipForBestJ(){
//static double bestJ = MAX_DOUBLE;
	//The test is made in execute
		  //if (this->idxAndStatistic[0]<this->indexs[0]->getPreivouslyLastValueAdded()){
			//bestJ = this->idxAndStatistic[0];
			for (int c=0;c<this->nClusters;c++){
				for (int i=0;i<this->nIndividuals;i++){
					this->bestMembershipJ[c][i] = this->membership[c][i];
					this->bestInitialMembershipJ[c][i] = this->tmpInitialMembershipJ[c][i];					
				}//for i
				for (int f=0;f<this->nFeatures;f++){
					this->centerClusterForBestJ[c][f] = this->centerCluster[c][f];
					this->centerClusterInitialForBestJ[c][f] = this->centerClusterTmp[c][f];
				}//for f

			}//for c

		  //}//if

};





void FuzzyPartitioner::makeOutPutReportHeader(string fileInName,int runs, int iterations){
	fprintf(this->out,"Results:\n");
    fprintf(this->out,"\n%28s = %s","Algorithm",this->algDescription.c_str());
	if (fileInName.size()>0)
	  fprintf(this->out,"\n%28s = %s","Input data",fileInName.c_str());
	fprintf(this->out,"\n%28s = %f","ParameterM",this->parameterM);
	fprintf(this->out,"\n%28s = %d","Clusters",this->nClusters);	
	fprintf(this->out,"\n%28s = %d\n%28s = %d","Run Number",runs,"Maximun Iteration Number",iterations);
	fflush(this->out);
    //put this information in the interpretation output file
	fprintf(this->outInte,"Results:\n");
    fprintf(this->outInte,"\n%28s = %s","Algorithm",this->algDescription.c_str());
	if (fileInName.size()>0)
	  fprintf(this->outInte,"\n%28s = %s","Input data",fileInName.c_str());
	fprintf(this->outInte,"\n%28s = %f","ParameterM",this->parameterM);
	fprintf(this->outInte,"\n%28s = %d","Clusters",this->nClusters);	
	fprintf(this->outInte,"\n%28s = %d\n%28s = %d","Run Number",runs,"Maximun Iteration Number",iterations);
	fflush(this->outInte);

};

void FuzzyPartitioner::makeOutPutCustomReport(){

	Partitioner::makeOutPutCustomReport();

	this->printMembershipMatrix(this->bestInitialMembershipJ,"Initial Membership matrix for the minimal J obtained",this->outCustomReport);

	this->printHardPartition(this->bestInitialMembershipJ,"Hard partition generated from Initial Membership matrix for the minimal J obtained",this->outCustomReport);

	this->printCentroids(this->centerClusterInitialForBestJ,"Center Cluster for the Initial Membership",this->outCustomReport);

	this->printIndividualDistributionMatrix(this->hardPartition,"Individuals Distribution into clusters according to Initial Membership Matrix",this->outCustomReport);

    //Information related to Best J
	this->printMembershipMatrix(this->bestMembershipJ,"Membership matrix for the minimal J obtained",this->outCustomReport);

	this->printHardPartition(this->bestMembershipJ,"Hard partition obtained for the minimal J obtained",this->outCustomReport);

	this->printCentroids(this->centerClusterForBestJ,"Center Cluster for the minimal J obtained",this->outCustomReport);

	this->printIndividualDistributionMatrix(this->hardPartition,"Individuals Distribution into clusters according to Final Membership Matrix",this->outCustomReport);
};


void FuzzyPartitioner::printMembershipMatrix(vector<tVecDouble> membership,string Header,FILE *out){

	int i,c;

	fprintf(out,"\n\n------------------- %s ------------------- ",Header.c_str());
    
	fprintf(out,"\n Individuali");
    
	fprintf(out,"\t%22s","Hard Cluster Obtained");
	
	if (this->labelFeature!=-1){
		fprintf(out,"\t%22s","A priori partition");
	}//if

	for (c=0;c<this->nClusters;c++){
      fprintf(out,"\t%10s%02d","Mi",c+1);
	}//for c


	for (i=0;i<this->nIndividuals;i++){
		fprintf(out,"\n%12d",i+1);
		double maxM=MIN_DOUBLE;
		int clustHard=0;
		//Print partition obtained
		for (c=0;c<this->nClusters;c++){			
			if (maxM<membership[c][i]){
				maxM = membership[c][i];
				clustHard = c;
			}
		}//for c
		fprintf(out,"\t%22d",clustHard+1);
		//Print label
		if (this->labelFeature!=-1){
			fprintf(out,"\t%22.0f",this->labelFeatureData[i]);
		}//if
		//Print Memberships
        for (c=0;c<this->nClusters;c++){
			fprintf(out,"\t%12.8f",membership[c][i]);		
		}//for c
	}//for i
};

void FuzzyPartitioner::printHardPartition(vector<tVecDouble> membership,string Header, FILE *out){

	int c;

	fprintf(out,"\n\n%s",Header.c_str());
	this->createHardPartition(membership,this->hardPartition);
    for ( c=0;c<this->hardPartition.size();c++){
        fprintf(out,"\nCluster %d -> ",c+1);
		int countIndLine=0;
		for (int j=0;j<this->hardPartition[c].size();j++,countIndLine++){
			fprintf(out,"%d,",this->hardPartition[c][j]+1);
			if (countIndLine==50){
			  fprintf(out,"\n");
			  countIndLine = 0;
			}
		}
		fprintf(out,"\nCluster size-> %d",this->hardPartition[c].size());
	}//for c
};

void FuzzyPartitioner::makeOutPutReportFooter(){
//creates the output file

	
	this->makeOutPutCustomReport();

	this->printMembershipMatrix(this->bestInitialMembershipJ,"Initial Membership matrix for the minimal J obtained",this->out);

	this->printMembershipMatrix(this->bestMembershipJ,"Membership matrix for the minimal J obtained",this->out);

	this->printHardPartition(this->bestMembershipJ,"Hard partition obtained for the minimal J obtained",this->out);
	
    
	Partitioner::makeOutPutReportFooter();

};//void FuzzyPartitioner::makeOutPutReport(){

void FuzzyPartitioner::workOutBTR(double &B,double &T, double &R,double &nonFuz,double &dPrime){
  double retR=0,retB=0,retT=0,retNonFuz=0,retDPrime=0;
/*  for (int i=0;i<this->nClusters;i++){
	  double difCCG=0;
	  for (int f=0;f<this->nFeatures;f++){
		  difCCG += pow(this->centerCluster[i][f]-this->globalCenter[f],2.0);
	  }//for f
	  double u=0,difCGI=0;
	  for (int k=0;k<this->nIndividuals;k++){
		  double auxU = pow(this->membership[i][k],this->parameterM);
		  u += auxU;
		  for (int p=0;p<this->nFeatures;p++){
			  difCGI += pow((*this->tabData)[k][p]-this->globalCenter[p],2);
		  }//for p
		  retT += auxU*difCGI;
	  }//for j
	  retB += u*difCCG;
  }//for i*/
  int i,k,p;
  
  //It calculates NonFuz, T
  for (i=0;i<this->nClusters;i++){
	  for (k=0;k<this->nIndividuals;k++){
		  double dist=0;
		  for (p=0;p<this->nFeatures;p++){
			  dist += pow((*this->tabData)[k][p]-this->globalCenter[0][p],2.0);
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

#ifndef NO_INTERPRETATION_INDICE
  //It calculates B
  for (i=0;i<this->nClusters;i++){
	  double dist=0;
	  for  (p=0;p<this->nFeatures;p++){
		  dist += pow(this->centerCluster[i][p]-this->globalCenter[0][p],2.0);
	  }//for p
	  double u=0;
      for (k=0;k<this->nIndividuals;k++){
		  u += this->membershipRaisedToM[i][k];
	  }//for k
      retB += dist*u;
  }//for i
#endif

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



void FuzzyPartitioner::upDateCenterClusterForAdaptiveFuzzy(){
//it updates the center cluster and the global center
//it's exactally a copy of void FuzzyCMeans::upDateCenterClusters()
	int distZero=0;
        int c;
	for (c=0;c<this->nClusters;c++){
		for (int f=0;f<this->nFeatures;f++){
		  double deno = 0;
		  double num = 0;
		  int i;
		  for (i=0;i<this->nIndividuals;i++){
			double memb = this->membership[c][i];
			double membershipM = this->membershipRaisedToM[c][i];
			deno += membershipM;
			num += membershipM*(*this->tabData)[i][f];
			double inf=numeric_limits<double>::infinity();
			double NaN=numeric_limits<double>::quiet_NaN();
			double SNaN = numeric_limits<double>::signaling_NaN();
			if ( (deno==INDF) || (num==INDF))
             cerr<<"\nERROR1 - Indefinition -> upDateCenterClusters()";
		  }//for i
		  if (deno==0)
			  cerr<<"\nERROR2 - Division by Zero -> upDateCenterClusters()";
		  this->centerCluster[c][f] = num/deno;		  
          
		}//for f
		
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
              double tmp = this->membershipRaisedToM[i][k];
			  tmpGCPNum += tmp*(*this->tabData)[k][f];
			  tmpGCPDeno += tmp;
			}//for k
		}//for i
		this->globalCenter[1][f] = tmpGCPNum/tmpGCPDeno;
	}//for f
#endif
#ifdef DEBUG_CLUSTER_CENTER
	cout<<endl<<"Cluster Centers";
	for (c=0;c<this->nClusters;c++){
		cout<<endl<<"Centro "<<(c+1)<<" -> ";
		for (int f=0;f<this->nFeatures;f++){
			cout<<"C"<<(c+1)<<(f+1)<<"="<<this->centerCluster[c][f]<<",";
		}
	}
#endif

//Updates the difference from individuals to centercluster[c][f]
  //cout<<endl;
  for (c=0;c<this->nClusters;c++){
	for (int f=0;f<this->nFeatures;f++){
      for (int i=0;i<this->nIndividuals;i++){
        this->difIndFromCenter[i][c][f] = pow((*this->tabData)[i][f]-this->centerCluster[c][f],2.0);
#ifndef NO_DIST_IND_CENTER_ZERO_CHECK        
		if (this->difIndFromCenter[i][c][f]==0){
	      distZero++;
          double featI = (*this->tabData)[i][f];
	      double featC = this->centerCluster[c][f];
          cout<<"d("<<(i+1)<<","<<(c+1)<<","<<(f+1)<<")="<<this->difIndFromCenter[i][c][f]<<"  ";

		}//if 
#endif
	  }//for i
	}//for f
  }//for c
  //cout<<endl<<"Dist equals to zero="<<distZero<<endl;
};


void FuzzyPartitioner::cleanUp(){

	Partitioner::cleanUp();

	this->membership.clear();
	this->bestMembershipJ.clear();
	this->bestInitialMembershipJ.clear();
	this->tmpInitialMembershipJ.clear();

    
	
};