#include "AdjustedRandIndex.h"

using namespace std;

#include <iostream>

#define DEBUGRC 0

AdjustedRandIndex::AdjustedRandIndex(vector<double> &partitionVariables,int nClasses):
  partitionVariables(partitionVariables),nClasses(nClasses){
  /*cout<<endl<<&this->partitionVariables<<endl;
  for (int i=0;i<this->partitionVariables.size();i++){
	  for (int j=0;j<partitionFeature.size();j++)
		  cout<<this->partitionVariables[i][j]<<",";
	  cout<<endl;
  }*/
  

}

long AdjustedRandIndex::factorial(int m){
	long ret = 1;
	for (int i=1;i<=m;i++)
		ret *=i;
	return ret;
};

double AdjustedRandIndex::mChosenk(int m,int k){
  double ret=1;
  for (int i=m;i>(m-k);i--)
	  ret *= i;
  return ((double)ret/(double)factorial(k));
};
/*  This method does the assumption of SODAS´s structure retrievals a field which corresponds to 
   the element class. This field is considered as being a vector of n elements where just the element which
   correspondes to the element class has a value diferent from zero
*/
void AdjustedRandIndex::makeContingencyTable(int clusters,vector<int> **&algorithmPartition){

	int i,j;

	this->contingencyTable.clear();
	this->contingencyTable.resize(clusters);
	
	/*cout<<"\nInitial Partition";
	for (int j=0;j<this->partitionFeature.size();j++){
	  cout<<"\nClass "<<(j+1)<<" | ";
	  for (int k=0;k<this->partitionVariables.size();k++){
        if (this->partitionVariables[k][j]>0){
          cout<<k<<",";
		}
	  }
	}

	cout<<"\nAlgorithm Partition";
	for (j=0;j<clusters;j++){
	  cout<<"\nClass "<<(j+1)<<" | ";
	  for (int k=0;k<algorithmPartition[j]->size();k++){
          cout<<(*algorithmPartition[j])[k]<<",";		
	  }
	}*/


/*cout<<endl<<&this->partitionVariables<<endl;
  for (i=0;i<this->partitionVariables.size();i++){
	  for (j=0;j<partitionFeature.size();j++)
		  cout<<this->partitionVariables[i][j]<<",";
	  cout<<endl;
  }*/
    DEBUGRC?(cout<<"\nContingencyTable",1):0;
	  int totalEle=0;
	for (i=0;i<clusters;i++){//for1
		DEBUGRC?cout<<endl,1:0;
	  int intersectionsSum = 0;
      this->contingencyTable[i].clear();
      this->contingencyTable[i].resize(this->nClasses+1);
	  for (int tmp=0;tmp<this->contingencyTable[i].size();tmp++)
		  this->contingencyTable[i][tmp]=0;
	  for (j=0;j<this->nClasses;j++){//for 2
		  int intersection = 0;
		  //for (int k=0;k<this->partitionVariables.size();k++){//for 3
		  for (int l=0;l<algorithmPartition[i]->size();l++){
			     //the classes are number starting from 1
				 if ( this->partitionVariables[(*algorithmPartition[i])[l]]==j+1){
				   intersection++;				  
				 }//if
		  }//for l
		  //}//for k
          intersectionsSum += intersection;		 
		  DEBUGRC?cout<<"  ("<<i<<","<<j<<") = "<<intersection,1:0;
          this->contingencyTable[i][j] = intersection;
		  if (intersectionsSum==algorithmPartition[i]->size())
			  break;
	  }//for j
      DEBUGRC?cout<<"   "<<intersectionsSum,1:0;
	  totalEle += intersectionsSum;
	  this->contingencyTable[i][this->nClasses] = intersectionsSum;
	}//for i
	DEBUGRC?cout<<endl<<"Total of elements: "<<totalEle<<endl,1:0;

};

/*  
     i indexes rows; j indexes columns
	 This method consider that if i equals -1  then the result should be the sum of all elements
   in the column j. If j equals -1 then the result should be the sum of all elements in the row
   i.
*/
long AdjustedRandIndex::getContingencyTableElement(int i,int j){
   if (j==-1)
     return this->contingencyTable[i][this->contingencyTable[i].size()-1];
   if (i==-1){
	 long ret =0;
	 for (int k=0;k<this->contingencyTable.size();k++)
       ret += this->contingencyTable[k][j];
     return ret;
   }
   return this->contingencyTable[i][j];

};


double AdjustedRandIndex::getA(){
	double a=0.0;
	for (int r=0;r<this->contingencyTable.size();r++)
		for (int c=0;c<this->contingencyTable[r].size()-1;c++)
			a += this->mChosenk(this->getContingencyTableElement(r,c));
	return a;

};

double AdjustedRandIndex::getAdjustedRandIndex(int clusters,vector<int> **&algorithmPartition){
 this->makeContingencyTable(clusters,algorithmPartition);
 double a = this->getA();
 double fator1=0.0;
 for (int r=0;r<this->contingencyTable.size();r++)
   fator1 += this->mChosenk(this->getContingencyTableElement(r,-1));

 double fator2=0.0;
 for (int c=0;c<this->nClasses;c++)
   fator2 += this->mChosenk(this->getContingencyTableElement(-1,c)); 
 double fator3 = this->mChosenk(this->partitionVariables.size());
 double fator4 = (fator1*fator2/fator3);
 //cout<<endl<<"a="<<a<<",fator1="<<fator1<<",fator2="<<fator2<<",fator3="<<fator3<<",fator4="<<fator4;
 return (a - fator4 ) / (0.5*(fator1+fator2)-fator4);

};

AdjustedRandIndex::~AdjustedRandIndex(void){
}
