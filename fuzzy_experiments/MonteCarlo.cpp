// MonteCarlo.cpp: implementation of the MonteCarlo class.
//
//////////////////////////////////////////////////////////////////////


#include "MonteCarlo.h"
#define DEBUGMC 1
#define NOT_USE_CREATESYMBOLICS


//////////////////////////////////////////////////////////////////////
// MonteCarloUserInterface Class
//////////////////////////////////////////////////////////////////////
MonteCarloUserInterface::MonteCarloUserInterface(MonteCarlo *mc):mc(mc){
};

void MonteCarloUserInterface::readParameters(){	
	//this->mc->setClusterNumber();
	this->readOutPutFileName();
	this->readVariableNumber();
	this->readClassParameters();
	this->readRepetitions();
	#ifndef NOT_USE_CREATESYMBOLICS
		this->readIntervals();
	#endif
};


void MonteCarloUserInterface::readOutPutFileName(){
  string fileName;
  cout<<endl<<"Type output file (name.sds): ";
  cin>>fileName;
  this->mc->setSDSFileName(fileName);
};

void MonteCarloUserInterface::readRepetitions(){
	int nRepetitions;
	do{
	  std::cout<<endl<<"Number of repetitions? ";
	  std::cin>>nRepetitions;
	}while (nRepetitions<1);
	this->mc->setMCRepetitions(nRepetitions);
};

void MonteCarloUserInterface::readVariableNumber(){
	int nVariables;
	std::cout<<endl<<"Number of variables? ";
	std::cin>>nVariables;
	this->mc->setVariableNumber(nVariables);
};

void MonteCarloUserInterface::readIntervals(){
	int nInterval;
	cout<<"\n	Type the number of intervals: ";
	cin>>nInterval;
	this->mc->setIntervalNumber(nInterval);	
	double lower,upper;
	for(int i=0;i<nInterval;i++){
		do{
		  cout<<"\nInterval:"<<i+1;
		  cout<<"\n	Lower Bound:";
		  cin>>lower;
		  cout<<"\n	Upper Bound:";
		  cin>>upper;
		  if( lower > upper )
			cout<<"\nERRO, repeat"; 

		}while( lower > upper );
		this->mc->setInterval(i,lower,upper);
	}//for
};//void MonteCarloUserInterface::readIntervals(){

void MonteCarloUserInterface::readClassParameters(){
	for (int i=0;i<this->mc->getNumClasses();i++){//for 0
	    cout<<"\n			 Parameter of the Class "<<(i+1)<<"\n";
		cout<<"\n			--------------------------";		
		double correlation;
		do{		  
		  cout<<"\nType correlation for Class "<<(i+1)<<":";
		  cin>>correlation;
		  if(correlation>1.0 || correlation<-1.0)
			cout<<" Invalid, repeat\n";
		}while(correlation>1.0 || correlation<-1.0);
		this->mc->setCorrelation(i,correlation);
		int nObj;
		cout<<"\nNumber of objects:";
		cin>>nObj;		
        this->mc->setIndividuals(i,nObj);
		double average,standardDeviation;
		for(int k=0;k<this->mc->getNumVariables();k++){			
			cout<<"\nFor variable "<<(1+k)<<":";
			cout<<"\n                Average:";
			cin>>average;
			cout<<"\n                Varriance:";
			cin>>standardDeviation;
			this->mc->setAverage(i,k,average);
			this->mc->setStandardDeviation(i,k,pow(standardDeviation,0.5));
		}
		cout<<endl;
	}//for 0
	this->mc->createPartitionVector();
};//void MonteCarloUserInterface::readClassParameters(){

//////////////////////////////////////////////////////////////////////
// MonteCarlo Class
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MonteCarlo::MonteCarlo(int &nClasses){
  this->nClasses = nClasses;
  this->setClassesNumber(nClasses);
  this->usrInt = new MonteCarloUserInterface(this);  
}

MonteCarlo::~MonteCarlo(){
  delete this->usrInt;
  this->freePartitionVector();
}

/*tbl &MonteCarlo::getTabData(){
	return this->tabData;
};

void MonteCarlo::setTabData(tbl &tabData){
	this->tabData = tabData;
};*/

vector<double *> MonteCarlo::getPartitionVariable(){
  return this->partitionVariables;
};

vector<char *> MonteCarlo::getPartitionFeature(){
	return this->partitionFeature;
};


int MonteCarlo::numBase=0;

string MonteCarlo::getBaseName(){
	string ret;
	int i=100;
	//ret = ret+(i);
	ret.append("ArtificialBaseGeneratedByMonteCarloExperiment_"+(++MonteCarlo::numBase));
	//MonteCarlo::numBase++;
	//ret.append(MonteCarlo::numBase);
	return ret;
};


void MonteCarlo::freePartitionVector(){
/*	for (int i=0;i<this->partitionFeature.size();i++){
		free(this->partitionFeature[i]);
	}
	for (i=0;i<this->partitionVariables.size();i++){		
		free(this->partitionVariables[i]);
	}

	this->partitionFeature.clear();
	this->partitionVariables.clear();
	*/
};

void MonteCarlo::createPartitionVector(){  
  this->freePartitionVector();
  int nInd = 0, i;
  for ( i=0;i<this->individuals.size();i++)
      nInd += this->individuals[i];	
  this->nIndividuals = nInd;

  int classCounter=0;
  int classNumber = 0;
  this->individualLabels.resize(this->nIndividuals);
  for (i=0;i<this->individualLabels.size();i++){
      this->individualLabels[i] = classNumber+1;
	  if (++classCounter>=this->individuals[classNumber]){
		  classCounter =0;
		  classNumber++;
	  }//if
  }//for i

 /* this->tabData.ni = nInd;
  this->partitionFeature.resize(this->tabData.nclasses,NULL);
  this->partitionVariables.resize(this->tabData.ni,NULL);
  for (i=0;i<this->partitionFeature.size();i++){	
	char *buff = (char *)malloc(sizeof(char)*20);
	itoa(i+1,buff,10);
	this->partitionFeature[i] = strdup(buff);
	free(buff);
  }


  int classCounter=0;
  int classNumber = 0;
  for (i=0;i<this->partitionVariables.size();i++){
	  this->partitionVariables[i] = (double *)malloc(sizeof(double)*this->partitionFeature.size());
	  for (int j=0;j<this->partitionFeature.size();j++)
        this->partitionVariables[i][j] = 0.0;
	  //memset(this->partitionVariables[i],0,sizeof(double *)*this->partitionFeature.size());
      this->partitionVariables[i][classNumber] = 1.0;
	  if (++classCounter>=this->individuals[classNumber]){
		  classCounter =0;
		  classNumber++;
	  }
  }
*/

};

void MonteCarlo::memoryAllocation(){
/*	for(i = 0; i < tbl->nclasses; i++){//for 0
		tabela[i] = (double  **) malloc(tbl->nv * sizeof(double  *));		
		for(j = 0; j < tbl->nv; j++){//for 1
			if(((tbl->v[j].type_var == Qualitative) || (tbl->v[j].type_var == Ordinale)) || (tbl->v[j].type_var == Modal)){
				tabela[i][j] = (double  *)malloc(tbl->v[j].nm * sizeof(double ));
				for(k = 0; k < tbl->v[j].nm; k++){
					tabela[i][j][k] = 0.0;
				}//for
			}//if
			else if(tbl->v[j].type_var == Quantitative){			
				tabela[i][j] = (double  *)malloc(1 * sizeof(double ));
				for(k = 0; k < 2; k++){
					tabela[i][j][k] = 0.0;
				}//for				
			}//else if
		}//for 1
	}//for 0*/
};

void MonteCarlo::initialize(){
	

	/*this->tabData.mat_poids = NULL;	
    this->tabData.titre = strdup("Artificial Base");
	this->tabData.ti = (char  **)malloc(this->tabData.ni * sizeof(char  *));	
	for (int n=0;n<this->tabData.ni;n++)
		this->tabData.ti[n] = strdup("EXE");
    this->tabData.v = (struct var  *)malloc(this->tabData.nv * sizeof(struct var ));
	for(int i=0; i < this->tabData.nv; i++){
	  this->tabData.v[i].type_var=3;	    
      this->tabData.v[i].macroClassification = QUANTITATIVE_FEATURE;
      this->tabData.v[i].interv = (struct interv *)malloc(sizeof(struct interv));
	  this->tabData.v[i].interv->inf = 0;
	  this->tabData.v[i].interv->sup = 0;
	  this->tabData.v[i].nm = 0;
	  this->tabData.v[i].tv = NULL;
	  this->tabData.v[i].mod = NULL;	  
	}
    
	this->tabData.mat = (double  ***)malloc(this->tabData.ni * sizeof(double  **));
	for(int m=0; m < this->tabData.ni; m++) {
		this->tabData.mat[m] = (double  **)malloc(this->tabData.nv * sizeof(double  *));
		for(int l = 0; l < this->tabData.nv;l++) {
			this->tabData.mat[m][l] = (double  *)malloc(2 * sizeof(double ));
			this->tabData.mat[m][l][0] = 0.0;
			this->tabData.mat[m][l][1] = 0.0;
		}//for		
	}//for */
int i;
	#ifndef NOT_USE_CREATESYMBOLICS
		this->tempValue.resize(this->nIndividuals);
                 
		for (i=0;i<this->nIndividuals;i++){//for 
			this->tempValue[i].resize(this->nVariables);
			for (int j=0;j<this->nVariables;j++){//for
			  this->tempValue[i][j].resize(2);
			  this->tempValue[i][j][0] = 0.0;
			  this->tempValue[i][j][1] = 0.0;
			}//for
		}//for		
	#else
		this->tabData.resize(this->nIndividuals);
		for (i=0;i<this->tabData.size();i++)
			this->tabData[i].resize(this->nVariables,0.0);
	#endif
    int pos;
    
	string fileOutName = this->fileOutPreffixName.substr(0, (pos = this->fileOutPreffixName.find_last_of(".",this->fileOutPreffixName.length()))==string::npos?this->fileOutPreffixName.length():pos);  
	this->mc = fopen((fileOutName+"_montecarlo_L2.wri").c_str(),"a");
	this->mcAUX = fopen(string(fileOutName+"_AUX_DATA.wri").c_str(),"w");
	this->table =  fopen((fileOutName+"_usualdata.wri").c_str(),"w");
	this->symbolicTable = fopen((fileOutName+"_simbolico.txt").c_str(),"a");
    this->sdsFile = fopen(this->fileSDSName.c_str(),"w");
	//the number 9 below refers to the others indexs observed according to the J criteria
	this->indexs = (StatisticalSummary **)malloc(sizeof(StatisticalSummary*)*(2+9));
	this->indexs[0] = new StatisticalSummary("RC",true);
	this->indexs[1] = new StatisticalSummary("Objective Function",false);
	for (i=2;i<11;i++){
		this->indexs[i] = NULL;
	}//for i
};//void MonteCarlo::initialize(){

void MonteCarlo::finalize(){
	for (int i=0;i<11;i++){
		if (this->indexs[i]!=NULL)
			delete this->indexs[i];
	}//for i 
	free(this->indexs);
	fclose(this->mc);
	fclose(this->mcAUX);
//	fclose(this->table);
	fclose(this->symbolicTable);
	fclose(this->sdsFile);
	//Application is responsible for liberating memory alocated on tbl structure
};

//int sodasFiles=0;
void MonteCarlo::createSODASFile(){
	int i, j,l;
	double min,max;		
	vector<MonteCarlo::vecDoubleType> limit;

	limit.resize(this->nVariables);	
	for(i=0;i<limit.size();i++){
		limit[i].resize(2);			
	}//for
	
	//look for the minimun and maximun values of each variable 
	for(j = 0; j < this->nVariables; j++){//for 0	
        min = HUGE_VAL;
		max = -HUGE_VAL;		
		//menor do limite inferior
		for(l=0; l < this->nIndividuals; l++){
			if (min>this->tabData[l][j])
				min = this->tabData[l][j];
			if (max<this->tabData[l][j])
				max = this->tabData[l][j];
		}//for			
		limit[j][0]=min;		
		limit[j][1]=max;	
	}//for 0

	//It Creates SODA File	
    int nind, nvar;
	
	nind = this->nIndividuals;
	nvar = this->nVariables;    
	//fprintf(this->sdsFile,"\n[%d]\n",sodasFiles+1);
	fprintf(this->sdsFile,"SODAS = (\n");
    fprintf(this->sdsFile,"CONTAINS = (\n");
    fprintf(this->sdsFile,"   FILES, HEADER, INDIVIDUALS, VARIABLES, RECTANGLE_MATRIX\n");
    fprintf(this->sdsFile,"),\n"); 
    fprintf(this->sdsFile,"FILE =  (\n"); 
    fprintf(this->sdsFile,"   procedure_name = \"DB2SO\" ,\n");
    fprintf(this->sdsFile,"   version = \"sans\" ,\n");
    fprintf(this->sdsFile,"   create_date = \"\"\n");
    fprintf(this->sdsFile,"),\n");
    fprintf(this->sdsFile,"HEADER =  (\n"); 
    fprintf(this->sdsFile,"   title = \"Titre\" ,\n");
    fprintf(this->sdsFile,"   sub_title = \"Sous-titre\" ,\n");
    fprintf(this->sdsFile,"   indiv_nb = %d,\n",nind);
	//It adds the label variable into the SDS
    fprintf(this->sdsFile,"   var_nb = %d ,\n",nvar+1);
    fprintf(this->sdsFile,"   rules_nb = 0 ,\n");
    fprintf(this->sdsFile,"   nb_var_set = 0 ,\n");
    fprintf(this->sdsFile,"   nb_indiv_set = 0 ,\n");
    fprintf(this->sdsFile,"   nb_var_nom = 0 ,\n");
    fprintf(this->sdsFile,"   nb_var_cont = 0 ,\n");
    fprintf(this->sdsFile,"   nb_var_text = 0 ,\n");
    fprintf(this->sdsFile,"   nb_var_cont_symb = %d ,\n",nvar);
    //It´s here due to put label variable
    fprintf(this->sdsFile,"   nb_var_nom_symb = 1 ,\n");
    fprintf(this->sdsFile,"   nb_var_nom_mod = 0 ,\n");
    fprintf(this->sdsFile,"   nb_na = 0 ,\n");
    fprintf(this->sdsFile,"   nb_null = 0 ,\n");
    fprintf(this->sdsFile,"   nb_nu = 0 ,\n");
    fprintf(this->sdsFile,"   nb_hierarchies = 0\n");
    fprintf(this->sdsFile,"),\n\n");
    fprintf(this->sdsFile,"INDIVIDUALS = (\n");
  
    //It prints the individuals 
	int cla=0,indAnt=0;
    for(i=0; i < (nind); i++) {
	  if (((i+1)-indAnt)>this->individuals[cla]){
		  indAnt += this->individuals[cla];
		  cla++;
	  }
	  fprintf(this->sdsFile,"  (%d,\"EXE\", \"%02d\" )%c\n",i,cla+1,(i==nind-1?'\0':','));
    }//for	
	fprintf(this->sdsFile,"),\n");

    //It prints the variables
    fprintf(this->sdsFile,"VARIABLES =  (\n");
    for(i=0; i < nvar; i++) {
     fprintf(this->sdsFile," (%d ,inter_cont ,\"\" ,\"IIII\" ,\"position  %d\" , 0, 0, %f, %f),\n",(i+1),(i+1),limit[i][0],limit[i][1]); 
    }//for
    fprintf(this->sdsFile," (%d ,nominal ,\"\" ,\"IIII\" ,\"Class\" , 0, 0, %d, (",(i+1),this->nClasses); 
	for (i=0;i<this->nClasses;i++){
	  fprintf(this->sdsFile,"\n(%02d,\"AA%2d\",\"CLA%02d\",0)%c",(i+1),(i+1),(i+1),(i==this->nClasses-1?')':','));
	}//for
    fprintf(this->sdsFile,"\n)\n),\n");

	//It prints the rectangle matrix
    fprintf(this->sdsFile,"RECTANGLE_MATRIX = (\n");
  	cla=0,indAnt=0;
    for(i=0; i < (nind); i++) {//for 0         
	  fprintf(this->sdsFile,"(");
	  for(j=0; j< (nvar); j++){		      
	    fprintf(this->sdsFile," ( %f : %f ),",this->tabData[i][j], this->tabData[i][j]); 
	  }//for          	  
	  if (((i+1)-indAnt)>this->individuals[cla]){
		  indAnt += this->individuals[cla];
		  cla++;
	  }
      fprintf(this->sdsFile,"%d)",cla+1);
      fprintf(this->sdsFile,"%c\n",(i==nind-1?')':','));  		
    }//for 0
 
	fprintf(this->sdsFile,")\n");
	fprintf(this->sdsFile,"  END\n\n");    
	
    fflush(this->sdsFile);
	
};//void MonteCarlo::createSODASFile(){

void MonteCarlo::createNormals(bool bMakeReport){
   #define REPORT 1	
   #define PI 3.1415926535
   int  k, valor, i, limite;
   double r, expo, ro_qua, media_cond, desvio_cond;
   double temp_a;
   double u1, u2, x1, x2;
   double z1, z2;   
   vector<int> limi_in;     


   /*if (DEBUGMC){
	   for(k=0; k<this->tabData.nclasses; k++){
		   cout<<endl<<"Class "<<k+1<<" Correlation="<<this->correlation[k];
		   for (int v=0;v<2;v++){
			 cout<<endl<<"Variable "<<v+1<<" : "<<endl;
		     cout<<"  -> Average="<<this->average[k][v]<<"StdDeviation="<<this->standardDeviation[k][v];
		   }//for v
	   }//for k
   }*/

   limi_in.resize(this->nClasses);
   for(i = 0; i < this->nClasses; i++)
		limi_in[i]=0;   
   
   for(i = 1; i < this->nClasses; i++)
	 limi_in[i]=limi_in[i-1] + this->individuals[i-1];  

   for(k=0; k<this->nClasses; k++){//for 0
		i=limi_in[k];
		limite= limi_in[k] + this->individuals[k];
		if(bMakeReport)	{
			//REPORT?fprintf(this->table,"\n Class: %d\n",k+1):0;			
		}
		while(i<limite)	{//while 0
			// gera uniforme
			valor = rand()%100;
			u1=(double)valor/100;
			//gera exponencial
			do{//do
                 valor = rand()%100;
			     u2=(double)valor/100;
            }while(u2 == 0.0);//do-while
			expo = -log(u2);
			r=expo*2.0;
			expo=sqrt(r);	   
			//REPORT z1 e z2
			//<REVER>
			// é para por Pi?
			temp_a=6.28*u1;	//falta pi
			z1=expo*cos(temp_a);
			z2=expo*sin(temp_a);			
			//<REVER>
			// é para para usar uma média e desvio fixos?
    		x2=this->average[k][1]+(this->standardDeviation[k][1]*z2);            
			media_cond=this->average[k][0]+(((this->correlation[k]*this->standardDeviation[k][0])/this->standardDeviation[k][1])*(x2-this->average[k][1]));            				
			ro_qua=pow(this->correlation[k],2);			
			ro_qua=1.0-ro_qua;
			desvio_cond=this->standardDeviation[k][0]*z1*pow(ro_qua,0.5);			
			x1=media_cond+desvio_cond;
			if(bMakeReport)	{
				REPORT?fprintf(this->table,"%d\t%f\t",k+1, x1):0;
			 	REPORT?fprintf(this->table,"%f\n", x2):0;
			}	
            #ifdef NOT_USE_CREATESYMBOLICS
               this->tabData[i][0] = x1;
			   this->tabData[i][1] = x2;
            #else
			  this->tempValue[i][0][0]=x1;
			  this->tempValue[i][0][1]=x1;
			  this->tempValue[i][1][0]=x2;
			  this->tempValue[i][1][1]=x2;
            #endif
			i++;	    
		}//while 0		
   }//for 0
   if(bMakeReport)	{
	   fflush(this->table);
   }
};//void MonteCarlo::createNormals(bool bMakeReport){

void MonteCarlo::createSymbolics(double lower, double upper,bool bMakeReport){
   #define REPORT_SYMBOLIC 1
   int j, k, valor,diferenca;
   float range;   

   diferenca=upper - lower;
   diferenca=diferenca*10 + 1;
	if(bMakeReport){
		REPORT_SYMBOLIC?fprintf(this->symbolicTable,"Number of individuals = %d\n",this->nIndividuals):0;
		REPORT_SYMBOLIC?fprintf(this->symbolicTable,"Number of variables = %d\n",this->nVariables):0;
	}
    //srand( (unsigned)time( NULL ) );
	for(j=0; j<this->nIndividuals; j++){//for 0
		for(k=0; k<this->nVariables; k++){
			valor = rand()%diferenca;
			range=(float)valor/10.0 + lower;
			range=range/2.0;
			if (random(2)==0)
			  this->tabData[j][k]=this->tempValue[j][k][0]-range;
			else
	          this->tabData[j][k]=this->tempValue[j][k][1]+range;
			if(bMakeReport){
				REPORT_SYMBOLIC?fprintf(this->symbolicTable,"%f %f ", this->tabData[j][k], this->tabData[j][k]):0;
			}
       	}//for
		if(bMakeReport){
			REPORT_SYMBOLIC?fprintf(this->symbolicTable,"\n"):0;	
		}
   }//for 0
   if(bMakeReport){
  	  fflush(this->symbolicTable);	
   }
};//void MonteCarlo::createSymbolics(double lower, double upper,bool bMakeReport){

void MonteCarlo::printSummaryTitle(){
  fprintf(this->mc,"\n");
int i;
  for (i=0;i<100;i++){
	fprintf(this->mc,"*");			
  }
  fflush(this->mc);
  
  
  fprintf(this->mcAUX,"\n");
  for ( i=0;i<100;i++){
	fprintf(this->mcAUX,"*");			
  }
  			
  fflush(this->mcAUX);
  //<FAZER>
  //colocar aqui informacoes q sumarizem o status do FSCM
};

void MonteCarlo::printOutInformation(string str){
	fprintf(this->mc,"\n%s\n",str.c_str());
	fprintf(this->mc,"\nNumber of replications %d",this->mcRepetitions);
	fflush(this->mc);
};


void MonteCarlo::printOutStatistics(double lower, double upper, double bestRC,double worstRC, double ofForBestRC,double averageRC, double stdDeviationRC, double bestOF, double bestOFRC, double averageOF, double stdDeviationOF,double &NumIte,vector<double>& bestJAssociatedValues,vector<string>& bestJAssociatedValueDescriptions){
  fprintf(this->mc,"\n");
  #ifndef    NOT_USE_CREATESYMBOLICS
    fprintf(this->mc," Range [%f,%f]\n",lower,upper);
  #endif
  fprintf(this->mc,"\nAvg number of iterations = %6.3f",NumIte);
  fprintf(this->mc,"\n\nInformation about Adjusted Rand Index (RC)\n");
  fprintf(this->mc,"  Best |    Worst    |       OF       | Average | Standard Deviation\n");
  fprintf(this->mc,"%7.5f|%13.7f|%16.7f|%10.7f|%19.7f\n", bestRC,worstRC,ofForBestRC,averageRC, stdDeviationRC);
  fprintf(this->mc,"Information about Objective Function values (OF)\n");
  fprintf(this->mc,"      Best      |      Worst     |    Average     | Standard Deviation\n");
  fprintf(this->mc,"%16.7f|%16.7f|%16.7f|%19.7f\n", bestOF,this->indexs[1]->getWorstValue(),averageOF, stdDeviationOF);
  fprintf(this->mc,"\n");
  fprintf(this->mc,"\n-------- Indexs values to the minimal J obtained over all replications-------- ");
  
  int i;
  
  for (i=0;i<bestJAssociatedValueDescriptions.size();i++){
	  fprintf(this->mc,"\n\t->%30s = %-16.8f",this->indexs[i+2]->getValueDescription().c_str(),this->indexs[1]->getBestValueAssociationsValue(i));      	  
  }
  fprintf(this->mc,"\n-------- Indexs related to the minimal J obtained in each replication-------- ");
  
  
  for (i=0;i<bestJAssociatedValueDescriptions.size();i++){
	  fprintf(this->mc,"\n\t->%-30s",this->indexs[i+2]->getValueDescription().c_str());
      fprintf(this->mc,"\n\t------------------------------------------------------------------------------------------------------------");
	  fprintf(this->mc,"\n\t\t    Avg           |  Std Deviation  |      Best         |     Worst");
	  fprintf(this->mc,"\n\t\t%16.8f   %16.8f   %16.8f   %16.8f",this->indexs[i+2]->getMean(),this->indexs[i+2]->getStdDeviation(),this->indexs[i+2]->getBestValue(),this->indexs[i+2]->getWorstValue());
	  
  }
  fflush(this->mc);
  fprintf(this->mcAUX,"\n%.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f", ofForBestRC,bestRC,worstRC,averageRC, stdDeviationRC,bestOF,this->indexs[1]->getWorstValue(),averageOF, stdDeviationOF);  
     // T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
    //indexs in best J associated values
   //   1,2,0 ,3,   4          ,  5         , 6 ,  7, 8
  for (i=0;i<bestJAssociatedValueDescriptions.size();i++){
    fprintf(this->mcAUX," %.6f %.6f %.6f %.6f %.6f",this->indexs[1]->getBestValueAssociationsValue(i), this->indexs[i+2]->getMean(),this->indexs[i+2]->getStdDeviation(),this->indexs[i+2]->getBestValue(),this->indexs[i+2]->getWorstValue());  
  }//for i
  for (i=0;i<bestJAssociatedValueDescriptions.size();i++){
	  fprintf(this->mcAUX," %.6f",bestJAssociatedValues[i]);
  }
  
  fprintf(this->mcAUX," %.6f",NumIte);
  #ifdef DEBUGTIME
    fprintf(this->mcAUX," %.6fs%dmil",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);
  #endif
  fflush(this->mcAUX);

};

void MonteCarlo::printOutPartitionerAlgorithmCfg(Partitioner &algorithm){
  vector<string> strTmp = algorithm.getConfiguration();
  for (vector<string>::iterator itr=strTmp.begin();itr!=strTmp.end();itr++){
	  fprintf(this->mc,itr->c_str());
  }

  fflush(this->mc);
};


int MonteCarlo::mcRepetitions = 100;

void MonteCarlo::execute(Partitioner &algorithm,string fileOutName,int runs,int iterations, double parameterM){
    this->fileOutPreffixName = fileOutName;
 	this->initialize();
	bool bStart = true;
	double NumIte=0;


#ifndef NOT_USE_CREATESYMBOLICS
    for (int inter=0;inter<this->intervals.size();inter++){//for 0
		double lower = this->intervals[inter].lowerBound,upper=this->intervals[inter].upperBound;
#endif
    
        #ifdef DEBUGTIME
		#ifndef _WIN32 
          ftime( &timebufferBeforeRun );           
		#else
		  _ftime( &timebufferBeforeRun );           
		#endif
        #endif

		for (int mcRep=0;mcRep<this->mcRepetitions;mcRep++){//for 1
            
			createNormals(bStart);
			#ifndef NOT_USE_CREATESYMBOLICS
			  createSymbolics(lower,upper,bStart);
			  DEBUGMC?cout<<endl<<"Interval ["<<lower<<","<<upper<<"] replica "<<mcRep<<endl,1:0;
			#endif

			cout<<endl<<"Replication "<<(mcRep+1)<<endl;
			 //It is necessary to create a different file about run information
			 //int pos;
	         //string fileOutName = this->fileOutPreffixName.substr(0, (pos = this->fileOutPreffixName.find_last_of(".",this->fileOutPreffixName.length()))==string::npos?this->fileOutPreffixName.length():pos);  
			 //fileOutName.append("RunsForReplication");
			 //char buf[10];
			 //itoa((mcRep+1),buf,10);
			 //fileOutName.append(buf);
			 //fileOutName.append(".wri");
			algorithm.execute(this->tabData,this->individualLabels,this->nClasses,-1,"","",runs,iterations,this->nClasses,parameterM,true,true);
			NumIte += algorithm.getAvgIte();
			if (bStart){
			  this->createSODASFile();
			  this->printSummaryTitle();
			  this->printOutPartitionerAlgorithmCfg(algorithm);
			}
			this->indexs[0]->addValue(algorithm.getBestRC());
		    this->indexs[0]->setBestValueAssociationsValue(0,algorithm.getJForBestRC());
			
			this->indexs[1]->addValue(algorithm.getBestJ());

			int k;
			

			if (bStart){
				
				vector<string> vecValDesc = algorithm.getBestJAssociatedValuesDescription();            
				for (k=0;k<vecValDesc.size();k++){
				  this->indexs[1]->setBestValueAssociationDescription(k,vecValDesc[k]);
				}//for k
				StatisticalSummary **staSum = algorithm.getStatisticalIndexs();
				int idxs[] = {3,1,2,4,5,6,7,8,9};
                for (k=2;k<11;k++){
					this->indexs[k] = new StatisticalSummary(staSum[idxs[k-2]]->getValueDescription(),staSum[idxs[k-2]]->getBestCriteria());
				}//for k

			}

			vector<double> vecVal = algorithm.getBestJAssociatedValues();            
			for (k=0;k<vecVal.size();k++){
			  this->indexs[1]->setBestValueAssociationsValue(k,vecVal[k]);
			  
			  this->indexs[k+2]->addValue(vecVal[k]);
			}//for k
			  fprintf(this->mc,"\n\nReplication: %d",mcRep+1);   
              fprintf(this->mc,"\n      J             |      T           |      R           |      B           |      RC          |   Non Fuzziness    |       Dk(U)'       |      Tp          |       Bp         |     Rp           |     Avg Iterations");
              fprintf(this->mc,"\n-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------");
	          
                 //0,1,2,3, 4, 5            , 6          , 7 ,  8, 9																
	 			// J,T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
			    // T,R,RC,B, Non Fuzziness, Dk(U) Prime, Tp, Bp, Rp
                //indexs in best J associated values
		        // 1,2,0 ,3,   4          ,  5         , 6 ,  7, 8
			  
              fprintf(this->mc,"\n   %16.8f   %16.8f   %16.8f   %16.8f   %16.8f    %16.8f      %16.8f   %16.8f   %16.8f   %16.8f  %16.8f",this->indexs[1]->getLastValueAdded(),vecVal[1],vecVal[2],vecVal[3],
	            vecVal[0],vecVal[4],vecVal[5],vecVal[6],vecVal[7],vecVal[8],algorithm.getAvgIte());
			  
              #ifdef DEBUGTIME
		        #ifndef _WIN32 
                  ftime( &timebufferAfterRun );           
		        #else
		          _ftime( &timebufferAfterRun );           
		        #endif
			    timebuffer.millitm = timebufferAfterRun.millitm - timebufferBeforeRun.millitm;
			    (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?timebuffer.millitm += 1000 : 0; 
			    fprintf(this->mc,"\n\n-Time used to execute the replication: %.3f s %d mil\n",difftime(timebufferAfterRun.time,timebufferBeforeRun.time)+( (timebufferAfterRun.millitm<timebufferBeforeRun.millitm)?-1:0 ),timebuffer.millitm);
              #endif
              fflush(this->mc);			
            algorithm.doCleanUp();       			
			bStart = false;
		}//for 1		
		vector<string> vecDescriptions=this->indexs[1]->getBestValueAssociationsDescription();
	    vector<double> vecValues=this->indexs[1]->getBestValueAssociationsValues();
		NumIte /= this->mcRepetitions;		
#ifndef NOT_USE_CREATESYMBOLICS
		this->printOutStatistics(lower,upper,this->indexs[0]->getBestValue(),this->indexs[0]->getWorstValue(),this->indexs[0]->getBestValueAssociationsValue(0),this->indexs[0]->getMean(),this->indexs[0]->getStdDeviation(),this->indexs[1]->getBestValue(),this->indexs[1]->getBestValueAssociationsValue(0),this->indexs[1]->getMean(),this->indexs[1]->getStdDeviation(),NumIte,vecValues,vecDescriptions);
	}//for 0
#else
 	    this->printOutStatistics(0,0,this->indexs[0]->getBestValue(),this->indexs[0]->getWorstValue(),this->indexs[0]->getBestValueAssociationsValue(0),this->indexs[0]->getMean(),this->indexs[0]->getStdDeviation(),this->indexs[1]->getBestValue(),this->indexs[1]->getBestValueAssociationsValue(0),this->indexs[1]->getMean(),this->indexs[1]->getStdDeviation(),NumIte,vecValues,vecDescriptions);
#endif
	
	this->finalize();
};

MonteCarloUserInterface *MonteCarlo::getUserInterface(){
	return this->usrInt;
};

void MonteCarlo::setSDSFileName(string sdsFileName){
	this->fileSDSName = sdsFileName;
};

void MonteCarlo::setClassesNumber(int n){
	this->correlation.resize(n);
	this->average.resize(n);
	this->standardDeviation.resize(n);
	this->individuals.resize(n);
	this->individualLabels.resize(n);
};

void MonteCarlo::setVariableNumber(int n){
	this->nVariables = n;
	for (int i=0;i<this->correlation.size();i++){
	  this->average[i].resize(n);
	  this->standardDeviation[i].resize(n);	
	}
};

void MonteCarlo::setCorrelation(int cluster,double correlation){
	this->correlation[cluster] = correlation;
};

void MonteCarlo::setAverage(int cluster,int variable,double average){
	this->average[cluster][variable] = average;
};

void MonteCarlo::setStandardDeviation(int cluster,int variable, double standardDeviation){
	this->standardDeviation[cluster][variable] = standardDeviation;
};


void MonteCarlo::setIndividuals(int cluster,int number){
	this->individuals[cluster] = number;
};

void MonteCarlo::setIntervalNumber(int n){
  this->intervals.resize(n);
};

void MonteCarlo::setInterval(int interval,double lower,double upper){
	this->intervals[interval].lowerBound = lower;
    this->intervals[interval].upperBound = upper;
};

double MonteCarlo::getNumClasses(){
	return this->nClasses;
};

double MonteCarlo::getNumVariables(){
	return this->nVariables;
};

void MonteCarlo::setMCRepetitions(int repetitions){
	MonteCarlo::mcRepetitions = repetitions;
};