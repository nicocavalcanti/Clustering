// Application.cpp: implementation of the Application class.
//
//////////////////////////////////////////////////////////////////////

#include "Application.h"
#include "ClusteringAlgorithm.h"
//Debug macros
//#define DEBUG_LOAD_DATA

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application::Application():bNeedParamM(true){
  this->tabData.resize(0);
  this->individuals.resize(0);
  this->features.resize(0);
}

Application::~Application(){

}

//////////////////////////////////////////////////////////////////////
//it should observe the dependecies among the other readX methods for example if a method X depends on B and C as pre-condition so call B and C before you can call X
//////////////////////////////////////////////////////////////////////
void Application::readParameters(){
	this->readApplicationKind();
	
	 switch (this->appKind){
	  case 1:
		  //Monte Carlo experiment
		  this->readFileOutName();
	      this->readAlgorithm();
		  //k-means
          if (this->algorithm==4){
			  this->bNeedParamM = false;
		  }
		  if (this->bNeedParamM){
	        this->readMParameter();
		  }	      
	      this->readRuns();
	      this->readIterations();
		  this->labelFeature = 1;
		  this->readNumClasses();
          this->mc = new MonteCarlo(this->nClasses);
		  this->usrInterface = this->mc->getUserInterface();
		  this->usrInterface->readParameters();
		  break;
	  default:
		  //Normal execution
	      this->readIndividuals();
	      this->readFeatures();
	      this->readLabelFeature();  
		  this->readNumClasses();
	      this->readClusterNumber();
	      this->readDataFileName();
	      this->readFileOutName();
	      this->readAlgorithm();
		  //k-means
          if (this->algorithm==4){
			  this->bNeedParamM = false;
		  }
		  if (this->bNeedParamM){
	        this->readMParameter();
		  }
	      this->readRuns();
	      this->readIterations();	  
		 break;
  };	
	
	
	
};

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void Application::execute(){
	this->readParameters();
	Partitioner * part;
	switch (this->appKind){
	  case 1:
        part = ClusteringAlgorithm::getAlgorithm(this->algorithm);
		this->mc->execute(*part,this->fileOutName,this->runs,this->iterations,this->parameterM);
		delete part;
        delete this->mc;  
		break;
	  default:
	    part = ClusteringAlgorithm::getAlgorithm(this->algorithm);
	    part->execute(this->tabData,this->individualLabels,this->nClasses,this->labelFeature,this->fileInName,this->fileOutName,this->runs,this->iterations,this->nClusters,this->parameterM,false);
	    delete part;
		break;
	};
};


void Application::readApplicationKind(){
do{
		std::cout<<endl<<"Would you like to execute Monte Carlo experiment? (0=No 1=Yes)"<<endl;
        std::cin>>this->appKind;
		//cout<<endl<<"Valor lido="<<this->applicationKind;
	}while (this->appKind<0 || this->appKind>1);
cin.ignore(numeric_limits<int>::max(),'\n');
};
//////////////////////////////////////////////////////////////////////
// Pre-condition: 
//   (1) the individual to be used have to be known
//   (2) the features to be read have to be known
//   (3) the feature label must be known               
// Coments:
//   It reads the file input name and load the data
//////////////////////////////////////////////////////////////////////
void Application::readDataFileName(){
  cout<<endl<<endl<<"Please, inform the name of a file contaning correct data?"<<endl;
  cin>>this->fileInName;
  this->loadData();
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1) The file must contain data in the form:
//        1-each line must correspond to a individual data
//        2-each line must have all the features possible with some value
//        3-the features in a line are distinguished by a single space
//        4-the features values are just numeric
//   (n)
// Post-conditions:
//   (1)
//   (n)
// Coments:
//    It tries to forecast the number of individuals the file will have by resizing the vector to the maximum element number that can be
//  conclued from the tokens
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::loadData(){
   //FILE *fp = fopen(this->fileInName.c_str(),"r");
	fstream file(this->fileInName.c_str(),ios::in);
	double tmpFeature;
	int ind=0;
	while ( (file.eof()==0) && (ind<this->individuals.size()) && (this->features.size()>0) ){
		tVecDouble vFeatures(0);
		double labelFea=INDF;
		bool labelFeaRead=false;
		int f;
		for (f=0;f<this->features.size();f++){
			file>>tmpFeature;
			if (f+1==this->labelFeature){
              labelFea=tmpFeature;
			  labelFeaRead = true;
			}
		    if (!this->features[f]) continue;
			//even the label feature can be used to clustering purpose
               vFeatures.push_back(tmpFeature);
		}//for f

		//prevent from having to have the label feature as the first
		while(!labelFeaRead){
			file>>tmpFeature;
			if (f+1==this->labelFeature){
              labelFea=tmpFeature;
			  labelFeaRead = true;
			}
		}
		
		if (vFeatures.size()>0 && this->individuals[ind]){
		  this->tabData.push_back(vFeatures);
		  if (this->labelFeature>0){//this indicates either the presence of the label on the dataset or the desire of using it to work out the RC
		    this->individualLabels.push_back(labelFea);
		  }
		}//if
		//Descarde features not to be used for the individual
        while ( file.get()!='\n' && !file.eof())
			;
		ind++;
	}//while 
    
   #ifdef DEBUG_LOAD_DATA
	  cout<<endl<<"Individuals and their features"<<endl;
	  for (int i=0;i<this->tabData.size();i++){
		  tVecDouble vFeatures = tabData.at(i);
		  cout<<(i+1)<<"-> ";
		  for (int f=0;f<vFeatures.size();f++){
			  cout<<vFeatures[f]<<" | ";
		  }//for f
		  cout<<endl;
	  }//for i
	  cout<<endl<<endl<<"Label Feature"<<endl;
	  cout<<"Individual | Class"<<endl;
	  for (int l=0;l<this->individualLabels.size();l++){
        cout<<(l+1)<<" | "<<this->individualLabels[l]<<endl;
	  }//for l
   #endif
};//method

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void Application::readTokens( string tokens,tFuncProcessToken func){
  string strDelimiters(", ");//the delimeters are 
  char *token;
  char *tmpStr = strdup(tokens.c_str());
  token=strtok(tmpStr,strDelimiters.c_str());
   while ( token!=NULL){
	  (this->*func)(token);
   	  token = strtok(NULL,strDelimiters.c_str());
   };
   free(tmpStr);
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1)
//   (n)
// Coments:
//    It tries to forecast the number of individuals the file will have by resizing the vector to the maximum element number that can be
//  conclued from the tokens
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::processIndividualTokens(char *token){
  //static char *lastToken=NULL;
  int begin,end;
  if ( sscanf(token,"%d-%d",&begin,&end)==2){
	  int max = (begin>end?begin:end);
      if (this->individuals.size()<max)
		  this->individuals.resize(max,false);
	  for (int i=begin-1;i<end;i++)
		  this->individuals[i] = true;
  }else if ( sscanf(token,"%d",&begin)==1 ){
     if (this->individuals.size()<begin)
	   this->individuals.resize(begin,false);
	 this->individuals[begin-1] = true;
  }else{
    cerr<<endl<<"Erro parsing individuals number to be read";
  }
};
	
//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void Application::readIndividuals(){
	//string auxStr;
	//TODO:
	//Not use pre-defined line size
	char  str[600];
	cout<<endl<<"Choose the individuals from the file:"<<endl;
	cin.getline(str,600);
	//long flags = cin.flags();
	//long nskipws = ~ios::skipws;
	//cin.setf(flags & nskipws);
    //cin>>auxStr;
	//cin.setf(flags);
	string auxStr(str);	
//tFuncProcessToken f = Application::processIndividualTokens;
	this->readTokens(auxStr,&Application::processIndividualTokens);
};

//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void Application::processFeatureTokens(char *token){
int begin,end;
  if ( sscanf(token,"%d-%d",&begin,&end)==2){
	  int max = (begin>end?begin:end);
      if (this->features.size()<max)
		  this->features.resize(max,false);
	  for (int i=begin-1;i<end;i++)
		  this->features[i] = true;
  }else if ( sscanf(token,"%d",&begin)==1 ){
     if (this->features.size()<begin)
	   this->features.resize(begin,false);
	 this->features[begin-1] = true;
  }else{
    cerr<<endl<<"Erro parsing features number to be read";
  }
};
//////////////////////////////////////////////////////////////////////
//
//////////////////////////////////////////////////////////////////////
void Application::readFeatures(){
	cout<<endl<<"Choose the features from the file:"<<endl;
    //string auxStr;
	//cin>>auxStr;
	char str[600];
	cin.getline(str,600);
	string auxStr(str);
	this->readTokens((string)auxStr,&Application::processFeatureTokens);
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1)
//   (n)
// Coments:
//    XXX
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::readAlgorithm(){
	short tmpAlg;
	do{
		cout<<endl<<"Which algorithm do you like to use?"<<endl;
		cout<<"[1] - Fuzzy c-means non adaptative"<<endl;
		cout<<"[2] - Fuzzy c-means adaptative1"<<endl;
        cout<<"[3] - Fuzzy c-means adaptative2"<<endl;
        cout<<"[4] - k-means"<<endl;
		cout<<"Type your option(1-3): ";
		cin>>tmpAlg;
	}while( tmpAlg<1 || tmpAlg>4);
	switch(tmpAlg){
	case 1:
		this->algorithmDesc = "Fuzzy c-means non adaptative";
		break;
	case 2:
		this->algorithmDesc = "Fuzzy c-means adaptative1";
		break;
	case 3:
		this->algorithmDesc = "Fuzzy c-means adaptative2";
		break;
	case 4:
		this->algorithmDesc = "k-means";
		break;
	};
	this->algorithm = tmpAlg;
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1)
//   (n)
// Coments:
//    XXX
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::readMParameter(){
	double tmpM;
	//Accepts 1 as input in order to simulate KMeans by a FCM implementation
	do{
		cout<<endl<<"What's the parameter M value? "<<endl;
		cin>>tmpM;
	}while(tmpM<1);
	this->parameterM = tmpM;
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1)
//   (n)
// Coments:
//    XXX
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::readRuns(){
  int tmpRuns;
  do{
	  cout<<endl<<"How many runs do you like?"<<endl;
	  cin>>tmpRuns;
  }while(tmpRuns<1);
  this->runs = tmpRuns;
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1)
//   (n)
// Coments:
//    XXX
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::readIterations(){
  int tmpIte;
  do{
	  cout<<endl<<"How many iterations per run do you like?"<<endl;
	  cin>>tmpIte;
  }while(tmpIte<1);
  this->iterations = tmpIte;
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1)
//   (n)
// Coments:
//    XXX
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::readFileOutName(){
  cout<<endl<<endl<<"Please, inform a name to the output file?"<<endl;
  string strOut;
  cin>>strOut;
  this->fileOutName = strOut;
};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1)
//   (n)
// Post-conditions:
//   (1) the field labelFeature must have a value different from -1 if the labelFeature is present on the dataset and will be used to work out the RC
//       labelFeature must be -1 otherwise
//   (n)
// Coments:
//    XXX
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::readLabelFeature(){
	int tmpLabel=-1;
	bool hasLabel=false;
	char yn;
	do{
		cout<<endl<<"Does the base have a label feature and you like to use it to work out the RC?(Y/N) "<<endl;
        cin>>yn;
	}while(yn!='Y' && yn!='N');
	hasLabel = yn=='Y';
	if (hasLabel){
	  do{
		cout<<endl<<"What is the label feature? "<<endl;
		cin>>tmpLabel;
	  }while(tmpLabel<1);
	}//if
	this->labelFeature = tmpLabel; 
};
 //pre-conditions: -the user has to be previously asked to want or not labelFeature
void Application::readNumClasses(){
  int tmpNClasses=-1;
  if (this->labelFeature<0)
	  return;
  do{
	cout<<endl<<"How many classes does have the data being used? "<<endl;
	cin>>tmpNClasses;
  }while(tmpNClasses<1);
  this->nClasses = tmpNClasses;

};

////////////////////////////////////////////////////////////////////////////////////////
// Pre-conditions:
//   (1) the elements number must be known
//   (n)
// Post-conditions:
//   (1) 
//   (n)
// Coments:
//    XXX
// Parameters:
//   -param1:
//       XXXXX
//   -paramN:
//       XXXXX
// Return:
//    XXXXX
////////////////////////////////////////////////////////////////////////////////////////
void Application::readClusterNumber(){
  int tmpClust;
  do{
	  cout<<endl<<"How many clusters does the dataset have? "<<endl;
	  cin>>tmpClust;
  }while(tmpClust<2 || tmpClust>=this->individuals.size());
  this->nClusters = tmpClust;
};