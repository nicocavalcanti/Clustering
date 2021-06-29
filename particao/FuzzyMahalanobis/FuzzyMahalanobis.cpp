// FuzzyMahalanobis.cpp: implementation of the FuzzyMahalanobis class.
//
//////////////////////////////////////////////////////////////////////

#include "FuzzyMahalanobis.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyMahalanobis::FuzzyMahalanobis()
{

}

FuzzyMahalanobis::~FuzzyMahalanobis()
{
	for(int k=0;k<totalPadroes;k++)
	{
		delete melhorParticao[k];
	}
	delete melhorParticao;

	delete z;

	for(int i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			delete melhorM[i][j];
		}
		delete melhorM[i];
	}
	delete melhorM;

	delete Ti;

	delete Ji;

	delete mi;

	delete Bi;

	delete relativeT;

	delete relativeJ;

	delete relativeB;
}

FuzzyMahalanobis::FuzzyMahalanobis(Pattern *p,char *nDados,char *nSaida) : Particao(p,nDados,nSaida)
{
	melhorParticao=new double*[totalPadroes];
	for(int k=0;k<totalPadroes;k++)
	{
		melhorParticao[k]=new double[numeroClasses];
	}

	z=new double[numeroVariaveis];

	melhorM=new double**[numeroClasses];
	for(int i=0;i<numeroClasses;i++)
	{
		melhorM[i]=new double*[numeroVariaveis];
		for(int j=0;j<numeroVariaveis;j++)
		{
			melhorM[i][j]=new double[numeroVariaveis];
		}
	}

	Ti=new double[numeroClasses];

	Ji=new double[numeroClasses];

	mi=new double[numeroClasses];

	Bi=new double[numeroClasses];

	relativeT=new double[numeroClasses];

	relativeJ=new double[numeroClasses];

	relativeB=new double[numeroClasses];
}

//////////////////////////////////////////////////////////////////////
//Esta função inicializa o algoritmo

void FuzzyMahalanobis::inicializacao()
{
	//Inicialização dos graus de pertinência
	for(int i=0;i<totalPadroes;i++)
	{
		double limite=1.0;

		for(int j=0;j<numeroClasses-1;j++)
		{
			double semente=(rand()%100)/100.0;

			double valor=semente*limite;

			pattern[i].setGrauPertinencia(j,valor);

			limite-=valor;
		}

		pattern[i].setGrauPertinencia(j,limite);

	}

/*	////////////////////////////////////////////
	//Monitoramento
	ofstream fout("Graus de Pertinencia Iniciais.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"\n\n";

	for(i=0;i<totalPadroes;i++)
	{
		double total=0.0;

		double *gp=pattern[i].getGrauPertinencia();

		for(int j=0;j<numeroClasses;j++)
		{
			fout<<setw(12)<<gp[j];
			total+=gp[j];
		}

		fout<<"  Total: "<<setw(3)<<total;
		fout<<"\n";
	}

	fout.close();

	////////////////////////////////////////////*/

/*	////////////////////////////////////////////
	//Monitoramento
	ofstream tout("Classes Iniciais.wri",ios::ate);
	tout<<"\n\n";
	for(int k=0;k<totalPadroes;k++)
	{
		double max=0.0;
		int c;

		for(int i=0;i<numeroClasses;i++)
		{
			if(pattern[k].getGrauPertinencia(i)>0.5)
			{
				c=i;
				break;
			}
			else
			{
				if(pattern[k].getGrauPertinencia(i)>max)
				{
					c=i;
					max=pattern[k].getGrauPertinencia(i);
				}
			}
		}
	
		tout<<"padrão "<<k<<": "<<c<<"\n";

	}
	////////////////////////////////////////////*/

}

//////////////////////////////////////////////////////////////////////
//Esta etapa atualiza os graus de pertinência de cada indivíduo

void FuzzyMahalanobis::alocacao()
{
	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			double numerador=distanciaMahalanobis(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getM());
				
			double acumulado=0.0;

			for(int h=0;h<numeroClasses;h++)
			{
				double denominador=distanciaMahalanobis(pattern[k].getX(),cluster[h].getPrototipo(),cluster[h].getM());

				acumulado+=pow((numerador/denominador),(1.0/(m-1)));
			}

			pattern[k].setGrauPertinencia(i,pow(acumulado,-1));
		}
	}

/*	////////////////////////////////////////////
	//Monitoramento
	ofstream fout("Graus de Pertinencia.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"\n---------- Iteração "<<iteracao<<" ---------\n";

	for(k=0;k<totalPadroes;k++)
	{
		double total=0.0;

		double *gp=pattern[k].getGrauPertinencia();

		for(int i=0;i<numeroClasses;i++)
		{
			fout<<setw(12)<<gp[i];
			total+=gp[i];
		}

		fout<<"  Total: "<<setw(3)<<total;
		fout<<"\n";
	}

	fout.close();
	////////////////////////////////////////////*/
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::atualizaMelhorParticao()
{
	for(int	k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			melhorParticao[k][i]=pattern[k].getGrauPertinencia(i);
		}
	}

	for(int i=0;i<numeroClasses;i++)
	{
		double **M=cluster[i].getM();
		
		for(int j=0;j<numeroVariaveis;j++)
		{
			for(int l=0;l<numeroVariaveis;l++)
			{
				melhorM[i][j][l]=M[j][l];
			}
		}
	}
}

//////////////////////////////////////////////////////////////////////
//Essa função calcula a quantidade final de padrões em cada classe

int * FuzzyMahalanobis::classesFinais()
{
	int *classesDosPadroes=new int[totalPadroes];

	for(int k=0;k<totalPadroes;k++)
	{
		double max=0.0;
		int c;

		for(int i=0;i<numeroClasses;i++)
		{
			if(melhorParticao[k][i]>0.5)
			{
				c=i;
				break;
			}
			else
			{
				if(melhorParticao[k][i]>max)
				{
					c=i;
					max=melhorParticao[k][i];
				}
			}
		}
	
		classesDosPadroes[k]=c;
		pattern[k].setClasse(c);

	}

	return classesDosPadroes;
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaZ()
{
	for(int j=0;j<numeroVariaveis;j++)
	{
		double numerador=0.0;
		double denominador=0.0;

		for(int i=0;i<numeroClasses;i++)
		{
			for(int k=0;k<totalPadroes;k++)
			{
				double u=pattern[k].getGrauPertinencia(i);
				double g=cluster[i].getPrototipo(j);
				
				numerador+=(pow(u,m))*g;
				denominador+=pow(u,m);
			}
		}

		z[j]=numerador/denominador;
	}

/*	/////////////////////////////
	//Monitoramento
	ofstream fout("Z.wri",ios::ate);

	for(j=0;j<numeroVariaveis;j++)
	{
		fout<<z[j]<<"\t";
	}
	fout<<"\n";

	/////////////////////////////*/
}


//////////////////////////////////////////////////////////////////////

double FuzzyMahalanobis::distanciaMahalanobis(double *x, double *g, double **M)
{
	double distancia=0.0;

	for(int i=0;i<numeroVariaveis;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			distancia+=(x[i]-g[i])*(x[j]-g[j])*M[i][j];
		}
	}

	return distancia;
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaDC()
{
	double Fc_numerador=0.0;

	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			double u=pattern[k].getGrauPertinencia(i);

			Fc_numerador+=pow(u,2);
		}
	}

	double Fc=Fc_numerador/totalPadroes;

	DC=(numeroClasses*Fc-1.0)/(numeroClasses-1.0);
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::setMelhoresIndices()
{
	melhorJ=J;
	melhorT=T;
	melhorB=B;
	melhorR=R;
	melhorCR=CR;
	melhorDC=DC;
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::setClassificacaoHard()
{
	for(int k=0;k<totalPadroes;k++)
	{
		double max=0.0;
		int c;
		double *u=pattern[k].getGrauPertinencia();

		for(int i=0;i<numeroClasses;i++)
		{
			if(u[i]>0.5)
			{
				c=i;
				break;
			}
			else
			{
				if(u[i]>max)
				{
					c=i;
					max=u[i];
				}
			}
		}
	
		pattern[k].setClasse(c);
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaIndices()
{
	calculaZ();
	calculaJ();
	calculaT();
	calculaMi();
	calculaB();
	calculaR();
	calculaCorrectedRandIndex();
	calculaDC();
}

//////////////////////////////////////////////////////////////////////

double FuzzyMahalanobis::getMelhorDC()
{
	return melhorDC;
}

//////////////////////////////////////////////////////////////////////

double FuzzyMahalanobis::getMelhorJ()
{
	return melhorJ;
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::setMelhorParticao(int *classes)
{
	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			pattern[k].setGrauPertinencia(i,melhorParticao[k][i]);
		}
	}

	for(int i=0;i<numeroClasses;i++)
	{
		cluster[i].setM(melhorM[i]);
	}

	for(i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoFuzzy(pattern);
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaT()
{
	T=0.0;

	for(int i=0;i<numeroClasses;i++)
	{
		for(int k=0;k<totalPadroes;k++)
		{
			double u=pattern[k].getGrauPertinencia(i);
			
			double distancia=distanciaMahalanobis(pattern[k].getX(),z,cluster[i].getM());

			T+=distancia*pow(u,m);
		}
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaTi()
{
	//Cálculo do class-specific
	for(int i=0;i<numeroClasses;i++)
	{
		Ti[i]=0.0;

		for(int k=0;k<totalPadroes;k++)
		{
			double u=pattern[k].getGrauPertinencia(i);
			
			double distancia=distanciaMahalanobis(pattern[k].getX(),z,cluster[i].getM());

			Ti[i]+=distancia*pow(u,m);
		}
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaJi()
{
	for(int i=0;i<numeroClasses;i++)
	{
		Ji[i]=0.0;
		for(int k=0;k<totalPadroes;k++)
		{
			double u=pattern[k].getGrauPertinencia(i);
			double distancia=distanciaMahalanobis(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getM());

			Ji[i]+=distancia*pow(u,m);
		}
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaMi()
{
	for(int i=0;i<numeroClasses;i++)
	{
		mi[i]=0.0;
		for(int k=0;k<totalPadroes;k++)
		{
			double u=pattern[k].getGrauPertinencia(i);

			mi[i]+=pow(u,m);
		}
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaB()
{
	B=0.0;
	for(int i=0;i<numeroClasses;i++)
	{
		double distancia=distanciaMahalanobis(z,cluster[i].getPrototipo(),cluster[i].getM());		

		B+=mi[i]*distancia;
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaBi()
{
	for(int i=0;i<numeroClasses;i++)
	{
		double distancia=distanciaMahalanobis(z,cluster[i].getPrototipo(),cluster[i].getM());		

		Bi[i]=mi[i]*distancia;
	}	
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaR()
{
	R=B/T;
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaRelativeT()
{
	for(int i=0;i<numeroClasses;i++)
	{
		relativeT[i]=Ti[i]/melhorT;
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaRelativeB()
{
	for(int i=0;i<numeroClasses;i++)
	{
		relativeB[i]=Bi[i]/melhorB;
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaRelativeJ()
{
	for(int i=0;i<numeroClasses;i++)
	{
		relativeJ[i]=Ji[i]/melhorJ;
	}
}

//////////////////////////////////////////////////////////////////////

double FuzzyMahalanobis::getMelhorT()
{
	return melhorT;
}

//////////////////////////////////////////////////////////////////////

double FuzzyMahalanobis::getMelhorR()
{
	return melhorR;
}

//////////////////////////////////////////////////////////////////////

double FuzzyMahalanobis::getMelhorB()
{
	return melhorB;
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::imprimeInterpretadores()
{
	ofstream saida;
	double soma;
	
	/////////////////////////////
	//Criação do arquivo para impressão dos indices interpretadores dos clusters
	char extensao2[35];
	strcpy(extensao2,"_ClusterResult_Interpretation.wri");
	char *arquivoVariaveis=new char[95];
	strcpy(arquivoVariaveis,nomeSaida);
	strncat(arquivoVariaveis,extensao2,34);
	saida.open(arquivoVariaveis);
	/////////////////////////////
	
//	saida<<setiosflags(ios::fixed);

	saida<<"Results: \n\n";
	saida<<setw(25)<<"Algorithm: "<<"Fuzzy c-means";
	saida<<"\n"<<setw(25)<<"Parameter m: "<<m;
	saida<<"\n"<<setw(25)<<"Clusters: "<<numeroClasses;
	saida<<"\n"<<setw(25)<<"Initialization Number: "<<numeroInicializacoes;
	saida<<"\n"<<setw(25)<<"Iteration Limit: "<<limiteIteracao;
	saida<<"\n"<<setw(25)<<"Variables Number: "<<numeroVariaveis;

	saida<<"\n\n/-------------------------------------  Interpretation -------------------------------------/";

	saida<<"\n\n/------------------------------  Global Inertia - T ------------------------------/";

	saida<<"\nData set global inertia = "<<melhorT;

	saida<<"\n\nThe global inertia for each cluster";
	saida<<"\n\n===========================\n";
	saida<<setw(7)<<"CLUSTER"<<setw(20)<<"VALUE";
	saida<<"\n---------------------------\n";
	soma=0.0;
	for(int i=0;i<numeroClasses;i++)
	{
		saida<<setw(7)<<i+1<<setw(20)<<Ti[i];
		saida<<"\n";
		soma+=Ti[i];
	}
	saida<<"===========================\n";
	saida<<setw(7)<<"Sum"<<setw(20)<<soma;


	saida<<"\n\n/------------------------------  Within-cluster Inertia - J ------------------------------/";
	saida<<"\nData set within-cluster inertia = "<<melhorJ;

	saida<<"\n\nThe within-cluster inertia for each cluster";
	saida<<"\n\n===========================\n";
	saida<<setw(7)<<"CLUSTER"<<setw(20)<<"VALUE";
	saida<<"\n---------------------------\n";
	soma=0.0;
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(7)<<i+1<<setw(20)<<Ji[i];
		saida<<"\n";
		soma+=Ji[i];
	}
	saida<<"===========================\n";
	saida<<setw(7)<<"Sum"<<setw(20)<<soma;

	
	saida<<"\n\n/------------------------------  Between-cluster Inertia - B ------------------------------/";
	saida<<"\nData set between-cluster inertia = "<<melhorB;

	saida<<"\n\nThe between-cluster inertia for each cluster";
	saida<<"\n\n===========================\n";
	saida<<setw(7)<<"CLUSTER"<<setw(20)<<"VALUE";
	saida<<"\n---------------------------\n";
	soma=0.0;
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(7)<<i+1<<setw(20)<<Bi[i];
		saida<<"\n";
		soma+=Bi[i];
	}
	saida<<"===========================\n";
	saida<<setw(7)<<"Sum"<<setw(20)<<soma;

	saida<<"\n\n/------------------------------  General Index - R ------------------------------/";
	saida<<"\nThe proportion of inertia explaines by the cluster = "<<melhorR;

	saida<<"\n\n/------------------------------  Cluster description ------------------------------/";
	
	saida<<"\n\nThe proportion of the global inertia explained by clusters - T(i)";
	saida<<"\n\n===========================\n";
	saida<<setw(7)<<"CLUSTER"<<setw(20)<<"VALUE";
	saida<<"\n---------------------------\n";
	soma=0.0;
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(7)<<i+1<<setw(20)<<relativeT[i];
		saida<<"\n";
		soma+=relativeT[i];
	}
	saida<<"===========================\n";
	saida<<setw(7)<<"Sum"<<setw(20)<<soma;


	saida<<"\n\nThe relative contribution of clusters to the between-cluster inertia - B(i)";
	saida<<"\n\n===========================\n";
	saida<<setw(7)<<"CLUSTER"<<setw(20)<<"VALUE";
	saida<<"\n---------------------------\n";
	soma=0.0;
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(7)<<i+1<<setw(20)<<relativeB[i];
		saida<<"\n";
		soma+=relativeB[i];
	}
	saida<<"===========================\n";
	saida<<setw(7)<<"Sum"<<setw(20)<<soma;


	saida<<"\n\nThe relative contribution of clusters to the within-cluster inertia - J(i)";
	saida<<"\n\n===========================\n";
	saida<<setw(7)<<"CLUSTER"<<setw(20)<<"VALUE";
	saida<<"\n---------------------------\n";
	soma=0.0;
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(7)<<i+1<<setw(20)<<relativeJ[i];
		saida<<"\n";
		soma+=relativeJ[i];
	}
	saida<<"===========================\n";
	saida<<setw(7)<<"Sum"<<setw(20)<<soma;

	saida<<"\n\n/------------------------------  Classification Error ------------------------------/";


	saida<<"\n\nThe classification error for each cluster (relative to the a priori classification)";
	saida<<"\n\n===========================\n";
	saida<<setw(7)<<"CLUSTER"<<setw(20)<<"VALUE";
	saida<<"\n---------------------------\n";
	soma=0.0;
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(7)<<i+1<<setw(20)<<erroClasse[i];
		saida<<"\n";
	}
	saida<<"===========================\n";
	saida<<setw(7)<<"Global"<<setw(20)<<erroGlobal;

}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaIndicesClasseVariavel()
{
	calculaZ();
	calculaMi();

	calculaJi();

	calculaTi();

	calculaBi();

	calculaRelativeB();
	calculaRelativeJ();
	calculaRelativeT();

	calculaErroClassificacao();

}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::imprimeVariaveisMonitoramento()
{
	ofstream fout;

	///////////////////////////////////////////////////////

	fout.open("Prototipos.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroVariaveis+3+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	fout<<setw(3)<<" ";
	for(int j=0;j<numeroVariaveis;j++)
	{
		fout<<setw(12)<<"Variavel "<<setw(3)<<j+1;
	}
	fout<<"\n";
	fout<<setw(15*numeroVariaveis+3+1)<<setfill('-')<<" "<<setfill(' ');
	fout<<"\n";

	double *prototipo;
	
	for(int i=0;i<numeroClasses;i++)
	{
		prototipo=cluster[i].getPrototipo();
	
		fout<<setw(3)<<i+1;

		for(int j=0;j<numeroVariaveis;j++)
		{
			fout<<setw(15)<<prototipo[j];
		}

		fout<<"\n";
	}
	
	fout<<setw(15*numeroVariaveis+3+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////

	fout.open("Graus de Pertinência.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroClasses+4+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	fout<<setw(4)<<"k";
	for(i=0;i<numeroClasses;i++)
	{
		fout<<setw(14)<<"u"<<i+1;
	}
	fout<<"\n";
	fout<<setw(15*numeroClasses+4+1)<<setfill('-')<<" "<<setfill(' ');
	fout<<"\n";

	for(int k=0;k<totalPadroes;k++)
	{
		fout<<setw(4)<<k+1;

		double *gp=pattern[k].getGrauPertinencia();

		for(int i=0;i<numeroClasses;i++)
		{
			fout<<setw(15)<<gp[i];
		}

		fout<<"\n";
	}

	fout<<setw(15*numeroClasses+4+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////

	fout.open("Z.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroVariaveis+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		fout<<setw(12)<<"Variavel"<<setw(3)<<j+1;
	}
	fout<<"\n";
	fout<<setw(15*numeroVariaveis+1)<<setfill('-')<<" "<<setfill(' ');
	fout<<"\n";
	
	for(j=0;j<numeroVariaveis;j++)
	{
		fout<<setw(15)<<z[j];
	}

	fout<<"\n";

	fout<<setw(15*numeroVariaveis+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////

	fout.open("Mi.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroClasses+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	for(i=0;i<numeroClasses;i++)
	{
		fout<<setw(14)<<"Classe"<<i+1;
	}
	fout<<"\n";
	fout<<setw(15*numeroClasses+1)<<setfill('-')<<" "<<setfill(' ');
	fout<<"\n";

	for(i=0;i<numeroClasses;i++)
	{
		fout<<setw(15)<<mi[i];
	}

	fout<<"\n";

	fout<<setw(15*numeroClasses+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////


	///////////////////////////////////////////////////////

	fout.open("Matriz M.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroVariaveis+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	fout<<setw(15*numeroVariaveis+1)<<setfill('-')<<" "<<setfill(' ');

	for(i=0;i<numeroClasses;i++)
	{
		fout<<"\nClasse "<<i+1;
		fout<<"\n";
		double **M=cluster[i].getM();
		for(int j=0;j<numeroVariaveis;j++)
		{
			for(int l=0;l<numeroVariaveis;l++)
			{
				fout<<setw(15)<<M[j][l];
			}
			fout<<"\n";
		}

		fout<<"\n";
	}

	fout<<setw(15*numeroVariaveis+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////

	fout.open("Distancias Padrão_Prototipo.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroClasses+4+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	fout<<setw(4)<<"k";
	for(i=0;i<numeroClasses;i++)
	{
		fout<<setw(14)<<"u"<<i+1;
	}
	fout<<"\n";
	fout<<setw(15*numeroClasses+4+1)<<setfill('-')<<" "<<setfill(' ');
	fout<<"\n";

	for(k=0;k<totalPadroes;k++)
	{
		fout<<setw(4)<<k+1;

		for(int i=0;i<numeroClasses;i++)
		{
			double distancia=distanciaMahalanobis(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getM());
			fout<<setw(15)<<distancia;
		}

		fout<<"\n";
	}

	fout<<setw(15*numeroClasses+4+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////

	///////////////////////////////////////////////////////

	fout.open("Distancias Prototipo_Z.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroClasses+4+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	fout<<setw(4)<<"k";
	for(i=0;i<numeroClasses;i++)
	{
		fout<<setw(14)<<"u"<<i+1;
	}
	fout<<"\n";
	fout<<setw(15*numeroClasses+4+1)<<setfill('-')<<" "<<setfill(' ');
	fout<<"\n";

	for(i=0;i<numeroClasses;i++)
	{
		fout<<setw(4)<<i+1;

		for(int r=0;r<numeroClasses;r++)
		{
			double distancia=distanciaMahalanobis(z,cluster[i].getPrototipo(),cluster[i].getM());
			fout<<setw(15)<<distancia;
		}

		fout<<"\n";
	}

	fout<<setw(15*numeroClasses+4+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////
}

//////////////////////////////////////////////////////////////////////