// FuzzyCMeans.cpp: implementation of the FuzzyCMeans class.
//
//////////////////////////////////////////////////////////////////////

#include "FuzzyCMeans.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyCMeans::FuzzyCMeans()
{

}

FuzzyCMeans::~FuzzyCMeans()
{
	for(int k=0;k<totalPadroes;k++)
	{
		delete melhorParticao[k];
	}
	delete melhorParticao;

	for(int i=0;i<numeroClasses;i++)
	{
		delete melhorWeightVector[i];
	}
	delete melhorWeightVector;

	delete z;

	for(i=0;i<numeroClasses;i++)
	{
		delete Bij[i];
	}

	delete Bij;

	delete Bj;

	delete Bi;
	delete Ji;
	delete Ti;

	delete CORj;

	for(i=0;i<numeroClasses;i++)
	{
		delete CORij[i];
	}
	delete CORij;

	delete CTRj;

	for(i=0;i<numeroClasses;i++)
	{
		delete CTRij[i];
	}
	delete CTRij;

	delete Tj;

	for(i=0;i<numeroClasses;i++)
	{
		delete CEij[i];
	}

	delete CEij;

	for(i=0;i<numeroClasses;i++)
	{
		delete Tij[i];
	}
	delete Tij;

	for(i=0;i<numeroClasses;i++)
	{
		delete Jij[i];
	}
	delete Jij;

	delete Jj;

	delete mi;

	delete relativeT;

	delete relativeJ;

	delete relativeB;
}

FuzzyCMeans::FuzzyCMeans(Pattern *p,char *nDados,char *nSaida) : Particao(p,nDados,nSaida)
{
	melhorParticao=new double*[totalPadroes];
	for(int k=0;k<totalPadroes;k++)
	{
		melhorParticao[k]=new double[numeroClasses];
	}

	melhorWeightVector=new double*[numeroClasses];
	for(int i=0;i<numeroClasses;i++)
	{
		melhorWeightVector[i]=new double[numeroVariaveis];
	}

	z=new double[numeroVariaveis];

	Ji=new double[numeroClasses];
	Bi=new double[numeroClasses];
	Ti=new double[numeroClasses];

	Bij=new double*[numeroClasses];
	for(i=0;i<numeroClasses;i++)
	{
		Bij[i]=new double[numeroVariaveis];
	}

	Bj=new double[numeroVariaveis];

	CORj=new double[numeroVariaveis];
	
	CORij=new double*[numeroClasses];
	for(i=0;i<numeroClasses;i++)
	{
		CORij[i]=new double[numeroVariaveis];
	}

	CTRj=new double[numeroVariaveis];
	
	CTRij=new double*[numeroClasses];
	for(i=0;i<numeroClasses;i++)
	{
		CTRij[i]=new double[numeroVariaveis];
	}

	Tj=new double[numeroVariaveis];

	CEij=new double*[numeroClasses];
	for(i=0;i<numeroClasses;i++)
	{
		CEij[i]=new double[numeroVariaveis];
	}

	Tij=new double*[numeroClasses];
	for(i=0;i<numeroClasses;i++)
	{
		Tij[i]=new double[numeroVariaveis];
	}

	Jij=new double*[numeroClasses];
	for(i=0;i<numeroClasses;i++)
	{
		Jij[i]=new double[numeroVariaveis];
	}

	Jj=new double[numeroVariaveis];

	mi=new double[numeroClasses];

	relativeT=new double[numeroClasses];

	relativeJ=new double[numeroClasses];

	relativeB=new double[numeroClasses];
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::atualizaMelhorParticao()
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
		double *wV=cluster[i].getWeightVector();

		for(int j=0;j<numeroVariaveis;j++)
		{
			melhorWeightVector[i][j]=wV[j];
		}
	}
}

//////////////////////////////////////////////////////////////////////

int * FuzzyCMeans::classesFinais()
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

void FuzzyCMeans::inicializacao()
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

void FuzzyCMeans::alocacao()
{
	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			double numerador=distanciaL2MinkowskyAdaptativo(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getWeightVector());
				
			double acumulado=0.0;

			for(int h=0;h<numeroClasses;h++)
			{
				double denominador=distanciaL2MinkowskyAdaptativo(pattern[k].getX(),cluster[h].getPrototipo(),cluster[h].getWeightVector());

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

double FuzzyCMeans::calculaW()
{
	double W=0.0;

	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			double grauPertinencia=pattern[k].getGrauPertinencia(i);
			double distancia=distanciaL2MinkowskyAdaptativo(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getWeightVector());

			W+=pow(grauPertinencia,m)*distancia;
		}
	}

	return W;
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaZ()
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

/*	////////////////////////////
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
//Overall sum of squares

void FuzzyCMeans::calculaT()
{
	///////////////////////////
	T=0.0;

	for(int i=0;i<numeroClasses;i++)
	{
		for(int k=0;k<totalPadroes;k++)
		{
			double u=pattern[k].getGrauPertinencia(i);
			
			double distancia=distanciaL2MinkowskyAdaptativo(pattern[k].getX(),z,cluster[i].getWeightVector());

			T+=distancia*pow(u,m);
		}
	}
	///////////////////////////
}

//////////////////////////////////////////////////////////////////////
//Within-class sum of squares

void FuzzyCMeans::calculaJ()
{
	J=0.0;

	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			double grauPertinencia=pattern[k].getGrauPertinencia(i);
			double distancia=distanciaL2MinkowskyAdaptativo(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getWeightVector());

			J+=pow(grauPertinencia,m)*distancia;
		}
	}
}

//////////////////////////////////////////////////////////////////////
//Between-class sum of squares

void FuzzyCMeans::calculaB()
{
	/////////////////////
	B=0.0;
	for(int i=0;i<numeroClasses;i++)
	{
		double distancia=distanciaL2MinkowskyAdaptativo(z,cluster[i].getPrototipo(),cluster[i].getWeightVector());		

		B+=mi[i]*distancia;
	}

	/////////////////////
}

//////////////////////////////////////////////////////////////////////

double FuzzyCMeans::distanciaL2MinkowskyAdaptativo(double *x, double *g, double *weightVector)
{
	double distancia=0.0;
	
	for(int j=0;j<numeroVariaveis;j++)
	{
		distancia+=weightVector[j]*pow((x[j]-g[j]),2);
	}

	return distancia;
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaR()
{
	R=B/T;
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaDC()
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

void FuzzyCMeans::calculaIndices()
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

void FuzzyCMeans::setMelhorParticao(int *classes)
{
	for(int i=0;i<numeroClasses;i++)
	{
		cluster[i].setWeightVector(melhorWeightVector[i]);
	}

	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			pattern[k].setGrauPertinencia(i,melhorParticao[k][i]);
		}
	}

	for(i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoFuzzy(pattern);
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::setMelhoresIndices()
{
	melhorJ=J;
	melhorT=T;
	melhorB=B;
	melhorR=R;
	melhorCR=CR;
	melhorDC=DC;
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaCORj()
{
	for(int j=0;j<numeroVariaveis;j++)
	{
		CORj[j]=Bj[j]/Tj[j];
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaCORij()
{
	for(int i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			CORij[i][j]=Bij[i][j]/Tj[j];
		}
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaCTRj()
{
	for(int j=0;j<numeroVariaveis;j++)
	{
		CTRj[j]=Bj[j]/melhorB;
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaCTRij()
{
	for(int i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			CTRij[i][j]=Bij[i][j]/Bi[i];
		}
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::setClassificacaoHard()
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

void FuzzyCMeans::calculaCEij()
{
	for(int i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			CEij[i][j]=Bij[i][j]/melhorB;
		}

	}
}

//////////////////////////////////////////////////////////////////////

double FuzzyCMeans::getMelhorJ()
{
	return melhorJ;
}

//////////////////////////////////////////////////////////////////////

double FuzzyCMeans::getMelhorT()
{
	return melhorT;
}

//////////////////////////////////////////////////////////////////////

double FuzzyCMeans::getMelhorR()
{
	return melhorR;		
}

//////////////////////////////////////////////////////////////////////

double FuzzyCMeans::getMelhorB()
{
	return melhorB;
}

//////////////////////////////////////////////////////////////////////

double FuzzyCMeans::getMelhorDC()
{
	return melhorDC;
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::imprimeInterpretadores()
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


	saida<<"\n\nThe global inertia for each variable";
	saida<<"\n\n============================\n";
	saida<<setw(8)<<"VARIABLE"<<setw(20)<<"VALUE";
	saida<<"\n----------------------------\n";
	soma=0.0;
	for(int j=0;j<numeroVariaveis;j++)
	{
		saida<<setw(8)<<j+1<<setw(20)<<Tj[j];
		saida<<"\n";
		soma+=Tj[j];
	}
	saida<<"============================\n";
	saida<<setw(8)<<"Sum"<<setw(20)<<soma;


	saida<<"\n\nThe global inertia for each variable and cluster";
	saida<<"\n\n"<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');
	saida<<"\n";
	saida<<setw(15)<<" ";
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(16)<<"CLUSTER"<<setw(4)<<i+1;
	}
	saida<<"\n"<<setw(20*numeroClasses+16)<<setfill('-')<<' '<<setfill(' ');
	saida<<"\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<"VARIABLE"<<setw(7)<<j+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(20)<<Tij[i][j];
		}
		saida<<"\n";
	}
	saida<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');



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

	saida<<"\n\nThe within-cluster inertia for each variable";
	saida<<"\n\n============================\n";
	saida<<setw(8)<<"VARIABLE"<<setw(20)<<"VALUE";
	saida<<"\n----------------------------\n";
	soma=0.0;
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<setw(8)<<j+1<<setw(20)<<Jj[j];
		saida<<"\n";
		soma+=Jj[j];
	}
	saida<<"============================\n";
	saida<<setw(8)<<"Sum"<<setw(20)<<soma;


	saida<<"\n\nThe within-cluster inertia for each variable and cluster";
	saida<<"\n\n"<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');
	saida<<"\n";
	saida<<setw(15)<<" ";
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(16)<<"CLUSTER"<<setw(4)<<i+1;
	}
	saida<<"\n"<<setw(20*numeroClasses+16)<<setfill('-')<<' '<<setfill(' ');
	saida<<"\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<"VARIABLE"<<setw(7)<<j+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(20)<<Jij[i][j];
		}
		saida<<"\n";
	}
	saida<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');

	
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

	saida<<"\n\nThe between-cluster inertia for each variable";
	saida<<"\n\n============================\n";
	saida<<setw(8)<<"VARIABLE"<<setw(20)<<"VALUE";
	saida<<"\n----------------------------\n";
	soma=0.0;
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<setw(8)<<j+1<<setw(20)<<Bj[j];
		saida<<"\n";
		soma+=Bj[j];
	}
	saida<<"============================\n";
	saida<<setw(8)<<"Sum"<<setw(20)<<soma;


	saida<<"\n\nThe between-cluster inertia for each variable and cluster";
	saida<<"\n\n"<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');
	saida<<"\n";
	saida<<setw(15)<<" ";
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(16)<<"CLUSTER"<<setw(4)<<i+1;
	}
	saida<<"\n"<<setw(20*numeroClasses+16)<<setfill('-')<<' '<<setfill(' ');
	saida<<"\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<"VARIABLE"<<setw(7)<<j+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(20)<<Bij[i][j];
		}
		saida<<"\n";
	}
	saida<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');

	saida<<"\n\n/------------------------------  General Index - R ------------------------------/";
	saida<<"\nThe proportion of inertia explaines by the cluster = "<<melhorR;

	saida<<"\n\n/------------------------------  Variables contribution ------------------------------/";

	saida<<"\n\nThe proportion of inertia of the variables taken into account by the cluster - COR(j)";
	saida<<"\n\n============================\n";
	saida<<setw(8)<<"VARIABLE"<<setw(20)<<"VALUE";
	saida<<"\n----------------------------\n";
	soma=0.0;
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<setw(8)<<j+1<<setw(20)<<CORj[j];
		saida<<"\n";
		soma+=CORj[j];
	}
	saida<<"============================\n";
	saida<<setw(8)<<"Sum"<<setw(20)<<soma;

	saida<<"\n\nThe relative contribution of the variables to the betwwen-class inertia - CTR(j)";
	saida<<"\n\n============================\n";
	saida<<setw(8)<<"VARIABLE"<<setw(20)<<"VALUE";
	saida<<"\n----------------------------\n";
	soma=0.0;
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<setw(8)<<j+1<<setw(20)<<CTRj[j];
		saida<<"\n";
		soma+=CTRj[j];
	}
	saida<<"============================\n";
	saida<<setw(8)<<"Sum"<<setw(20)<<soma;

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

	saida<<"\n\n/------------------------------  Cluster description by variables------------------------------/";
	
	saida<<"\n\nThe proportion of the discriminant power of variables taken into account by clusters - COR(j,i)";
	saida<<"\n\n"<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');
	saida<<"\n";
	saida<<setw(15)<<" ";
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(16)<<"CLUSTER"<<setw(4)<<i+1;
	}
	saida<<"\n"<<setw(20*numeroClasses+16)<<setfill('-')<<' '<<setfill(' ');
	saida<<"\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<"VARIABLE"<<setw(7)<<j+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(20)<<CORij[i][j];
		}
		saida<<"\n";
	}
	saida<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');


	saida<<"\n\nThe relative contribution of variables to the between-cluster inertia explained by clusters - CTR(j,i)";
	saida<<"\n\n"<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');
	saida<<"\n";
	saida<<setw(15)<<" ";
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(16)<<"CLUSTER"<<setw(4)<<i+1;
	}
	saida<<"\n"<<setw(20*numeroClasses+16)<<setfill('-')<<' '<<setfill(' ');
	saida<<"\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<"VARIABLE"<<setw(7)<<j+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(20)<<CTRij[i][j];
		}
		saida<<"\n";
	}
	saida<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');

	saida<<"\n\nThe relative contribution of variables and clusters to the between-cluster inertia - CE(j,i)";
	saida<<"\n\n"<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');
	saida<<"\n";
	saida<<setw(15)<<" ";
	for(i=0;i<numeroClasses;i++)
	{
		saida<<setw(16)<<"CLUSTER"<<setw(4)<<i+1;
	}
	saida<<"\n"<<setw(20*numeroClasses+16)<<setfill('-')<<' '<<setfill(' ');
	saida<<"\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		saida<<"VARIABLE"<<setw(7)<<j+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(20)<<CEij[i][j];
		}
		saida<<"\n";
	}
	saida<<setw(20*numeroClasses+16)<<setfill('=')<<' '<<setfill(' ');


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

void FuzzyCMeans::calculaBi()
{
	/////////////////////
	//Cálculo do class-specific

	for(int i=0;i<numeroClasses;i++)
	{
		double distancia=distanciaL2MinkowskyAdaptativo(z,cluster[i].getPrototipo(),cluster[i].getWeightVector());		

		Bi[i]=mi[i]*distancia;
	}
	/////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaBj()
{
	/////////////////////
	//Cálculo do variable-especific
	for(int j=0;j<numeroVariaveis;j++)
	{
		Bj[j]=0.0;

		for(int i=0;i<numeroClasses;i++)
		{
			Bj[j]+=Bij[i][j];
		}
	}
	/////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaBij()
{
	/////////////////////
	//Cálculo da dissimilarity
	for(int i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			double g=cluster[i].getPrototipo(j);
			double weightVector=cluster[i].getWeightVector(j);

			Bij[i][j]=mi[i]*weightVector*(pow(g-z[j],2));
		}
	}

	/////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaMi()
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

void FuzzyCMeans::calculaJi()
{
	//Cálculo do class-specific
	for(int i=0;i<numeroClasses;i++)
	{
		Ji[i]=0.0;
		
		for(int j=0;j<numeroVariaveis;j++)
		{
			Ji[i]+=Jij[i][j];
		}
	}
}

//////////////////////////////////////////////////////////////////////


void FuzzyCMeans::calculaJj()
{
	///////////////////////
	//Cálculo do variable-specific

	for(int j=0;j<numeroVariaveis;j++)
	{
		Jj[j]=0.0;
		for(int i=0;i<numeroClasses;i++)
		{
			for(int k=0;k<totalPadroes;k++)
			{
				double u=pattern[k].getGrauPertinencia(i);
				double x=pattern[k].getX(j);
				double g=cluster[i].getPrototipo(j);
				double weightVector=cluster[i].getWeightVector(j);

				Jj[j]+=weightVector*(pow(u,m))*(pow(x-g,2));
			}
		}
	}
	///////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaJij()
{
	///////////////////////
	//Cálculo do within-cluster

	for(int i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			Jij[i][j]=0.0;

			for(int k=0;k<totalPadroes;k++)
			{
				double u=pattern[k].getGrauPertinencia(i);
				double x=pattern[k].getX(j);
				double g=cluster[i].getPrototipo(j);
				double weightVector=cluster[i].getWeightVector(j);

				Jij[i][j]+=weightVector*(pow(u,m))*(pow(x-g,2));
			}
		}
	}
	///////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaTi()
{
	///////////////////////////
	//Cálculo do class-specific
	for(int i=0;i<numeroClasses;i++)
	{
		Ti[i]=0.0;

		for(int k=0;k<totalPadroes;k++)
		{
			double u=pattern[k].getGrauPertinencia(i);
			
			double distancia=distanciaL2MinkowskyAdaptativo(pattern[k].getX(),z,cluster[i].getWeightVector());

			Ti[i]+=distancia*pow(u,m);
		}
	}
	///////////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaTj()
{
	///////////////////////////
	//Cálculo do variable-specific

	for(int j=0;j<numeroVariaveis;j++)
	{
		Tj[j]=0.0;

		for(int i=0;i<numeroClasses;i++)
		{
			Tj[j]+=Tij[i][j];
		}
	}
	///////////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaTij()
{
	///////////////////////////
	//Calculo do partial

	for(int i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			Tij[i][j]=0.0;

			for(int k=0;k<totalPadroes;k++)
			{
				double u=pattern[k].getGrauPertinencia(i);

				double x=pattern[k].getX(j);

				double weightVector=cluster[i].getWeightVector(j);

				Tij[i][j]+=weightVector*(pow(u,m))*(pow(x-z[j],2));
			}
		}
	}
	///////////////////////////
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaIndicesClasseVariavel()
{
	calculaZ();
	calculaMi();

	calculaJij();
	calculaJi();
	calculaJj();

	calculaTij();
	calculaTi();
	calculaTj();

	calculaBij();
	calculaBi();
	calculaBj();
	
	calculaCORj();
	calculaCTRj();
	calculaCORij();
	calculaCTRij();
	calculaCEij();

	calculaRelativeB();
	calculaRelativeJ();
	calculaRelativeT();

	calculaErroClassificacao();
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaRelativeT()
{
	for(int i=0;i<numeroClasses;i++)
	{
		relativeT[i]=Ti[i]/melhorT;
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaRelativeJ()
{
	for(int i=0;i<numeroClasses;i++)
	{
		relativeJ[i]=Ji[i]/melhorJ;
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::calculaRelativeB()
{
	for(int i=0;i<numeroClasses;i++)
	{
		relativeB[i]=Bi[i]/melhorB;
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeans::imprimeVariaveisMonitoramento()
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

	fout.open("Vetor de Pesos.wri",ios::ate);

	fout<<"\n";
	fout<<setw(15*numeroVariaveis+3+1)<<setfill('=')<<" "<<setfill(' ');	
	fout<<"\nInicialização: "<<inicializacaoAtual+1<<"\n";
	fout<<"\nIteração: "<<iteracao+1<<"\n";
	fout<<"\n";
	fout<<setw(3)<<" ";
	for(j=0;j<numeroVariaveis;j++)
	{
		fout<<setw(12)<<"Variavel "<<setw(3)<<j+1;
	}
	fout<<"\n";
	fout<<setw(15*numeroVariaveis+3+1)<<setfill('-')<<" "<<setfill(' ');
	fout<<"\n";

	double *wV;
	
	for(i=0;i<numeroClasses;i++)
	{
		wV=cluster[i].getWeightVector();
	
		fout<<setw(3)<<i+1;

		for(int j=0;j<numeroVariaveis;j++)
		{
			fout<<setw(15)<<wV[j];
		}

		fout<<"\n";
	}
	
	fout<<setw(15*numeroVariaveis+3+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

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
			double distancia=distanciaL2MinkowskyAdaptativo(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getWeightVector());
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
			double distancia=distanciaL2MinkowskyAdaptativo(z,cluster[i].getPrototipo(),cluster[i].getWeightVector());
			fout<<setw(15)<<distancia;
		}

		fout<<"\n";
	}

	fout<<setw(15*numeroClasses+4+1)<<setfill('-')<<" "<<setfill(' ')<<"\n";

	fout.close();

	///////////////////////////////////////////////////////

}

//////////////////////////////////////////////////////////////////////
