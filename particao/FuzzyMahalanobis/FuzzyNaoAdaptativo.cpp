// FuzzyNaoAdaptativo.cpp: implementation of the FuzzyNaoAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#include "FuzzyNaoAdaptativo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyNaoAdaptativo::FuzzyNaoAdaptativo()
{

}

FuzzyNaoAdaptativo::~FuzzyNaoAdaptativo()
{

}

FuzzyNaoAdaptativo::FuzzyNaoAdaptativo(Pattern *p,char *nDados,char *nSaida) : FuzzyMahalanobis(p,nDados,nSaida)
{

}

//////////////////////////////////////////////////////////////////////

double FuzzyNaoAdaptativo::varianciaCovariancia(int var1,int var2)
{
	double varianceCovariance=0.0;

	for(int k=0;k<totalPadroes;k++)
	{
		varianceCovariance+=(pattern[k].getX(var1)-prototipoGeral[var1])*(pattern[k].getX(var2)-prototipoGeral[var2]);
	}

	//Precisa verificar se divide pelo total de padr�es ou pelo numero de padr�es menos 1
	varianceCovariance/=totalPadroes-1;

	return varianceCovariance;
}

//////////////////////////////////////////////////////////////////////

void FuzzyNaoAdaptativo::calculaMatrizVarianciaCovariancia()
{
	//Calcula o prototipo geral
	prototipoGeral=new double[numeroVariaveis];
	
	for(int j=0;j<numeroVariaveis;j++)
	{
		prototipoGeral[j]=0.0;

		for(int k=0;k<totalPadroes;k++)
		{
			prototipoGeral[j]+=pattern[k].getX(j);
		}
		
		prototipoGeral[j]/=totalPadroes;
	}


	//Calcula a matriz Q
	for(j=0;j<numeroVariaveis;j++)
	{
		for(int k=0;k<numeroVariaveis;k++)
		{
			if(j<=k)
			{
				Q[j][k]=varianciaCovariancia(j,k);
			}
			else
			{
				Q[j][k]=Q[k][j];
			}
		}
	}

/*	//////////////////////////////////////
	//Monitoramento
	ofstream fout("Q.wri",ios::ate);
//	fout<<setiosflags(ios::fixed);
	fout<<"Itera��o: "<<iteracao<<"\n";

	fout<<"Protipo geral:\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		fout<<prototipoGeral[j]<<"\t";
	}

	fout<<"\n\nMatriz Q:\n";
	for(j=0;j<numeroVariaveis;j++)
	{
		for(int  k=0;k<numeroVariaveis;k++)
		{
			fout<<setw(15)<<Q[j][k];
		}
		fout<<"\n";
	}

	fout.close();

	//////////////////////////////////////*/

	//Desaloca o prot�tipo geral
	delete prototipoGeral;
}

/////////////////////////////////////////////////////////////////////

void FuzzyNaoAdaptativo::representacao()
{
	for(int i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoFuzzy(pattern);
	}

/*	////////////////////////////
	//Monitoramento
	ofstream fout("Prototipos.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"Itera��o: "<<iteracao<<"\n";

	double *prototipo;
	for(i=0;i<numeroClasses;i++)
	{
		prototipo=cluster[i].getPrototipo();

		fout<<"\nPrototipo da Classe "<<i+1<<": ";
		for(int j=0;j<numeroVariaveis;j++)
		{
			fout<<setw(6)<<prototipo[j];
		}

		fout<<"\n";
	}

	fout.close();
	/////////////////////////////*/
}

//////////////////////////////////////////////////////////////////////

void FuzzyMahalanobis::calculaJ()
{
	J=0.0;

	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			double grauPertinencia=pattern[k].getGrauPertinencia(i);
			double distancia=distanciaMahalanobis(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getM());

			J+=pow(grauPertinencia,m)*distancia;
		}
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyNaoAdaptativo::calculaM()
{
	//Aloca��o de mem�ria para Q
	Q=new double*[numeroVariaveis];
	for(int i=0;i<numeroVariaveis;i++)
	{
		Q[i]=new double[numeroVariaveis];
	}
	///////////////////////////////////////

	calculaMatrizVarianciaCovariancia();

	double det=determinante(Q,numeroVariaveis);

	double **inv=inversa(Q,numeroVariaveis);

	//Desaloca��o de Q
	for(i=0;i<numeroVariaveis;i++)
	{
		delete Q[i];
	}
	delete Q;
	//////////////////////////////////////

	M=new double*[numeroVariaveis];
	for(i=0;i<numeroVariaveis;i++)
	{
		M[i]=new double[numeroVariaveis];
	}

	double multiplicador=pow(det,(1.0/numeroVariaveis));

	for(i=0;i<numeroVariaveis;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			M[i][j]=multiplicador*inv[i][j];
		}
	}

/*	//////////////////////////////////////
	//Monitoramento
	ofstream fout("M.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"\n";

	fout<<"\n\nMatriz M:\n";
	for(int j=0;j<numeroVariaveis;j++)
	{
		for(int  k=0;k<numeroVariaveis;k++)
		{
			fout<<setw(15)<<M[j][k];
		}
		fout<<"\n";
	}

	fout.close();
	//////////////////////////////////////*/


	//Desaloca��o de inv
	for(i=0;i<numeroVariaveis;i++)
	{
		delete inv[i];
	}
	delete inv;
	/////////////////////////////////
	
	for(i=0;i<numeroClasses;i++)
	{
		cluster[i].setM(M);
	}
}

//////////////////////////////////////////////////////////////////////

void FuzzyNaoAdaptativo::executarAplicacao()
{
	calculaM();

	bool primeiroValor=true;

	/////////////////////////////
	//Cria��o do Arquivo de Sa�da
	ofstream saida;
	char extensao[5];
	strcpy(extensao,".wri");
	char *arquivoResultado=new char[65];
	strcpy(arquivoResultado,nomeSaida);
	strncat(arquivoResultado,extensao,4);
	saida.open(arquivoResultado);
	/////////////////////////////

	////////////////////////////
	//Impress�o
//	saida<<setiosflags(ios::fixed);
	saida<<"Results: \n\n";
	saida<<setw(25)<<"Algorithm: "<<"Fuzzy Mahalanobis N�o Adaptativo";
	saida<<"\n"<<setw(25)<<"Input Data: "<<nomeDados;
	saida<<"\n"<<setw(25)<<"Parameter m: "<<m;
	saida<<"\n"<<setw(25)<<"Clusters: "<<numeroClasses;
	saida<<"\n"<<setw(25)<<"Initialization Number: "<<numeroInicializacoes;
	saida<<"\n"<<setw(25)<<"Iteration Limit: "<<limiteIteracao;

	////////////////////////////
	
	for(inicializacaoAtual=0;inicializacaoAtual<numeroInicializacoes;inicializacaoAtual++)
	{
		////////////////////////////
		//Impress�o
		saida<<"\n\n"<<setw(1+9+20*6)<<setfill('=')<<' '<<setfill(' ');
		saida<<"\nInitialization: "<<inicializacaoAtual+1;
		saida<<"\n"<<setw(9)<<"Iteration"<<setw(20)<<"J"<<setw(20)<<"T"<<setw(20)<<"R"<<setw(20)<<"B"<<setw(20)<<"CR"<<setw(20)<<"DC";
		saida<<"\n"<<setw(1+9+20*6)<<setfill('-')<<' '<<setfill(' ');
		////////////////////////////

		inicializacao();

		double Janterior;
		double Jatual;

		for(iteracao=0;iteracao<limiteIteracao;iteracao++)		
		{
			representacao();
	
			setClassificacaoHard();
			calculaIndices();

			Jatual=J;

			if(primeiroValor)
			{
				setMelhoresIndices();

				atualizaMelhorParticao();

				primeiroValor=false;
			}

			if(melhorJ>Jatual)
			{
				setMelhoresIndices();

				atualizaMelhorParticao();
			}

			//////////////////////////
			//Impress�o
			saida<<"\n"<<setw(9)<<iteracao+1;
			saida<<setw(20)<<J<<setw(20)<<T<<setw(20)<<R<<setw(20)<<B<<setw(20)<<CR<<setw(20)<<DC;
			//////////////////////////

			if(iteracao)
			{
				if(fabs(Jatual-Janterior)<=epsilon)
				{
					break;
				}
			}

			alocacao();
		
			Janterior=Jatual;
		}
	}

	int *classes=classesFinais();
	setMelhorParticao(classes);

	/////////////////////////
	//Impress�o
	saida<<"\n\n------------------- Membership matrix for the minimal J obtained -------------------\n";
	saida<<"Pattern";
	for(int i=0;i<numeroClasses;i++)
	{
		saida<<setw(14)<<"u"<<i+1;
	}
	saida<<setw(15)<<"Hard Cluster"<<"\n";

	for(int k=0;k<totalPadroes;k++)
	{
		saida<<setw(7)<<k+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(15)<<pattern[k].getGrauPertinencia(i);
		}
		saida<<setw(15)<<classes[k]+1<<"\n";
	}

	saida<<"\nHard partition obtained for the minimal J obtained\n";
	
	for(i=0;i<numeroClasses;i++)
	{
		saida<<"Cluster "<<i+1<<" -> ";
		for(int k=0;k<totalPadroes;k++)
		{
			if(classes[k]==i)
			{
				saida<<k+1<<",";
			}
		}
		saida<<"\n";
	}

	saida<<"\n";
	
	saida<<"About J\n";
	saida<<"\t\tBest = "<<melhorJ;
	saida<<"\n\n";
	saida<<"\t\t---- The values of the others indexs for the minimal J ----\n";
	saida<<"\t\t"<<setw(30)<<"CR = "<<melhorCR;
	saida<<"\n";
	saida<<"\t\t"<<setw(30)<<"DC = "<<melhorDC;
	saida<<"\n";
	saida<<"\t\t"<<setw(30)<<"T = "<<melhorT;
	saida<<"\n";
	saida<<"\t\t"<<setw(30)<<"R = "<<melhorR;
	saida<<"\n";
	saida<<"\t\t"<<setw(30)<<"B = "<<melhorB;

	saida<<"\n";

	saida.close();
	delete arquivoResultado;

	/////////////////////////

	calculaIndicesClasseVariavel();
	imprimeInterpretadores();

	for(int j=0;j<numeroVariaveis;j++)
	{
		delete M[j];
	}
	delete M;

	delete classes;

}

//////////////////////////////////////////////////////////////////////

void FuzzyNaoAdaptativo::executarSimulacao()
{
	calculaM();

	bool primeiroValor=true;

	int totalIteracoes=0;

	for(int i=0;i<numeroInicializacoes;i++)
	{
		inicializacao();

		double Janterior;
		double Jatual;

		for(iteracao=0;iteracao<limiteIteracao;iteracao++)		
		{
			representacao();

			calculaJ();
			Jatual=J;

			if(primeiroValor)
			{
				melhorJ=J;
				atualizaMelhorParticao();

				primeiroValor=false;
			}

			if(melhorJ>Jatual)
			{
				melhorJ=J;
				atualizaMelhorParticao();
			}

			if(iteracao)
			{
				if(fabs(Jatual-Janterior)<=epsilon)
				{
					break;
				}
			}

			alocacao();

			Janterior=Jatual;
		}

		totalIteracoes+=iteracao;
	}

	int *classes=classesFinais();
	setMelhorParticao(classes);

	calculaIndices();

	setMelhoresIndices();

	calculaErroClassificacao();

	numeroIteracaoMedio=(double)totalIteracoes/numeroInicializacoes;

	delete classes;

	for(int j=0;j<numeroVariaveis;j++)
	{
		delete M[j];
	}
	delete M;
}

//////////////////////////////////////////////////////////////////////