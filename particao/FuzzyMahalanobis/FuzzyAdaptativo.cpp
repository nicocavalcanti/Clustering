// FuzzyAdaptativo.cpp: implementation of the FuzzyAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#include "FuzzyAdaptativo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyAdaptativo::FuzzyAdaptativo()
{

}

FuzzyAdaptativo::~FuzzyAdaptativo() 
{

}

FuzzyAdaptativo::FuzzyAdaptativo(Pattern *p,char *nDados,char *nSaida) : FuzzyMahalanobis(p,nDados,nSaida)
{

}

//////////////////////////////////////////////////////////////////////

void FuzzyAdaptativo::representacao()
{
	/////////////////////////////
	//Estágio 1
	for(int i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoFuzzy(pattern);
	}
	/////////////////////////////

	/////////////////////////////
	//Estágio 2
	sinalInvertivel=true;

	for(i=0;i<numeroClasses;i++)
	{
		double **Q=calculaMatrizVarianciaCovarianciaFuzzy(i);

		double det=determinante(Q,numeroVariaveis);

		if(det<0.0000000001)
		{
			sinalInvertivel=false;

			/////////////////////////////
			//Criação do Arquivo de Saída
			ofstream saida;
			char extensao[15];
			strcpy(extensao,"_ALERTA.wri");
			char *arquivoResultado=new char[75];
			strcpy(arquivoResultado,nomeSaida);
			strncat(arquivoResultado,extensao,11);
			saida.open(arquivoResultado,ios::ate);
			/////////////////////////////

			///////////////////////////
			//Impressão
			saida<<"\n--------------------------\n";
			saida<<"\nInicialização: "<<inicializacaoAtual+1;
			saida<<"\nIteração: "<<iteracao+1;
			saida<<"\nMatriz de determinante nulo:\n";
			for(int l=0;l<numeroVariaveis;l++)
			{
				for(int i=0;i<numeroVariaveis;i++)
				{
					saida<<setw(15)<<Q[l][i];
				}
				saida<<"\n";
			}
			saida<<"\n--------------------------\n";

			saida.close();
			///////////////////////////



			/////////////////////
			//libera memória alocada e que não é mais necessária
			
			//Desaloca Q
			for(int j=0;j<numeroVariaveis;j++)
			{
				delete Q[j];
			}
			delete Q;
			
			//Desaloca os M's antes calculados
			for(j=0;j<i;j++)
			{
				cluster[j].desalocaM();
			}

			break;
			/////////////////////
		}
		else
		{
			double **M=inversa(Q,numeroVariaveis);
			
/*			/////////////////////////////
			//Monitoramento
			ofstream tout("Matriz Q.wri",ios::ate);
			tout<<setprecision(7);
			tout<<"----------   Iteração: "<<iteracao<<"----------\n";

			tout<<"\nMatriz Q da Classe "<<i+1<<":\n";
			for(int k=0;k<numeroVariaveis;k++)
			{
				for(int j=0;j<numeroVariaveis;j++)
				{
					tout<<Q[k][j]<<'\t';
				}
				tout<<"\n";
			}
			tout.close();
			/////////////////////////////*/

			/////////////////////////////////
			//Desaloca Q
			for(int j=0;j<numeroVariaveis;j++)
			{
				delete Q[j];
			}
			delete Q;
			////////////////////////////////

			double multiplicador=pow(det,(1.0/numeroVariaveis));

			for(j=0;j<numeroVariaveis;j++)
			{
				for(int k=0;k<numeroVariaveis;k++)
				{
					M[j][k]*=multiplicador;
				}
			}


/*			/////////////////////////////
			//Monitoramento
			ofstream fout("Matriz M.wri",ios::ate);
			fout<<setprecision(7);
			fout<<"----------   Iteração: "<<iteracao<<"----------\n";

			fout<<"\nMatriz Q da Classe "<<i+1<<":\n";
			for(k=0;k<numeroVariaveis;k++)
			{
				for(int j=0;j<numeroVariaveis;j++)
				{
					fout<<M[k][j]<<'\t';
				}
				fout<<"\n";
			}
			fout.close();
			/////////////////////////////*/

			cluster[i].setM(M);
		}

	}
	/////////////////////////////
	
/*	/////////////////////////////
	//Monitoramento
	ofstream fout("Prototipos.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"Iteração: "<<iteracao<<"\n";

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


double ** FuzzyAdaptativo::calculaMatrizVarianciaCovarianciaFuzzy(int i)
{
	//Aloca memória para Q;
	double **Q=new double*[numeroVariaveis];
	for(int j=0;j<numeroVariaveis;j++)
	{
		Q[j]=new double[numeroVariaveis];
	}

	//Calcula a matriz Q
	for(j=0;j<numeroVariaveis;j++)
	{
		for(int k=0;k<numeroVariaveis;k++)
		{
			if(j<=k)
			{
				Q[j][k]=varianciaCovarianciaFuzzy(j,k,i);
			}
			else
			{
				Q[j][k]=Q[k][j];
			}
		}
	}

	return Q;
}

//////////////////////////////////////////////////////////////////////

double FuzzyAdaptativo::varianciaCovarianciaFuzzy(int var1, int var2, int i)
{
	double numerador=0.0;
//	double denominador=0.0;

	double *prototipo=cluster[i].getPrototipo();
	
	for(int k=0;k<totalPadroes;k++)
	{
		numerador+=pow(pattern[k].getGrauPertinencia(i),m)*(pattern[k].getX(var1)-prototipo[var1])*(pattern[k].getX(var2)-prototipo[var2]);
//		denominador+=pow(pattern[k].getGrauPertinencia(i),m);
	}
	
	double varCovFuzzy=numerador;

	return varCovFuzzy;
}

//////////////////////////////////////////////////////////////////////

void FuzzyAdaptativo::executarAplicacao()
{
	bool primeiroValor=true;

	/////////////////////////////
	//Criação do Arquivo de Saída
	ofstream saida;
	char extensao[5];
	strcpy(extensao,".wri");
	char *arquivoResultado=new char[65];
	strcpy(arquivoResultado,nomeSaida);
	strncat(arquivoResultado,extensao,4);
	saida.open(arquivoResultado);
	/////////////////////////////

	////////////////////////////
	//Impressão
//	saida<<setiosflags(ios::fixed);
	saida<<"Results: \n\n";
	saida<<setw(25)<<"Algorithm: "<<"Fuzzy Mahalanobis Adaptativo";
	saida<<"\n"<<setw(25)<<"Input Data: "<<nomeDados;
	saida<<"\n"<<setw(25)<<"Parameter m: "<<m;
	saida<<"\n"<<setw(25)<<"Clusters: "<<numeroClasses;
	saida<<"\n"<<setw(25)<<"Initialization Number: "<<numeroInicializacoes;
	saida<<"\n"<<setw(25)<<"Iteration Limit: "<<limiteIteracao;

	////////////////////////////

	for(inicializacaoAtual=0;inicializacaoAtual<numeroInicializacoes;inicializacaoAtual++)
	{
		////////////////////////////
		//Impressão
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

			if(!sinalInvertivel)
			{
				//Caso seja não invertível
				cout<<"Inicializacao "<<inicializacaoAtual;
				cout<<"\nIteracao "<<iteracao;
				cout<<" - Encontrada matriz nao inversivel\n\n";
				cout.flush();
				break;
			}

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
			//Impressão
			saida<<"\n"<<setw(9)<<iteracao+1;
			saida<<setw(20)<<J<<setw(20)<<T<<setw(20)<<R<<setw(20)<<B<<setw(20)<<CR<<setw(20)<<DC;
			//////////////////////////

			if(iteracao)
			{
				if(fabs(Jatual-Janterior)<=epsilon)
				{
					for(int l=0;l<numeroClasses;l++)
					{
						cluster[l].desalocaM();
					}

					break;
				}
			}

			alocacao();

			for(int l=0;l<numeroClasses;l++)
			{
				cluster[l].desalocaM();
			}

			Janterior=Jatual;

			cout<<"Inicializacao "<<inicializacaoAtual<<"\nIteracao "<<iteracao<<": OK\n";
			cout.flush();
		}
	}

	if(primeiroValor)
	{
		return;
	}

	int *classes=classesFinais();
	setMelhorParticao(classes);

	/////////////////////////
	//Impressão
	saida<<"\n\n------------------- Membership matrix for the minimal J obtained -------------------\n";
	saida<<"Pattern";
	for(int i=0;i<numeroClasses;i++)
	{
		saida<<setw(14)<<"u"<<i+1;
	}
	saida<<setw(15)<<"Hard Cluster"<<setw(20)<<"A Priori Cluster"<<"\n";

	for(int k=0;k<totalPadroes;k++)
	{
		saida<<setw(7)<<k+1;
		for(int i=0;i<numeroClasses;i++)
		{
			saida<<setw(15)<<pattern[k].getGrauPertinencia(i);
		}
		saida<<setw(15)<<classes[k]+1<<setw(20)<<pattern[k].getClassePriori()+1<<"\n";

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

	delete classes;


}

//////////////////////////////////////////////////////////////////////

void FuzzyAdaptativo::executarSimulacao()
{
	bool primeiroValor=true;

	int totalIteracoes=0;

	for(inicializacaoAtual=0;inicializacaoAtual<numeroInicializacoes;inicializacaoAtual++)
	{
		inicializacao();

		double Janterior;
		double Jatual;

		for(iteracao=0;iteracao<limiteIteracao;iteracao++)	
		{
			representacao();

			if(!sinalInvertivel)
			{
				//Caso seja não invertível
				cout<<"Inicializacao "<<inicializacaoAtual+1;
				cout<<"\nIteracao "<<iteracao+1;
				cout<<" - Encontrada matriz nao inversivel\n\n";
				cout.flush();
				break;
			}

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
					for(int l=0;l<numeroClasses;l++)
					{
						cluster[l].desalocaM();
					}

					break;
				}
			}

			alocacao();

			for(int l=0;l<numeroClasses;l++)
			{
				cluster[l].desalocaM();
			}

			Janterior=Jatual;
		}
	
		totalIteracoes+=iteracao;
	}

	if(primeiroValor)
	{
		return;
	}

	int *classes=classesFinais();

	setMelhorParticao(classes);

	calculaIndices();

	setMelhoresIndices();

	calculaErroClassificacao();

	numeroIteracaoMedio=(double)totalIteracoes/numeroInicializacoes;

	delete classes;
	
}

//////////////////////////////////////////////////////////////////////