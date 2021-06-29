// FuzzyCMeansAdaptativoUnico.cpp: implementation of the FuzzyCMeansAdaptativoUnico class.
//
//////////////////////////////////////////////////////////////////////

#include "FuzzyCMeansAdaptativoUnico.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

FuzzyCMeansAdaptativoUnico::FuzzyCMeansAdaptativoUnico()
{

}

FuzzyCMeansAdaptativoUnico::~FuzzyCMeansAdaptativoUnico()
{

}

FuzzyCMeansAdaptativoUnico::FuzzyCMeansAdaptativoUnico(Pattern *p,char *nDados,char *nSaida) : FuzzyCMeans(p,nDados,nSaida)
{

}

//////////////////////////////////////////////////////////////////////

void FuzzyCMeansAdaptativoUnico::representacao()
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
	
	weightVector=new double[numeroVariaveis];

	for(int j=0;j<numeroVariaveis;j++)
	{
		//Calculo do numerador
		double produtorio=1.0;
		
		for(int h=0;h<numeroVariaveis;h++)
		{
			double somatorioNumerador=0.0;
			
			for(int k=0;k<totalPadroes;k++)
			{
				for(int i=0;i<numeroClasses;i++)
				{
					double X=pattern[k].getX(h);
					double grauPertinencia=pattern[k].getGrauPertinencia(i);
					double g=cluster[i].getPrototipo(h);

					somatorioNumerador+=pow(grauPertinencia,m)*pow((X-g),2);
				}
			}

			produtorio*=somatorioNumerador;
		}

		double numerador=pow(produtorio,(1.0/numeroVariaveis));

		//Calculo do denominador

		double denominador=0.0;

		for(int k=0;k<totalPadroes;k++)
		{
			for(int i=0;i<numeroClasses;i++)
			{
				double X=pattern[k].getX(j);
				double grauPertinencia=pattern[k].getGrauPertinencia(i);
				double g=cluster[i].getPrototipo(j);

				denominador+=pow(grauPertinencia,m)*pow((X-g),2);
			}
		}

		weightVector[j]=numerador/denominador;
	}

	for(i=0;i<numeroClasses;i++)
	{
		cluster[i].setWeightVector(weightVector);
	}

/*	/////////////////////////////
	//Monitoramento
	ofstream tout("Weight Vectors.wri",ios::ate);
	tout<<setprecision(4);
	tout<<"\n\n--------- Iteração "<<iteracao<<" ----------\n";

	double total=1.0;

	for(j=0;j<numeroVariaveis;j++)
	{
		tout<<setw(12)<<weightVector[j];
		total*=weightVector[j];
	}
	
	tout<<"  Total: "<<setw(3)<<total;

	tout.close();
	/////////////////////////////*/


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

void FuzzyCMeansAdaptativoUnico::executarAplicacao()
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
	saida<<setw(25)<<"Algorithm: "<<"Fuzzy C-Means Adaptativo Único";
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

//			imprimeVariaveisMonitoramento();

			if(iteracao)
			{
				if(fabs(Jatual-Janterior)<=epsilon)
				{
					////////////////////////
					//Apaga weightVector
					delete weightVector;
					////////////////////////

					break;
				}
			}

 			alocacao();

			////////////////////////
			//Apaga weightVector
			delete weightVector;
			////////////////////////

			Janterior=Jatual;

			cout<<"Inicializacao "<<inicializacaoAtual<<"\nIteracao "<<iteracao<<": OK\n";
			cout.flush();

		}

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

void FuzzyCMeansAdaptativoUnico::executarSimulacao()
{
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
					////////////////////////
					//Apaga weightVector
					delete weightVector;
					////////////////////////

					break;
				}
			}

			alocacao();

			////////////////////////
			//Apaga weightVector
			delete weightVector;
			////////////////////////

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
	
}

//////////////////////////////////////////////////////////////////////
