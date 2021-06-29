// HardNaoAdaptativo.cpp: implementation of the HardNaoAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#include "HardNaoAdaptativo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

HardNaoAdaptativo::HardNaoAdaptativo()
{

}

HardNaoAdaptativo::~HardNaoAdaptativo()
{

}

HardNaoAdaptativo::HardNaoAdaptativo(Pattern *p,char *nSaida) : HardMahalanobis(p,nSaida)
{

}

//////////////////////////////////////////////////////////////////////

void HardNaoAdaptativo::calculaM()
{
	///////////////////////////////////////
	//Alocação de memória para Q
	Q=new double*[numeroVariaveis];
	for(int i=0;i<numeroVariaveis;i++)
	{
		Q[i]=new double[numeroVariaveis];
	}
	///////////////////////////////////////

	calculaMatrizVarianciaCovariancia();

	double det=determinante(Q,numeroVariaveis);

	double **inv=inversa(Q,numeroVariaveis);

	//////////////////////////////////////
	//Desalocação de Q
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

	//////////////////////////////////////
	//Monitoramento
	ofstream fout("M.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"\n";

	fout<<"\n\nMatriz M:\n";
	for(int j=0;j<numeroVariaveis;j++)
	{
		for(int  k=0;k<numeroVariaveis;k++)
		{
			fout<<M[j][k]<<"\t";
		}
		fout<<"\n";
	}

	fout.close();
	//////////////////////////////////////*/


	//Desalocação de inv
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

void HardNaoAdaptativo::calculaMatrizVarianciaCovariancia()
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

	//////////////////////////////////////
	//Monitoramento
	ofstream fout("Q.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"Iteração: "<<iteracao<<"\n";

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
			fout<<Q[j][k]<<"\t";
		}
		fout<<"\n";
	}

	fout.close();

	//////////////////////////////////////*/

	//Desaloca o protótipo geral
	delete prototipoGeral;
}

//////////////////////////////////////////////////////////////////////

double HardNaoAdaptativo::varianciaCovariancia(int var1, int var2)
{
	double varianceCovariance=0.0;

	for(int k=0;k<totalPadroes;k++)
	{
		varianceCovariance+=(pattern[k].getX(var1)-prototipoGeral[var1])*(pattern[k].getX(var2)-prototipoGeral[var2]);
	}

	//Precisa verificar se divide pelo total de padrões ou pelo numero de padrões menos 1
	varianceCovariance/=totalPadroes-1;

	return varianceCovariance;

}

//////////////////////////////////////////////////////////////////////

void HardNaoAdaptativo::representacao()
{
	for(int i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoHard(pattern);
	}

/*	/////////////////////////////
	//Monitoramento
	ofstream fout("PrototiposHard.wri",ios::ate);
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

int * HardNaoAdaptativo::executar()
{
	calculaM();

	bool primeiroValor=true;

	for(int i=0;i<numeroInicializacoes;i++)
	{
		inicializacao();

		bool maisDeUmValor=false;

		double Wanterior;
		double Watual;

		for(iteracao=0;iteracao<limiteIteracao;iteracao++)	
		{
			representacao();

			Watual=calculaW();

			if(iteracao)
			{
				if(Watual==Wanterior)
				{
					break;
				}
			}

			if(primeiroValor)
			{
				Wminimo=Watual;
				atualizaMelhorParticao();

				primeiroValor=false;
			}

			if(Wminimo>Watual)
			{
				Wminimo=Watual;
				atualizaMelhorParticao();
			}

			alocacao();
		
			Wanterior=Watual;
		}
	}

	int *classes=new int[totalPadroes];
	
	for(int k=0;k<totalPadroes;k++)
	{
		classes[k]=melhorParticao[k];
	}

	for(i=0;i<numeroVariaveis;i++)
	{
		delete M[i];
	}
	delete M;

	setMelhorParticao(classes);

	return classes;
}

//////////////////////////////////////////////////////////////////////
