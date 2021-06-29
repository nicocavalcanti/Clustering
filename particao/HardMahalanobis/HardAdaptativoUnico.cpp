// HardAdaptativoUnico.cpp: implementation of the HardAdaptativoUnico class.
//
//////////////////////////////////////////////////////////////////////

#include "HardAdaptativoUnico.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

HardAdaptativoUnico::HardAdaptativoUnico()
{

}

HardAdaptativoUnico::~HardAdaptativoUnico()
{

}

HardAdaptativoUnico::HardAdaptativoUnico(Pattern *p,char *nSaida) : HardMahalanobis(p,nSaida)
{

}

//////////////////////////////////////////////////////////////////////

int * HardAdaptativoUnico::executar()
{
	bool primeiroValor=true;

	for(int i=0;i<numeroInicializacoes;i++)
	{
		inicializacao();
		
		double Wanterior;
		double Watual;

		for(iteracao=0;iteracao<limiteIteracao;iteracao++)	
		{
			representacao();
			
			if(!sinalInvertivel)
			{
				//Caso seja não invertível
				break;
			}

			Watual=calculaW();

 			if(iteracao)
			{
				if(fabs(Watual-Wanterior)<=epsilon)
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

			for(int l=0;l<numeroVariaveis;l++)
			{
				delete M[l];
			}
			delete M;

			Wanterior=Watual;
		}
	}

	int *classes=new int[totalPadroes];
	for(int k=0;k<totalPadroes;k++)
	{
		classes[k]=melhorParticao[k];
	}

	setMelhorParticao(classes);

	return classes;
}

//////////////////////////////////////////////////////////////////////

double ** HardAdaptativoUnico::calculaMatrizVarianciaCovarianciaHard(int i)
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
				Q[j][k]=varianciaCovarianciaHardAlterada(j,k,i);
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

double HardAdaptativoUnico::varianciaCovarianciaHardAlterada(int var1, int var2, int i)
{
	double varCov=0.0;

	double *prototipo=cluster[i].getPrototipo();
	
	for(int k=0;k<totalPadroes;k++)
	{
		if(i==pattern[k].getClasse())
		{
			varCov+=(pattern[k].getX(var1)-prototipo[var1])*(pattern[k].getX(var2)-prototipo[var2]);
		}
	}
	
	return varCov;
}

//////////////////////////////////////////////////////////////////////

void HardAdaptativoUnico::representacao()
{
	/////////////////////////////
	//Estágio 1
	for(int i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoHard(pattern);
	}
	/////////////////////////////

	/////////////////////////////
	//Estágio 2
	sinalInvertivel=true;

	double ***Q=new double**[numeroClasses];
	for(i=0;i<numeroClasses;i++)
	{
		Q[i]=calculaMatrizVarianciaCovarianciaHard(i);
	}
	
	double **Qpooled=new double*[numeroVariaveis];
	for(int j=0;j<numeroVariaveis;j++)
	{
		Qpooled[j]=new double[numeroVariaveis];
	}

	for(j=0;j<numeroVariaveis;j++)
	{
		for(int k=0;k<numeroVariaveis;k++)
		{
			Qpooled[j][k]=0.0;

			for(int i=0;i<numeroClasses;i++)
			{
				Qpooled[j][k]+=Q[i][j][k]/(totalPadroes-numeroClasses);
			}
		}
	}

/*	/////////////////////////////
	//Monitoramento
	ofstream tout("Matriz Q.wri",ios::ate);
	tout<<setprecision(7);
	tout<<"----------   Iteração: "<<iteracao<<"----------\n";

	for(i=0;i<numeroClasses;i++)
	{
		tout<<"\nMatriz Q da Classe "<<i+1<<":\n";
		for(int j=0;j<numeroVariaveis;j++)
		{
			for(int k=0;k<numeroVariaveis;k++)
			{
				tout<<Q[i][j][k]<<'\t';
			}
			tout<<"\n";
		}

		tout<<"\n";
	}

	tout<<"\nMatriz Qpooled\n";
	for(i=0;i<numeroVariaveis;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			tout<<Qpooled[i][j]<<'\t';
		}
		tout<<"\n";
	}
	tout.close();
	/////////////////////////////*/


	/////////////////
	//Apaga	Q
	for(i=0;i<numeroClasses;i++)
	{
		for(int j=0;j<numeroVariaveis;j++)
		{
			delete Q[i][j];
		}
		delete Q[i];
	}

	delete Q;
	/////////////////

	double det=determinante(Qpooled,numeroVariaveis);

	if(det<0.0000000001)
	{
		sinalInvertivel=false;

		/////////////////////
		//libera memória alocada e que não é mais necessária

		//Apaga Qpooled
		for(j=0;j<numeroVariaveis;j++)
		{
			delete Qpooled[j];
		}
		delete Qpooled;
		/////////////////////
	}
	else
	{
		M=inversa(Qpooled,numeroVariaveis);

		/////////////////
		//Apaga Qpooled
		for(j=0;j<numeroVariaveis;j++)
		{
			delete Qpooled[j];
		}
		delete Qpooled;
		/////////////////

		double multiplicador=pow(det,(1.0/numeroVariaveis));

		for(j=0;j<numeroVariaveis;j++)
		{
			for(int k=0;k<numeroVariaveis;k++)
			{
				M[j][k]*=multiplicador;
			}
		}

		for(i=0;i<numeroClasses;i++)
		{
			cluster[i].setM(M);
		}
	}


/*	/////////////////////////////
	//Monitoramento
	ofstream fout("Prototipos.wri",ios::ate);
	fout<<setprecision(4);
	fout<<"-----------  Iteração: "<<iteracao<<" -------------\n";

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
