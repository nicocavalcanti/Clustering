// HardAdaptativo.cpp: implementation of the HardAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#include "HardAdaptativo.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

HardAdaptativo::HardAdaptativo()
{

}

HardAdaptativo::~HardAdaptativo()
{

}

HardAdaptativo::HardAdaptativo(Pattern *p,char *nSaida) : HardMahalanobis(p,nSaida)
{

}

//////////////////////////////////////////////////////////////////////

int * HardAdaptativo::executar()
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

			for(int l=0;l<numeroClasses;l++)
			{
				cluster[l].desalocaM();
			}

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

double ** HardAdaptativo::calculaMatrizVarianciaCovarianciaHard(int i)
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
				Q[j][k]=varianciaCovarianciaHard(j,k,i);
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

double HardAdaptativo::varianciaCovarianciaHard(int var1, int var2, int i)
{
	double numerador=0.0;
	int denominador=0;

	double *prototipo=cluster[i].getPrototipo();
	
	for(int k=0;k<totalPadroes;k++)
	{
		if(i==pattern[k].getClasse())
		{
			numerador+=(pattern[k].getX(var1)-prototipo[var1])*(pattern[k].getX(var2)-prototipo[var2]);
			denominador++;
		}		
	}
	
	double varCov=numerador/(denominador-1);

	return varCov;
}

//////////////////////////////////////////////////////////////////////

void HardAdaptativo::representacao()
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

	for(i=0;i<numeroClasses;i++)
	{
		double **Q=calculaMatrizVarianciaCovarianciaHard(i);

		double det=determinante(Q,numeroVariaveis);

		if(det<0.0000000001)
		{
			sinalInvertivel=false;

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
			
	/*		/////////////////////////////
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
