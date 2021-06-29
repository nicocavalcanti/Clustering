// SCAD.cpp: implementation of the SCAD class.
//
//////////////////////////////////////////////////////////////////////

#include "SCAD.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SCAD::SCAD()
{

}

SCAD::~SCAD()
{

}

SCAD::SCAD(Pattern *p,char *nSaida) : Particao(p,nSaida)
{
	melhorParticao=new double*[totalPadroes];
	for(int k=0;k<totalPadroes;k++)
	{
		melhorParticao[k]=new double[numeroClasses];
	}

	delta=new double[numeroClasses];

	K=2;
}

//////////////////////////////////////////////////////////////////////

void SCAD::inicializacao()
{
	//Inicialização dos graus de pertinência
	for(int k=0;k<totalPadroes;k++)
	{
		double limite=1.0;

		for(int i=0;i<numeroClasses-1;i++)
		{
			double semente=(rand()%100)/100.0;

			double valor=semente*limite;

			pattern[k].setGrauPertinencia(i,valor);

			limite-=valor;
		}

		pattern[k].setGrauPertinencia(i,limite);

	}

	//Inicializa Feature Weight
	for(int i=0;i<numeroClasses;i++)
	{
		double *featureWeight=new double[numeroVariaveis];
		
		for(int j=0;j<numeroVariaveis;j++)
		{
			featureWeight[j]=1.0/numeroVariaveis;
		}

		cluster[i].setFeatureWeight(featureWeight);
	}

	//Inicializa Prototipos
	for(i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoFuzzy(pattern);
	}
	

	////////////////////////////////////////////
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

void SCAD::atualizaDelta()
{
	for(int i=0;i<numeroClasses;i++)
	{
		/////////////////////
		//Cálculo do numerador

		double numerador=0.0;

		for(int k=0;k<totalPadroes;k++)
		{
			double soma=0.0;

			for(int j=0;j<numeroVariaveis;j++)
			{
				double v=cluster[i].getFeatureWeight(j);
				double x=pattern[k].getX(j);
				double g=cluster[i].getPrototipo(j);

				soma+=v*pow(x-g,2.0);
			}

			double u=pattern[k].getGrauPertinencia(i);

			numerador+=soma*pow(u,m);
		}


		////////////////////
		//Cálculo do denominador

		double denominador=0.0;

		for(int j=0;j<numeroVariaveis;j++)
		{
			double v=cluster[i].getFeatureWeight(j);

			denominador+=pow(v,2.0);
		}

		delta[i]=K*numerador/denominador;

	}
}

//////////////////////////////////////////////////////////////////////

void SCAD::atualizaPrototipos()
{
	for(int i=0;i<numeroClasses;i++)
	{
		cluster[i].calculaPrototipoFuzzy(pattern);
		
		for(int j=0;j<numeroVariaveis;j++)
		{
			double featureWeight=cluster[i].getFeatureWeight(j);

			if(featureWeight==0.0)
			{
				double *g=cluster[i].getPrototipo();

				g[j]=0.0;
			}
		}
	}

	/////////////////////////////
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

int * SCAD::executar()
{
	bool primeiroValor=true;

	for(int i=0;i<numeroInicializacoes;i++)
	{
		inicializacao();

		atualizaDelta();

		for(int i=0;i<numeroClasses;i++)
		{
			cluster[i].desalocaFeatureWeight();
		}

		bool maisDeUmValor=false;

		double Wanterior;
		double Watual;

		//Inicializa o número de iterações
		iteracao=0;

		do
		{
			atualizaFeatureWeight();

			atualizaGrauPertinencia();

			atualizaPrototipos();

			Watual=calculaW();

/*			////////////////////////
			//Monitoramento
			ofstream tout("W.wri",ios::ate);
			tout<<setprecision(15);
			tout<<"\n"<<Watual;

			tout.close();
			/////////////////////////*/

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

			atualizaDelta();

			for(int l=0;l<numeroClasses;l++)
			{
				cluster[l].desalocaFeatureWeight();
			}

			if(maisDeUmValor)
			{
				if(fabs(Watual-Wanterior)<=epsilon)
				{
					break;
				}
			}

			Wanterior=Watual;

			iteracao++;

			maisDeUmValor=true;

		}
		while(iteracao<limiteIteracao);
	}

/*	////////////////////////
	//Monitoramento
	ofstream fout("melhorW.wri",ios::ate);
	fout<<"\n"<<Wminimo;

	fout.close();
	/////////////////////////*/

//	calculaDC();

	int *classes=classesFinais();

	return classes;

}

//////////////////////////////////////////////////////////////////////

void SCAD::atualizaFeatureWeight()
{
	double *temp=new double[numeroVariaveis];
	for(int j=0;j<numeroVariaveis;j++)
	{
		temp[j]=1.0;
	}

	for(int i=0;i<numeroClasses;i++)
	{
		double *featureWeight=new double[numeroVariaveis];

		for(int j=0;j<numeroVariaveis;j++)
		{
			double soma=0.0;

			for(int k=0;k<totalPadroes;k++)
			{
				double *x=pattern[k].getX();
				double *g=cluster[i].getPrototipo();

				double distancia=distanciaEuclidianaAdaptativa(x,g,temp);

				double u=pattern[k].getGrauPertinencia(i);

				soma+=pow(u,m)*((distancia/numeroVariaveis)-pow(x[j]-g[j],2));
			}

			double valor=(1.0/numeroVariaveis)+(soma/(2.0*delta[i]));
			
			
			if(valor<0.0)
			{
				featureWeight[j]=0.0;
			}
			else
			{
				if(valor>1.0)
				{
					featureWeight[j]=1.0;
				}
				else
				{
					featureWeight[j]=valor;
				}
			}
		}

		cluster[i].setFeatureWeight(featureWeight);
	}

	delete temp;

	//////////////////////////////////
	//Monitoramento
	ofstream fout("Feature Weight.wri",ios::ate);
	
	fout<<"\n ------------- Iteração "<<iteracao<<" ------------------\n";

	for(i=0;i<numeroClasses;i++)
	{
		double *v=cluster[i].getFeatureWeight();
		
		double soma=0.0;

		for(int j=0;j<numeroVariaveis;j++)
		{
			fout<<v[j]<<"\t";
			soma+=v[j];
		}

		fout<<soma<<"\n";
	}

	///////////////////////////////////*/
}

//////////////////////////////////////////////////////////////////////

double SCAD::distanciaEuclidianaAdaptativa(double *x, double *g,double *v)
{
	double distancia=0.0;
	
	for(int j=0;j<numeroVariaveis;j++)
	{
		distancia+=v[j]*pow((x[j]-g[j]),2);
	}

	return distancia;
}

//////////////////////////////////////////////////////////////////////

void SCAD::atualizaGrauPertinencia()
{
	for(int k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			double numerador=distanciaEuclidianaAdaptativa(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getFeatureWeight());
				
			double acumulado=0.0;

			for(int h=0;h<numeroClasses;h++)
			{
				double denominador=distanciaEuclidianaAdaptativa(pattern[k].getX(),cluster[h].getPrototipo(),cluster[h].getFeatureWeight());

				acumulado+=pow((numerador/denominador),(1.0/(m-1)));
			}

			pattern[k].setGrauPertinencia(i,pow(acumulado,-1));
		}
	}

	////////////////////////////////////////////
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

double SCAD::calculaW()
{
	double W=0.0;

	for(int i=0;i<numeroClasses;i++)
	{
		//Parte 1
		double soma1=0.0;

		for(int k=0;k<totalPadroes;k++)
		{
			double distancia=distanciaEuclidianaAdaptativa(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getFeatureWeight());

			double u=pattern[k].getGrauPertinencia(i);
			
			soma1+=pow(u,m)*distancia;
		}

		//Parte 2
		double soma2=0.0;

		for(int j=0;j<numeroVariaveis;j++)
		{
			double v=cluster[i].getFeatureWeight(j);

			soma2+=pow(v,2.0);
		}

		W+=soma1+delta[i]*soma2;
	}

	return W;
}

//////////////////////////////////////////////////////////////////////

void SCAD::atualizaMelhorParticao()
{
	for(int	k=0;k<totalPadroes;k++)
	{
		for(int i=0;i<numeroClasses;i++)
		{
			melhorParticao[k][i]=pattern[k].getGrauPertinencia(i);
		}
	}
}

//////////////////////////////////////////////////////////////////////

int * SCAD::classesFinais()
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

	}

	return classesDosPadroes;
}

//////////////////////////////////////////////////////////////////////