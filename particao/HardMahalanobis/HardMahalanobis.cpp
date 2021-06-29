// HardMahalanobis.cpp: implementation of the HardMahalanobis class.
//
//////////////////////////////////////////////////////////////////////

#include "HardMahalanobis.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

HardMahalanobis::HardMahalanobis()
{

}

HardMahalanobis::~HardMahalanobis()
{
	delete melhorParticao;
}

HardMahalanobis::HardMahalanobis(Pattern *p,char *nSaida) : Particao(p,nSaida)
{
	melhorParticao=new int[totalPadroes];
}

//////////////////////////////////////////////////////////////////////

void HardMahalanobis::alocacao()
{
	for(int k=0;k<totalPadroes;k++)
	{
		double distanciaMinima=distanciaMahalanobis(pattern[k].getX(),cluster[0].getPrototipo(),cluster[0].getM());
		int c=0;
		for(int i=1;i<numeroClasses;i++)
		{
			double distancia=distanciaMahalanobis(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getM());
			
			if(distancia<distanciaMinima)
			{
				distanciaMinima=distancia;
				c=i;
			}
		}

		pattern[k].setClasse(c);
	}

/*	////////////////////////////////////////////
	//Monitoramento
	ofstream tout("Classes Atuais.wri",ios::ate);
	tout<<"\n\n------------- Iteração "<<iteracao<<" -----------------\n";

	for(k=0;k<totalPadroes;k++)
	{
		tout<<"padrão "<<k<<": "<<pattern[k].getClasse()<<"\n";
	}
	////////////////////////////////////////////*/
}

//////////////////////////////////////////////////////////////////////

double HardMahalanobis::calculaW()
{
	double W=0.0;

	for(int k=0;k<totalPadroes;k++)
	{
		int i=pattern[k].getClasse();

		double distancia=distanciaMahalanobis(pattern[k].getX(),cluster[i].getPrototipo(),cluster[i].getM());

		W+=distancia;

/*		//////////////////
		//Monitoramento
		ofstream fout("distancias.wri",ios::ate);

		fout<<"\n---------- Iteração "<<iteracao<<" ----------\n";
		fout<<"Padrão "<<k<<": "<<distancia;
		fout.close();
		/////////////////*/

	}

	return W;
}

//////////////////////////////////////////////////////////////////////

void HardMahalanobis::inicializacao()
{
	//Inicialização dos graus de pertinência
	for(int k=0;k<totalPadroes;k++)
	{
		int i=rand()%numeroClasses;

		pattern[k].setClasse(i);
	}

/*	////////////////////////////////////////////
	//Monitoramento
	ofstream tout("Classes Iniciais.wri",ios::ate);
	tout<<"\n\n";
	for(k=0;k<totalPadroes;k++)
	{
		tout<<"padrão "<<k<<": "<<pattern[k].getClasse()<<"\n";
	}
	////////////////////////////////////////////*/

}

//////////////////////////////////////////////////////////////////////

void HardMahalanobis::atualizaMelhorParticao()
{
	for(int k=0;k<totalPadroes;k++)
	{
		melhorParticao[k]=pattern[k].getClasse();
	}

/*	/////////////////////
	//Monitoramento
	ofstream fout("Melhor Partição.wri",ios::ate);
	fout<<"\n\n---------- Iteração "<<iteracao<<" -----------\n";
	
	for(k=0;k<totalPadroes;k++)
	{
		fout<<"\nPadrão "<<k<<": "<<melhorParticao[k];
	}

	fout.close();
	/////////////////////*/
}

//////////////////////////////////////////////////////////////////////

double HardMahalanobis::distanciaMahalanobis(double *x, double *g, double **M)
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