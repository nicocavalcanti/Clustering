// NaoAdaptativo.h: interface for the NaoAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYNAOADAPTATIVO_H__FC5671BB_015B_45F1_90A2_48C40F309E81__INCLUDED_)
#define AFX_FUZZYNAOADAPTATIVO_H__FC5671BB_015B_45F1_90A2_48C40F309E81__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FuzzyMahalanobis.h"

class FuzzyNaoAdaptativo : public FuzzyMahalanobis  
{
private:
	double **M;
	double **Q;					//Matriz de variâncias e covariâncias
	double *prototipoGeral;		//Este vetor representa um protótipo geral: a média dos dados

public:
	void executarSimulacao();
	void executarAplicacao();
	FuzzyNaoAdaptativo();
	FuzzyNaoAdaptativo(Pattern *p,char *nDados,char *nSaida);
	virtual ~FuzzyNaoAdaptativo();
protected:
	void calculaM();
	void representacao();
	void calculaMatrizVarianciaCovariancia();
	double varianciaCovariancia(int var1,int var2);
};

#endif // !defined(AFX_FUZZYNAOADAPTATIVO_H__FC5671BB_015B_45F1_90A2_48C40F309E81__INCLUDED_)
