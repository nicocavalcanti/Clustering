// FuzzyAdaptativo.h: interface for the FuzzyAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYADAPTATIVO_H__57B07C3B_AA0A_4EEB_AB36_26EF8E25333A__INCLUDED_)
#define AFX_FUZZYADAPTATIVO_H__57B07C3B_AA0A_4EEB_AB36_26EF8E25333A__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FuzzyMahalanobis.h"

class FuzzyAdaptativo : public FuzzyMahalanobis  
{
public:
	void executarSimulacao();
	void executarAplicacao();
	FuzzyAdaptativo(Pattern *p,char *nDados,char *nSaida);
	FuzzyAdaptativo();
	virtual ~FuzzyAdaptativo();

protected:
	double varianciaCovarianciaFuzzy(int var1,int var2,int i);
	void representacao();
	double ** calculaMatrizVarianciaCovarianciaFuzzy(int i);
};

#endif // !defined(AFX_FUZZYADAPTATIVO_H__57B07C3B_AA0A_4EEB_AB36_26EF8E25333A__INCLUDED_)
