// AdaptativoUnico.h: interface for the AdaptativoUnico class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYADAPTATIVOUNICO_H__BBA8A07E_3F6E_43EC_944F_36333FABA85B__INCLUDED_)
#define AFX_FUZZYADAPTATIVOUNICO_H__BBA8A07E_3F6E_43EC_944F_36333FABA85B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FuzzyMahalanobis.h"

class FuzzyAdaptativoUnico : public FuzzyMahalanobis  
{
public:
	void executarSimulacao();
	void executarAplicacao();
	FuzzyAdaptativoUnico(Pattern *p,char *nDados,char *nSaida);
	FuzzyAdaptativoUnico();
	virtual ~FuzzyAdaptativoUnico();

protected:
	double **M;
	double varianciaCovarianciaFuzzyAlterada(int var1,int var2,int i);
	double ** calculaMatrizVarianciaCovarianciaFuzzy(int i);
	void representacao();
};

#endif // !defined(AFX_FUZZYADAPTATIVOUNICO_H__BBA8A07E_3F6E_43EC_944F_36333FABA85B__INCLUDED_)
