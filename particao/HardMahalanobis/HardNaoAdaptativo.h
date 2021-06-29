// HardNaoAdaptativo.h: interface for the HardNaoAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_HARDNAOADAPTATIVO_H__3B159F8C_C040_40AD_8D3A_FEBC603770D1__INCLUDED_)
#define AFX_HARDNAOADAPTATIVO_H__3B159F8C_C040_40AD_8D3A_FEBC603770D1__INCLUDED_

#include "HardMahalanobis.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class HardNaoAdaptativo  : public HardMahalanobis
{
public:
	HardNaoAdaptativo(Pattern *p,char *nSaida);
	HardNaoAdaptativo();
	virtual ~HardNaoAdaptativo();
	int * executar();

protected:

	void representacao();
	double varianciaCovariancia(int var1,int var2);
	void calculaMatrizVarianciaCovariancia();
	void calculaM();
private:
	double **Q;
	double *prototipoGeral;
	double **M;
};

#endif // !defined(AFX_HARDNAOADAPTATIVO_H__3B159F8C_C040_40AD_8D3A_FEBC603770D1__INCLUDED_)
