// HardAdaptativo.h: interface for the HardAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_HARDADAPTATIVO_H__1F34615E_5332_4342_BE90_3FBFAAF6D7D2__INCLUDED_)
#define AFX_HARDADAPTATIVO_H__1F34615E_5332_4342_BE90_3FBFAAF6D7D2__INCLUDED_

#include "HardMahalanobis.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class HardAdaptativo  : public HardMahalanobis
{
public:
	int * executar();
	HardAdaptativo(Pattern *p,char *nSaida);
	HardAdaptativo();
	virtual ~HardAdaptativo();

protected:
	double varianciaCovarianciaHard(int var1,int var2,int i);
	void representacao();
	double ** calculaMatrizVarianciaCovarianciaHard(int i);
};

#endif // !defined(AFX_HARDADAPTATIVO_H__1F34615E_5332_4342_BE90_3FBFAAF6D7D2__INCLUDED_)
