// HardAdaptativoUnico.h: interface for the HardAdaptativoUnico class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_HARDADAPTATIVOUNICO_H__68DA4EB2_32FC_4EB8_862C_EB9337842B53__INCLUDED_)
#define AFX_HARDADAPTATIVOUNICO_H__68DA4EB2_32FC_4EB8_862C_EB9337842B53__INCLUDED_

#include "HardMahalanobis.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class HardAdaptativoUnico  : public HardMahalanobis
{
public:
	int * executar();
	HardAdaptativoUnico(Pattern *p,char *nSaida);
	HardAdaptativoUnico();
	virtual ~HardAdaptativoUnico();

private:
	double **M;
protected:
	void representacao();
	double varianciaCovarianciaHardAlterada(int var1,int var2,int i);
	double ** calculaMatrizVarianciaCovarianciaHard(int i);
};

#endif // !defined(AFX_HARDADAPTATIVOUNICO_H__68DA4EB2_32FC_4EB8_862C_EB9337842B53__INCLUDED_)
