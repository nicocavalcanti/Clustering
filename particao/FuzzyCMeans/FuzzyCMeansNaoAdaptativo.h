// FuzzyCMeansNaoAdaptativo.h: interface for the FuzzyCMeansNaoAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYCMEANSNAOADAPTATIVO_H__59513D47_1EB4_446E_A4A5_9FF38BB15CA5__INCLUDED_)
#define AFX_FUZZYCMEANSNAOADAPTATIVO_H__59513D47_1EB4_446E_A4A5_9FF38BB15CA5__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FuzzyCMeans.h"

class FuzzyCMeansNaoAdaptativo : public FuzzyCMeans  
{
public:
	void executarAplicacao();
	void executarSimulacao();
	FuzzyCMeansNaoAdaptativo(Pattern *p,char *nDados,char *nSaida);
	FuzzyCMeansNaoAdaptativo();
	virtual ~FuzzyCMeansNaoAdaptativo();

protected:
	void representacao();
};

#endif // !defined(AFX_FUZZYCMEANSNAOADAPTATIVO_H__59513D47_1EB4_446E_A4A5_9FF38BB15CA5__INCLUDED_)
