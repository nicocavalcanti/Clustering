// FuzzyCMeansAdaptativo.h: interface for the FuzzyCMeansAdaptativo class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYCMEANSADAPTATIVO_H__52E3A27B_016E_405C_8E16_4B8155B49DA1__INCLUDED_)
#define AFX_FUZZYCMEANSADAPTATIVO_H__52E3A27B_016E_405C_8E16_4B8155B49DA1__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FuzzyCMeans.h"

class FuzzyCMeansAdaptativo : public FuzzyCMeans  
{
public:
	void executarSimulacao();
	void executarAplicacao();
	FuzzyCMeansAdaptativo(Pattern *p,char *nDados,char *nSaida);
	FuzzyCMeansAdaptativo();
	virtual ~FuzzyCMeansAdaptativo();

protected:
	void representacao();
};

#endif // !defined(AFX_FUZZYCMEANSADAPTATIVO_H__52E3A27B_016E_405C_8E16_4B8155B49DA1__INCLUDED_)
