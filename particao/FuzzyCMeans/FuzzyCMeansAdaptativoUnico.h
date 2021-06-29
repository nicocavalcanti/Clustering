// FuzzyCMeansAdaptativoUnico.h: interface for the FuzzyCMeansAdaptativoUnico class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYCMEANSADAPTATIVOUNICO_H__87E2C089_71CC_400F_B96D_EC2B7AE1A554__INCLUDED_)
#define AFX_FUZZYCMEANSADAPTATIVOUNICO_H__87E2C089_71CC_400F_B96D_EC2B7AE1A554__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "FuzzyCMeans.h"

class FuzzyCMeansAdaptativoUnico : public FuzzyCMeans  
{
public:
	FuzzyCMeansAdaptativoUnico(Pattern *p,char *nDados,char *nSaida);
	FuzzyCMeansAdaptativoUnico();
	virtual ~FuzzyCMeansAdaptativoUnico();
	void executarAplicacao();
	void executarSimulacao();
protected:
	void representacao();
	double * weightVector;
};

#endif // !defined(AFX_FUZZYCMEANSADAPTATIVOUNICO_H__87E2C089_71CC_400F_B96D_EC2B7AE1A554__INCLUDED_)
