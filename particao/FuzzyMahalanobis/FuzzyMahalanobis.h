// FuzzyMahalanobis.h: interface for the FuzzyMahalanobis class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYMAHALANOBIS_H__C6725CD1_36CC_40CE_8D71_CB0B56091FFF__INCLUDED_)
#define AFX_FUZZYMAHALANOBIS_H__C6725CD1_36CC_40CE_8D71_CB0B56091FFF__INCLUDED_

#include "fstream.h"
#include "..\Particao.h"	// Added by ClassView

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class FuzzyMahalanobis : public Particao  
{
protected:
	void imprimeVariaveisMonitoramento();
	void calculaIndicesClasseVariavel();
	void imprimeInterpretadores();
	void calculaRelativeJ();
	void calculaRelativeB();
	void calculaRelativeT();
	double *relativeB;
	double *relativeJ;
	double *relativeT;
	void calculaR();
	double melhorR;
	double R;
	void calculaBi();
	double *Bi;
	void calculaB();
	double melhorB;
	double B;
	void calculaMi();
	double *mi;
	void calculaJi();
	double *Ji;
	double melhorT;
	void calculaTi();
	double *Ti;
	double T;
	void calculaT();
	double ***melhorM;
	void setMelhorParticao(int *classes);
	void calculaIndices();
	void setClassificacaoHard();
	void setMelhoresIndices();
	double melhorDC;
	double melhorJ;
	double DC;
	double J;
	void calculaDC();
	double distanciaMahalanobis(double *x,double *g,double **M);
	void calculaZ();
	double *z;					//Overall prototype
	int * classesFinais();
	void atualizaMelhorParticao();
	double **melhorParticao;
	void calculaJ();
	void alocacao();
	void inicializacao();
public:
	double getMelhorB();
	double getMelhorR();
	double getMelhorT();
	double getMelhorJ();
	double getMelhorDC();
	FuzzyMahalanobis(Pattern *p,char *nDados,char *nSaida);
	FuzzyMahalanobis();
	virtual ~FuzzyMahalanobis();
};

#endif // !defined(AFX_FUZZYMAHALANOBIS_H__C6725CD1_36CC_40CE_8D71_CB0B56091FFF__INCLUDED_)
