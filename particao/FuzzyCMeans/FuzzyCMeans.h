// FuzzyCMeans.h: interface for the FuzzyCMeans class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_FUZZYCMEANS_H__D026EB14_F53D_48AE_9759_458AADD2A3C0__INCLUDED_)
#define AFX_FUZZYCMEANS_H__D026EB14_F53D_48AE_9759_458AADD2A3C0__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "..\Particao.h"

class FuzzyCMeans : public Particao  
{
public:
	double getMelhorDC();
	double getMelhorB();
	double getMelhorR();
	double getMelhorT();
	double getMelhorJ();
	FuzzyCMeans(Pattern *p,char *nDados,char *nSaida);
	FuzzyCMeans();
	virtual ~FuzzyCMeans();

protected:
	void imprimeVariaveisMonitoramento();
	void calculaRelativeB();
	void calculaRelativeJ();
	void calculaRelativeT();
	double *relativeB;
	double *relativeJ;
	double *relativeT;
	void calculaIndicesClasseVariavel();
	void calculaTij();
	void calculaTj();
	void calculaTi();
	void calculaJij();
	void calculaJj();
	void calculaJi();
	void calculaMi();
	double *mi;
	void calculaBij();
	void calculaBj();
	void calculaBi();
	double **melhorParticao;
	double *Jj;
	double **Tij;
	double **Jij;
	void imprimeInterpretadores();
	void calculaCEij();
	double **CEij;
	double melhorDC;
	void setClassificacaoHard();
	double *Tj;
	void calculaCTRij();
	void calculaCTRj();
	void calculaCORij();
	void calculaCORj();
	double **CTRij;
	double *CTRj;
	double **CORij;
	double *CORj;
	void setMelhoresIndices();
	double *Bj;
	double **Bij;
	double *Ti;
	double *Bi;
	double *Ji;
	double melhorR;
	double melhorB;
	double melhorT;
	double melhorJ;
	double B;
	double J;
	double R;
	double T;
	double DC;
	void setMelhorParticao(int *classes);
	double **melhorWeightVector;
	void calculaIndices();
	void calculaDC();
	void calculaR();
	double distanciaL2MinkowskyAdaptativo(double *x,double *g,double *weightVector);
	void calculaB();
	void calculaJ();
	void calculaT();
	void calculaZ();
	double *z;					//Overall prototype
	double calculaW();
	void alocacao();
	void inicializacao();
	int * classesFinais();
	void atualizaMelhorParticao();
};

#endif // !defined(AFX_FUZZYCMEANS_H__D026EB14_F53D_48AE_9759_458AADD2A3C0__INCLUDED_)
