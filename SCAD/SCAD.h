// SCAD.h: interface for the SCAD class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_SCAD_H__642CC1E6_89CD_4C90_A47A_5091D46D6A8F__INCLUDED_)
#define AFX_SCAD_H__642CC1E6_89CD_4C90_A47A_5091D46D6A8F__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "..\Particao.h"

class SCAD : public Particao  
{
public:
	int * executar();
	SCAD(Pattern *p,char *nSaida);
	SCAD();
	virtual ~SCAD();

private:
	double K;
	double *delta;
	double **melhorParticao;
protected:
	int * classesFinais();
	void atualizaMelhorParticao();
	double calculaW();
	void atualizaGrauPertinencia();
	double distanciaEuclidianaAdaptativa(double *x, double *g,double *v);
	void atualizaFeatureWeight();
	void atualizaPrototipos();
	void atualizaDelta();
	void inicializacao();
};

#endif // !defined(AFX_SCAD_H__642CC1E6_89CD_4C90_A47A_5091D46D6A8F__INCLUDED_)
