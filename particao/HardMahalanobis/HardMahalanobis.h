// HardMahalanobis.h: interface for the HardMahalanobis class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_HARDMAHALANOBIS_H__C7188FBD_291A_4FD3_930E_7881C77068D5__INCLUDED_)
#define AFX_HARDMAHALANOBIS_H__C7188FBD_291A_4FD3_930E_7881C77068D5__INCLUDED_

#include "..\Particao.h"	// Added by ClassView

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class HardMahalanobis  : public Particao
{
public:
	HardMahalanobis(Pattern *p,char *nSaida);
	HardMahalanobis();
	virtual ~HardMahalanobis();

protected:
	double distanciaMahalanobis(double *x, double *g, double **M);
	void atualizaMelhorParticao();
	void inicializacao();
	double calculaW();
	void alocacao();
	int *melhorParticao;
};

#endif // !defined(AFX_HARDMAHALANOBIS_H__C7188FBD_291A_4FD3_930E_7881C77068D5__INCLUDED_)
