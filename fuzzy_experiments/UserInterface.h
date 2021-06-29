// UserInterface.h: interface for the UserInterface class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_USERINTERFACE_H__E3738AB4_2D32_4D5C_95EC_541A16218BDE__INCLUDED_)
#define AFX_USERINTERFACE_H__E3738AB4_2D32_4D5C_95EC_541A16218BDE__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000


class UserInterface  
{
protected:
public:
	UserInterface();
	virtual void readParameters()=0;
	virtual ~UserInterface();

};

#endif // !defined(AFX_USERINTERFACE_H__E3738AB4_2D32_4D5C_95EC_541A16218BDE__INCLUDED_)
