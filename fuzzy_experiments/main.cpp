#include "Application.h"

/*#include <stdio.h>
#include <stdlib.h>



void main( void )
{
string buf;
buf.assign("teste.txt");
buf.replace(5,8,"_AUX.wri");
   char buffer[50];
   double source = -3.1415;
   _gcvt( source, 10, buffer );
   printf( "source: %f  buffer: '%s'\n", source, buffer );
   _gcvt( source, 10, buffer );
   printf( "source: %e  buffer: '%s'\n", source, buffer );
}*/





int main(int argc, char **argv){
	//char c='M',b='F',d='I';
	//printf("%d %d %d",c,b,d);
	Application app;
	app.execute();
	return 0;
}
