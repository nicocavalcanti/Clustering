#include "GlobalDefinitions.h"

double INDF = numeric_limits<double>::signaling_NaN();

  int random(int num){
	return (int)( ((long)rand()*num)/(RAND_MAX+1) );
  };

  void randomize(void){
    //srand((unsigned)time(NULL));
  };