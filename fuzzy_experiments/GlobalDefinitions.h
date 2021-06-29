#ifndef FUZZYEXPERIMENTS_GLOBAL_DEFINITIONS
  
  #define FUZZYEXPERIMENTS_GLOBAL_DEFINITIONS
  

  #include <limits>
  #include <vector>
  #include <time.h>
  using namespace std;
  //define types that will be used by all experiments
  typedef vector<double> tVecDouble;
  typedef vector<tVecDouble> tVecDouble2;
  typedef vector<tVecDouble2> tVecDouble3;
  typedef vector<tVecDouble3> tVecDouble4;

  typedef vector<int> tVecInt;
  typedef vector<tVecInt> tVecInt2;
  typedef vector<tVecInt2> tVecInt3;

  //#define INDF numeric_limits<double>::quiet_NaN()
  #define MAX_DOUBLE 1.7E+308
  #define MIN_DOUBLE -MAX_DOUBLE

  //#define DEBUG_CREATE_HARD_PARTITION
//#define DEBUG_MEMBERSHIP
//#define DEBUG_CLUSTER_CENTER

#define DEBUGTIME
#define NO_RUN_REPORT
//#define DEBUG_MEMBERSHIP_INTIALIZATION
//The following defines are intended to avoid some operations being made
#define NO_INDF_CHECK
#define NOT_CHECK_DIVISION_BY_ZERO
#define NO_INTERPRETATION_REPORT
#define NO_OSCILATION_CHECK
#define NO_DIST_IND_CENTER_ZERO_CHECK 
#define NO_INTERPRETATION_INDICE

#define NUM_OF_GLOBAL_PROTOTYPE 2

  extern double INDF;

  int random(int num);

  void randomize(void);
#endif//FUZZYEXPERIMENTS_GLOBAL_DEFINITIONS