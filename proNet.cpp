#include <err.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <time.h>
#include <pthread.h>
#include <stdio.h>
#include <unistd.h>  // usleep()
#include <math.h>
#include <algorithm>
#include <cmath>

//#ifndef KEI
//#include <regex>
//#endif

using namespace std;

// output file stream
// Pactual 
ofstream ofs;

// pvalue
ofstream pvfs;

// global variable
// defauld 10 step
int simulationStep=10;

// in files

// io line (firm base)
string ioLineFile;

// c vector (firm base)
string cVectorFile;

// pini vector (firm base)
string piniVectorFile;

// optional
// delta vector
// i delta
// (0 if there is no input)
string deltaFile="";

// out files

// Pactual file name
string pactOutputFile;

// value added file name
string pValueAddedOutputFile;

// Final consumer file name
// and file stream
string finalConsumeOutputFile="";
ofstream  fcOfs;

// Stats file name
// and file stream
string statsOutputFile="";
ofstream  stOfs;

// Adaptation file name
// and file stream
string adaptationOutputFile="";
ofstream  adapOfs;

// debug
// iniとactの比
// string  piniPactRatioFile="piniPactRatio.txt";
// ofstream  piniPactRatioFS;

// debug
// damage amount
// string  damagedFirmFile="damagedFirm.txt";
// ofstream  damageFS;
// int damagedNum;

// debug
// full debug option (output for step by step)
// t, i, pini, delta, pcap, d, pmax, pact, d*, s_j,...
int iFullDebugFrag=0;
string  fullDebugFile="fullDebug.csv";
ofstream  debugFS;

// rationing algorithm switch
int iRationingTypeFrag=0;

// order algorithm switch
int iOrderTypeFrag=0;

// internal structure

// A table 
// direction is the same as Hallegatte's paper
// it means the direction is the opposite to IOtables
// supplier -> client volume
// forward: client -> supplier volume
// backward: supplier -> client volume
unordered_map<int, unordered_map<int, double> > fATableHoH;
unordered_map<int, unordered_map<int, double> > bATableHoH;

// C vector 
unordered_map<int, double> cVectorH;

// Pini vector
// total amount output at initial
unordered_map<int, double> piniVectorH;

// input amount
unordered_map<int, double> piniInputVectorH;

// valueAdded ratio
unordered_map<int, double> valueAddedVectorH;

// firm affi hash
unordered_map<int, int> firmAffiH;

// firm name hash
unordered_map<int, int> firmH;

// order between firm client -order-> supplier
// this tells the total input of client
unordered_map<int, unordered_map<int, double> > fOFirmFirmHoH;

// adjusted order between firm client -order-> supplier
unordered_map<int, unordered_map<int, double> > fOAdjustedFirmFirmHoH;

// backward order between firm supplier <-order- client
unordered_map<int, unordered_map<int, double> > bOFirmFirmHoH;

// backward order between firm supplier <-order- client at start
unordered_map<int, unordered_map<int, double> > bOStartFirmFirmHoH;

// backward adjusted order between firm supplier <-order- client
unordered_map<int, unordered_map<int, double> > bOAdjustedFirmFirmHoH;

// demand on firm
unordered_map<int, double> DFirmVectorH;

// adjusted demand on firm
unordered_map<int, double> DAdjustedFirmVectorH;

// adjusted c on firm
unordered_map<int, double> cAdjustedVectorH;

// Pini Cap vector
unordered_map<int, double> piniCapVectorH;

// S on firm for firm (client <-- supplier)
unordered_map<int, unordered_map<int, double> > SFirmFirmHoH;

// S on firm for sector (total)
unordered_map<int, unordered_map<int, double> > SFirmSectorHoH;

// A on firm for sector (total)
unordered_map<int, unordered_map<int, double> > AFirmSectorHoH;

// Pprop on firm for sector
// PsProp
// Ratio based inventor for each sector(product)
// (However, here we use the ratio multiplied by Pini
unordered_map<int, unordered_map<int, double> > PpropFirmSectorHoH;

// Pini Max vector
// (not ratio, ammount)
unordered_map<int, double> PiniMaxVectorH;

// Pactual vector
unordered_map<int, double> PactualVectorH;

// Damage vector 0<=delta<=1
unordered_map<int, double> DeltaVectorH;

// file stream
ifstream ifsDelta;

// data for Adaptation

// adaptation type
int iAdapType=-1;

// adapatation max counter
unordered_map<int, int> iAdaptH;

// adapatation done counter
unordered_map<int, int> iAdaptDoneH;

// secFirmHoL
unordered_map<int, vector<int> > secFirmHoL;

// link that have taken adaptation
unordered_map<int, vector<int> > adaptHoL;

// inventory size(day)
double inventoryN=15;

// adjustment ratio for inventory
double tau=6;

// for output stats
double lastTotalPactual;
double lastTotalValueAdded;

// get snapshot pact option
int iPrintSnapFrag=0;

// snapshot step
unordered_map<int,int> hPrintSnapStep;
string sPrintSnapStep;

// snapshot file
string sPrintSnapFileBase="";

// functions

#ifndef Pi 
#define Pi 3.141592653589793238462643 
#endif 

//#if defined(INTLC) || defined(KEI)
#if defined(INTLC)
inline bool isEqual(double x, double y)
{
  const double epsilon = 1e-5 /* some small number such as 1e-5 */;
  return std::abs(x - y) <= epsilon * std::abs(x);
  // see Knuth section 4.2.2 pages 217-218
}

#else

template<class T>
typename std::enable_if<!std::numeric_limits<T>::is_integer, bool>::type
    almost_equal(T x, T y, int ulp)
{
    // the machine epsilon has to be scaled to the magnitude of the values used
    // and multiplied by the desired precision in ULPs (units in the last place)
    return std::abs(x-y) < std::numeric_limits<T>::epsilon() * std::abs(x+y) * ulp
    // unless the result is subnormal
           || std::abs(x-y) < std::numeric_limits<T>::min();
}

#endif

double cnd_manual(double x)
{
  double L, K, w ;
  /* constants */
  double const a1 = 0.31938153, a2 = -0.356563782, a3 = 1.781477937;
  double const a4 = -1.821255978, a5 = 1.330274429;

  L = fabs(x);
  K = 1.0 / (1.0 + 0.2316419 * L);
  w = 1.0 - 1.0 / sqrt(2 * Pi) * exp(-L *L / 2) * (a1 * K + a2 * K *K + a3 * pow(K,3) + a4 * pow(K,4) + a5 * pow(K,5));

  if (x < 0 ){
    w= 1.0 - w;
  }
  return w;
}

double adapProbFunc(double x){

  double ret;

  // standart sigmoid

  // debug
//   cerr << "x: " << x << endl;
//   cerr << "exp(-(x-0.5)) " << exp(-(x-0.5)) << endl;
//   cerr << "(1.0+exp(-(x-0.5)) " << (1.0+exp(-(x-0.5))) << endl;
//   cerr << "(1.0/(1.0+exp(-(x-0.5))) " << (1.0/(1.0+exp(-(x-0.5)))) << endl;

//   ret=1.0-(1.0/(1.0+exp(-(x-0.5))));

//   // debug
//   cerr << "ret: " << ret << endl;

  // normal cumulative gaussiun
  
  // linear
  // ret=x;

  // square
  ret=pow(x, 2.0);

  // debug
  //   cerr << "x: " << x << endl;
  //   cerr << "exp(-(x-0.5)) " << exp(-(x-0.5)) << endl;
  //   cerr << "(1.0+exp(-(x-0.5)) " << (1.0+exp(-(x-0.5))) << endl;
  //   cerr << "(1.0/(1.0+exp(-(x-0.5))) " << (1.0/(1.0+exp(-(x-0.5)))) << endl;

  // alpha
  // double alpha=10;

  //  ret=(1.0/(1.0+exp(-alpha*(x-0.5))));

  return ret;
}

void printHash1LayerIntDouble(char const *name, unordered_map<int, double> *hash){

  cout << name << " is" << endl;
  map<int, double> sortedVectorH(hash->begin(), hash->end());
  for(auto itr = sortedVectorH.begin(); itr != sortedVectorH.end(); ++itr) {
    // show keys values
    cout << "key:" << itr->first << " val:" << itr->second << endl;
  }

}

void printHash1LayerIntInt(char const *name, unordered_map<int, int> *hash){

  cout << name << " is" << endl;
  map<int, int> sortedVectorH(hash->begin(), hash->end());
  for(auto itr = sortedVectorH.begin(); itr != sortedVectorH.end(); ++itr) {
    // show keys values
    cout << "key:" << itr->first << " val:" << itr->second << endl;
  }

}

void printVector1LayerString(char const *name, vector<string> *sv){

  cout << name << " is" << endl;
  int i=0;
  for(auto itr = sv->begin(); itr != sv->end(); ++itr) {
    // show values
    cout << "index:" << i << " val:" << *itr << endl;
    i++;
  }

}

void printHash2LayerIntIntDouble(char const *name, unordered_map<int, unordered_map<int, double> > *hash){

  cout << name << " is" << endl;

  for(auto itr = hash->begin(); itr != hash->end(); ++itr) {
    // show keys
    cout << "from: " << itr->first << endl;
    map<int, double> sortedVector2H((itr->second).begin(), (itr->second).end());

    for(auto itr2 = sortedVector2H.begin(); itr2 != sortedVector2H.end(); ++itr2) {
      // show keys
      cout << "    to: " << itr2->first << endl;
      // show values
      cout << "      val = " << itr2->second << "\n";
    }
  }
}

void printHashVectorIntInt(char const *name, unordered_map<int, vector<int> > *hash){

  cout << name << " is" << endl;

  for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

    // show keys
    cout << "sector: " << itr->first << endl;

    cout << "      val =";
    for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {
      cout << " " << *itr2;
    }
    cout << endl;
  }
}


static void usage(char *argv){
  cerr << "Usage: " << endl;
  cerr << "  " << argv << " [option] aTableLineFile cVectorFile piniVectorFile firmAffiFile pactOutputFile" << endl;
  cerr << "Option: " << endl;
  fprintf(stderr,"  -a (int): adaptation (currently, 0: infty, 1: only once)\n");
  fprintf(stderr,"  -A (str): adaptation counter file iff -a is on\n");
  fprintf(stderr,"  -d (str): delta file (default null)\n");
  fprintf(stderr,"  -f (str): final consumption output file (default null)\n");
  fprintf(stderr,"  -F : fulldebug\n");
  fprintf(stderr,"  -o (int): order algorithm type. 0: normal, 1: keep initial demand, 2: ignore negative inv adjustment\n");
  fprintf(stderr,"  -p (file:int:int...): get snapshot of production at indicated step(s). for more than one snapshot, int should be separated by :\n");
  fprintf(stderr,"  -r (int): random seed (default 331)\n");
  fprintf(stderr,"  -R (int): rationing algorithm 0. proportional, 1. FC low priority, 2. lower has high priority, 3. 1*2 (default proportional)\n");
  fprintf(stderr,"  -s (str): stats output file (default null)\n");
  fprintf(stderr,"  -t (int): simulation step (default 10)\n");
  fprintf(stderr,"  -v (str): value added output file (default null)\n");
  exit(1);
}

// full debug
void fullDebug(int t){

  unordered_map<int, int> *hashD;
  hashD=&firmH;
  map<int, int> sortedVectorH(hashD->begin(), hashD->end());
  for(auto itr = sortedVectorH.begin(); itr != sortedVectorH.end(); ++itr) {
    // show keys values
    // i, pini, pact, ratio

    int i = itr->first;
    double pini = piniVectorH[i];
    double delta = DeltaVectorH[i];
    double pcap = piniCapVectorH[i];
    double demand = DFirmVectorH[i];
    double pmax = PiniMaxVectorH[i];
    double pact = PactualVectorH[i];
    double demandA = DAdjustedFirmVectorH[i];
    unordered_map<int, double> *hashS=&SFirmFirmHoH[i];
    unordered_map<int, double> *hashO=&fOFirmFirmHoH[i];

    debugFS << t << ", " << i << ", " << pini << ", " << delta << ", " << pcap << ", " << demand << ", " << pmax << ", " << pact << ", " << demandA ;

    //       // inventory s
    map<int, int> sortedVector2H(hashS->begin(), hashS->end());
    for(auto itr = sortedVector2H.begin(); itr != sortedVector2H.end(); ++itr) {

      debugFS << ", s(" << i << "," << itr->first << ")_" << itr->second;

    }

    // order
    map<int, int> sortedVector3H(hashO->begin(), hashO->end());
    for(auto itr = sortedVector3H.begin(); itr != sortedVector3H.end(); ++itr) {

      debugFS << ", O(" << i << "," << itr->first << ")_" << itr->second;

    }
    debugFS << endl;
      
  }

}


void printAllInternal(){
  cerr << "----- All Internal -----" << endl;

  // debug
  printHash2LayerIntIntDouble("SFirmFirmHoH", &SFirmFirmHoH);

  // debug
  printHash2LayerIntIntDouble("fOFirmFirmHoH", &fOFirmFirmHoH);

  // debug
  printHash2LayerIntIntDouble("bOFirmFirmHoH", &bOFirmFirmHoH);

  // debug
  printHash1LayerIntDouble("DFirmVectorH", &DFirmVectorH);

  // debug
  printHash2LayerIntIntDouble("AFirmSectorHoH", &AFirmSectorHoH);

  // debug
  printHash2LayerIntIntDouble("SFirmSectorHoH", &SFirmSectorHoH);

  // debug
  printHash2LayerIntIntDouble("PpropFirmSectorHoH", &PpropFirmSectorHoH);

  // debug
  printHash1LayerIntDouble("PactualVectorH", &PactualVectorH);

  // debug
  printHash1LayerIntDouble("piniVectorH", &piniVectorH);

  // debug
  printHash1LayerIntDouble("piniInputVectorH", &piniInputVectorH);
  
  // debug
  printHash2LayerIntIntDouble("fOAdjustedFirmFirmHoH", &fOAdjustedFirmFirmHoH);

  // debug
  printHash2LayerIntIntDouble("bOAdjustedFirmFirmHoH", &bOAdjustedFirmFirmHoH);

  // debug
  printHash1LayerIntDouble("DAdjustedFirmVectorH", &DAdjustedFirmVectorH);

}


void initializeModel(void){

  // create initial (t=0)
  
  // S(i,j)
  // SFirmFirmHoH

  // order
  // O(i,j)
  // fOFirmFirmHoH
  // O*(i,j)
  // fOAdjustedFirmFirmHoH

  // opposite direction of order
  // bOFirmFirmHoH
  // bOAdjustedFirmFirmHoH

  // using fAtable
  for(auto itr = fATableHoH.begin(); itr != fATableHoH.end(); ++itr) {
    // row direction: client
    for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {
      
      // column direction: supplier

      // does client exists in hash?

      // forward creating hash
      unordered_map<int, unordered_map<int, double> >::iterator findI;
      findI = SFirmFirmHoH.find(itr->first);

      if(findI == SFirmFirmHoH.end()){
	// no

	// debug
	//      cout << "not found" <<endl;

	unordered_map<int, double> SFirmH;
	SFirmFirmHoH[itr->first]=SFirmH;

	unordered_map<int, double> OFirmH;
	fOFirmFirmHoH[itr->first]=OFirmH;

	unordered_map<int, double> OAdjustedFirmH;
	fOAdjustedFirmFirmHoH[itr->first]=OAdjustedFirmH;

      }else{
	// yes
	// debug
	//     cout << "found" << endl;
      }

      (SFirmFirmHoH[itr->first])[itr2->first]=itr2->second * inventoryN;
      (fOFirmFirmHoH[itr->first])[itr2->first]=itr2->second;
      (fOAdjustedFirmFirmHoH[itr->first])[itr2->first]=itr2->second;

      // input vector
      piniInputVectorH[itr->first]=piniInputVectorH[itr->first]+itr2->second;

      // backward creating hash
      unordered_map<int, unordered_map<int, double> >::iterator findI2;
      findI2 = bOFirmFirmHoH.find(itr2->first);

      if(findI2 == bOFirmFirmHoH.end()){
	// no

	// debug
	//      cout << "not found" <<endl;

	unordered_map<int, double> bOFirmH;
	bOFirmFirmHoH[itr2->first]=bOFirmH;

	unordered_map<int, double> bOAdjustedFirmH;
	bOAdjustedFirmFirmHoH[itr2->first]=bOAdjustedFirmH;

      }else{
	// yes
	// debug
	//     cout << "found" << endl;
      }

      (bOFirmFirmHoH[itr2->first])[itr->first]=itr2->second;

      // debug
//       if(isnan(itr2->second)){
// 	cerr << itr2->first << " " << itr->first << endl;
// 	exit(1);
//       }

      (bOAdjustedFirmFirmHoH[itr2->first])[itr->first]=itr2->second;

    }
  }

  // initialize start hash
  bOStartFirmFirmHoH=bOFirmFirmHoH;

  // calculate value added
  for(auto itr = piniVectorH.begin(); itr != piniVectorH.end(); ++itr) {
    int firm=itr->first;
    double output=itr->second;

    //debug
    //    cout << "firm " << firm << " output " << output << endl;

    if(output<=0){
      cout << "firm " << firm << " output(Pini) is <=0" << endl;
    }

    double input=piniInputVectorH[firm];

    if(output<=0){
      valueAddedVectorH[itr->first]=0;

      // debug
      cout << "firm " << firm << " output is <=0" << endl;
      
    }else if((output-input)<=0){

      // debug
      cout << "firm " << firm << " value added is <=0" << endl;

      //      valueAddedVectorH[itr->first]=0;

      // allow it
      valueAddedVectorH[itr->first]=(output-input)/output;
	    
    }else{
      
      // output-input / output
      valueAddedVectorH[itr->first]=(output-input)/output;
    }

    // debug
    //    cerr << "firm " << firm << endl;
    //    cerr << "input " << input << endl;
    //    cerr << "output " << output << endl;
    //    cerr << "value added " << valueAddedVectorH[firm] << endl;
  }
  

  // demand

  // initialize
  // bOFirmFirmHoH
  // cAdjusted

  // cVectorH is defined for every firm

  unordered_map<int, int> *hashII;
  hashII=&firmH;
    
  for(auto itr = hashII->begin(); itr != hashII->end(); ++itr) {
    DFirmVectorH[itr->first]=cVectorH[itr->first];
  }

  // initilize cAdjusted
  cAdjustedVectorH=cVectorH;

  for(auto itr = bOFirmFirmHoH.begin(); itr != bOFirmFirmHoH.end(); ++itr) {
    // fOrder: client -> supplier
    // bOrder: supplier -> client
    // demand means orders toward the firm

    double totalDemand=0;
    for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {
      totalDemand=totalDemand+itr2->second;
    }

    DFirmVectorH[itr->first]=DFirmVectorH[itr->first]+totalDemand;

//     // debug
//     //    if(itr->first==400227878){521036828
//     if(itr->first==521036828){
//       cerr <<     DAdjustedFirmVectorH[itr->first] << endl;
//       exit(1);
//     }

  }

  DAdjustedFirmVectorH=DFirmVectorH;

  // copy Actual to D

  PactualVectorH=DFirmVectorH;
  PiniMaxVectorH=DFirmVectorH;
  piniCapVectorH=DFirmVectorH;

  // initialize adaptation

  if(iAdapType!=-1){

    int initNum=iAdapType;

    if(iAdapType==0){
      initNum=-1;
    }

    unordered_map<int, int> *hash;
    hash=&firmH;

    for(auto itr = hash->begin(); itr != hash->end(); ++itr) {
      iAdaptH[itr->first]=initNum;
      iAdaptDoneH[itr->first]=0;
    }
  }

  // output

  // pact
   double totalPactual=0;
   for(auto itr = PactualVectorH.begin(); itr != PactualVectorH.end(); ++itr) {
     totalPactual=totalPactual+itr->second;
   }
  // disaster happes at 0, so as a pre-disaster, t=-1 is output
  ofs << -1 << " " << totalPactual << endl;  

  double totalVAdded=0;
  for(auto itr = PactualVectorH.begin(); itr != PactualVectorH.end(); ++itr) {
    totalVAdded=totalVAdded+valueAddedVectorH[itr->first]*itr->second;
  }

  if(pvfs){
    pvfs << -1 << " " << totalVAdded << endl;  
  }

  // create A total (i,s) (Fixed)
  // create S total (i,s)

  unordered_map<int, unordered_map<int, double> > *hash;

  hash=&fATableHoH;

  for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

    // itr->first is i
    int i=itr->first;

    unordered_map<int, unordered_map<int,double> >::iterator findI;
    findI = AFirmSectorHoH.find(i);

    if(findI == AFirmSectorHoH.end()){
      // not found

      // debug
      //      cout << "not found" <<endl;

      unordered_map<int, double> ASectorH;
      AFirmSectorHoH[i]=ASectorH;

      unordered_map<int, double> SSectorH;
      SFirmSectorHoH[i]=SSectorH;

    }else{
      // yes
      // debug
      //     cout << "found" << endl;
    }

    for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {

      // itr2->first is j
      int j=itr2->first;

      // itr2->second is A(i,j)
      double a=itr2->second;

      // get sector
      int sj=firmAffiH[j];

      unordered_map<int,double>::iterator findI2;
      findI2 = AFirmSectorHoH[i].find(sj);

      if(findI2 == AFirmSectorHoH[i].end()){
	// not found

	// debug
	//      cout << "not found" <<endl;

	(AFirmSectorHoH[i])[sj]=0;

	(SFirmSectorHoH[i])[sj]=0;

      }else{
	// yes
	// debug
	//     cout << "found" << endl;
      }

      (AFirmSectorHoH[i])[sj]=(AFirmSectorHoH[i])[sj]+a;

      (SFirmSectorHoH[i])[sj]=(SFirmSectorHoH[i])[sj]+inventoryN*a;

    }
  }

  // final consume file open
  if(finalConsumeOutputFile.size()>0){
    fcOfs.open(finalConsumeOutputFile.c_str());
    if(!fcOfs){
      cerr << "Unable to open " << finalConsumeOutputFile << " as finalConsumeOutputFile" << endl;
      exit(1);
    }
  }

  // stats file open
  if(statsOutputFile.size()>0){
    stOfs.open(statsOutputFile.c_str());
    if(!stOfs){
      cerr << "Unable to open " << statsOutputFile << " as statsOutputFile" << endl;
      exit(1);
    }
  }

  // output to stat file
  if(stOfs){
        stOfs << "totalPactutalStart " << totalPactual << endl;
        stOfs << "totalValueAddedStart " << totalVAdded << endl;
  }

  // output
  // final consumption
  if(fcOfs){
    double totalFC=0;
    for(auto itr = cAdjustedVectorH.begin(); itr != cAdjustedVectorH.end(); ++itr) {
      totalFC=totalFC+itr->second;
    }
    fcOfs << "-1 " << totalFC << endl;
  }

  // debug
  //    printAllInternal();
  //    exit(1);
  
}

void hazard(){

  // initialize deltaVector with 0
  for(auto itr = piniVectorH.begin(); itr != piniVectorH.end(); ++itr) {
    DeltaVectorH[itr->first]=0;
  }

  // deltaFile
  if(deltaFile.size()>0){

    string str;
    while(getline(ifsDelta,str)){
      string tmp;
      istringstream stream(str);
      vector<string> result;
      while(getline(stream,tmp,' ')){

	// debug
	//	  cout<< tmp << endl;

	result.push_back(tmp);
      }

      // does source exists in hash?
      unordered_map<int, double>::iterator findI;

      DeltaVectorH[atoi(result[0].c_str())]=atof(result[1].c_str());

    }

  }

  // total direct damage
  double totalDamage=0;

  // total firms with direct damage
  double damagedNum=0;

  // input damage
  for(auto itr = piniVectorH.begin(); itr != piniVectorH.end(); ++itr) {
    piniCapVectorH[itr->first]=piniVectorH[itr->first]*(1-DeltaVectorH[itr->first]);


    if(piniCapVectorH[itr->first]<piniVectorH[itr->first]){
      damagedNum++;

    // debug
      //      cerr << itr->first << " " << piniVectorH[itr->first] << " " << piniCapVectorH[itr->first] << " " << piniCapVectorH[itr->first]/piniVectorH[itr->first] << endl;

    }

    // stats
    totalDamage=totalDamage+piniVectorH[itr->first]*DeltaVectorH[itr->first];
  }

  // debug
  //  cerr << "totalDamage " << totalDamage << endl;

  if(stOfs){
    stOfs << "totalDamage " << totalDamage << endl;
    stOfs << "totalDamageFirmNum " << damagedNum << endl;
  }

  // debug
  //  printHash1LayerIntDouble("piniCapVectorH", &piniCapVectorH);

}

// adapt
void adapt(){

  // NB: there are several types of adaptation functions
  
  // clear adapt
  adaptHoL.clear();

  // for each firm
  // probability for adapatation is calculated by using Ai,j and O*i,j
  // choose the most redundant and better situation
  // probability is caluculated in adapProbFunc

  // debug
  //  cerr << "---------- adapt" << endl;

  // check all firms
  // (order should be random in the future implementation)

  unordered_map<int, int> *hash;
  hash=&firmH;
  for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

    int i=itr->first;

    if(iAdaptH[i]==0){
      continue;
    }

    // debug
    //    cerr << "---------- checking firm " << i << endl;

    // For all Ai,j
    
    // Find small O* 

    // (Ai,j-O*i,j)/Ai,j
    // ordered
    unordered_map<int,double> orderUnFullfillH;

    unordered_map<int, double> *hash2;      
    try{
      //      hash2=&(fOAdjustedFirmFirmHoH.at(i));
      hash2=&(fOAdjustedFirmFirmHoH[i]);
      }
    catch(const std::out_of_range& oor){

      // it is possible that there is no Oi,j
      // terminal node
      continue;
    }

    // Oi,j unullfill
    // hash key vector
    vector<int> keyL;

    // max
    int maxIndex=-1;
    double maxVolume=-1;

    for(auto itr2 = hash2->begin(); itr2 != hash2->end(); ++itr2) {
      
      int j=itr2->first;

      // debug
      //      cerr << "firm j: " << j << endl;

      keyL.push_back(j);

      // calculate Oi,j unfullfillment
      // (A - adjusted order) / A
      //      orderUnFullfillH[j]=( (fATableHoH[i])[j]-itr2->second )/ (fATableHoH[i])[j] ;

      // debug
      //      cerr << "(fATableHoH[" << i << "])[" << j << "] " << (fATableHoH[i])[j] << endl;

      // (order - adjusted order) / order
      orderUnFullfillH[j]=( (fOFirmFirmHoH[i])[j]-itr2->second )/ (fOFirmFirmHoH[i])[j] ;

      // debug
      //      cerr << "(fOFirmFirmHoH[" << i << "])[" << j << "] " << (fOFirmFirmHoH[i])[j] << endl;

      // debug
      //      cerr << "(fOAdjustedFirmFirmHoH[" << i << "])[" << j << "] " << itr2->second << endl;

      // debug
      //      cerr << "orderUnFullfillH[" << j << "] " << orderUnFullfillH[j] << endl;

      if(maxVolume<orderUnFullfillH[j]){
	maxVolume=orderUnFullfillH[j];
	maxIndex=j;

	// debug
	//	cerr << "maxIndex: " << maxIndex << endl;
	//	cerr << "maxVolume: " << maxVolume << endl;
      }

    }

    // debug
//      unordered_map<int,double> *dHash;
//      dHash=&(orderUnFullfillH);
//      for(auto itrd = dHash->begin(); itrd != dHash->end(); ++itrd) {
//        int j=itrd->first;
//        cerr << "i: " << i << " j: " << j << " unfullfill : " << itrd->second << endl;
//      }

    // debug
//     vector<int> *dv;
//     dv=&(keyL);
//     cerr << "keyL:" << endl;
//     for(auto itrd = dv->begin(); itrd != dv->end(); ++itrd) {
//       cerr << " " << *itrd ;
//     }
//     cerr << endl;
    
    // debug
    //    cerr << "shuffle" << endl;
    
    // shuffle key list
    //    random_shuffle(keyL.begin(), keyL.end());

    // debug
    //      vector<int> *dv;
    //      dv=&(keyL);
    //      cerr << "keyL:" << endl;
    //      for(auto itrd = dv->begin(); itrd != dv->end(); ++itrd) {
      
    //        cerr << " " << *itrd ;

    //      }
    //      cerr << endl;

    // use maximum unfulfillment
    // use only max
    // (however a lot of other algorithm is possible)

    if(maxIndex==-1){
      // does not exist

      // debug
      //      cerr << "no unfullfill candidate" << endl;

      continue;
    }else{
      // exist
      
    }

    // max=j
    int j=maxIndex;

    // adap prob
    double adapProb=adapProbFunc(orderUnFullfillH[j]);

    // debug
    //    cerr << "i (" << i << ") -> j (" << j << ") unfullfill: " << orderUnFullfillH[j] << " adap prob: " << adapProb << endl;

    // randm
    // use 7 digit
    double rNum=((double)(rand()%1000000))/1000000;

    // debug
    //    cerr << "rNum " << rNum << endl;
    
    // rand
    if(rNum>adapProb){
      // do nothing

      // debug
      //      cerr << "no swap by prob" << endl;

      continue;

    }else{

      // debug
      //      cerr << "swap activated" << endl;

    }

    // this i,j is the target of adaptation
    // get j's sector
    int jSec=firmAffiH[j];

    // debug
    //    cerr << "j(" << j << ")'s sector is " << jSec << endl;

    // get redundancy

    vector<int> *pv;
    pv=&(secFirmHoL[jSec]);

// 	int maxIndex=-1;
// 	double max=-1;

    // debug
    //    cerr << "redundant list: " << endl;

    // search in the sector
    int maxIndexK=-1;
    double maxVolumeK=-1;

    // hash2 is set to the link of i's target
    hash2=&(fOAdjustedFirmFirmHoH[i]);

    for(auto itr = pv->begin(); itr != pv->end(); ++itr) {

      int k =*itr;
      double volume;

      // discard unnecessary candidate

      // debug
      //      cerr << "k: " << k << endl;

      // i itself
      if(i==k){

	// debug
	//	cerr << "j=k" << endl;

	continue;
      }

      // the link i already has
      unordered_map<int, double>::iterator findI;
      findI = hash2->find(k);
      if(findI == hash2->end()){
	// no

	// debug
	//	cout << "k is not found in link of j" <<endl;


      }else{
	// yes

	// debug
	//	cout << "found k in link of j" << endl;

	// ignore it
	continue;

      }

      // redundancy（k's supplier, from i's viewpoint, it is supplier's supplier）
      
      // check order and adjusted order of k's supplier
      
      // calculate redundancy

      // consider current redundancy
      unordered_map<int, double> *hashOK;
      int iCatchF=0;
      try{
	hashOK=&(fOFirmFirmHoH[k]);
      }
      catch(const std::out_of_range& oor){
	// no supplier
	// terminal
	iCatchF=1;
      }

      int iLessOrderF=0;
      
      if(iCatchF==0){

	// there are k's suppliers
	for(auto itrOK = hashOK->begin(); itrOK != hashOK->end(); ++itrOK) {
	  int l=itrOK->first;
	  double order=itrOK->second;
	  double aOrder=(fOAdjustedFirmFirmHoH[k])[l];
	  if(order>aOrder){
	    iLessOrderF=1;
	    break;
	  }
	}

      }

      if(iLessOrderF){
	// k's adjustedOrder is less than order
	continue;
      }

      // possible productivity based on capacity and inventory
      volume=PiniMaxVectorH[k]-DAdjustedFirmVectorH[k];

      // debug
//       cerr << "PiniMaxVectorH[k] " << PiniMaxVectorH[k] << endl;
//       cerr << "DAdjustedFirmVectorH[k] " << DAdjustedFirmVectorH[k] << endl;
//       cerr << "redundant volume " << volume << endl;

      // debug
      //      cerr << "firm " << k << " volume " << volume << endl;
	    
      // 	  firmAndVolume data;
      // 	  data.setVolume(volume);
      // 	  data.setFirm(j);
      // 	  redOrderedL.push_back(data);

      if(volume>maxVolumeK){

	maxIndexK=k;
	maxVolumeK=volume;

	// debug
// 	cerr << "max is recorded" << endl;
// 	cerr << "maxIndexK " << k << endl;
// 	cerr << "maxVolumeK " << maxVolumeK << endl;

      }

    }

    if(maxIndexK!=-1){

    }else{
      // not found

      // debug
      //      cerr << "no candidate" << endl;

      continue;
    }

    // k is fixed
    int k=maxIndexK;
    double kvolume=maxVolumeK;

    // is the volume enough

    // debug
//     cerr << "kvolume " << kvolume << endl;
//     cerr << "orderUnFullfillH[j] " << orderUnFullfillH[j] << endl;
//     cerr << "(fATableHoH[i])[j] " << (fATableHoH[i])[j] << endl;

    double intrinsicOrderIJ;
    intrinsicOrderIJ=(fOFirmFirmHoH[i])[j];

    // debug
    //    cerr << "Oij's adjusted order is " << intrinsicOrderIJ << endl;

    // redundancy is more than order?

    // this equation follows original paper
//    if(kvolume<( orderUnFullfillH[j] * (fATableHoH[i])[j] )){

    // kvolumeとadjustedOrderIJの差も本来判断に入るはず

    // this is wrong
    // maxVolumeは
    //    if(kvolume <= adjustedOrderIJ){

    // Redundancy is greater than order (not adjusted order)

    if(kvolume <= intrinsicOrderIJ){

      // no
      // less than order

      // debug
      //      cerr << "substitute volume short" << endl;
      
      continue;

    }else{

      // debug
      //      cerr << "substitute volume fullfill" << endl;

    }

    // debug
//     cerr << "i: " << i << " j: " << j << " k: " << k << endl;
//     cerr << "redandant volume of k " << kvolume << endl;
//     cerr << "intrinsicOrder of i j " << intrinsicOrderIJ << endl;

    // swith j to k

    vector<int> adaptV;
    adaptV.push_back(j);
    adaptV.push_back(k);
    adaptHoL[i]=adaptV;

    // fATableHoH: move
    // bATableHoH: move
    // cVectorH: no move
    // piniVectorH: move（partially）
    // firmAffiH: no move
    // firmH: no move
    // fOFirmFirmHoH: move
    // fOAdjustedFirmFirmHoH: move
    // bOFirmFirmHoH: move
    // bOStartFirmFirmHoH: move
    // bOAdjustedFirmFirmHoH: move
    // DFirmVectorH: move（partially）
    // DAdjustedFirmVectorH: move（partially）
    // cAdjustedVectorH: no move
    // piniCapVectorH: move（and modification）
    // SFirmFirmHoH: move
    // SFirmSectorHoH: move
    // AFirmSectorHoH: move
    // PpropFirmSectorHoH: move
    // PiniMaxVectorH: move（and modification）
    // PactualVectorH: move（and modification）
    // DeltaVectorH: no move

    // fATableHoH: move
    double copyFAij;
    copyFAij=(fATableHoH[i])[j];
    (fATableHoH[i])[k]=copyFAij;
    (fATableHoH[i]).erase(j);

    // bATableHoH: move
    double copyBA;
    copyBA=(bATableHoH[j])[i];
    (bATableHoH[k])[i]=copyBA;
    (bATableHoH[j]).erase(i);

    // cVectorH: no move

    // piniVectorH: no move
    //piniVectorH[j]=piniVectorH[j]-copyFAij;
    //piniVectorH[k]=piniVectorH[k]+copyFAij;

    // firmAffiH: no move

    // firmH: no move

    // fOFirmFirmHoH: move
    // but not necessary
    double copyFOij;
    copyFOij=(fOFirmFirmHoH[i])[j];
    (fOFirmFirmHoH[i])[k]=copyFOij;
    (fOFirmFirmHoH[i]).erase(j);

    // fOAdjustedFirmFirmHoH: move
    double copyFOAij;
    copyFOAij=(fOAdjustedFirmFirmHoH[i])[j];
    (fOAdjustedFirmFirmHoH[i])[k]=copyFOAij;
    (fOAdjustedFirmFirmHoH[i]).erase(j);

    // bOFirmFirmHoH: move
    double copyBOji;
    copyBOji=(bOFirmFirmHoH[j])[i];
    (bOFirmFirmHoH[k])[i]=copyBOji;
    (bOFirmFirmHoH[j]).erase(i);

    // bOStartFirmFirmHoH: move
    double copyBOSji;
    copyBOSji=(bOStartFirmFirmHoH[j])[i];
    (bOStartFirmFirmHoH[k])[i]=copyBOSji;
    (bOStartFirmFirmHoH[j]).erase(i);

    // bOAdjustedFirmFirmHoH: move
    double copyBOAji;
    copyBOAji=(bOAdjustedFirmFirmHoH[j])[i];
    (bOAdjustedFirmFirmHoH[k])[i]=copyBOAji;
    (bOAdjustedFirmFirmHoH[j]).erase(i);

    // DFirmVectorH: no move
    // DAdjustedFirmVectorH: move（partially）
    // cAdjustedVectorH: no move

    // piniCapVectorH: move（and modification）
    //piniCapVectorH[j]=piniVectorH[j]*(1-DeltaVectorH[j]);
    //piniCapVectorH[k]=piniVectorH[k]*(1-DeltaVectorH[k]);

    // SFirmFirmHoH: move
    double copySij;
    copySij=(SFirmFirmHoH[i])[j];
    (SFirmFirmHoH[i]).erase(j); // no inventory for j
    (SFirmFirmHoH[i])[k]=copySij; // move inventory

    // SFirmSectorHoH: no move（but calculated）

    // AFirmSectorHoH: no move

    // PpropFirmSectorHoH: no move（but calculated）
    // PiniMaxVectorH: no move（but calculated）
    // PactualVectorH: no move（but calculated）
    // DeltaVectorH: no move

    // adaptation count minus one
    if(iAdaptH[i]==-1){

    }else{
      iAdaptH[i]--;
      iAdaptDoneH[i]++;
    }
    
  }

  // debug
  //  cerr << "adapt end" << endl;

}


void simulation(){
  // simulation start
  // all variables are already initialized

  for(int t=0;t<simulationStep;t++){

    // debug
    //    cout << "*************************************************************** time " << t << endl;
  
    // process 3 ==============================================================
    // renewal of inventory
    // SFirmFirmHoH

    // debug
    //    cerr << "inv renewal" << endl;

    // original paper uses j,i but here uses i,j
    unordered_map<int, unordered_map<int, double> > *hash;

    hash=&SFirmFirmHoH;

    for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

      // itr->first is i

      for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {

	// itr2->first is j
	// itr2->second is s(i,j)

	// if piniVector=0, renewal equation is different

	double newInventory;
	
	if(piniVectorH[itr->first]==0){
	  // avoid zero division
	  // no change
	  newInventory=itr2->second;
	}else{
	  newInventory=
	    itr2->second // S(i,j)
	    + (fOAdjustedFirmFirmHoH[itr->first])[itr2->first] // O*(i,j)
	    - (fATableHoH[itr->first])[itr2->first] * PactualVectorH[itr->first] / piniVectorH[itr->first] // A(i,j)*Pia(t-1)/Piini
	    ;
	}

	// debug
// 	if(isnan(newInventory)){
// 	  cerr << "itr->first " << itr->first << endl;
// 	  cerr << "itr2->first " << itr2->first << endl;
// 	  cerr << "itr2->second " << itr2->second << endl;
// 	  cerr << "(fOAdjustedFirmFirmHoH[itr->first])[itr2->first] " << (fOAdjustedFirmFirmHoH[itr->first])[itr2->first] << endl;
// 	  cerr << "(fATableHoH[itr->first])[itr2->first] " << (fATableHoH[itr->first])[itr2->first] << endl;
// 	  cerr << "PactualVectorH[itr->first] " << PactualVectorH[itr->first] << endl;
// 	  cerr << "piniVectorH[itr->first] " << piniVectorH[itr->first] << endl;
// 	  exit(1);
// 	}

	// debug

// 	  cerr << "itr->first " << itr->first << endl;
// 	  cerr << "itr2->first " << itr2->first << endl;
// 	  cerr << "itr2->second " << itr2->second << endl;
// 	  cerr << "(fOAdjustedFirmFirmHoH[itr->first])[itr2->first] " << (fOAdjustedFirmFirmHoH[itr->first])[itr2->first] << endl;
// 	  cerr << "(fATableHoH[itr->first])[itr2->first] " << (fATableHoH[itr->first])[itr2->first] << endl;
// 	  cerr << "PactualVectorH[itr->first] " << PactualVectorH[itr->first] << endl;
// 	  cerr << "piniVectorH[itr->first] " << piniVectorH[itr->first] << endl;

	// here minus inventory is possible
	// this is firm level inventory
	// and production is based on sector
	// therefore, at consumption, there is possible that inventory lower than zero
	if(newInventory<0){
	  newInventory=0;
	}

	itr2->second = newInventory;
      }
    }

    // debug
    //    printHash2LayerIntIntDouble("SFirmFirmHoH", &SFirmFirmHoH);


    // check adaptation

    // adapt
    if(iAdapType!=-1){
      adapt();
    }

    // process 4 ==============================================================
    // renewal of Order
    // fOFirmFirmHoH, bOFirmFirmHoH
    // Adjusted is later
    
    // NB: 'b' means backward and the opposite of 'f'
    
    // use hash again
    hash=&fOFirmFirmHoH;

    // debug
    //    cerr << "order renewal" << endl;

    for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

      // itr->first is i

      // debug
      //      cerr << "  i: " << itr->first << endl;

      int i=itr->first;

      for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {

	// itr2->first is j
	// itr2->second is o(i,j)
	int j=itr2->first;
      
	// use Dadjusted

	// debug
	//	cerr << "    j: " << itr2->first << endl;

	// debug
	//	cerr << "      pre: " << itr2->second << endl;

	double newOrder;

	// accept pini=0

	// debug
	//	cerr << "iOrderTypeFrag " << iOrderTypeFrag << endl;

	if(piniVectorH[itr->first]==0){

	  // debug
	  //	  cerr << "0 pass" << endl;
	  //	  cerr << piniVectorH[itr->first] << endl;

	  // order=0 is not accepted
	  // newOrder!=0

	  // order the amount based on fATable
	  newOrder=(fATableHoH[itr->first])[itr2->first];
	  // Si,j is constant
	  // consumption is always equalt to fATable

	}else{

	  if(iOrderTypeFrag==0){

	    // normal order
	    newOrder=
	      // original
	      (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] // 1st term
	      + 1 / tau
	      * ( inventoryN * (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] // 2nd term target inventory
		  - (SFirmFirmHoH[itr->first])[itr2->first]); // 3rd term

	  }else if(iOrderTypeFrag==1){

	    // inventory ini keep
	    newOrder=
	      (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] // 1st term
	      + 1 / tau
	      * ( inventoryN * (fATableHoH[itr->first])[itr2->first] // 2nd term （not depend on D. try to keep first inventory）
		  - (SFirmFirmHoH[itr->first])[itr2->first]); // 3rd term

	  }else if(iOrderTypeFrag==2){

	    // ignore negative inventory

	    if((inventoryN * (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] // 2nd term. target inventory
		- (SFirmFirmHoH[itr->first])[itr2->first])<=0) // 3rd term
	      {

		newOrder=
		  (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first]; // 1st term
	      }else{
	    
	      newOrder=
		// original
		(fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first]

		// 1st term
		
		+ 1 / tau
		* ( inventoryN * (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first]
		    // 2nd term. target inventory

		    - (SFirmFirmHoH[itr->first])[itr2->first]);
	      // 3rd term
	    }

	  }else{
	    exit(1);
	  }

	  // order boost verstion (not depend on demand）
	  /*
	  if((adaptHoL.count(i)!=0) && ((adaptHoL[i])[1]==j)){
	    // adaptation target: reset it

	    newOrder=
	      // inventory ini keep
	      (fATableHoH[itr->first])[itr2->first]  // 1st term D/Pini=1
	      + 1 / tau
	      * ( inventoryN * (fATableHoH[itr->first])[itr2->first] // 2nd term （not depend on D. keep initial inventory）
		  - (SFirmFirmHoH[itr->first])[itr2->first]); // 3rd term
	    
	  }else{

	    // nomarl
	    newOrder=


	      // original
//	      (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] // 1st term
//	      + 1 / tau
//	      * ( inventoryN * (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] // 2nd term
//	      - (SFirmFirmHoH[itr->first])[itr2->first]); // 3rd term

	      // inventory ini keep
	      (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] // 1st term
	      + 1 / tau
	      * ( inventoryN * (fATableHoH[itr->first])[itr2->first] // 2nd term (not depend on D. keep initial inventory)
		  - (SFirmFirmHoH[itr->first])[itr2->first]); // 3rd term

	  }
	  // order boost version end
	  
	  */

	}

	// if order is minus, set it to zero
	if(newOrder<0){
	  newOrder=0;
	}

	// renewal of forder
	itr2->second = newOrder;

	// debug
	//	cerr << "newOrder " << newOrder << endl;

	// debug
// 	if(isnan(newOrder)){
// 	  cerr << "nan " << itr2->first << " " << itr->first << endl;

//   	  cerr << "border i: " << itr2->first << " j: " << itr->first << " border " << (bOFirmFirmHoH[itr2->first])[itr->first] << " newOrder " << newOrder <<endl;
// // // 	  // print all
//   	  cerr << (fATableHoH[itr->first])[itr2->first] << endl;
//   	  cerr << DAdjustedFirmVectorH[itr->first] << endl;
//   	  cerr << piniVectorH[itr->first] << endl;
//   	  cerr << inventoryN / tau << endl;
//   	  cerr << (SFirmFirmHoH[itr->first])[itr2->first] << endl;
//  	  cerr << cVectorH[itr->first] << endl;

// 	  exit(1);
// 	}

	// debug
// 	if((bOFirmFirmHoH[itr2->first])[itr->first]>newOrder){
//  	  cerr << "border i: " << itr2->first << " j: " << itr->first << " border " << (bOFirmFirmHoH[itr2->first])[itr->first] << " newOrder " << newOrder <<endl;
// // 	  // print all
//  	  cerr << (fATableHoH[itr->first])[itr2->first] << endl;
//  	  cerr << DAdjustedFirmVectorH[itr->first] << endl;
//  	  cerr << piniVectorH[itr->first] << endl;
//  	  cerr << inventoryN / tau << endl;
//  	  cerr << (SFirmFirmHoH[itr->first])[itr2->first] << endl;
// 	  cerr << cVectorH[itr->first] << endl;
//  	}


	(bOFirmFirmHoH[itr2->first])[itr->first]=newOrder;
	

	// debug
// 	cerr << " 1st term:" << endl;
// 	cerr << (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] << endl;
// 	cerr << " 2nd term:" << endl;
// 	cerr << (fATableHoH[itr->first])[itr2->first] * DAdjustedFirmVectorH[itr->first] / piniVectorH[itr->first] * inventoryN / tau << endl;
// 	cerr << " 3rd term:" << endl;
// 	cerr << " - " << (SFirmFirmHoH[itr->first])[itr2->first] / tau << endl;

// 	// debug
// 	cerr << "      post: " << itr2->second << endl;

      }
    }

  // debug
    //  printHash2LayerIntIntDouble("fOFirmFirmHoH", &fOFirmFirmHoH);

  // debug
    //  printHash2LayerIntIntDouble("bOFirmFirmHoH", &bOFirmFirmHoH);

    // process 5 ==============================================================
    // renewal of Demand
    // renewal of DFirmVectorH

    // debug
    //    cerr << "demand renewal" << endl;

    // demand is reset by cvector out of consideration of terminal nodes
    // do not reset by zero
    // use hash
    unordered_map<int, int> *hashII;
    hashII=&firmH;
    
    for(auto itr = hashII->begin(); itr != hashII->end(); ++itr) {
      DFirmVectorH[itr->first]=cVectorH[itr->first];
    }
    
    // use hash
    hash=&bOFirmFirmHoH;

    // debug
    //    double diff=0;

    for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

      // debug
//       if((itr->second).size()==0){
// 	cerr << "i: " << itr->first << " doesnot have client but border recorded" << endl;
//       }

      // debug
      //      cerr << "i: " << itr->first << endl;

      // fOrder: client -> supplier
      // bOrder: supplier -> client
      // demand is order toward the firm
      // use bOrder and total client order
      double totalDemand=0;

      for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {

	// debug
	//	cerr << "j: " << itr2->first << endl;
	//	cerr << "volume: " << itr2->second << endl;

	totalDemand=totalDemand+itr2->second;
      }

      // debug
      //      cerr << totalDemand << endl;

      // add C
      // Order is changed backward, forward
      // but C is never affected by disaster
      totalDemand=totalDemand+cVectorH[itr->first];

      // debug
      //      cerr << totalDemand << endl;

      // debug
//       if(DFirmVectorH[itr->first]<totalDemand){
// 	diff = diff + (totalDemand-DFirmVectorH[itr->first]);
//       }

      DFirmVectorH[itr->first]=totalDemand;

      // debug
      // 	cerr << "i: " << itr->first << " D " << DFirmVectorH[itr->first] << " totalDemand " << totalDemand << endl;
      //      if(isnan(totalDemand)){
      //	exit(1);
      //      }

      // debug
      //            cerr << "i " << itr->first << endl;
      //            cerr << "DFirmVector[i] " << totalDemand << endl;

    }

    // debug
    //    cerr << "diff " << diff << endl;

    // process 6 ==============================================================
    // calculate of Pactual

    // Pcap is decided

    // first renew Stotal i,s

    // process 7 is included because it is simple

    // debug
    //    cerr << "pactual renewal" << endl;

    // clear sector inventory with 0
    // use hash
    hash=&SFirmSectorHoH;

    for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

      // itr->first is i
      int i=itr->first;

      for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {

	// itr2->first is sj
	//	int sj=itr2->first;
	//	(SFirmSectorHoH[i])[sj]=0;

	itr2->second=0;

      }
    }
    
    // use hash
    hash=&SFirmFirmHoH;

    // Stotal i,s does not require reset because it is set by Atotal i,s

    // create setor inventory
    // and create PpropFirmSectorHoH(PsProp)

    for(auto itr = hash->begin(); itr != hash->end(); ++itr) {

      // itr->first is i
      int i=itr->first;

      for(auto itr2 = (itr->second).begin(); itr2 != (itr->second).end(); ++itr2) {

	// itr2->first is j
	int j=itr2->first;

	// itr2->second is S(i,sj)
	double s=itr2->second;

	// get sector
	int sj=firmAffiH[j];

	(SFirmSectorHoH[i])[sj]=(SFirmSectorHoH[i])[sj]+s;

      }

      // prop can be calculated because i is know at this moment
      unordered_map<int, double> *hash2;

      // get i's SFirmSectorHoH
      hash2=&(SFirmSectorHoH[i]);

      // debug
      //      cerr << "SFirmSectorHoH[" << i << "]" << endl;
      //      printHash1LayerIntDouble("SFirmSectorHoH[i]", hash2);

      // for all sector
      for(auto itr2 = hash2->begin(); itr2 != hash2->end(); ++itr2) {
	int sj=itr2->first;

	if((AFirmSectorHoH[i])[sj]==0){
	  cerr << "0 div (AFirmSectorHoH[i])[sj]: " << (AFirmSectorHoH[i])[sj] << endl;
	  exit(1);
	}

	// debug
	//	cerr << "(SFirmSectorHoH[i])[sj] " << (SFirmSectorHoH[i])[sj] << endl;
	//	cerr << "(AFirmSectorHoH[i])[sj] " << (AFirmSectorHoH[i])[sj] << endl;
	//	cerr << "piniVectorH[i] " << piniVectorH[i] << endl;

	(PpropFirmSectorHoH[i])[sj]=(SFirmSectorHoH[i])[sj]/(AFirmSectorHoH[i])[sj]*piniVectorH[i];

	// debug
	//	cerr << "i: " << i << " sj: " << sj << " PpropFirmSectorHoH[i])[sj]: " << (PpropFirmSectorHoH[i])[sj] << endl;

      }
      
    }
    // PpropFirmSectorHoH is calculated at this moment

    // debug
    //    printHash2LayerIntIntDouble("PropFirmSectorHoH", &PpropFirmSectorHoH);
    

    // create Pini Max Vector
    // use hash
    for(auto itr = firmH.begin(); itr != firmH.end(); ++itr) {

      // debug
      //            cerr << "i " << itr->first << endl;

      int i=itr->first;

      unordered_map<int, double> *hash2;

      // get i's SFirmSectorHoH
      // if there is none of them, that is beginning node
      try{
	hash2=&(PpropFirmSectorHoH[i]);
      }
      catch(const std::out_of_range& oor){
	// beginning node
	PiniMaxVectorH[i]=piniCapVectorH[i];

	// debug
	//	cerr << "caught" << endl;
	
	continue;
      }

      // hash is possibly empty
      if(hash2->empty()){

	// debug
	//	cerr << "empty" << endl;

	PiniMaxVectorH[i]=piniCapVectorH[i];
	continue;

      }

      // debug
      //      cerr << "no caught" << endl;
      //      printHash1LayerIntDouble("propfirmsectorhoh[i]", hash2);
      

      // for all i's sector
      double minProp;
      int first=1;
      for(auto itr2 = hash2->begin(); itr2 != hash2->end(); ++itr2) {
	if(first==1){
	  minProp=itr2->second;
	  first=0;
	}else{
	  if(minProp>itr2->second){
	    minProp=itr2->second;
	  }
	}
      }

      // at this moment, minProp is the mimimum of Pprop
      
      if(piniCapVectorH[i]>minProp){
	PiniMaxVectorH[i]=minProp;
      }else{
	PiniMaxVectorH[i]=piniCapVectorH[i];
      }
    }

    // create PactualVectorH
    // demand is the newest D

    for(auto itr = firmH.begin(); itr != firmH.end(); ++itr) {

      int i = itr->first;

      // debug
      //      cerr << "i " << i << " PiniMaxVectorH[i] " << PiniMaxVectorH[i] << " DFirmVectorH[i] " << DFirmVectorH[i] << endl;
      //            cerr << "old PactualVectorH[i] " << PactualVectorH[i] << endl;

      // compare D and set Actual
      if(PiniMaxVectorH[i]>DFirmVectorH[i]){
	PactualVectorH[i]=DFirmVectorH[i];
      }else{
	PactualVectorH[i]=PiniMaxVectorH[i];
      }

      // debug
      //            cerr << "new PactualVectorH[i] " << PactualVectorH[i] << endl;

      // set Dadjusted
      DAdjustedFirmVectorH[i]=PactualVectorH[i];

      // debug
       // if(PactualVectorH[i] < (piniVectorH[i]/2)){
       // 	 cerr << "firm " << i << "==========" << endl;
       // 	 printHash1LayerIntDouble("Inventory",&SFirmFirmHoH[i]);
       // 	 printHash1LayerIntDouble("Order",&fOFirmFirmHoH[i]);
       // 	 printHash1LayerIntDouble("BOrder(Demand)",&bOFirmFirmHoH[i]);
       // 	 printHash1LayerIntDouble("Adjusted Order",&fOAdjustedFirmFirmHoH[i]);
       // 	 cerr << "D " << DFirmVectorH[i] << endl;
       // 	 cerr << "DAdjusted " << DAdjustedFirmVectorH[i] << endl;
       // 	 cerr << "pini " << piniVectorH[i] << endl;
       // 	 cerr << "pact " << PactualVectorH[i] << endl;

       // 	 // int j=350204446;
       // 	 // cerr << "firm " << j << "==========" << endl;
       // 	 // printHash1LayerIntDouble("Inventory",&SFirmFirmHoH[j]);
       // 	 // printHash1LayerIntDouble("Order",&fOFirmFirmHoH[j]);
       // 	 // printHash1LayerIntDouble("BOrder(Demand)",&bOFirmFirmHoH[j]);
       // 	 // cerr << "D " << DFirmVectorH[j] << endl;
       // 	 // cerr << "pini " << piniVectorH[j] << endl;
       // 	 // cerr << "pact " << PactualVectorH[j] << endl;
      //       }


      // rationing begins from here
       
      // get OrderRatio
      // mostly, that is 1
      double adjustRatio;

      // OrderとFCのrationing

      if(iRationingTypeFrag==0){
      
      // this is original paper's rationing

      // 0 possibly exists
      if(DFirmVectorH[i]==0){

	// debug
	//	cerr << "DFirmVectorH[i]==0 0 div" << endl;
	//	cerr << "i :" << i << endl;
	//	exit(1);

	adjustRatio=0;
      }else{
	adjustRatio=PactualVectorH[i]/DFirmVectorH[i];
      }

      // renew C* first
      // cadjusted
      cAdjustedVectorH[i]=adjustRatio*cVectorH[i];

      // get i's fOFirmFirmHoH and bOFirmFirmHoH
      unordered_map<int, double> *hash2;      
      hash2=&(bOFirmFirmHoH[i]);

      // for all i's supplier
      for(auto itr2 = hash2->begin(); itr2 != hash2->end(); ++itr2) {

	int j=itr2->first;

	// adjusted f
	(bOAdjustedFirmFirmHoH[i])[j]=itr2->second*adjustRatio;

	// adjusted b
	(fOAdjustedFirmFirmHoH[j])[i]=itr2->second*adjustRatio;

      }

      // original rationing ends here
      
      }else if(iRationingTypeFrag==1){

	// rationing but FC is always considered last

	double shortD=DFirmVectorH[i]-PactualVectorH[i];
	double FC=cVectorH[i];
	if(shortD<0){
	  shortD=0;
	}
      
	if(FC>=shortD){
	  // FC is enough to fulfill

	  // cadjusted
	  cAdjustedVectorH[i]=FC-shortD;

	  adjustRatio=1;

	}else{
	  // FC is no enough
	  // rationing occurs

	  // cadjusted
	  cAdjustedVectorH[i]=0;

	  // total amount of order
	  double order=DFirmVectorH[i]-FC;

	  // there possibly exists 0
	  if(DFirmVectorH[i]==0){

	    // debug
	    //	cerr << "DFirmVectorH[i]==0 0 div" << endl;
	    //	cerr << "i :" << i << endl;
	    //	exit(1);

	    adjustRatio=0;
	  }else{
	    adjustRatio=PactualVectorH[i]/order;
	  }
	}

	// get i's fOFirmFirmHoH and bOFirmFirmHoH
	unordered_map<int, double> *hash2;      
	hash2=&(bOFirmFirmHoH[i]);

	// for all i's supplier
	for(auto itr2 = hash2->begin(); itr2 != hash2->end(); ++itr2) {

	  int j=itr2->first;

	  // adjusted f
	  (bOAdjustedFirmFirmHoH[i])[j]=itr2->second*adjustRatio;

	  // adjusted b
	  (fOAdjustedFirmFirmHoH[j])[i]=itr2->second*adjustRatio;

	}

	// rationing: fc is last, ends here

      }else if(iRationingTypeFrag==2){

	// soothing algorithm
	// inoue & todo 17 uses this algorithm
	// FC's ratio=1

	// get i's bOFirmFirmHoH
	unordered_map<int, double> *bhash;      
	bhash=&(bOFirmFirmHoH[i]);

	//	bhash->operator[](-1)=cVectorH[i];

	// get i's bOStartFirmFirmHoH
	unordered_map<int, double> *bshash;      
	bshash=&(bOStartFirmFirmHoH[i]);

	//	bshash->operator[](-1)=cVectorH[i];

	// bOrder ratio
	unordered_map<int, double> bRatH;

	// bAdjust ratio
	unordered_map<int, double> bAdjustH;

	for(auto itrb = bhash->begin(); itrb != bhash->end(); ++itrb) {

	  unordered_map<int, double>::iterator findI=bshash->find(itrb->first);

	  // debug
	  // if(findI->second<=0){
	  //   cerr << "back order <=0" << endl;
	  //   exit;
	  // }
	  
	  bRatH[itrb->first]=itrb->second/findI->second;
	  bAdjustH[itrb->first]=0;
	}

	double ratFC=1;
	double adjustFC=0;

	// rationing

	// remaining act
	double remainingAct=PactualVectorH[i];

	// debug
	//	cerr << "firm " << i << " rationing" << endl;
	//	cerr << "PactualVectorH[i] " << remainingAct << endl;

	// debug
	// demand ratio
	//	cerr << "ratFC: " << ratFC << endl;
	//	printHash1LayerIntDouble("bRatH", &bRatH);

	// debug
	//	int counter=0;
	//	int bOnum=bshash->size();

	// bRatH, find minimum bOrder's ratio

	while(1){

	  // debug
	  //	  cerr << "remainingAct: " << remainingAct << endl;
	  //	  counter++;
// 	  if((bOnum+2)<counter){
// 	    cerr << "counter over" << endl;
// 	    cerr << "backorder num " << bOnum << endl;
// 	    cerr << "counter " << counter << endl;
// 	    cerr << "i " << i << endl;
// 	    exit(1);
// 	  }

//	  int minI; 
	  double minD=-1;

	  for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {

	    if(
#ifdef INTLC
	       isEqual(itrb->second,(double)0)
#else
	       almost_equal(itrb->second,(double)0,2)
#endif
	       ){
	      
	      // skip
	      // it is 0
	    }else if(minD<0){

	      //	      minI=itrb->first;
	      minD=itrb->second;
		
	    }else if(minD>itrb->second){
	      //	      minI=itrb->first;
	      minD=itrb->second;
	    }
          }

	  if(minD==-1){

	    if(
#ifdef INTLC
	       isEqual(ratFC,(double)0)
#else
	       almost_equal(ratFC,(double)0,2)
#endif
	       ){

	      break;

	    }

	    //	    minI=-1;
	    minD=ratFC;

	  }

	  // if demand is fulfilled, end
	  // consider 10 digit for ratio base=0
	  if(minD<(double)0.00000000001){
	    break;
	  }

	  // consider minimum

	  // how much is required for the fulfillment of minimum ratio
	  double totalMin=0;
	  for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
	    unordered_map<int, double>::iterator findI = bshash->find(itrb->first);

	    totalMin=totalMin+minD*findI->second;
	  }
	  // fc
	  totalMin=totalMin+minD*cVectorH[i];

	  // debug
	  //	  cerr << "totalMin " << totalMin << endl;
	  //	  cerr << "remainingAct: " << remainingAct << endl;

	  // is it possible to fulfill the minimum
	  if(remainingAct>totalMin){
	    // possible

	    // debug
	    //	    cerr << "remainingAct > totalMin" << endl;

	    // add minD to the achieved
	    // and delete from bRatH
	    for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
	      if(itrb->second-minD>=0){
		itrb->second=itrb->second-minD;
		bAdjustH[itrb->first]=bAdjustH[itrb->first]+minD;
	      }
	    }

	    // fc
	    ratFC=ratFC-minD;
	    adjustFC=adjustFC+minD;

	    remainingAct=remainingAct-totalMin;




	  }else{
	    //  not possible

	    // debug
	    //	    cerr << "remainingAct <= totalMin" << endl;

	    // find the ratio for all remaining demand that can be fulfilled

	    double vTotal=0;
	    for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
	      if(itrb->second>0){
		// bshash
		unordered_map<int, double>::iterator findI = bshash->find(itrb->first);
		vTotal=vTotal+findI->second;
	      }
	    }

	    // fc
	    if(ratFC>0){
	      vTotal=vTotal+cVectorH[i];
	    }

	    double ratio;

	    if(vTotal<=0){

	      // debug
	      //cerr << "vTotal <=0" << endl;
	      //exit(1);

	      ratio=0;

	    }else{
	      ratio=remainingAct/vTotal;
	    }

	    // add to fulfilled ratio
	    for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
	      if(itrb->second>0){
		bAdjustH[itrb->first]=bAdjustH[itrb->first]+ratio;
	      }
	    }

	    if(ratFC>0){
	      adjustFC=adjustFC+ratio;
	    }

	    break;
	  }
	}

	// debug
	// fulfilled ratio
	//	cerr << "fc adjust " << adjustFC << endl;
	//	printHash1LayerIntDouble("bAdjustH", &bAdjustH);

	// debug
	//	double realizeVolume=0;

	// activate fulfilled ratio

	// fc
	cAdjustedVectorH[i]=adjustFC*cVectorH[i];

	// debug
	//	realizeVolume=realizeVolume+cAdjustedVectorH[i];
	//	cerr << "realized FC " << cAdjustedVectorH[i] << endl;
	
	for(auto itrb = bAdjustH.begin(); itrb != bAdjustH.end(); ++itrb) {

	  int j=itrb->first;
	  
	  // inter firm
	    
	  // adjusted f
	  (bOAdjustedFirmFirmHoH[i])[j]=itrb->second*(bOStartFirmFirmHoH[i])[j];

	  // adjusted b
	  (fOAdjustedFirmFirmHoH[j])[i]=itrb->second*(bOStartFirmFirmHoH[i])[j];

	  // debug
	  //	  realizeVolume=realizeVolume+itrb->second*(bOStartFirmFirmHoH[i])[j];
	  //	  cerr << "realized volume to " << j << " : " << itrb->second*(bOStartFirmFirmHoH[i])[j] << endl;

	}

	// debug
	//	cerr << "realizeVolume " << realizeVolume << endl;
	
	// soothing algorithm ends here

      }else if(iRationingTypeFrag==3){
	
	// soothing algorithm and considers FC last

	// get i's fOFirmFirmHoH and bOFirmFirmHoH
	unordered_map<int, double> *bhash;      
	bhash=&(bOFirmFirmHoH[i]);

	// get i's bOStartFirmFirmHoH
	unordered_map<int, double> *bshash;      
	bshash=&(bOStartFirmFirmHoH[i]);

	// bOrder
	unordered_map<int, double> bRatH;

	// actual fulfill ratio
	unordered_map<int, double> bAdjustH;

	for(auto itrb = bhash->begin(); itrb != bhash->end(); ++itrb) {

	  unordered_map<int, double>::iterator findI=bshash->find(itrb->first);
	  
	  bRatH[itrb->first]=itrb->second/findI->second;
	  bAdjustH[itrb->first]=0;
	}	  

	
	double shortD=DFirmVectorH[i]-PactualVectorH[i];
	double FC=cVectorH[i];
	if(shortD<0){
	  shortD=0;
	}

	// FC is enough to fulfill demand
	if(FC>=shortD){
	  // enough

	  // cadjusted
	  cAdjustedVectorH[i]=FC-shortD;

	  // activate
	  bAdjustH=bRatH;

	}else{
	  // not enough
	  // rationing

	  // cadjusted
	  cAdjustedVectorH[i]=0;

	  // total amount of order
	  //double order=DFirmVectorH[i]-FC;

	  // remaining act
	  double remainingAct=PactualVectorH[i];

	  // debug
	  //	   int counter=0;
	  //	    int bOnum=bshash->size();

	  while(1){

	    // debug
// 	    counter++;
// 	     if(bOnum<counter){
// 	       cerr << "counter over" << endl;
// 	       cerr << "backorder num " << bOnum << endl;
// 	       cerr << "counter " << counter << endl;
// 	       cerr << "i " << i << endl;
// 	     }

//	    int minI; 
	    double minD=-1;

	    //debug
// 	    if(i==842021116){
// 	      	    cerr << "remainingAct " << remainingAct << endl;
// 		    printHash1LayerIntDouble("bRatH", &bRatH);
// 	    }

	    for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {

#ifdef INTLC
	      if(isEqual(itrb->second,(double)0)){
#else
	      if(almost_equal(itrb->second,(double)0,2)){
#endif
		// skip
		// it is 0
	      }else if(minD<0){

		//		minI=itrb->first;
		minD=itrb->second;
		
	      }else if(minD>itrb->second){
		//	minI=itrb->first;
		minD=itrb->second;
	      }
	    }

	    if(minD==-1){
	      break;
	    }

	    // minimum is obtained

	    // how much is necessary to fulfill the minimum ratio
	    double totalMin=0;
	    for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
	      unordered_map<int, double>::iterator findI = bshash->find(itrb->first);

	      totalMin=totalMin+minD*findI->second;
	    }

	    // is it possible to fulfill the minimum
	    if(remainingAct>totalMin){
	      // possible

	      // add minD as activation
	      // remove from bRatH
	      for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
		if(itrb->second-minD>=0){
		  itrb->second=itrb->second-minD;
		  bAdjustH[itrb->first]=bAdjustH[itrb->first]+minD;
		}
	      }

	      remainingAct=remainingAct-totalMin;
	      
	    }else{
	      // no possible

	      // find the ratio for all remaining demand that can be fulfilled

	      double vTotal=0;
	      for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
		if(itrb->second>0){
		// bshash
		  unordered_map<int, double>::iterator findI = bshash->find(itrb->first);
		  vTotal=vTotal+findI->second;
		}
	      }

	      double ratio=remainingAct/vTotal;

	      // add to fulfilled ratio
	      for(auto itrb = bRatH.begin(); itrb != bRatH.end(); ++itrb) {
		if(itrb->second>0){
		  bAdjustH[itrb->first]=bAdjustH[itrb->first]+ratio;
		}
	      }

	      break;
	    }
	      
	  }

	}

	// debug
	// fulfilled ratio
	//	printHash1LayerIntDouble("bAdjustH", &bAdjustH);

	// activate fulfilled ratio
	for(auto itrb = bAdjustH.begin(); itrb != bAdjustH.end(); ++itrb) {
	  int j=itrb->first;

	  // adjusted f
	  (bOAdjustedFirmFirmHoH[i])[j]=itrb->second*(bOStartFirmFirmHoH[i])[j];

	  // adjusted b
	  (fOAdjustedFirmFirmHoH[j])[i]=itrb->second*(bOStartFirmFirmHoH[i])[j];
	}
	
	// soothing algorithm with FC last


	}else{
	exit(1);
      }

    }

    // debug
    //    printAllInternal();
    // for(int i=1;i<=3;i++){
    //   cerr << "firm " << i << "==========" << endl;
    //   printHash1LayerIntDouble("Inventory",&SFirmFirmHoH[i]);
    //   printHash1LayerIntDouble("Order",&fOFirmFirmHoH[i]);
    //   printHash1LayerIntDouble("Adjusted Order",&fOAdjustedFirmFirmHoH[i]);
    //   printHash1LayerIntDouble("BOrder(Demand)",&bOFirmFirmHoH[i]);
    //   printHash1LayerIntDouble("Adjusted BOrder(Demand)",&bOAdjustedFirmFirmHoH[i]);
    //   cerr << "D " << DFirmVectorH[i] << endl;
    //   cerr << "DAdjusted " << DAdjustedFirmVectorH[i] << endl;
    //   cerr << "pini " << piniVectorH[i] << endl;
    //   cerr << "pact " << PactualVectorH[i] << endl;
    // }

    // output 
    // for every turn
//     unordered_map<int, double> *hashO;
//     hashO=&PactualVectorH;
//     map<int, int> sortedVectorH(hashO->begin(), hashO->end());
//     for(auto itr = sortedVectorH.begin(); itr != sortedVectorH.end(); ++itr) {
//       // show keys values
//       // i, actual P
//       ofs << t << " " << itr->first << " " << itr->second << endl;
//     }

    // output pactual
    // for every turn
    // double totalPactual=0;
    // for(auto itr = PactualVectorH.begin(); itr != PactualVectorH.end(); ++itr) {
    //   totalPactual=totalPactual+itr->second;
    // }
    // lastTotalPactual=totalPactual;
    // ofs << t << " " << totalPactual << endl;

    // output value added
    double totalValueAdded=0;
    double totalPactual=0;
    for(auto itr = PactualVectorH.begin(); itr != PactualVectorH.end(); ++itr) {
      totalPactual=totalPactual+itr->second;
      totalValueAdded=totalValueAdded+valueAddedVectorH[itr->first]*itr->second;
    }
    lastTotalPactual=totalPactual;
    lastTotalValueAdded=totalValueAdded;

    // value added
    ofs << t << " " << totalPactual << endl;
    pvfs << t << " " << totalValueAdded << endl;

    // output final consumption
    if(fcOfs){
      double totalFC=0;
      for(auto itr = cAdjustedVectorH.begin(); itr != cAdjustedVectorH.end(); ++itr) {
	totalFC=totalFC+itr->second;
      }
      fcOfs << t << " " << totalFC << endl;
    }

    // debug
    if(iFullDebugFrag){
      fullDebug(t);
    }

    // output snapshot
    unordered_map<int, int>::iterator findI;
    findI = hPrintSnapStep.find(t);


    if(findI == hPrintSnapStep.end()){
      // not found
    }else{

      // found
      string sFile=sPrintSnapFileBase;

#if defined(INTLC) || defined(KEI)
      char buf[256];
      sprintf(buf, "%d", t);
      sFile=sFile+buf;
#else
      sFile=sFile+to_string(t);
#endif

      ofstream snapOfs;

      snapOfs.open(sFile.c_str());

      if(!snapOfs){
	cerr << "Unable to open " << sFile << " as snap output file at step " << t << endl;
	exit(1);
      }


      unordered_map<int, double> *hashO;
      hashO=&PactualVectorH;
      map<int, double> sortedVectorH(hashO->begin(), hashO->end());
      
      for(auto itr = sortedVectorH.begin(); itr != sortedVectorH.end(); ++itr) {
	snapOfs << itr->first << " " << itr->second << endl;;
      }

    }


  } // for t, simulation step 

  // output stat
  if(stOfs){
        stOfs << "totalPactutalEnd " << lastTotalPactual << endl;
        stOfs << "totalValueAddedEnd " << lastTotalValueAdded << endl;
  }

  // output adaptation
  if(adapOfs){
    for(auto itr =iAdaptDoneH.begin(); itr != iAdaptDoneH.end(); ++itr) {
      adapOfs << itr->first << " " << itr->second << endl;
    }
  }
  
}

vector<string> split(const string &s, char delim) {
  vector<string> elems;
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    if (!item.empty()) {
       elems.push_back(item);
    }
  }
  return elems;
}

int main(int argc, char *argv[]){

  // rand seed
  int randSeed=331;

  // debug
  //  cout << relFileS << endl;

  time_t btime;
  time_t etime;

  // start time
  time(&btime);

  // options
  int		ch;
  extern char	*optarg;
  extern int	optind, opterr;

  char* proName=argv[0];
  // minimum input files
  int argNum=5; 

  while ((ch = getopt(argc, argv, "r:t:d:f:s:o:p:a:FR:v:A:")) != -1){
    switch (ch){
    case 'r':
      randSeed=atoi(optarg);
      break;
    case 't':
      simulationStep=atoi(optarg);
      break;
    case 'd':
      deltaFile=optarg;
      break;
    case 'f':
      finalConsumeOutputFile=optarg;
      break;
    case 'o':
      iOrderTypeFrag=atoi(optarg);
      break;
    case 's':
     statsOutputFile=optarg;
      break;
    case 'a':
      iAdapType=atoi(optarg);
      break;
    case 'F':
      iFullDebugFrag=1;
      break;
    case 'R':
      iRationingTypeFrag=atoi(optarg);
      break;
    case 'p':
      iPrintSnapFrag=1;
      sPrintSnapStep=optarg;
      break;
    case 'v':
      pValueAddedOutputFile=optarg;
      break;
    case 'A':
      adaptationOutputFile=optarg;
      break;
    default:
      usage(proName);
    }
  }
  argc -= optind;
  argv += optind;
  
  // initialize random

  // debug
  //  cerr << "randSeed " << randSeed << endl;

  srand(randSeed);

  // debug
  //  cerr << "argc " << argc << endl;

  // argc
  if(argc!=argNum){
    usage(proName);
    return 1;
  }

  // ioLineFile
  char *ioLineFile=argv[argc-5];
  ifstream ifs(ioLineFile);

  if(!ifs){
    if(argc>=argNum){
      cerr << "Unable to open " << ioLineFile << " as aTableLineFile" << endl;
    }else{
      usage(proName);
    }
    return 1;
  }

  // debug
  //  cerr << "ioLineFile " << ioLineFile << endl;

  // cVectorFile

  char *cVectorFile=argv[argc-4];
  ifstream ifs2(cVectorFile);
  if(!ifs2){
    if(argc>=argNum){
      cerr << "Unable to open " << cVectorFile << " as cVectorFile" << endl;
    }else{
      usage(proName);
    }
    return 1;
  }

  // debug
  //  cerr << "cVectorFile " << cVectorFile << endl;

  // piniVectorFile

  char *piniVectorFile=argv[argc-3];
  ifstream ifs3(piniVectorFile);
  if(!ifs3){
    if(argc>=argNum){
      cerr << "Unable to open " << piniVectorFile << " as piniVectorFile" << endl;
    }else{
      usage(proName);
    }
    return 1;
  }

  // debug
  //  cerr << "piniVectorFile " << piniVectorFile << endl;

  // firmAffiFile

  char *firmAffiFile=argv[argc-2];
  ifstream ifs4(firmAffiFile);
  if(!ifs4){
    if(argc>=argNum){
      cerr << "Unable to open " << firmAffiFile << " as firmAffiFile" << endl;
    }else{
      usage(proName);
    }
    return 1;
  }

  //deltafile
  if(deltaFile.size()>0){
    ifsDelta.open(deltaFile.c_str());
    if(!ifsDelta){
      if(argc>=argNum){
	cerr << "Unable to open " << deltaFile << " as deltaFile" << endl;
      }else{
	usage(proName);
      }
      return 1;
    }
  }

  // pactOutputFile

  char *pactOutputFile=argv[argc-1];
  ofs.open(pactOutputFile);
  if(!ofs){

    // debug
    cerr << "argc: " << argc << endl;
    cerr << "argNum: " << argNum << endl;
    
    if(argc>=argNum){
      cerr << "Unable to open " << pactOutputFile << " as pactOutputFile" << endl;
    }else{
      usage(proName);
    }
    return 1;
  }

  // pvalueAddedFile
  if(pValueAddedOutputFile.size()>0){
    pvfs.open(pValueAddedOutputFile.c_str());
    if(!pvfs){
      if(argc>=argNum){
	cerr << "Unable to open " << pValueAddedOutputFile << " as pValueAddedFile" << endl;
      }else{
	usage(proName);
      }
      return 1;
    }
  }

  // adaptationOutputFile
  if(adaptationOutputFile.size()>0){
    adapOfs.open(adaptationOutputFile.c_str());
    if(!adapOfs){
      if(argc>=argNum){
	cerr << "Unable to open " << adaptationOutputFile << " as adaptationOutputFile" << endl;
      }else{
	usage(proName);
      }
      return 1;
    }
  }

  // get snapshot info
  if(iPrintSnapFrag){

    vector<string> vS;
    vS=split(sPrintSnapStep, ':');

    if(vS.size()<=1){
      cerr << "Invalid arg for option p (should be file:step:...)" << endl;
      usage(proName);
      return(1);
    }

    sPrintSnapFileBase=vS[0];

    auto itr = vS.begin();
    itr++;
    for(; itr != vS.end(); ++itr) {
      hPrintSnapStep[atoi((*itr).c_str())]=1;
    }

  }

  // debug
  //  cout << "sPrintSnapFileBase is " << sPrintSnapFileBase << endl;
  //  printHash1LayerIntInt("hPrintSnapStep", &hPrintSnapStep);
  //  return(1);

  // debug
  // output firms with damage
  // count the number
//   damageFS.open(damagedFirmFile);
//   if(!damageFS){
//     cerr << "Unable to open" << endl;
//     return 1;
//   }

// debug
// fulldebug

  if(iFullDebugFrag){
    debugFS.open(fullDebugFile.c_str());
    if(!debugFS){
      cerr << "Unable to open" << endl;
      return 1;
    }
    debugFS << "t, i, pini, delta, pcap, d, pmax, pact, dadjusted, s..., o..., " << endl;
  }

  // debug
  //  cerr << "ioLineFile " << ioLineFile << endl;

  // read ATableLine
  // column: 0, 1, 2 means
  // 0: supplier, 1: client, 2: volume
  // however, hallegatte paper says
  // row client - product -> column supplierなので
  // forward should be initialized by 1, 0, 2
  // backward should be initialized by 0, 1, 2

  string str;
  while(getline(ifs,str)){
    string tmp;
    istringstream stream(str);
    vector<string> result;
    while(getline(stream,tmp,' ')){

      // debug
      //	  cout<< tmp << endl;

      result.push_back(tmp);
    }

    if(result[2].c_str()<=0){
      continue;
    }

    // does source exists in hash?
    unordered_map<int, unordered_map<int, double> >::iterator findI;

    findI = fATableHoH.find (atoi(result[1].c_str()));
    if(findI == fATableHoH.end()){
      // no

      // debug
      //      cout << "not found" <<endl;

      unordered_map<int, double> fATableH;
      fATableHoH[atoi(result[1].c_str())]=fATableH;

    }else{
      // yes

      // debug
      //     cout << "found" << endl;

    }

    (fATableHoH[atoi(result[1].c_str())])[atoi(result[0].c_str())]=atof(result[2].c_str());

    // opposite
    findI = bATableHoH.find (atoi(result[0].c_str()));
    if(findI == bATableHoH.end()){
      // no

      // debug
      //      cout << "not found" <<endl;

      unordered_map<int, double> bATableH;
      bATableHoH[atoi(result[0].c_str())]=bATableH;

    }else{
      // yes

      // debug
      //     cout << "found" << endl;

    }

    (bATableHoH[atoi(result[0].c_str())])[atoi(result[1].c_str())]=atof(result[2].c_str());

    // debug
//     if(isnan((bATableHoH[atoi(result[0].c_str())])[atoi(result[1].c_str())])){
//       cerr << atoi(result[0].c_str()) << " " << atoi(result[1].c_str()) << endl;
//       exit(1);
//     }

    // for check
    firmH[atoi(result[0].c_str())]=1;
    firmH[atoi(result[1].c_str())]=1;

  }

  if(fATableHoH.size()==0){
    cerr << "ATableLineFile is empty" << endl;
  }

  if(bATableHoH.size()==0){
    cerr << "ATableLineFile is empty" << endl;
  }

  // input cVector

  // initialize with 0
  unordered_map<int, int> *hashII;
  hashII=&firmH;
    
  for(auto itr = hashII->begin(); itr != hashII->end(); ++itr) {
    cVectorH[itr->first]=0;
  }

  while(getline(ifs2,str)){
    string tmp;
    istringstream stream(str);
    vector<string> result;
    while(getline(stream,tmp,' ')){

      // debug
      //	  cout<< tmp << endl;

      result.push_back(tmp);
    }

    // does source exists in hash?
    unordered_map<int, int>::iterator findI;

    findI = firmH.find (atoi(result[0].c_str()));
    if(findI == firmH.end()){
      // no
      // do not record it
      // cerr << "There is a firm C vector " << atoi(result[0].c_str()) << " without element in A table" << endl;

    }else{
      // yes
      // record it
      cVectorH[atoi(result[0].c_str())]=atof(result[1].c_str());
    }

  }

  if(cVectorH.size()==0){
    cerr << "cVectorFile is empty" << endl;
  }

  // input piniVector

  // initialize with 0
  hashII=&firmH;
    
  for(auto itr = hashII->begin(); itr != hashII->end(); ++itr) {
    piniVectorH[itr->first]=0;
    piniInputVectorH[itr->first]=0;
  }

  while(getline(ifs3,str)){
    string tmp;
    istringstream stream(str);
    vector<string> result;
    while(getline(stream,tmp,' ')){

      // debug
      //	  cout<< tmp << endl;

      result.push_back(tmp);
    }

    // does source exists in hash?
    unordered_map<int, int>::iterator findI;

    findI = firmH.find (atoi(result[0].c_str()));
    if(findI == firmH.end()){
      // no
      //      cerr << "There is a firm  of Pini " << atoi(result[0].c_str()) << " vector without element in A table" << endl;;
      // do not record it

    }else{
      // yes
      // record it
      piniVectorH[atoi(result[0].c_str())]=atof(result[1].c_str());

    }

  }

  if(piniVectorH.size()==0){
    cerr << "piniVectorFile is empty" << endl;
  }


  // read firm affiliation

  while(getline(ifs4,str)){
    string tmp;
    istringstream stream(str);
    vector<string> result;
    while(getline(stream,tmp,' ')){

      // debug
      //	  cout<< tmp << endl;

      result.push_back(tmp);
    }

    // does source exists in hash?
    //    unordered_map<int, double>::iterator findI;

    int firm=atoi(result[0].c_str());
    int sector=atoi(result[1].c_str());
    firmAffiH[firm]=sector;
    
    unordered_map<int, vector<int> >::iterator findI;

    findI = secFirmHoL.find (sector);
    if(findI == secFirmHoL.end()){
      // no
      vector<int> v;
      secFirmHoL[sector]=v;

    }else{
      // yes

      // debug
      //     cout << "found" << endl;

    }

    secFirmHoL[sector].push_back(firm);

  }

  // debug
  //  printHashVectorIntInt("secFirmHoL", &secFirmHoL);


  // inputfile check start

  // debug
  //  printHash2LayerIntIntDouble("fATableHoH", &fATableHoH);

  // debug
  //  printHash2LayerIntIntDouble("bATableHoH", &bATableHoH);

  // debug
  //  printHash1LayerIntDouble("piniVectorH", &piniVectorH);


  // debug
  //  printHash1LayerIntDouble("cVectorH", &cVectorH);

  // debug
  //  printHash1LayerIntInt("firmAffiH", &firmAffiH);

  // inputfile check end

  // initialize model

  initializeModel();

  // hazard t=0
  hazard();

  if(iFullDebugFrag){
    fullDebug(-1);
  }

  // debug
  //  printAllInternal();

  // simulation
  simulation();

  // if output at last
//   unordered_map<int, double> *hashO;
//   hashO=&PactualVectorH;
//   map<int, int> sortedVectorH(hashO->begin(), hashO->end());
//   for(auto itr = sortedVectorH.begin(); itr != sortedVectorH.end(); ++itr) {
//     // show keys values
//     // i, actual P
//     ofs << itr->first << " " << itr->second << endl;
//   }

  // debug pini pact ratio
//   piniPactRatioFS.open(piniPactRatioFile);
//   if(!piniPactRatioFS){
//     cerr << "Unable to open" << endl;
//     return 1;
//   }
//   unordered_map<int, double> *hashD;
//   hashD=&PactualVectorH;
//   map<int, int> sortedVectorH(hashD->begin(), hashD->end());
//   for(auto itr = sortedVectorH.begin(); itr != sortedVectorH.end(); ++itr) {
//     // show keys values
//     // i, pini, pact, ratio
//     piniPactRatioFS << itr->first << " " << piniVectorH[itr->first] << " " << PactualVectorH[itr->first] << " " << (PactualVectorH[itr->first]/piniVectorH[itr->first]) << endl;
//   }


  // end time
  time(&etime);
  int diffTime=difftime(etime, btime);
  printf("%d sec. elapsed\n",diffTime);

  return 0;

}
