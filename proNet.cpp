#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <map>
#include <unordered_map>
#include <vector>
#include <ctime>
#include <cstdio>
#include <unistd.h>
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include <random>
#include <chrono>
#include <boost/container/flat_map.hpp>
#include <omp.h>
#include <chrono>

using namespace std;

void sim_predefinedDelta(int t);

// predefined table
class Predefined{
public:
  int start;
  int end;
  double delta;
  Predefined(int s, int e, double d){start=s;end=e;delta=d;}
  int GetS(){return start;}
  int GetE(){return end;}
  double GetD(){return delta;}
};
unordered_map<int, std::vector<Predefined>> predefinedHoL;

// debug
// full debug option (output for step by step)
// t, i, pini, delta, pcap, d, pmax, pact, d*, s_j,...
int iFullDebugFlag = 0;

// mute flag
int iMuteFlag = 0;

// pact suppredd flag
int iPActSupFlag = 0;

// internal structure

// A table
// direction is the same as Hallegatte's paper
// it means the direction is the opposite to IOtables
// supplier -> client volume
// forward: client -> supplier volume
// backward: supplier -> client volume
size_t NUM_NODES = 0;
std::vector<boost::container::flat_map<int, double>> fATableHoH;

// C vector
std::vector<double> cVectorH;

// Pini vector
// total amount output at initial
std::vector<double> piniVectorH;

// input amount
std::vector<double> piniInputVectorH;

// valueAdded ratio
std::vector<double> valueAddedVectorH;

// firm affi hash
std::vector<int> firmAffiH;

// firm name hash
std::vector<int> firmH;

// order between firm client -order-> supplier
// this tells the total input of client
std::vector<boost::container::flat_map<int, double> > fOFirmFirmHoH;

// adjusted order between firm client -order-> supplier
std::vector<boost::container::flat_map<int, double> > fOAdjustedFirmFirmHoH;

// backward order between firm supplier <-order- client
std::vector<boost::container::flat_map<int, double> > bOFirmFirmHoH;

// backward order between firm supplier <-order- client at start
std::vector<boost::container::flat_map<int, double> > bOIniFirmFirmHoH;

// backward adjusted order between firm supplier <-order- client
std::vector<boost::container::flat_map<int, double> > bOAdjustedFirmFirmHoH;

// demand on firm
std::vector<double> DFirmVectorH;

// demand on firm (previous)
std::vector<double> preDFirmVectorH;

// adjusted demand on firm
std::vector<double> DAdjustedFirmVectorH;

// adjusted c on firm
std::vector<double> cAdjustedVectorH;

// Pini Cap vector
std::vector<double> piniCapVectorH;

// S on firm for firm (client <-- supplier)
std::vector<boost::container::flat_map<int, double> > SFirmFirmHoH;

// S on firm for sector (total)
std::vector<boost::container::flat_map<int, double> > SFirmSectorHoH;

// A on firm for sector (total)
std::vector<boost::container::flat_map<int, double> > AFirmSectorHoH;

// Pprop on firm for sector
// PsProp
// How much can produce if we consider only the inventory of the sector
// (sfirm)/(afirm)*pini
std::vector<boost::container::flat_map<int, double> > PpropFirmSectorHoH;

// Pini Max vector
// (not ratio, ammount)
std::vector<double> PiniMaxVectorH;

// Pprop (but amount) (inventory based minimum amount of production)
std::vector<double> MinPropAmountVectorH;

// Pactual vector
std::vector<double> PactualVectorH;

// Damage vector 0<=delta<=1
std::vector<double> DeltaVectorH;

// Damage reserve vector
std::vector<double> DeltaReserveVectorH;
// Damage Start vector 0<=day
std::vector<int> DeltaStartVectorH;


// secFirmHoL
unordered_map<int, vector<int> > secFirmHoL;

// inventory size(day)
// each firm has its own target inventory size that has poisson dist.
// lambda is inventoryN
int inventoryN = 15;

// inventory distribution type
int iInvDistType = 0;

// minimal inventory
int iInventoryMin = 2;

// target inventory vector (day based)
std::vector<double> tgtInvVectorH;

// rationing algorithm switch
int iRationingTypeFlag = 0;

// order algorithm switch
int iOrderTypeFlag = 0;


// adjustment ratio for inventory
double tau = 6;

// recovery day
int iRecoveryDay = 0;

// recovery type
int iRecoveryType = -1;


// forcibly stop day
int iForceStopDay = 0;

// for output stats
double lastTotalPactual;
double lastTotalValueAdded;

// get snapshot pact option
int iPrintSnapFlag = 0;

// recovery step
std::vector<double> recoveryStepH;

// hidden delta
std::vector<double> hiddenDeltaVectorH;

// snapshot step
unordered_map<int, int> hPrintSnapStep;
string sPrintSnapStep;

// snapshot file
string sPrintSnapFileBase = "";

typedef std::unordered_map<std::string, int> firm_index_t;

typedef std::unordered_map<int, std::string> rev_firm_index_t;
rev_firm_index_t rev_firm_index;

// selective output flag
int iSelectiveOutput = 0;
string sSelectiveOutput;
vector<map<int, int> > selectiveOutputLoH;

// Daily CVector
int iDailyCVector = 0;
string sDailyCVectorFile;
string sDailyCVectorFileBase;
std::unordered_map<int, int> dailyCVectorDayH;
std::vector<double> cVectorOrgH;

// delta flag
int iDeltaFile=0;

// predefined delta flag
int iPredefinedDeltaFile=0;

// functions

#ifndef Pi
#define Pi 3.141592653589793238462643
#endif

void muteSwitchedCerr(const std::string s) {
  if (!iMuteFlag) {
    std::cerr << s;
  }
}

template<class F>
void MeasureElapsed(std::string tag, F f) {
  auto start = std::chrono::system_clock::now();
  f();
  auto end = std::chrono::system_clock::now();
  auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::stringstream msg;
  msg << "elapsed time for " << tag << ": " << dt / 1000.0 << std::endl;
  muteSwitchedCerr(msg.str());
}

const double epsilon = 1e-20 /* some small number such as 1e-5 */;

inline bool almost_equal(double x, double y) {
  return std::abs(x - y) <= epsilon;
}

void printSnapshot(int t) {

  if (hPrintSnapStep.find(t) != hPrintSnapStep.end()) {
    std::ostringstream oss;
    oss << sPrintSnapFileBase << t;
    std::ofstream snapOfs(oss.str());

    if (!snapOfs) {
      cerr << "Unable to open " << oss.str() << " as snap output file at step " << t << endl;
      exit(1);
    }

    map<int, double> sortedVectorH;
    for (size_t i = 0; i < PactualVectorH.size(); i++) {
      sortedVectorH[i] = PactualVectorH[i];
    }

    snapOfs << "firm" << " " << "pini" << " " << "pact" << " " << "demand" << " " << "supply" << " " << "delta"
            << endl;;

    for (size_t i = 0; i < PactualVectorH.size(); i++) {
      double pini = piniVectorH[i];
      double pact = PactualVectorH[i];
      double demand = DFirmVectorH[i];
      double supply = MinPropAmountVectorH[i];
      double delta = DeltaVectorH[i];

      string id = rev_firm_index[i];
      
      // simple
      //	snapOfs << itr->first << " " << setprecision(17) << itr->second << endl;;

      // ini, actual, supply, demand
      snapOfs << std::setfill('0') << std::right << std::setw(9) << id << " " << setprecision(6) << pini << " " << pact
              << " " << demand << " " << supply << " " << delta << endl;;
    }

    std::ostringstream ostar;
    ostar << sPrintSnapFileBase << "OStar" << t;

    ofstream snapOfs2(ostar.str());
    if (!snapOfs2) {
      cerr << "Unable to open " << ostar.str() << " as snap output file at step " << t << endl;
      exit(1);
    }

    // header
    snapOfs2 << "client" << " " << "supplier" << " " << "Aij" << " " << "AdjustedOrder" << endl;

    // order star
    for (size_t i = 0; i < fATableHoH.size(); i++) {
      for (auto pi : fATableHoH[i]) {
        int j = pi.first;

        double dAdjustedOrder = fOAdjustedFirmFirmHoH[i][j];
        double dA = pi.second;

	string idi = rev_firm_index[i];
        string idj = rev_firm_index[j];
	
        snapOfs2 << idi << " " << idj << " " << dA << " " << dAdjustedOrder << endl;
      }
    }
  }
}

static void usage(char *argv) {
  cerr << "Usage: " << endl;
  cerr << "  " << argv << " [option] aTableLineFile cVectorFile piniVectorFile firmAffiFile" << endl;
  cerr << "Option: " << endl;
  fprintf(stderr,
          "  -b (int): recovery step type (0: proportional, 1: uniform (1/recday regardless of delta), 2: spike (immeidate), 3: renewal equation (cannot use without -c or -C), 4: exponential, 5: damped renewal (default: 0)\n");
  fprintf(stderr, "  -c (int): proportional recovery day (default: off. greater than 0) (cannot use with -C)\n");
  fprintf(stderr, "  -d (str): delta file (default null) (col.1 firm id, col.2 delta, (optional) col.3 day to be given (default 0)\n");
  fprintf(stderr,
          "  -D (int): distribution of inventory and target (0: target and ini inventory uniform, 1: target uniform and ini inventory Poisson, 2: target poisson and ini inventory filled.) (default 0: uniform)\n");
  fprintf(stderr, "  -e : minimal inventory size (default 2) \n");
  fprintf(stderr, "  -f (str): final consumption output file (default null)\n");
  fprintf(stderr, "  -F : fulldebug\n");
  fprintf(stderr, "  -i : inventory size (possible used as average of -D) \n");
  fprintf(stderr, "  -M (str1):(str2):... daily CVector input file(s) (str1(base file name)str2(day).txt) (default null) designate id with FC different from CVector. others have CVector value as FC.\n");
  fprintf(stderr, "  -m : mute\n");
  fprintf(stderr,
          "  -o (int): order algorithm type. 0: normal, 1: keep initial demand, 2: ignore negative inv adjustment\n");
  fprintf(stderr,
          "  -p (file:int:int...): get snapshot of production at indicated step(s). for more than one snapshot, int should be separated by :\n");
  fprintf(stderr, "  -r (int): random seed (default 331)\n");
  fprintf(stderr,
          "  -R (int): rationing algorithm 0. proportional, 1. FC low priority, 2. lower has high priority, 3. 1*2 (default proportional)\n");
  fprintf(stderr, "  -s (str): stats output file (default null)\n");
  fprintf(stderr, "  -S (int): forcibly stop day (recovery occur AFTER the day) (delay of recovery) (default 0)\n");
  fprintf(stderr, "  -t (int): simulation step (default 10)\n");
  fprintf(stderr, "  -T (int): Tau (inventory adjustment ratio) (default 6)\n");
  fprintf(stderr, "  -u : supress PAct output\n");
  fprintf(stderr, "  -v (str): value added output file (default VA.txt)\n");
  fprintf(stderr, "  -x (str): output total direct production loss (default null)\n");
  fprintf(stderr, "  -X (str): PAct output file (default PAct.txt)\n");
  fprintf(stderr, "  -z (file:file...): selective output files (default: None. Only one column to give IDs.)\n");
  fprintf(stderr, "  -Z (str): predefined delta table (default null)\n");
  exit(1);
}

// print total direct production loss
void printTotalDirectProductionLoss(int t, std::ofstream &dlossOfs) {

  double dTotalCap = 0;
  double dTotalValueAdded = 0;
  for (size_t i = 0; i < DeltaVectorH.size(); i++) {
    double dCap = DeltaVectorH[i] * piniVectorH[i];
    double dVA = dCap * valueAddedVectorH[i];

    dTotalCap += dCap;
    dTotalValueAdded += dVA;
  }

  dlossOfs << t << " " << dTotalCap << " " << dTotalValueAdded << endl;
}

// full debug
void fullDebug(int t, std::ofstream &debugOfs) {

  for (size_t i = 0; i < firmH.size(); i++) {
    double pini = piniVectorH[i];
    double delta = DeltaVectorH[i];
    double pcap = piniCapVectorH[i];
    double demand = DFirmVectorH[i];
    double pmax = PiniMaxVectorH[i];
    double pact = PactualVectorH[i];
    double demandA = DAdjustedFirmVectorH[i];

    debugOfs << t << ", " << i << ", " << pini << ", " << delta << ", " << pcap << ", " << demand << ", " << pmax
             << ", "
             << pact << ", " << demandA;

    //       // inventory s
    const map<int, int> sortedVector2H(SFirmFirmHoH[i].begin(), SFirmFirmHoH[i].end());
    for (const auto &itr : sortedVector2H) {
      debugOfs << ", s(" << i << "," << itr.first << ")_" << itr.second;
    }

    // order
    const map<int, int> sortedVector3H(fOFirmFirmHoH[i].begin(), fOFirmFirmHoH[i].end());
    for (const auto &itr: sortedVector3H) {
      debugOfs << ", O(" << i << "," << itr.first << ")_" << itr.second;
    }
    debugOfs << endl;
  }
}

typedef long unsigned int luint;
luint poisson(luint lambda) {
  double L = exp(-double(lambda));
  luint k = 0;
  double p = 1;
  do {
    k++;
    double rNum = ((double) (rand() / (double) RAND_MAX));  // [TODO] stop using rand
    p *= rNum;
  } while (p > L);
  return (k - 1);
}

void initializeModel(uint64_t seed,
                     std::ofstream &pactOfs,
                     std::ofstream &pvalueOfs,
                     std::ofstream &fcOfs,
                     std::ofstream &statsOfs,
                     std::vector<std::ofstream> &selectOutPActL,
                     std::vector<std::ofstream> &selectOutVAL) {

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

  // decide inventory coefficient
  // initial inventory size
  std::vector<int> iniInventoryH(NUM_NODES, 0);
  std::mt19937 rng(seed);
  std::poisson_distribution<int> dist(inventoryN);

  for (size_t i = 0; i < firmH.size(); i++) {
    // target inventory vector
    // initialize tgtInvVectorH
    if (iInvDistType == 0 || iInvDistType == 1) {
      // uniform tgt and uniform ini
      // uniform tgt and poisson ini
      tgtInvVectorH[i] = inventoryN;
    } else if (iInvDistType == 2) {
      int tgt = 0;
      while (tgt < iInventoryMin) {
        tgt = dist(rng);
      }
      tgtInvVectorH[i] = tgt;
    } else {
      throw std::runtime_error("invalid option for -D");
    }

    if (iInvDistType == 0 || iInvDistType == 2) {
      // uniform tgt and uniform ini
      // poisson tgt and poisson ini
      iniInventoryH[i] = tgtInvVectorH[i];
    } else if (iInvDistType == 1) {
      // uniform tgt and poisson ini
      iniInventoryH[i] = dist(rng);
    } else {
      throw std::runtime_error("invalid option for -D");
    }

    // initialize min prop amount vector
    MinPropAmountVectorH[i] = piniVectorH[i] * iniInventoryH[i];

  }

  // using fAtable
  fOFirmFirmHoH = fATableHoH;
  fOAdjustedFirmFirmHoH = fATableHoH;
  // set inventory size
  SFirmFirmHoH = fATableHoH;
  for (size_t i = 0; i < SFirmFirmHoH.size(); i++) {
    for (auto &sij : SFirmFirmHoH[i]) {
      sij.second *= iniInventoryH[i];
    }
  }

  std::vector<std::vector<std::pair<int, double>>> vBLinks(NUM_NODES);
  for (size_t i = 0; i < fATableHoH.size(); i++) {
    // for each client
    // i: client
    // j: suppliers
    // suppliers <- order -- client
    for (const auto &ai : fATableHoH[i]) {
      const size_t j = ai.first;
      const double aij = ai.second;
      // input vector
      piniInputVectorH[i] += aij;
      // backward creating hash
      vBLinks[j].push_back(std::make_pair(i, aij));
    }
  }
  for (size_t j = 0; j < vBLinks.size(); j++) {
    bOFirmFirmHoH[j].insert(vBLinks[j].cbegin(), vBLinks[j].cend());
  }
  // initialize start hash
  bOAdjustedFirmFirmHoH = bOFirmFirmHoH;
  bOIniFirmFirmHoH = bOFirmFirmHoH;


  // demand

  // initialize
  // bOFirmFirmHoH
  // cAdjusted
  // cVectorH is defined for every firm

  DFirmVectorH = cVectorH;
  cAdjustedVectorH = cVectorH;

  for (size_t i = 0; i < bOFirmFirmHoH.size(); i++) {
    // fOrder: client -> supplier
    // bOrder: supplier -> client
    // demand means orders toward the firm
    double totalDemand = 0.0;
    for (const auto &it: bOFirmFirmHoH[i]) {
      totalDemand += it.second;
    }
    DFirmVectorH[i] += totalDemand;
  }

  // copy Actual to D
  DAdjustedFirmVectorH = DFirmVectorH;
  piniCapVectorH = DFirmVectorH;
  PiniMaxVectorH = DFirmVectorH;
  PactualVectorH = DFirmVectorH;

  // correct pini
  for (size_t i = 0; i < piniVectorH.size(); i++) {
    if(!almost_equal(DFirmVectorH[i],piniVectorH[i])){
    }
    piniVectorH[i] = DFirmVectorH[i];
  }
  
  // calculate value added
  for (size_t i = 0; i < piniVectorH.size(); i++) {
    double output = piniVectorH[i];
    double input = piniInputVectorH[i];
    valueAddedVectorH[i] = (output <= 0.0) ? 0.0 : (output - input) / output;
  }

  

  // output

  // pact
  double totalPactual = 0.0;
  for (double x : PactualVectorH) {
    totalPactual += x;
  }
  // disaster happens at 0, so as a pre-disaster, t=-1 is output

  if(!iPActSupFlag){
    pactOfs << -1 << " " << setprecision(17) << totalPactual << endl;
  }

  // selective
  if (iSelectiveOutput) {
    for (int i = 0; i < selectOutPActL.size(); i++) {

      double totalPactual = 0.0;
      for (auto j : selectiveOutputLoH[i]) {
        int id = j.first;
        totalPactual += PactualVectorH[id];
      }

      // disaster happens at 0, so as a pre-disaster, t=-1 is output
      selectOutPActL[i] << -1 << " " << setprecision(17) << totalPactual << endl;

    }
  }

  // va
  double totalVAdded = 0.0;
  for (size_t i = 0; i < PactualVectorH.size(); i++) {
    totalVAdded += valueAddedVectorH[i] * PactualVectorH[i];
  }
  if (pvalueOfs) {
    pvalueOfs << -1 << " " << setprecision(17) << totalVAdded << endl;
  }

  // selective
  if (iSelectiveOutput) {
    for (int i = 0; i < selectOutVAL.size(); i++) {

      double totalVA = 0.0;
      for (auto j : selectiveOutputLoH[i]) {
        int id = j.first;
        totalVA += valueAddedVectorH[id] * PactualVectorH[id];
      }

      // disaster happens at 0, so as a pre-disaster, t=-1 is output
      selectOutVAL[i] << -1 << " " << setprecision(17) << totalVA << endl;

    }
  }


  // create A total (i,s) (Fixed)
  // A is used for target inventory size
  // create S total (i,s)
  for (size_t i = 0; i < fATableHoH.size(); i++) {

    // check j <- i order
    for (const auto &it: fATableHoH[i]) {
      const size_t j = it.first;
      const double a = it.second;

      // get sector
      int sj = firmAffiH[j];
      if (AFirmSectorHoH[i].find(sj) == AFirmSectorHoH[i].end()) {
        AFirmSectorHoH[i][sj] = 0.0;
        SFirmSectorHoH[i][sj] = 0.0;
      }
      AFirmSectorHoH[i][sj] += a;
      SFirmSectorHoH[i][sj] += iniInventoryH[i] * fATableHoH[i][j];
    }
  }



  // output to stat file
  if (statsOfs) {
    statsOfs << "totalPactutalStart " << totalPactual << endl;
    statsOfs << "totalValueAddedStart " << totalVAdded << endl;
  }

  // output
  // final consumption
  if (fcOfs) {
    double totalFC = 0.0;
    for (double x: cAdjustedVectorH) {
      totalFC += x;
    }
    fcOfs << "-1 " << totalFC << endl;
  }
  
}

void hazard(const std::string &deltaFile, std::ofstream &statsOfs) {

  hiddenDeltaVectorH = DeltaReserveVectorH;

  // put recovery step
  if (iRecoveryDay > 0) {
    for (size_t i = 0; i < DeltaReserveVectorH.size(); i++) {
      recoveryStepH[i] = DeltaReserveVectorH[i] / iRecoveryDay;
    }
  }

}


void sim_inventoryRenewal() {

  // process 3 ==============================================================
  // renewal of inventory
  // SFirmFirmHoH

  // original paper uses j,i but here uses i,j

  //#pragma omp parallel for
  for (size_t i = 0; i < SFirmFirmHoH.size(); i++) {
    const auto row = &SFirmFirmHoH[i];

    // Inventory is reduced by equal proportion
    // obtain required amount (in ratio) of consumption for each sector
    auto RequiredConsumptionRatioFirmSectorH = AFirmSectorHoH[i];  // [TODO] copying. could be slow.

    for (auto &itrA : RequiredConsumptionRatioFirmSectorH) {
      if (almost_equal(SFirmSectorHoH[i][itrA.first], 0)) {
        itrA.second = 0;
      } else {
        itrA.second = (PactualVectorH[i] / piniVectorH[i]) * (itrA.second / SFirmSectorHoH[i][itrA.first]);
      }
    }

    for (auto &itr2 : SFirmFirmHoH[i]) {
      // itr2->first is j
      // itr2->second is s(i,j)
      // if piniVector=0, renewal equation is different
      const int j = itr2.first;
      const int sj = firmAffiH[j];
      const double sjConsRatio = RequiredConsumptionRatioFirmSectorH[sj];

      double newInventory;
      if (piniVectorH[i] == 0) {
        // avoid zero division
        // no change
        newInventory = itr2.second;
      } else {
        newInventory =
            itr2.second // S(i,j)
                + fOAdjustedFirmFirmHoH[i][j] // O*(i,j)
                - itr2.second * sjConsRatio;
        //	    - (fATableHoH[i])[j] * PactualVectorH[i] / piniVectorH[i] // A(i,j)*Pia(t-1)/Piini
      }

      // here minus inventory is possible
      // this is firm level inventory
      // and production is based on sector
      // therefore, at consumption, there is possible that inventory lower than zero
      if (newInventory < 0.0 || almost_equal(newInventory, 0.0)) {
        newInventory = 0.0;
      }
      itr2.second = newInventory;
    }
  }
}


void sim_order() {
  // process 4 ==============================================================
  // renewal of Order
  // fOFirmFirmHoH, bOFirmFirmHoH
  // Adjusted is later

  // NB: 'b' means backward and the opposite of 'f'

  // use hash again
  // original paper uses j,i but here uses i,j

#pragma omp parallel for
  for (size_t i = 0; i < fOFirmFirmHoH.size(); i++) {
    for (auto &itr2 : fOFirmFirmHoH[i]) {
      // itr2->first is j
      // itr2->second is o(i,j)
      const int j = itr2.first;

      // use Dadjusted
      double newOrder;
      if (piniVectorH[i] == 0) {
        // order=0 is not accepted
        // newOrder!=0
        // order the amount based on fATable
        newOrder = fATableHoH[i][j];
        // Si,j is constant
        // consumption is always equalt to fATable
      } else {
        if (iOrderTypeFlag == 0) {
          // normal order
          newOrder =
              // original
              fATableHoH[i][j] * DAdjustedFirmVectorH[i] / piniVectorH[i] // 1st term
                  + 1 / tau
                      * (tgtInvVectorH[i] * fATableHoH[i][j] * DAdjustedFirmVectorH[i]
                          / piniVectorH[i] // 2nd term target inventory
                          - SFirmFirmHoH[i][j]); // 3rd term

        } else if (iOrderTypeFlag == 1) {
          // inventory ini keep
          newOrder =
              fATableHoH[i][j] * DAdjustedFirmVectorH[i] / piniVectorH[i] // 1st term
                  + 1 / tau
                      * (tgtInvVectorH[i] * fATableHoH[i][j] // 2nd term i not depend on D. try to keep first inventory
                          - SFirmFirmHoH[i][j]); // 3rd term
        } else if (iOrderTypeFlag == 2) {
          // ignore negative inventory
          if ((tgtInvVectorH[i] * fATableHoH[i][j] * DAdjustedFirmVectorH[i]
              / piniVectorH[i] // 2nd term. target inventory
              - SFirmFirmHoH[i][j]) <= 0) // 3rd term
          {
            // target inventory is negative
            newOrder =
                fATableHoH[i][j] * DAdjustedFirmVectorH[i] / piniVectorH[i]; // 1st term
          } else {
            // target inventory is positive
            newOrder =
                // original
                fATableHoH[i][j] * DAdjustedFirmVectorH[i] / piniVectorH[i]
                    // 1st term
                    + 1 / tau
                        * (tgtInvVectorH[i] * fATableHoH[i][j] * DAdjustedFirmVectorH[i] / piniVectorH[i]
                            // 2nd term. target inventory
                            - SFirmFirmHoH[i][j]); // 3rd term
          }
        } else {
          exit(1);
        }
      }

      // if order is minus, set it to zero
      if (newOrder < 0.0) {
        newOrder = 0.0;
      }

      // renewal of order
      itr2.second = newOrder;
      bOFirmFirmHoH[j][i] = newOrder;
    }
  }
}

void sim_demand() {
  // process 5 ==============================================================
  // renewal of Demand
  // renewal of DFirmVectorH
  preDFirmVectorH = DFirmVectorH;

  // demand is reset by cvector out of consideration of terminal nodes
  // do not reset by zero
  // use hash
  for (size_t i = 0; i < firmH.size(); i++) {
    DFirmVectorH[i] = cVectorH[i];
  }

  // original paper uses j,i but here uses i,j
  for (size_t i = 0; i < bOFirmFirmHoH.size(); i++) {
    // fOrder: client -> supplier
    // bOrder: supplier -> client
    // demand is order toward the firm
    // use bOrder and total client order
    double totalDemand = 0;
    for (const auto &itr2 : bOFirmFirmHoH[i]) {
      totalDemand += itr2.second;
    }

    DFirmVectorH[i] += totalDemand;

    // ignore the error of the doubleing point
    // without the error causes fluctuation
    // or introducing gmp (not done)
    if (abs(preDFirmVectorH[i] - DFirmVectorH[i]) < epsilon) {
      DFirmVectorH[i] = preDFirmVectorH[i];
    }
  }
}

void sim_PAct() {
  // process 6 ==============================================================
  // calculate of Pactual

  // Pcap is decided
  // first renew Stotal i,s
  // process 7 is included because it is simple

  // clear sector inventory with 0
  // original paper uses j,i but here uses i,j
  for (size_t i = 0; i < SFirmSectorHoH.size(); i++) {
    for (auto &itr2: SFirmSectorHoH[i]) {
      itr2.second = 0.0;
    }
  }

  // Stotal i,s does not require reset because it is set by Atotal i,s
  // create setor inventory
  // and create PpropFirmSectorHoH(PsProp)
#pragma omp parallel for
  for (size_t i = 0; i < SFirmFirmHoH.size(); i++) {
    for (const auto &itr2: SFirmFirmHoH[i]) {
      int j = itr2.first;
      // itr2->second is S(i,sj)
      double s = itr2.second;
      // get sector
      int sj = firmAffiH[j];
      SFirmSectorHoH[i][sj] += s;
    }

    // prop can be calculated because i is know at this moment
    // for all sector
    for (const auto &itr2 : SFirmSectorHoH[i]) {
      const int sj = itr2.first;
      if (AFirmSectorHoH[i][sj] == 0.0) {
        cerr << "0 div (AFirmSectorHoH[i])[sj]: " << AFirmSectorHoH[i][sj] << endl;
        exit(1);
      }
      PpropFirmSectorHoH[i][sj] = SFirmSectorHoH[i][sj] / AFirmSectorHoH[i][sj] * piniVectorH[i];
    }
  }

#pragma omp parallel for
  for (size_t i = 0; i < firmH.size(); i++) {
    const auto hash2 = &PpropFirmSectorHoH[i];

    // get i's SFirmSectorHoH
    // if there is none of them, that is beginning node
    // hash is possibly empty
    if (hash2->empty()) {
      PiniMaxVectorH[i] = piniCapVectorH[i];
      continue;
    }

    // the firm i has supplier

    // for all i's sector
    double minProp = 0.0;
    int first = 1;
    for (const auto &itr2 : *hash2) {
      if (first == 1) {
        minProp = itr2.second;
        first = 0;
      } else {
        if (minProp > itr2.second) {
          minProp = itr2.second;
        }
      }
    }

    // at this moment, minProp is the mimimum of Pprop
    MinPropAmountVectorH[i] = minProp;
    PiniMaxVectorH[i] = (piniCapVectorH[i] > minProp) ? minProp : piniCapVectorH[i];
  }

  // create PactualVectorH
  // demand is the newest D
#pragma omp parallel for
  for (size_t i = 0; i < firmH.size(); i++) {
    // compare D and set Actual
    PactualVectorH[i] = (PiniMaxVectorH[i] > DFirmVectorH[i]) ? DFirmVectorH[i] : PiniMaxVectorH[i];

    // set Dadjusted
    DAdjustedFirmVectorH[i] = PactualVectorH[i];

    // ***** rationing begins from here *****
    // get OrderRatio
    // mostly, that is 1
    double adjustRatio;  // [TODO] initialize??

    // Order  FC  rationing
    if (iRationingTypeFlag == 0) {
      // this is original paper's rationing
      // 0 possibly exists
      adjustRatio = (DFirmVectorH[i] == 0.0) ? 0.0 : PactualVectorH[i] / DFirmVectorH[i];

      // renew C* first
      // cadjusted
      cAdjustedVectorH[i] = adjustRatio * cVectorH[i];

      // get i's fOFirmFirmHoH and bOFirmFirmHoH
      // for all i's supplier
      for (const auto &itr2 : bOFirmFirmHoH[i]) {
        const int j = itr2.first;
        // adjusted f
        bOAdjustedFirmFirmHoH[i][j] = itr2.second * adjustRatio;
        // adjusted b
        fOAdjustedFirmFirmHoH[j][i] = itr2.second * adjustRatio;
      }
      // original rationing ends here
    } else if (iRationingTypeFlag == 1) {
      // rationing but FC is always considered last
      double shortD = DFirmVectorH[i] - PactualVectorH[i];
      double FC = cVectorH[i];
      if (shortD < 0.0) {
        shortD = 0.0;
      }

      if (FC >= shortD) {
        // FC is enough to fulfill
        // cadjusted
        cAdjustedVectorH[i] = FC - shortD;
        adjustRatio = 1.0;
      } else {
        // FC is no enough
        // rationing occurs

        // cadjusted
        cAdjustedVectorH[i] = 0.0;

        // total amount of order
        double order = DFirmVectorH[i] - FC;

        // there possibly exists 0
        adjustRatio = (DFirmVectorH[i] == 0) ? 0.0 : PactualVectorH[i] / order;
      }

      // get i's fOFirmFirmHoH and bOFirmFirmHoH
      // for all i's supplier
      for (const auto &itr2 : bOFirmFirmHoH[i]) {
        const int j = itr2.first;
        // adjusted f
        bOAdjustedFirmFirmHoH[i][j] = itr2.second * adjustRatio;
        // adjusted b
        fOAdjustedFirmFirmHoH[j][i] = itr2.second * adjustRatio;
      }
      // rationing: fc is last, ends here
    } else if (iRationingTypeFlag == 2) {
      // soothing algorithm
      // inoue & todo 17 uses this algorithm
      // FC's ratio is not always 1 but cvec / cvecorg

      // get i's bOFirmFirmHoH
      const auto bhash = &bOFirmFirmHoH[i];

      // get i's bOIniFirmFirmHoH
      const auto bshash = &bOIniFirmFirmHoH[i];

      // bOrder ratio
      // bAdjust ratio
      std::vector<std::pair<int, double>> ini_rat;
      std::vector<std::pair<int, double>> ini_adjust;
      for (const auto &itrb : *bhash) {
        const int j = itrb.first;
        ini_rat.emplace_back(std::make_pair(j, itrb.second / bshash->find(j)->second));
        ini_adjust.emplace_back(std::make_pair(j, 0.0));
      }
      boost::container::flat_map<int, double> bRatH(ini_rat.cbegin(), ini_rat.cend());
      boost::container::flat_map<int, double> bAdjustH(ini_adjust.cbegin(), ini_adjust.cend());

      double ratFC;
      if(iDailyCVector){
	ratFC= cVectorH[i]/cVectorOrgH[i];
      }else{
	ratFC= 1;
      }
      double adjustFC = 0.0;

      // rationing
      // remaining act
      double remainingAct = PactualVectorH[i];

      // bRatH, find minimum bOrder's ratio
      while (1) {
        double minD = -1.0;

        for (const auto &itrb : bRatH) {
          if (almost_equal(itrb.second, 0.0)) {
            // skip
            // it is 0
          } else if (minD < 0) {
            minD = itrb.second;
          } else if (minD > itrb.second) {
            minD = itrb.second;
          }
        }

        if (almost_equal(ratFC, 0.0)) {
        } else if (ratFC <= 0.0) {
        } else if (minD < 0) {
          minD = ratFC;
        } else if (minD > ratFC) {
          minD = ratFC;
        }

        // consider minimum
        if (minD < 0.0 || almost_equal(minD, 0.0)) {
          break;
        }

        // how much is required for the fulfillment of minimum ratio
        double totalMin = 0;
        for (const auto &itrb : bRatH) {
          totalMin += minD * bshash->find(itrb.first)->second;
        }
        // fc
        totalMin += minD * cVectorOrgH[i];

        // is it possible to fulfill the minimum
        if (remainingAct > totalMin) {
          // add minD to the achieved
          // and delete from bRatH

          // record to delete
          std::vector<int> deleteBRatH;
          for (auto &itrb : bRatH) {
            itrb.second -= minD;
            bAdjustH[itrb.first] += minD;
            if (itrb.second < 0.0 || almost_equal(itrb.second, 0.0)) {
              deleteBRatH.emplace_back(itrb.first);
            }
          }

          for (const auto &i : deleteBRatH) {
            bRatH.erase(i);
          }

          // fc
          if (ratFC > 0) {
            ratFC -= minD;
            adjustFC += minD;

            if (ratFC < 0.0 || almost_equal(ratFC, 0.0)) {
              ratFC = 0.0;
            }
          }

          remainingAct -= totalMin;

          if (bRatH.empty() && almost_equal(ratFC, 0.0)) {
            break;
          }
          if (remainingAct < 0 || almost_equal(remainingAct, 0.0)) {
            break;
          }
        } else {
          // find the ratio for all remaining demand that can be fulfilled
          double vTotal = 0.0;
          for (const auto &itrb : bRatH) {
            if (itrb.second > 0.0) {
              vTotal += bshash->find(itrb.first)->second;  // [TODO] what if findI is end??
            }
          }

          // fc
          if (ratFC > 0) {
            vTotal += cVectorOrgH[i];
          }

          const double ratio = (vTotal <= 0) ? 0.0 : remainingAct / vTotal;

          // add to fulfilled ratio
          for (const auto &itrb : bRatH) {
            if (itrb.second > 0.0) {
              bAdjustH[itrb.first] += ratio;
            }
          }
          if (ratFC > 0.0) {
            adjustFC += ratio;
          }
          break;
        }
      }

      cAdjustedVectorH[i] = adjustFC * cVectorOrgH[i];
      for (const auto &itrb : bAdjustH) {
        const int j = itrb.first;
        // adjusted f
        bOAdjustedFirmFirmHoH[i][j] = itrb.second * bOIniFirmFirmHoH[i][j];
        // adjusted b
        fOAdjustedFirmFirmHoH[j][i] = itrb.second * bOIniFirmFirmHoH[i][j];
      }
    } else if (iRationingTypeFlag == 3) {
      // soothing algorithm and considers FC last

      // get i's bOFirmFirmHoH
      auto *bhash = &bOFirmFirmHoH[i];
      // get i's bOIniFirmFirmHoH
      auto *bshash = &bOIniFirmFirmHoH[i];
      // bOrder ratio
      unordered_map<int, double> bRatH;
      // bAdjust ratio
      unordered_map<int, double> bAdjustH;

      for (const auto &itrb : *bhash) {
        bRatH[itrb.first] = itrb.second / bshash->find(itrb.first)->second;
        bAdjustH[itrb.first] = 0.0;
      }

      // consider fc here
      double FC = cVectorH[i];
      double dDwithoutFC = DFirmVectorH[i] - FC;
      double dMarginActual = PactualVectorH[i] - dDwithoutFC;

      if (dMarginActual <= 0.0) {
        cAdjustedVectorH[i] = 0.0;
      } else if (dMarginActual > 0.0 && dMarginActual < FC) {
        cAdjustedVectorH[i] = dMarginActual;
      } else {
        cAdjustedVectorH[i] = FC;
      }

      double ratFC = 0.0;
      double adjustFC = 0.0;

      // rationing
      // remaining act
      double remainingAct = PactualVectorH[i] - cAdjustedVectorH[i];

      // bRatH, find minimum bOrder's ratio

      while (1) {
        //	  int minI;
        double minD = -1.0;
        for (const auto &itrb : bRatH) {
          if (almost_equal(itrb.second, 0.0)) {
          } else if (minD < 0.0) {
            minD = itrb.second;
          } else if (minD > itrb.second) {
            minD = itrb.second;
          }
        }
        // consider minimum
        if (minD < 0.0 || almost_equal(minD, 0.0)) {
          break;
        }

        // how much is required for the fulfillment of minimum ratio
        double totalMin = 0.0;
        for (const auto &itrb : bRatH) {
          totalMin += minD * bOIniFirmFirmHoH[i].find(itrb.first)->second;
        }

        // is it possible to fulfill the minimum
        if (remainingAct > totalMin) {
          // possible
          // add minD to the achieved
          // and delete from bRatH

          // record to delete
          unordered_map<int, int> deleteBRatH;
          for (auto &itrb : bRatH) {
            itrb.second -= minD;
            bAdjustH[itrb.first] += minD;
            if (itrb.second < 0.0 || almost_equal(itrb.second, 0.0)) {
              deleteBRatH[itrb.first] = 1;
            }
          }

          for (const auto &itrb : deleteBRatH) {
            bRatH.erase(itrb.first);
          }

          remainingAct -= totalMin;

          if (bRatH.empty()) {
            break;
          }
          if (remainingAct < 0 || almost_equal(remainingAct, 0.0)) {
            break;
          }
        } else {
          //  not possible
          //  or equal

          // find the ratio for all remaining demand that can be fulfilled
          double vTotal = 0.0;
          for (const auto &itrb : bRatH) {
            if (itrb.second > 0.0) {
              vTotal += bOIniFirmFirmHoH[i].find(itrb.first)->second;
            }
          }

          const double ratio = (vTotal <= 0) ? 0.0 : remainingAct / vTotal;

          // add to fulfilled ratio
          for (const auto &itrb : bRatH) {
            if (itrb.second > 0) {
              bAdjustH[itrb.first] += ratio;
            }
          }

          break;
        }
      }

      for (const auto &itrb : bAdjustH) {
        const int j = itrb.first;
        // adjusted f
        bOAdjustedFirmFirmHoH[i][j] = itrb.second * bOIniFirmFirmHoH[i][j];
        // adjusted b
        fOAdjustedFirmFirmHoH[j][i] = itrb.second * bOIniFirmFirmHoH[i][j];
      }
    } else {
      exit(1);
    }
  }
}



void sim_predefinedDelta(int t) {

  // predefine check
  // note that predefine is ignored in the recovery

  for (auto &itr : predefinedHoL) {

    int id = itr.first;
    auto &predefineL = itr.second;

    std::vector<int> eraseIndex;

    if(predefineL.size()==0){

      DeltaVectorH[id]=0;
      piniCapVectorH[id] = piniVectorH[id];
      
    }else{
      for (int i=0;i<predefineL.size();i++){
	int sdate=predefineL[i].GetS();
	int edate=predefineL[i].GetE();
	double delta=predefineL[i].GetD();

	if(sdate<=t && t<=edate){
	  DeltaVectorH[id] = delta;
	  piniCapVectorH[id] = (1.0 - DeltaVectorH[id]) * piniVectorH[id];
	}

	if(t>=edate){
	  eraseIndex.push_back(i);
	}
      }

      // erase reversely
      for(int i=(eraseIndex.size()-1);i>=0;i--){
	predefineL.erase(predefineL.begin()+eraseIndex[i]);
      }
    }
  }

}

void sim_dailyCVector(int t, const firm_index_t &firm_index){

  // daily cvector
  // if it is empty, it is not considered
  // only when it has value, the values rewrite the cvector TENTATIVELY

  if(iDailyCVector==0){
    // no daily delta
    return;
  }else{
    // daily delta on

      auto found = dailyCVectorDayH.find(t);
      if (found == dailyCVectorDayH.end()) {

	// t not found
	cVectorH=cVectorOrgH;
	return;

      }

      // t found
      std::string ts = std::to_string( t );
      string dailyCVectorFileT = sDailyCVectorFileBase + ts + ".txt";

      cVectorH=cVectorOrgH;

      // debug
      //      cerr << dailyCVectorFileT << endl;
      
      //dailycvectorfile
      std::ifstream ifsCVec(dailyCVectorFileT);
      if (!ifsCVec) { throw std::runtime_error("failed to open dailyCVectorFile"); }
      while (!ifsCVec.eof()) {
	std::string sx0;
	double x1;
	ifsCVec >> sx0 >> x1;
	if (sx0.empty()) continue;
	auto found = firm_index.find(sx0);
	if (found != firm_index.end()) {
	  const int x0 = found->second;
	  cVectorH[x0] = x1;

	  //	  cerr << sx0 << " found as " << x0 << endl;
	  
	}else{
	  //	  cerr << sx0 << " not found" << endl;
	}
      }

  }
}

void sim_recovery(int t) {
  if (iRecoveryType == -1 || iRecoveryType == 0 || iRecoveryType == 1
      || iRecoveryType == 2 || iRecoveryType == 3 || iRecoveryType == 4 || iRecoveryType == 5) {

    for (size_t i = 0; i < DeltaVectorH.size(); i++) {

      if (DeltaVectorH[i] > 0) {

	// forcibly stop day
	if ( (t - DeltaStartVectorH[i]) >= iForceStopDay) {
	    
	  if (iRecoveryType == -1 || iRecoveryType == 0) {
	    // proportional
	    DeltaVectorH[i] -= recoveryStepH[i];
	  } else if (iRecoveryType == 1) {
	    // uniform
	    DeltaVectorH[i] -= 1.0 / iRecoveryDay;
	  } else if (iRecoveryType == 2) {
	    // suddenly recover
	    DeltaVectorH[i] = 0.0;
	    piniCapVectorH[i] = (1.0 - DeltaVectorH[i]) * piniVectorH[i];
	  } else if (iRecoveryType == 3) {
	    // renewal equation
	    double beta = 1.0 / recoveryStepH[i];
	    DeltaVectorH[i] *= (1.0 - beta);
	  } else if (iRecoveryType == 4) {
	    // exponential
	    DeltaVectorH[i] = hiddenDeltaVectorH[i] * exp(-(double) t / iRecoveryDay);
	  } else if (iRecoveryType == 5) {
	    // damped renewal equation
	    // forward
	    int degree = 0;
	    double deltaSum = 0.0;
	    degree += fOFirmFirmHoH[i].size();
	    for (const auto &itr : fOFirmFirmHoH[i]) {
	      deltaSum += DeltaVectorH[itr.first];
	    }
	    // backward
	    degree += bOFirmFirmHoH[i].size();
	    for (const auto &itr : bOFirmFirmHoH[i]) {
	      deltaSum += DeltaVectorH[itr.first];
	    }
	    const double h = (degree == 0) ? 1.0 : 1.0 - deltaSum / degree;
	    if (h != 0) {
	      // damping factor
	      double beta = 1.0 / iRecoveryDay;
	      beta *= h;
	      DeltaVectorH[i] *= (1.0 - beta);
	    }
	  }

	  if (DeltaVectorH[i] < 0.0) {
	    DeltaVectorH[i] = 0.0;
	  }
	  piniCapVectorH[i] = (1.0 - DeltaVectorH[i]) * piniVectorH[i];
	}
      }
    }
  }
}

void sim_PrintOut(int t,
                  std::ofstream &pactOfs,
                  std::ofstream &pvalueOfs,
                  std::ofstream &fcOfs,
                  std::ofstream &dlossOfs,
                  std::ofstream &debugOfs,
                  std::vector<std::ofstream> &selectOutPActL,
                  std::vector<std::ofstream> &selectOutVAL) {
  // output value added
  double totalValueAdded = 0.0;
  double totalPactual = 0.0;
  for (size_t i = 0; i < PactualVectorH.size(); i++) {
    totalPactual += PactualVectorH[i];
    totalValueAdded += valueAddedVectorH[i] * PactualVectorH[i];
  }
  lastTotalPactual = totalPactual;
  lastTotalValueAdded = totalValueAdded;

  // pact
  if(!iPActSupFlag){
    pactOfs << t << " " << totalPactual << endl;
  }
  // value added
  pvalueOfs << t << " " << setprecision(17) << totalValueAdded << endl;

  // selective
  if (iSelectiveOutput) {
    for (int i = 0; i < selectOutPActL.size(); i++) {

      double totalPactual = 0.0;
      double totalVA = 0.0;
      for (auto j : selectiveOutputLoH[i]) {
        int id = j.first;
        totalPactual += PactualVectorH[id];
        totalVA += valueAddedVectorH[id] * PactualVectorH[id];
      }

      // disaster happens at 0, so as a pre-disaster, t=-1 is output
      selectOutPActL[i] << t << " " << setprecision(17) << totalPactual << endl;
      selectOutVAL[i] << t << " " << setprecision(17) << totalVA << endl;

    }
  }


  // output final consumption
  if (fcOfs) {
    double totalFC = 0.0;
    for (double x: cAdjustedVectorH) {
      totalFC += x;
    }
    fcOfs << t << " " << setprecision(17) << totalFC << endl;
  }

  // for debug
  //  if (iFullDebugFlag) {
  //    fullDebug(t, debugOfs);
  //  }

  // total direct production loss
  if (dlossOfs) {
    printTotalDirectProductionLoss(t, dlossOfs);
  }
}

void sim_setDelta(int t) {
  
  // input damage
  for (size_t i = 0; i < piniVectorH.size(); i++) {
    if (DeltaReserveVectorH[i]>=0) {
      // defined

      if(DeltaStartVectorH[i]==t){
	// at t

	DeltaVectorH[i]=DeltaReserveVectorH[i];
	
	// pinicapvector is used instead of pinivector because piniCapVector might have over production
	piniCapVectorH[i] = piniVectorH[i] * (1.0 - DeltaVectorH[i]);


      }
    }

  }
}

void simulation(int simulationStep,
                std::ofstream &pactOfs,
                std::ofstream &pvalueOfs,
                std::ofstream &fcOfs,
                std::ofstream &statsOfs,
                std::ofstream &dlossOfs,
                std::ofstream &debugOfs,
                std::vector<std::ofstream> &selectOutPActL,
                std::vector<std::ofstream> &selectOutVAL,
		const firm_index_t &firm_index){
  // simulation start
  // all variables are already initialized
  for (int t = 0; t < simulationStep; t++) {
    std::stringstream ss;

    ss << "time: " << t << std::endl;
    muteSwitchedCerr(ss.str());

    MeasureElapsed("sim_predefinedDelta", [&] {
      sim_predefinedDelta(t);
    });

    MeasureElapsed("sim_setDelta", [&] {
      sim_setDelta(t);
    });

    MeasureElapsed("sim_dailyCVector", [&] {
      sim_dailyCVector(t, firm_index);
    });
    
    MeasureElapsed("sim_inventoryRenewal", [&] {
      sim_inventoryRenewal();
    });

    MeasureElapsed("sim_order", [&] {
      sim_order();
    });

    MeasureElapsed("sim_demand", [&] {
      sim_demand();
    });

    MeasureElapsed("sim_PAct", [&] {
      sim_PAct();
    });

    MeasureElapsed("sim_PrintOut", [&] {
      sim_PrintOut(t, pactOfs, pvalueOfs, fcOfs, dlossOfs, debugOfs, selectOutPActL, selectOutVAL);
    });

    MeasureElapsed("sim_recovery", [&] {
      sim_recovery(t);
    });

    MeasureElapsed("printSnapshot", [&] {
      printSnapshot(t);
    });
  } // for t, simulation step

  // output stat
  if (statsOfs) {
    statsOfs << "totalPactutalEnd " << setprecision(17) << lastTotalPactual
             << endl;  // [TODO] stop writing duplicate values
    statsOfs << "totalValueAddedEnd " << setprecision(17) << lastTotalValueAdded << endl;
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

void OpenOutputStreams(const std::string &pactOutputFile,
                       const std::string &pValueAddedOutputFile,
                       const std::string &finalConsumerOutputFile,
                       const std::string &statsOutputFile,
                       const std::string &directLossOutputFile,
                       const std::vector<std::string> &selectOutInStrL,
                       std::ofstream &pactOfs,
                       std::ofstream &pvOfs,
                       std::ofstream &fcOfs,
                       std::ofstream &statsOfs,
                       std::ofstream &dlossOfs,
                       std::vector<std::ofstream> &selectOutPActL,
                       std::vector<std::ofstream> &selectOutVAL
) {

  if(!iPActSupFlag){
    if (!pactOutputFile.empty()) {
      pactOfs.open(pactOutputFile);
      if (!pactOfs) {
	throw std::runtime_error("unable to open pactOutputFile");
      }
    }
  }

  if (!pValueAddedOutputFile.empty()) {
    pvOfs.open(pValueAddedOutputFile);
    if (!pvOfs) {
      throw std::runtime_error("unable to open pValueAddedOutputFile");
    }
  }

  if (!finalConsumerOutputFile.empty()) {
    fcOfs.open(finalConsumerOutputFile);
    if (!fcOfs) {
      throw std::runtime_error("unable to open finalConsumerOutputFile");
    }
  }

  if (!statsOutputFile.empty()) {
    statsOfs.open(statsOutputFile);
    if (!statsOfs) {
      throw std::runtime_error("unable to open statsOutputFile");
    }
  }

  if (!directLossOutputFile.empty()) {
    dlossOfs.open(directLossOutputFile);
    if (!dlossOfs) {
      throw std::runtime_error("unable to open directLossOutputFile");
    }
  }

  std::vector<string> selectOutPActStrL;
  std::vector<string> selectOutVAStrL;

  selectOutPActL.resize(selectOutInStrL.size());
  selectOutVAL.resize(selectOutInStrL.size());

  // rename input file name to output file name and open output file stream
  for (int i = 0; i < selectOutInStrL.size(); i++) {
    char pact_name[255] = {'\0'};
    sprintf(pact_name, "Selective_PAct_%02d.txt", i);
    selectOutPActL[i].open(pact_name);
    if (!selectOutPActL[i]) {
      cerr << "Unable to open " << pact_name << " as selective output file" << endl;
      exit(0);
    }

    char va_name[255] = {'\0'};
    sprintf(va_name, "Selective_VA_%02d.txt", i);
    selectOutVAL[i].open(va_name);
    if (!selectOutVAL[i]) {
      cerr << "Unable to open " << va_name << " as selective output file" << endl;
      exit(0);
    }
  }
}

void ResizeAll() {
  assert(NUM_NODES > 0);
  // resize all
  fATableHoH.resize(NUM_NODES);
  cVectorH.resize(NUM_NODES, 0.0);
  piniVectorH.resize(NUM_NODES, 0.0);
  piniInputVectorH.resize(NUM_NODES, 0.0);
  valueAddedVectorH.resize(NUM_NODES, 0.0);
  firmAffiH.resize(NUM_NODES, 0);
  firmH.resize(NUM_NODES, 0);
  fOFirmFirmHoH.resize(NUM_NODES);
  fOAdjustedFirmFirmHoH.resize(NUM_NODES);
  bOFirmFirmHoH.resize(NUM_NODES);
  bOIniFirmFirmHoH.resize(NUM_NODES);
  bOAdjustedFirmFirmHoH.resize(NUM_NODES);
  DFirmVectorH.resize(NUM_NODES, 0.0);
  preDFirmVectorH.resize(NUM_NODES, 0.0);
  DAdjustedFirmVectorH.resize(NUM_NODES, 0.0);
  cAdjustedVectorH.resize(NUM_NODES, 0.0);
  piniCapVectorH.resize(NUM_NODES, 0.0);
  SFirmFirmHoH.resize(NUM_NODES);
  SFirmSectorHoH.resize(NUM_NODES);
  AFirmSectorHoH.resize(NUM_NODES);
  PpropFirmSectorHoH.resize(NUM_NODES);
  PiniMaxVectorH.resize(NUM_NODES, 0.0);
  MinPropAmountVectorH.resize(NUM_NODES, 0.0);
  PactualVectorH.resize(NUM_NODES, 0.0);
  DeltaVectorH.resize(NUM_NODES, 0.0);
  DeltaReserveVectorH.resize(NUM_NODES, -1);  
  DeltaStartVectorH.resize(NUM_NODES, 0.0);
  tgtInvVectorH.resize(NUM_NODES, 0.0);
  recoveryStepH.resize(NUM_NODES, 0.0);
  hiddenDeltaVectorH.resize(NUM_NODES, 0.0);

  // now it is meaningless
  firmH.resize(NUM_NODES);
  for (int i = 0; i < firmH.size(); i++) {
    firmH[i] = i;
  }
}

void ReadAtable(const std::string &ioLineFile,  const firm_index_t &firm_index) {

  ifstream ifs(ioLineFile);
  if (!ifs) {
    throw std::runtime_error(std::string("unable to open ") + ioLineFile);
  }

  std::unordered_map<int, std::vector<std::pair<int, double>>> vLinks;
  while (!ifs.eof()) {
    std::string x0, x1;
    double w;
    ifs >> x0 >> x1 >> w;
    if (x0.empty()) continue;
    if (firm_index.find(x0) == firm_index.end()) { continue; }
    int i0 = firm_index.at(x0);
    if (firm_index.find(x1) == firm_index.end()) { continue; }
    int i1 = firm_index.at(x1);

    vLinks[i1].push_back(std::make_pair(i0, w));
  }

  for (const auto &vi : vLinks) {
    const size_t i = vi.first;
    fATableHoH[i].insert(vi.second.cbegin(), vi.second.cend());
  }

}

void ReadCVector(const std::string &cVectorFile, const firm_index_t &firm_index) {
  ifstream ifs2(cVectorFile);
  if (!ifs2) {
    throw std::runtime_error("unable to open cVectorFile");
  }
  while (!ifs2.eof()) {
    std::string s0;
    double c;
    ifs2 >> s0 >> c;
    if (s0.empty()) continue;
    int x0 = firm_index.at(s0);
    cVectorH[x0] = c;
  }
}

void ReadPiniVector(const std::string &PiniVectorFile, const firm_index_t &firm_index) {

  // debug
  //  cerr << PiniVectorFile << endl;

  ifstream ifs3(PiniVectorFile);
  if (!ifs3) {
    throw std::runtime_error("unable to open PiniVectorFile");
  }
  while (!ifs3.eof()) {
    std::string s0;
    double pini;
    ifs3 >> s0 >> pini;
    if (s0.empty()) continue;

    // debug
    //    cerr << s0 << endl;
    
    int x0 = firm_index.at(s0);
    piniVectorH[x0] = pini;
  }
}

firm_index_t ReadFirmAffiliationFile(const std::string &firmAffiFile) {
  
  ifstream ifs4(firmAffiFile);
  if (!ifs4) {
    throw std::runtime_error("unable to open firmAffiFile");
  }

  int idx = 0;
  firm_index_t firm_index;
  unordered_map<int, int> tempIndexFirmH;
  while (!ifs4.eof()) {
    std::string sfirm;
    int sector;
    string x2;
    ifs4 >> sfirm >> sector >> x2;

    // debug
    //    cout << "sfirm " << sfirm << " sector " << sector << " sales " << x2 << endl;

    if (sfirm.empty()) continue;
    const auto found0 = firm_index.insert(std::make_pair(sfirm, idx));
    if (found0.second) { idx++; }

    int i0 = firm_index[sfirm];
    rev_firm_index[i0]=sfirm;

    tempIndexFirmH[i0]=sector;

    if (secFirmHoL.find(sector) == secFirmHoL.end()) {
      secFirmHoL[sector].clear();
    }
    secFirmHoL[sector].push_back(i0);
  }
  
  NUM_NODES = idx;

  ResizeAll();

  for(int i=0;i<NUM_NODES;i++){
    int sector=tempIndexFirmH[i];
    firmAffiH[i] = sector;
  }
  
  return std::move(firm_index);
}

void ReadSelectiveOutputFile(const std::vector<string> &selectOutInStrL, const firm_index_t &firm_index) {

  selectiveOutputLoH.resize(selectOutInStrL.size());

  // open input file stream and read seletive firm
  for (int i = 0; i < selectOutInStrL.size(); i++) {

    string s = selectOutInStrL[i];

    ifstream ifs(s);
    if (!ifs) {
      cerr << "unable to open " << s << " as selective output's input file" << endl;
    }
    while (!ifs.eof()) {
      std::string sfirm;
      ifs >> sfirm;
      if (sfirm.empty()) continue;
      if (firm_index.find(sfirm) == firm_index.end()) { continue; }
      int firm = firm_index.at(sfirm);

      //      selectedIDH[firm]=1;

      selectiveOutputLoH[i][firm] = 1;
    }

    //    selectiveOutputLoH.push_back(selectedIDH);
  }
}

void ReadPredefinedDeltaFile(const std::string &predefinedDeltaFile, const firm_index_t &firm_index) {

  if (predefinedDeltaFile.size() > 0) {


    // input
    ifstream pdIfs;
    pdIfs.open(predefinedDeltaFile);
    if (!pdIfs) {
      throw std::runtime_error("unable to open predefinedDeltaFile");
    }
    // deploy predefined input
    string str;
    while (getline(pdIfs, str)) {
      string tmp;
      istringstream stream(str);
      vector<string> result;
      while (getline(stream, tmp, ' ')) {
        result.push_back(tmp);
      }
      if (result.size() != 4) {
        throw std::runtime_error("illegal format in PredefinedDeltaFile");
      }

      if (firm_index.find(result[0]) == firm_index.end()) { continue; }
      
      int id = firm_index.at(result[0]);
      double delta = atof(result[1].c_str());
      int sdate = atoi(result[2].c_str());
      int edate = atoi(result[3].c_str());

      // debug
      //      cerr << id << " " << delta << " " << sdate << " " << edate << endl;

      if (predefinedHoL.find(id) == predefinedHoL.end()) {
        predefinedHoL[id].clear();
      }
      
      predefinedHoL[id].emplace_back(sdate,edate,delta);
      
    }
  }

  
}

int main(int argc, char *argv[]) {

  // rand seed
  int simulationStep = 10;
  int randSeed = 331;

  time_t btime;
  time_t etime;

  // start time
  time(&btime);

  // options
  int ch;
  extern char *optarg;
  extern int optind, opterr;

  char *proName = argv[0];
  // minimum input files
  int argNum = 4;
  std::string pValueAddedOutputFile = "VA.txt";
  std::string pactOutputFile = "PAct.txt";
  std::string deltaFile;  // [optional] delta vector
  std::string finalConsumeOutputFile;
  std::string statsOutputFile;
  std::string directLossOutputFile;
  std::string sPredefinedDeltaFile;

  while ((ch = getopt(argc, argv, "b:c:d:D:e:f:Fi:mM:o:O:p:P:r:R:s:S:t:T:uv:x:X:z:Z:")) != -1) {
    switch (ch) {
      case 'b':iRecoveryType = atoi(optarg);
        break;
      case 'c':iRecoveryDay = atoi(optarg);
        break;
      case 'd':deltaFile = optarg;
	iDeltaFile=1;
        break;
      case 'D':iInvDistType = atoi(optarg);
        break;
      case 'e':iInventoryMin = atoi(optarg);
        break;
      case 'f':finalConsumeOutputFile = optarg;
        break;
      case 'F':iFullDebugFlag = 1;
        break;
      case 'r':randSeed = atoi(optarg);
        break;
      case 't':simulationStep = atoi(optarg);
        break;
      case 'i':inventoryN = atoi(optarg);
        break;
      case 'm':iMuteFlag = 1;
        break;
      case 'M':iDailyCVector = 1;
        sDailyCVectorFile = optarg;
        break;
      case 'o':iOrderTypeFlag = atoi(optarg);
        break;
      case 'p':iPrintSnapFlag = 1;
        sPrintSnapStep = optarg;
        break;
      case 'R':iRationingTypeFlag = atoi(optarg);
        break;
      case 's':statsOutputFile = optarg;
        break;
      case 'S':iForceStopDay = atoi(optarg);
        break;
      case 'T':tau = atoi(optarg);
        break;
      case 'u':iPActSupFlag = 1;
        break;
      case 'v':pValueAddedOutputFile = optarg;
        break;
      case 'x':directLossOutputFile = optarg;
        break;
      case 'X':pactOutputFile = optarg;
        break;
      case 'z':iSelectiveOutput = 1;
        sSelectiveOutput = optarg;
        break;
      case 'Z':sPredefinedDeltaFile = optarg;
	iPredefinedDeltaFile=1;
        break;
      default:usage(proName);
    }
  }
  argc -= optind;
  argv += optind;

  if (iRecoveryType != -1) {
    if (!(iRecoveryDay > 0)) {
      cerr << "-b should be used with -c (Positive)" << endl;
      exit(1);
    }
  }

  // initialize random
  srand(randSeed);

  // argc
  if (argc != argNum) {
    cerr << "Invalid arg num" << endl;
    usage(proName);
    return 1;
  }

  const std::string ATableFile = argv[argc - 4];
  const std::string cVectorFile = argv[argc - 3];
  const std::string piniVectorFile = argv[argc - 2];
  const std::string firmAffiFile = argv[argc - 1];

  std::ofstream pactOfs;  // output stream for Pactual
  std::ofstream pvalueOfs;  // output stream for Pvalue
  std::ofstream fcOfs; // output stream for final consumer
  std::ofstream statsOfs;
  std::ofstream dlossOfs;  // output stream for direct loss
  std::ofstream debugOfs;

  std::vector<string> selectOutInStrL;
  std::vector<std::ofstream> selectOutPActL;
  std::vector<std::ofstream> selectOutVAL;

  // get selective output file
  if (iSelectiveOutput == 1) {

    vector<string> vS;
    vS = split(sSelectiveOutput, ':');

    if (vS.size() < 1) {
      cerr << "Invalid arg for option z (should be file:file...)" << endl;
      usage(proName);
      exit(1);
    }

    for (auto itr = vS.begin(); itr != vS.end(); itr++) {
      selectOutInStrL.push_back(*itr);
    }

  }

  OpenOutputStreams(pactOutputFile,
                    pValueAddedOutputFile,
                    finalConsumeOutputFile,
                    statsOutputFile,
                    directLossOutputFile,
                    selectOutInStrL,
                    pactOfs,
                    pvalueOfs,
                    fcOfs,
                    statsOfs,
                    dlossOfs,
                    selectOutPActL,
                    selectOutVAL
  );


  // get snapshot info
  if (iPrintSnapFlag) {
    vector<string> vS;
    vS = split(sPrintSnapStep, ':');
    if (vS.size() <= 1) {
      cerr << "Invalid arg for option p (should be file:step:...)" << endl;
      usage(proName);
      return (1);
    }
    sPrintSnapFileBase = vS[0];
    auto itr = vS.begin();
    itr++;
    for (; itr != vS.end(); ++itr) {
      hPrintSnapStep[atoi((*itr).c_str())] = 1;
    }
  }

  // recvoery day should be greater than 0
  if (iRecoveryDay < 0) {
    cerr << "recovery day should be greater than 0" << endl;
    return (1);
  }
  if (iFullDebugFlag) {
    debugOfs.open("FullDebugLog.csv");
    if (!debugOfs) {
      throw std::runtime_error("cannot open FullDebug file");
    }
    debugOfs << "t, i, pini, delta, pcap, d, pmax, pact, dadjusted, s..., o..., " << endl;
  }

  // read firm affiliation
  stringstream ss;

  ss.str("");
  ss.clear(stringstream::goodbit);
  ss << "reading firmAffiliation. File: " << firmAffiFile << std::endl;
  firm_index_t firm_index;
  muteSwitchedCerr(ss.str());
  MeasureElapsed("read firmAffiliation", [&] {
    firm_index = ReadFirmAffiliationFile(firmAffiFile);
  });

  // read ATable
  ss.str("");
  ss << "reading ATableLine. File: " << ATableFile << std::endl;
  muteSwitchedCerr(ss.str());
  MeasureElapsed("read ATableLine", [&] {
    ReadAtable(ATableFile, firm_index);
  });

  // read cVector
  ss.str("");
  ss.clear(stringstream::goodbit);
  ss << "reading cVector. File: " << cVectorFile << std::endl;
  muteSwitchedCerr(ss.str());
  MeasureElapsed("read cVector", [&] {
    ReadCVector(cVectorFile, firm_index);
  });

  // read pini
  ss.str("");
  ss.clear(stringstream::goodbit);
  ss << "reading Pini. File: " << piniVectorFile << std::endl;
  muteSwitchedCerr(ss.str());
  MeasureElapsed("read PiniVector", [&] {
    ReadPiniVector(piniVectorFile, firm_index);
  });
  
  // read selective output file
  if (iSelectiveOutput == 1) {
    ss.str("");
    ss.clear(stringstream::goodbit);
    ss << "reading Selective Output Files" << std::endl;
    muteSwitchedCerr(ss.str());
    ReadSelectiveOutputFile(selectOutInStrL, firm_index);
  }

  // get predefined delta table
  ss.str("");
  ss.clear(stringstream::goodbit);
  ss << "reading PredefinedDelta" << std::endl;
  muteSwitchedCerr(ss.str());
  ReadPredefinedDeltaFile(sPredefinedDeltaFile, firm_index);


  //deltafile
  if (!deltaFile.empty()) {
    std::ifstream ifsDelta(deltaFile);
    if (!ifsDelta) { throw std::runtime_error("failed to open deltaFile"); }
    string str;
    while (getline(ifsDelta,str)){
      std::vector<std::string> result; 
      std::istringstream iss(str); 
      for(std::string s; iss >> s; ){
	result.push_back(s);
      }

      if (result.empty()) continue;
      
      std::string x0=result[0];
      double x1=std::stod(result[1]);
      int x2=0;

      // allow column 3 is empty, which means start at 0
      if(result.size()==2){
      }else{
	x2=std::stoi(result[2]);
      }

      auto found = firm_index.find(x0);
      if (found != firm_index.end()) {
        const int x0 = found->second;
        DeltaReserveVectorH[x0] = x1;
	DeltaStartVectorH[x0] = x2;


      }

      
    }
  }

  // daily CVector
  if (iDailyCVector == 1) {

    vector<string> vS;
    vS = split(sDailyCVectorFile, ':');

    if (vS.size() < 1) {
      cerr << "Invalid arg for option M (should be file:day...)" << endl;
      usage(proName);
      exit(1);
    }

    auto itr = vS.begin();
    sDailyCVectorFileBase=*itr;

    for (auto itr = vS.begin()+1; itr != vS.end(); itr++) {
      dailyCVectorDayH[atoi((*itr).c_str())]=1;
    }

  }

  // check consistency of options

  if((iPredefinedDeltaFile == 1) && ((iRecoveryType != -1) || (iRecoveryDay != 0) || (iDeltaFile==1))){
      cerr << "Illegal option: -Z and (-b or -c or -d) cannot be used at the same time" << endl;
      exit(1);
  }

  
  // set cvec org
  cVectorOrgH=cVectorH;
  
  ss.str("");
  ss.clear(stringstream::goodbit);
  ss << "initialize model" << std::endl;
  muteSwitchedCerr(ss.str());
  MeasureElapsed("initializeModel", [&] {
    initializeModel(randSeed, pactOfs, pvalueOfs, fcOfs, statsOfs, selectOutPActL, selectOutVAL);
  });

  // hazard t=0
  hazard(deltaFile, statsOfs);
  
  if (iFullDebugFlag) {
    fullDebug(-1, debugOfs);
  }
  
  // simulation
  MeasureElapsed("simulation", [&] {
    simulation(simulationStep,
               pactOfs,
               pvalueOfs,
               fcOfs,
               statsOfs,
               dlossOfs,
               debugOfs,
               selectOutPActL,
               selectOutVAL,
	       firm_index);
  });

  // end time
  time(&etime);
  int diffTime = difftime(etime, btime);

  if (!iMuteFlag) {
    printf("%d sec. elapsed\n", diffTime);
  }
  return 0;
}
