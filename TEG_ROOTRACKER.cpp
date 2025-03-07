#ifndef TEG_ROOTRACKER
#define TEG_ROOTRACKER 1

#define CPLUS_INCLUDE_PATH $ROOTSYS
#include <fstream>
#include <string>
#include <vector>

#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "CLI11.hpp"

std::string InputFile = "";
std::string OutputFile = "TEG_rootracker.gtrac.root";

class Particle: public TObject {
  public:
    int    fEvent;
    int    fPDG;
    int    fState;
    double fpx;
    double fpy;
    double fpz;
    double fE;
    double fm;
};


std::vector<Particle> GetData(std::string InputFile) {
  
  std::vector<Particle> part;
  double PDG, in_out_state, px, py, pz, E, m;
  int EventCount = 0;

  std::string line;
  std::ifstream inFile;
  inFile.open(InputFile.c_str());

  for(std::string line; getline(inFile, line); ) { 
    //if line = <event>, increase event no. by 1, read in next 4 lines
    if (line == "<event>") {
      EventCount += 1;
      std::string lines[4];

      for (int i = 0; i < 4; i++) {
        getline(inFile, lines[i]);
	Particle c;
        
	std::stringstream ss(lines[i]);
        while(ss >> PDG >> in_out_state >> px >> py >> pz >> E >> m){
          c.fEvent     = EventCount;
          c.fPDG       = PDG;
          c.fState     = in_out_state;
          c.fpx        = px;
          c.fpy        = py;
          c.fpz        = pz;
          c.fE         = E;
          c.fm         = m;
          part.push_back(c);
        }
      }
    } //event loop
  }

  return part;
}


void ConvertToRooTracker(std::vector<Particle> data) {

  //find no. of events
  int nmax = data.size();
  int nmax_evt = nmax / 4;
  std::cout << "Converting " << nmax << " particles in " << nmax_evt << " events to rootracker format." << std::endl;

  int kNPmax = 250;                    

  int    brEvtNum;                     //Event number
  double brEvtXSec;                    //Cross section for selected event
  double brEvtDXSec;                   //Cross section for selected event kinematics
  UInt_t brEvtKPS;                     //Kinematic phase space variables as in KinePhaseSpace_t
  double brEvtWght;                    //Event weight
  double brEvtProb;                    //Probability for that event (given cross section, path lengths, etc)
  double brEvtVtx[4];                  //Event vertex position in detector co-ordinate system
  int    brStdHepN;                    //number of particles in particle array

  int    brStdHepPdg    [kNPmax];      //PDG codes
  int    brStdHepStatus [kNPmax];      //Generator-specific status code
  int    brStdHepRescat [kNPmax];      //Hadron transport model - specific rescattering code
  double brStdHepX4     [kNPmax][4];   //4-x (x,y,z,t) of particle in hit nucleus frame (fm)
  double brStdHepP4     [kNPmax][4];   //4-p (px,py,pz,E) of particle in LAB frame (GeV)
  double brStdHepPolz   [kNPmax][3];   //polarisation vector
  int    brStdHepFd     [kNPmax] ;     //first daughter
  int    brStdHepLd     [kNPmax] ;     //last daughter
  int    brStdHepFm     [kNPmax] ;     //first mother
  int    brStdHepLm     [kNPmax] ;     //last mother

  //open the output root file & tree
  TFile fout(OutputFile.c_str(), "RECREATE");
  TTree * rootracker_tree = new TTree("gRooTracker", "GENIE event tree rootracker format");

  //create output root tree branches
  rootracker_tree->Branch("EvtNum",      &brEvtNum,       "EvtNum/I");
  rootracker_tree->Branch("EvtXSec",     &brEvtXSec,      "EvtXSec/D");
  rootracker_tree->Branch("EvtDXSec",    &brEvtDXSec,     "EvtDXSec/D");
  rootracker_tree->Branch("EvtKPS",      &brEvtKPS,       "EvtKPS/i");
  rootracker_tree->Branch("EvtWght",     &brEvtWght,      "EvtWght/D");
  rootracker_tree->Branch("EvtProb",     &brEvtProb,      "EvtProb/D");
  rootracker_tree->Branch("EvtVtx",       brEvtVtx,       "EvtVtx[4]/D");
  rootracker_tree->Branch("StdHepN",     &brStdHepN,      "StdHepN/I");
  rootracker_tree->Branch("StdHepPdg",    brStdHepPdg,    "StdHepPdg[StdHepN]/I");
  rootracker_tree->Branch("StdHepStatus", brStdHepStatus, "StdHepStatus[StdHepN]/I");
  rootracker_tree->Branch("StdHepRescat", brStdHepRescat, "StdHepRescat[StdHepN]/I");
  rootracker_tree->Branch("StdHepX4",     brStdHepX4,     "StdHepX4[StdHepN][4]/D");
  rootracker_tree->Branch("StdHepP4",     brStdHepP4,     "StdHepP4[StdHepN][4]/D");
  rootracker_tree->Branch("StdHelPolz",   brStdHepPolz,   "StdHepPolz[StdHepN][3]/D");
  rootracker_tree->Branch("StdHepFd",     brStdHepFd,     "StdHepFd[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepLd",     brStdHepLd,     "StdHepLd[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepFm",     brStdHepFm,     "StdHepFm[StdHepN]/I"); 
  rootracker_tree->Branch("StdHepLm",     brStdHepLm,     "StdHepLm[StdHepN]/I"); 

  //loop through particles
  int PartCount = 0;
  for (std::vector<Particle>::iterator t=data.begin(); t!=data.end(); ++t) {
   
    Particle c = *t;
    PartCount += 1;

    //clear output tree branches
    brEvtNum   = 0;
    brEvtXSec  = 0;
    brEvtDXSec = 0;
    brEvtKPS   = 0;
    brEvtWght  = 0;

    for(int k; k<4; k++) {
      brEvtVtx[k] = 0;
    }

    brStdHepN = 0;

    for(int i=0; i<kNPmax; i++) {
      brStdHepPdg    [i]  = 0;
      brStdHepStatus [i] = -1;
      brStdHepRescat [i] = -1;
      
      for(int k=0; k<4; k++) {
        brStdHepX4 [i][k] = 0;
        brStdHepP4 [i][k] = 0;
      }

      for(int k=0; k<3; k++) {
        brStdHepPolz [i][k] = 0;
      }

      brStdHepFd    [i] = 0;
      brStdHepLd    [i] = 0;
      brStdHepFm    [i] = 0;
      brStdHepLm    [i] = 0;
    }


    //copy current event info to output tree
    brEvtNum    = c.fEvent;
    brStdHepPdg    [PartCount] = c.fPDG;
    brStdHepStatus [PartCount] = c.fState;
    brStdHepP4[PartCount][0] = c.fpx;
    brStdHepP4[PartCount][1] = c.fpy;
    brStdHepP4[PartCount][2] = c.fpz;
    brStdHepP4[PartCount][3] = c.fE;
    
    brStdHepN = PartCount;

    //fill tree
    rootracker_tree->Fill();
  } //event loop

  fout.Write();
  fout.Close();
}


int main(int argc, char** argv) {

  CLI::App app{"A program to convert the output of trident event generator into rootracker format"};
  app.add_option("-i, --input",     InputFile,   "Input file name (text file)")->required();
  app.add_option("-o, --output",    OutputFile,  "Output file name (gtrac.root)");
  CLI11_PARSE(app, argc, argv);
   
  std::vector<Particle> data = GetData(InputFile);
  ConvertToRooTracker(data);

}

#endif
