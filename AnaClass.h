//////////////////////////////////////////////////////////
// Simple Analysis Class, Lanzhou (2020)
// J.C. Zamora, University of Sao Paulo
// cardona@if.usp.br
//////////////////////////////////////////////////////////

#ifndef AnaClass_h
#define AnaClass_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


#include <TCanvas.h>
#include <TChain.h>
#include <TCutG.h>
#include "TF1.h"
#include <TF2.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TMultiGraph.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TVector3.h>
#include "TEnv.h"
#include "TSpectrum.h"
#include "TVirtualFitter.h"
#include "TSpectrumFit.h"
#include "TCutG.h"
#include "TApplication.h"


#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>      // std::iota
#include <algorithm>    // std::sort


using namespace std;



#define pi 3.14159265
///=============Mass definition============
//see for inctance https://wwwndc.jaea.go.jp/NuC/
#define aum 931.494043    // in MeV
#define e_mass  548.57990945e-6  // in aum

#define n_amu 1.0086654
#define n_ch 0
#define n_mass n_amu*aum //in MeV

#define H1  1.00782503223 // in aum
#define H1_ch 1
#define H1_mass (H1 -  H1_ch*e_mass )*aum //in MeV

#define He4  4.00260325413 // in aum
#define He4_ch 2
#define He4_mass (He4 -  He4_ch*e_mass )*aum //in MeV


#define Li6  6.01512288741 // in aum
#define Li6_ch 3
#define Li6_mass (Li6 -  Li6_ch*e_mass )*aum //in MeV

#define C11  11.011433611 // in aum
#define C11_ch 6
#define C11_mass (C11 -  C11_ch*e_mass )*aum //in MeV

#define C12  12.0000000000 // in aum
#define C12_ch 6
#define C12_mass (C12 -  C12_ch*e_mass )*aum //in MeV

#define O13  13.024814683 // in aum
#define O13_ch 8
#define O13_mass (O13 -  O13_ch*e_mass )*aum //in MeV


#define Au197  196.966568812 // in aum
#define Au197_ch 79
#define Au197_mass (Au197 -  Au197_ch*e_mass )*aum //in MeV


#define Pb208  207.976651189 // in aum
#define Pb208_ch 82
#define Pb208_mass (Pb208 -  Pb208_ch*e_mass )*aum //in MeV


// Header file for the classes stored in the TTree if any.

class AnaClass {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   void ReadConfigFile(const TString& filename);
   double Correct7BeTotEner(double Etot, int TelNumb);

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           DSSD1X[32];
   Int_t           DSSD1Y[32];
   Int_t           DSSD2X[32];
   Int_t           DSSD2Y[32];
   Int_t           DSSD3X[32];
   Int_t           DSSD3Y[32];
   Int_t           DSSD4X[32];
   Int_t           DSSD4Y[32];
   Int_t           DSSD5X[16];
   Int_t           DSSD5Y[16];
   Int_t           DSSD6X[16];
   Int_t           DSSD6Y[16];
   Float_t         DSSD1X_eng;
   Float_t         DSSD1Y_eng;
   Float_t         DSSD2X_eng;
   Float_t         DSSD2Y_eng;
   Float_t         DSSD3X_eng;
   Float_t         DSSD3Y_eng;
   Float_t         DSSD4X_eng;
   Float_t         DSSD4Y_eng;
   Float_t         DSSD5X_eng;
   Float_t         DSSD5Y_eng;
   Float_t         DSSD6X_eng;
   Float_t         DSSD6Y_eng;
   Float_t         DSSD1i;
   Float_t         DSSD1j;
   Float_t         DSSD2i;
   Float_t         DSSD2j;
   Float_t         DSSD3i;
   Float_t         DSSD3j;
   Float_t         DSSD4i;
   Float_t         DSSD4j;
   Float_t         DSSD5i;
   Float_t         DSSD5j;
   Float_t         DSSD6i;
   Float_t         DSSD6j;
   Float_t         DSSD1X_mul;
   Float_t         DSSD1Y_mul;
   Float_t         DSSD2X_mul;
   Float_t         DSSD2Y_mul;
   Float_t         DSSD3X_mul;
   Float_t         DSSD3Y_mul;
   Float_t         DSSD4X_mul;
   Float_t         DSSD4Y_mul;
   Float_t         DSSD5X_mul;
   Float_t         DSSD5Y_mul;
   Float_t         DSSD6X_mul;
   Float_t         DSSD6Y_mul;
   Float_t         DSSD1X_eng_MeV;
   Float_t         DSSD1Y_eng_MeV;
   Float_t         DSSD2X_eng_MeV;
   Float_t         DSSD2Y_eng_MeV;
   Float_t         DSSD3X_eng_MeV;
   Float_t         DSSD3Y_eng_MeV;
   Float_t         DSSD4X_eng_MeV;
   Float_t         DSSD4Y_eng_MeV;
   Float_t         DSSD5X_eng_MeV;
   Float_t         DSSD5Y_eng_MeV;
   Float_t         DSSD6X_eng_MeV;
   Float_t         DSSD6Y_eng_MeV;
   Int_t           SSD1;
   Int_t           SSD2;
   Int_t           SSD3;
   Int_t           SSD4;
   Float_t         SSD1_MeV;
   Float_t         SSD2_MeV;
   Float_t         SSD3_MeV;
   Float_t         SSD4_MeV;
   Float_t         tel1X_MeV;
   Float_t         tel2X_MeV;
   Float_t         tel3X_MeV;
   Float_t         tel4X_MeV;
   Float_t         tel1Y_MeV;
   Float_t         tel2Y_MeV;
   Float_t         tel3Y_MeV;
   Float_t         tel4Y_MeV;
   Float_t         tof1;
   Float_t         tof2;
   Float_t         tof;
   Float_t         position_1x;
   Float_t         position_1y;
   Float_t         position_1z;
   Float_t         position_2x;
   Float_t         position_2y;
   Float_t         position_2z;
   Float_t         position_3x;
   Float_t         position_3y;
   Float_t         position_3z;
   Float_t         position_4x;
   Float_t         position_4y;
   Float_t         position_4z;
   Float_t         position_5x;
   Float_t         position_5y;
   Float_t         position_5z;
   Float_t         position_6x;
   Float_t         position_6y;
   Float_t         position_6z;
   Float_t         target_x;
   Float_t         target_y;
   Float_t         target_z;

   // List of branches
   TBranch        *b_DSSD1X;   //!
   TBranch        *b_DSSD1Y;   //!
   TBranch        *b_DSSD2X;   //!
   TBranch        *b_DSSD2Y;   //!
   TBranch        *b_DSSD3X;   //!
   TBranch        *b_DSSD3Y;   //!
   TBranch        *b_DSSD4X;   //!
   TBranch        *b_DSSD4Y;   //!
   TBranch        *b_DSSD5X;   //!
   TBranch        *b_DSSD5Y;   //!
   TBranch        *b_DSSD6X;   //!
   TBranch        *b_DSSD6Y;   //!
   TBranch        *b_DSSD1X_eng;   //!
   TBranch        *b_DSSD1Y_eng;   //!
   TBranch        *b_DSSD2X_eng;   //!
   TBranch        *b_DSSD2Y_eng;   //!
   TBranch        *b_DSSD3X_eng;   //!
   TBranch        *b_DSSD3Y_eng;   //!
   TBranch        *b_DSSD4X_eng;   //!
   TBranch        *b_DSSD4Y_eng;   //!
   TBranch        *b_DSSD5X_eng;   //!
   TBranch        *b_DSSD5Y_eng;   //!
   TBranch        *b_DSSD6X_eng;   //!
   TBranch        *b_DSSD6Y_eng;   //!
   TBranch        *b_DSSD1i;   //!
   TBranch        *b_DSSD1j;   //!
   TBranch        *b_DSSD2i;   //!
   TBranch        *b_DSSD2j;   //!
   TBranch        *b_DSSD3i;   //!
   TBranch        *b_DSSD3j;   //!
   TBranch        *b_DSSD4i;   //!
   TBranch        *b_DSSD4j;   //!
   TBranch        *b_DSSD5i;   //!
   TBranch        *b_DSSD5j;   //!
   TBranch        *b_DSSD6i;   //!
   TBranch        *b_DSSD6j;   //!
   TBranch        *b_DSSD1X_mul;   //!
   TBranch        *b_DSSD1Y_mul;   //!
   TBranch        *b_DSSD2X_mul;   //!
   TBranch        *b_DSSD2Y_mul;   //!
   TBranch        *b_DSSD3X_mul;   //!
   TBranch        *b_DSSD3Y_mul;   //!
   TBranch        *b_DSSD4X_mul;   //!
   TBranch        *b_DSSD4Y_mul;   //!
   TBranch        *b_DSSD5X_mul;   //!
   TBranch        *b_DSSD5Y_mul;   //!
   TBranch        *b_DSSD6X_mul;   //!
   TBranch        *b_DSSD6Y_mul;   //!
   TBranch        *b_DSSD1X_eng_MeV;   //!
   TBranch        *b_DSSD1Y_eng_MeV;   //!
   TBranch        *b_DSSD2X_eng_MeV;   //!
   TBranch        *b_DSSD2Y_eng_MeV;   //!
   TBranch        *b_DSSD3X_eng_MeV;   //!
   TBranch        *b_DSSD3Y_eng_MeV;   //!
   TBranch        *b_DSSD4X_eng_MeV;   //!
   TBranch        *b_DSSD4Y_eng_MeV;   //!
   TBranch        *b_DSSD5X_eng_MeV;   //!
   TBranch        *b_DSSD5Y_eng_MeV;   //!
   TBranch        *b_DSSD6X_eng_MeV;   //!
   TBranch        *b_DSSD6Y_eng_MeV;   //!
   TBranch        *b_SSD1;   //!
   TBranch        *b_SSD2;   //!
   TBranch        *b_SSD3;   //!
   TBranch        *b_SSD4;   //!
   TBranch        *b_SSD1_MeV;   //!
   TBranch        *b_SSD2_MeV;   //!
   TBranch        *b_SSD3_MeV;   //!
   TBranch        *b_SSD4_MeV;   //!
   TBranch        *b_tel1X_MeV;   //!
   TBranch        *b_tel2X_MeV;   //!
   TBranch        *b_tel3X_MeV;   //!
   TBranch        *b_tel4X_MeV;   //!
   TBranch        *b_tel1Y_MeV;   //!
   TBranch        *b_tel2Y_MeV;   //!
   TBranch        *b_tel3Y_MeV;   //!
   TBranch        *b_tel4Y_MeV;   //!
   TBranch        *b_tof1;   //!
   TBranch        *b_tof2;   //!
   TBranch        *b_tof;   //!
   TBranch        *b_position_1x;   //!
   TBranch        *b_position_1y;   //!
   TBranch        *b_position_1z;   //!
   TBranch        *b_position_2x;   //!
   TBranch        *b_position_2y;   //!
   TBranch        *b_position_2z;   //!
   TBranch        *b_position_3x;   //!
   TBranch        *b_position_3y;   //!
   TBranch        *b_position_3z;   //!
   TBranch        *b_position_4x;   //!
   TBranch        *b_position_4y;   //!
   TBranch        *b_position_4z;   //!
   TBranch        *b_position_5x;   //!
   TBranch        *b_position_5y;   //!
   TBranch        *b_position_5z;   //!
   TBranch        *b_position_6x;   //!
   TBranch        *b_position_6y;   //!
   TBranch        *b_position_6z;   //!
   TBranch        *b_target_x;   //!
   TBranch        *b_target_y;   //!
   TBranch        *b_target_z;   //!



   TFile *file;
   TFile *cutfile;
   TCutG *CUTTEST;
   TH1F* hDSSD1X[32];
   TH1F* hDSSD1Y[32];
   TH1F* hDSSD2X[32];
   TH1F* hDSSD2Y[32];
   TH1F* hDSSD3X[32];
   TH1F* hDSSD3Y[32];
   TH1F* hDSSD4X[32];
   TH1F* hDSSD4Y[32];
   TH1F* hSiAX[16];
   TH1F* hSiAY[16];
   TH1F* hSiBX[16];
   TH1F* hSiBY[16];
   TH1F* hDSSD1X_Cal[32];
   TH1F* hDSSD1Y_Cal[32];
   TH1F* hDSSD2X_Cal[32];
   TH1F* hDSSD2Y_Cal[32];
   TH1F* hDSSD3X_Cal[32];
   TH1F* hDSSD3Y_Cal[32];
   TH1F* hDSSD4X_Cal[32];
   TH1F* hDSSD4Y_Cal[32];
   TH1F* hSSD1_Cal;
   TH1F* hSSD2_Cal;
   TH1F* hSSD3_Cal;
   TH1F* hSSD4_Cal;
   TH1F* htel1X_Cal;
   TH1F* htel1Y_Cal;
   TH1F* htel2X_Cal;
   TH1F* htel2Y_Cal;
   TH1F* htel3X_Cal;
   TH1F* htel3Y_Cal;
   TH1F* htel4X_Cal;
   TH1F* htel4Y_Cal;
   TH2F* hPID1;
   TH2F* hPID2;
   TH2F* hPID3;
   TH2F* hPID4;
   TVector3 vDSSD[4];
   TVector3 vSi[2];
   TVector3 vTarg;
   TH2F* posDSSD[4];
   TH2F* posSi[2];
   TH1F* angDSSD[4];
   TH1F* angBeam;
   TH1F* angCM[4];
   TH1F* ExRec[4];

   AnaClass(TTree *tree=0);
   virtual ~AnaClass();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(const TString& st);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   virtual double omega(double x, double y, double z);
   double* kine_2b(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject);
};

#endif

#ifdef AnaClass_cxx
AnaClass::AnaClass(TTree *tree) : fChain(0)
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("c11exp",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("c11exp","");
      chain->Add("11C.root/c11exp");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

AnaClass::~AnaClass()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AnaClass::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AnaClass::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void AnaClass::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("DSSD1X", DSSD1X, &b_DSSD1X);
   fChain->SetBranchAddress("DSSD1Y", DSSD1Y, &b_DSSD1Y);
   fChain->SetBranchAddress("DSSD2X", DSSD2X, &b_DSSD2X);
   fChain->SetBranchAddress("DSSD2Y", DSSD2Y, &b_DSSD2Y);
   fChain->SetBranchAddress("DSSD3X", DSSD3X, &b_DSSD3X);
   fChain->SetBranchAddress("DSSD3Y", DSSD3Y, &b_DSSD3Y);
   fChain->SetBranchAddress("DSSD4X", DSSD4X, &b_DSSD4X);
   fChain->SetBranchAddress("DSSD4Y", DSSD4Y, &b_DSSD4Y);
   fChain->SetBranchAddress("DSSD5X", DSSD5X, &b_DSSD5X);
   fChain->SetBranchAddress("DSSD5Y", DSSD5Y, &b_DSSD5Y);
   fChain->SetBranchAddress("DSSD6X", DSSD6X, &b_DSSD6X);
   fChain->SetBranchAddress("DSSD6Y", DSSD6Y, &b_DSSD6Y);
   fChain->SetBranchAddress("DSSD1X_eng", &DSSD1X_eng, &b_DSSD1X_eng);
   fChain->SetBranchAddress("DSSD1Y_eng", &DSSD1Y_eng, &b_DSSD1Y_eng);
   fChain->SetBranchAddress("DSSD2X_eng", &DSSD2X_eng, &b_DSSD2X_eng);
   fChain->SetBranchAddress("DSSD2Y_eng", &DSSD2Y_eng, &b_DSSD2Y_eng);
   fChain->SetBranchAddress("DSSD3X_eng", &DSSD3X_eng, &b_DSSD3X_eng);
   fChain->SetBranchAddress("DSSD3Y_eng", &DSSD3Y_eng, &b_DSSD3Y_eng);
   fChain->SetBranchAddress("DSSD4X_eng", &DSSD4X_eng, &b_DSSD4X_eng);
   fChain->SetBranchAddress("DSSD4Y_eng", &DSSD4Y_eng, &b_DSSD4Y_eng);
   fChain->SetBranchAddress("DSSD5X_eng", &DSSD5X_eng, &b_DSSD5X_eng);
   fChain->SetBranchAddress("DSSD5Y_eng", &DSSD5Y_eng, &b_DSSD5Y_eng);
   fChain->SetBranchAddress("DSSD6X_eng", &DSSD6X_eng, &b_DSSD6X_eng);
   fChain->SetBranchAddress("DSSD6Y_eng", &DSSD6Y_eng, &b_DSSD6Y_eng);
   fChain->SetBranchAddress("DSSD1i", &DSSD1i, &b_DSSD1i);
   fChain->SetBranchAddress("DSSD1j", &DSSD1j, &b_DSSD1j);
   fChain->SetBranchAddress("DSSD2i", &DSSD2i, &b_DSSD2i);
   fChain->SetBranchAddress("DSSD2j", &DSSD2j, &b_DSSD2j);
   fChain->SetBranchAddress("DSSD3i", &DSSD3i, &b_DSSD3i);
   fChain->SetBranchAddress("DSSD3j", &DSSD3j, &b_DSSD3j);
   fChain->SetBranchAddress("DSSD4i", &DSSD4i, &b_DSSD4i);
   fChain->SetBranchAddress("DSSD4j", &DSSD4j, &b_DSSD4j);
   fChain->SetBranchAddress("DSSD5i", &DSSD5i, &b_DSSD5i);
   fChain->SetBranchAddress("DSSD5j", &DSSD5j, &b_DSSD5j);
   fChain->SetBranchAddress("DSSD6i", &DSSD6i, &b_DSSD6i);
   fChain->SetBranchAddress("DSSD6j", &DSSD6j, &b_DSSD6j);
   fChain->SetBranchAddress("DSSD1X_mul", &DSSD1X_mul, &b_DSSD1X_mul);
   fChain->SetBranchAddress("DSSD1Y_mul", &DSSD1Y_mul, &b_DSSD1Y_mul);
   fChain->SetBranchAddress("DSSD2X_mul", &DSSD2X_mul, &b_DSSD2X_mul);
   fChain->SetBranchAddress("DSSD2Y_mul", &DSSD2Y_mul, &b_DSSD2Y_mul);
   fChain->SetBranchAddress("DSSD3X_mul", &DSSD3X_mul, &b_DSSD3X_mul);
   fChain->SetBranchAddress("DSSD3Y_mul", &DSSD3Y_mul, &b_DSSD3Y_mul);
   fChain->SetBranchAddress("DSSD4X_mul", &DSSD4X_mul, &b_DSSD4X_mul);
   fChain->SetBranchAddress("DSSD4Y_mul", &DSSD4Y_mul, &b_DSSD4Y_mul);
   fChain->SetBranchAddress("DSSD5X_mul", &DSSD5X_mul, &b_DSSD5X_mul);
   fChain->SetBranchAddress("DSSD5Y_mul", &DSSD5Y_mul, &b_DSSD5Y_mul);
   fChain->SetBranchAddress("DSSD6X_mul", &DSSD6X_mul, &b_DSSD6X_mul);
   fChain->SetBranchAddress("DSSD6Y_mul", &DSSD6Y_mul, &b_DSSD6Y_mul);
   fChain->SetBranchAddress("DSSD1X_eng_MeV", &DSSD1X_eng_MeV, &b_DSSD1X_eng_MeV);
   fChain->SetBranchAddress("DSSD1Y_eng_MeV", &DSSD1Y_eng_MeV, &b_DSSD1Y_eng_MeV);
   fChain->SetBranchAddress("DSSD2X_eng_MeV", &DSSD2X_eng_MeV, &b_DSSD2X_eng_MeV);
   fChain->SetBranchAddress("DSSD2Y_eng_MeV", &DSSD2Y_eng_MeV, &b_DSSD2Y_eng_MeV);
   fChain->SetBranchAddress("DSSD3X_eng_MeV", &DSSD3X_eng_MeV, &b_DSSD3X_eng_MeV);
   fChain->SetBranchAddress("DSSD3Y_eng_MeV", &DSSD3Y_eng_MeV, &b_DSSD3Y_eng_MeV);
   fChain->SetBranchAddress("DSSD4X_eng_MeV", &DSSD4X_eng_MeV, &b_DSSD4X_eng_MeV);
   fChain->SetBranchAddress("DSSD4Y_eng_MeV", &DSSD4Y_eng_MeV, &b_DSSD4Y_eng_MeV);
   fChain->SetBranchAddress("DSSD5X_eng_MeV", &DSSD5X_eng_MeV, &b_DSSD5X_eng_MeV);
   fChain->SetBranchAddress("DSSD5Y_eng_MeV", &DSSD5Y_eng_MeV, &b_DSSD5Y_eng_MeV);
   fChain->SetBranchAddress("DSSD6X_eng_MeV", &DSSD6X_eng_MeV, &b_DSSD6X_eng_MeV);
   fChain->SetBranchAddress("DSSD6Y_eng_MeV", &DSSD6Y_eng_MeV, &b_DSSD6Y_eng_MeV);
   fChain->SetBranchAddress("SSD1", &SSD1, &b_SSD1);
   fChain->SetBranchAddress("SSD2", &SSD2, &b_SSD2);
   fChain->SetBranchAddress("SSD3", &SSD3, &b_SSD3);
   fChain->SetBranchAddress("SSD4", &SSD4, &b_SSD4);
   fChain->SetBranchAddress("SSD1_MeV", &SSD1_MeV, &b_SSD1_MeV);
   fChain->SetBranchAddress("SSD2_MeV", &SSD2_MeV, &b_SSD2_MeV);
   fChain->SetBranchAddress("SSD3_MeV", &SSD3_MeV, &b_SSD3_MeV);
   fChain->SetBranchAddress("SSD4_MeV", &SSD4_MeV, &b_SSD4_MeV);
   fChain->SetBranchAddress("tel1X_MeV", &tel1X_MeV, &b_tel1X_MeV);
   fChain->SetBranchAddress("tel2X_MeV", &tel2X_MeV, &b_tel2X_MeV);
   fChain->SetBranchAddress("tel3X_MeV", &tel3X_MeV, &b_tel3X_MeV);
   fChain->SetBranchAddress("tel4X_MeV", &tel4X_MeV, &b_tel4X_MeV);
   fChain->SetBranchAddress("tel1Y_MeV", &tel1Y_MeV, &b_tel1Y_MeV);
   fChain->SetBranchAddress("tel2Y_MeV", &tel2Y_MeV, &b_tel2Y_MeV);
   fChain->SetBranchAddress("tel3Y_MeV", &tel3Y_MeV, &b_tel3Y_MeV);
   fChain->SetBranchAddress("tel4Y_MeV", &tel4Y_MeV, &b_tel4Y_MeV);
   fChain->SetBranchAddress("tof1", &tof1, &b_tof1);
   fChain->SetBranchAddress("tof2", &tof2, &b_tof2);
   fChain->SetBranchAddress("tof", &tof, &b_tof);
   fChain->SetBranchAddress("position_1x", &position_1x, &b_position_1x);
   fChain->SetBranchAddress("position_1y", &position_1y, &b_position_1y);
   fChain->SetBranchAddress("position_1z", &position_1z, &b_position_1z);
   fChain->SetBranchAddress("position_2x", &position_2x, &b_position_2x);
   fChain->SetBranchAddress("position_2y", &position_2y, &b_position_2y);
   fChain->SetBranchAddress("position_2z", &position_2z, &b_position_2z);
   fChain->SetBranchAddress("position_3x", &position_3x, &b_position_3x);
   fChain->SetBranchAddress("position_3y", &position_3y, &b_position_3y);
   fChain->SetBranchAddress("position_3z", &position_3z, &b_position_3z);
   fChain->SetBranchAddress("position_4x", &position_4x, &b_position_4x);
   fChain->SetBranchAddress("position_4y", &position_4y, &b_position_4y);
   fChain->SetBranchAddress("position_4z", &position_4z, &b_position_4z);
   fChain->SetBranchAddress("position_5x", &position_5x, &b_position_5x);
   fChain->SetBranchAddress("position_5y", &position_5y, &b_position_5y);
   fChain->SetBranchAddress("position_5z", &position_5z, &b_position_5z);
   fChain->SetBranchAddress("position_6x", &position_6x, &b_position_6x);
   fChain->SetBranchAddress("position_6y", &position_6y, &b_position_6y);
   fChain->SetBranchAddress("position_6z", &position_6z, &b_position_6z);
   fChain->SetBranchAddress("target_x", &target_x, &b_target_x);
   fChain->SetBranchAddress("target_y", &target_y, &b_target_y);
   fChain->SetBranchAddress("target_z", &target_z, &b_target_z);
   Notify();
}

Bool_t AnaClass::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AnaClass::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AnaClass::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


///=============Two Body Kinematics===========
//Kinematics based on the Mandelstam invariant variables
double AnaClass::omega(double x, double y, double z){
  return sqrt(x*x + y*y + z*z -2*x*y -2*y*z -2*x*z);
}
double* AnaClass::kine_2b(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject){

  //in this definition: m1(projectile); m2(target); m3(ejectile); and m4(recoil);
 double Et1 = K_proj + m1;
 double Et2 = m2;
 double Et3 = K_eject + m3;
 double Et4  = Et1 + Et2 - Et3;
 double m4_ex, Ex, theta_cm;
 double s,t,u; //---Mandelstam variables
 double p1, p3;
 double J_LtoCM; //jacobian Lab to CM

 s = pow(m1,2) + pow(m2,2) +2*m2*Et1;
 u = pow(m2,2) + pow(m3,2) - 2*m2*Et3;

 m4_ex = sqrt(  (cos(thetalab) * omega(s,pow(m1,2),pow(m2,2)) * omega(u,pow(m2,2),pow(m3,2)) - (s - pow(m1,2) - pow(m2,2))*(pow(m2,2) + pow(m3,2) - u) )/(2*pow(m2,2)) + s + u - pow(m2,2)  );
 Ex = m4_ex - m4;

 t =   pow(m2,2) + pow(m4_ex,2) - 2*m2*Et4;

 //for normal kinematics
 theta_cm = acos( ( pow(s,2) +s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) + (pow(m1,2) - pow(m2,2))*(pow(m3,2) - pow(m4_ex,2)) )/( omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4_ex,2))) ) ;

 //for inverse kinematics Note: this angle corresponds to the recoil
 //theta_cm = TMath::Pi() - acos( ( pow(s,2) +s*(2*t - pow(m1,2) - pow(m2,2) - pow(m3,2) - pow(m4_ex,2)) + (pow(m1,2) - pow(m2,2))*(pow(m3,2) - pow(m4_ex,2)) )/( omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4_ex,2))) ) ;

 p1= sqrt(pow(Et1,2)-pow(m1,2));
 p3 = sqrt(pow(Et3,2)-pow(m3,2));

 J_LtoCM = abs( ((omega(s,pow(m1,2),pow(m2,2))*omega(s,pow(m3,2),pow(m4,2)))/(4*s*p1*p3))*(1.+Et1/m2 - cos(thetalab)*(Et3*p1)/(m2*p3)) );


 static double output[3];
 output[0]= theta_cm;
 output[1]= Ex;
 output[2]= J_LtoCM;
 return output;

}
///==================================

#endif // #ifdef AnaClass_cxx
