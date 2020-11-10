#include "TFile.h"
#include "TCutG.h"

void savecut(void){

   TFile *fcoup=new TFile("cuttest.root","recreate");
   fcoup->cd();
   gROOT->FindObject("cuttest")->Write();
   
   delete fcoup; 

}

 
