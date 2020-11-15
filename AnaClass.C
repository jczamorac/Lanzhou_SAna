#define AnaClass_cxx
#include "AnaClass.h"



TChain* MakeChain(const TString& str2);

int main(int argc, char** argv) {


   TApplication *a = new TApplication("a", 0, 0);

  TString str(argv[1]);

  TChain *chain = MakeChain( str);

  AnaClass t(chain);

  t.ReadConfigFile(str);

  t.Loop(str);

  return 0;
}


TChain* MakeChain(const TString& str2) {
  auto *chain = new TChain("c11exp");



   //TString PathToFiles = "/home/juan/proyectos/tandar_6li_19f/Li6F19_Tandar_2019/2019-07-6Li/spm/";
  //chain->Add(PathToFiles + str2);
   chain->Add(str2);
   cout<<str2<<endl;



  return chain;
}





void AnaClass::ReadConfigFile(const TString& filename){



   cout<<"Reading Configuration File"<<endl;
   TEnv* CFile;
   //CFile = new TEnv("config/par.txt");


   //theta1Z = CFile->GetValue("Theta1Z",-1000.0);
   //theta2Z = CFile->GetValue("Theta2Z",-1000.0);

   delete CFile;

}

double AnaClass::Correct7BeTotEner(double Etot, int TelNumb){

  double A[4]={650.954,	664.522,	655.307,	673.395};
  double B[4]={-8.43807,	-8.65518,	-8.50244,	-8.79981};
  double C[4]={0.0401052,	0.0414934,	0.0404924,	0.0424319};
  double D[4]={-7.38e-05,	-7.70e-05,	-7.47e-05,	-7.91e-05};

  if(Etot>190 || Etot<80 || TelNumb>3) return 0;
  else{

    double Etotcorr = A[TelNumb] + B[TelNumb]*pow(Etot,2.0) + C[TelNumb]*pow(Etot,3.0) + D[TelNumb]*pow(Etot,4.0);
    Etotcorr += Etot;
    return Etotcorr;
  }

}


void AnaClass::Loop(const TString& st)
{

 //----------------------create a new outfile for the respective run number
  string infile(st);
  //string str2("Teflon");
  //string str3("HistosTeflon") ;
  //infile.replace(infile.find(str2),str2.length(),str3);
  string outfile("Histo" + infile);

  //cout<<infile<<endl;


//------------------------------------------


//========Loading the cut

	//cutfile = new TFile("cuts/cuttest.root");		   // open file
	//CUTTEST = (TCutG *)cutfile->Get("cuttest");	   // read CUTEC

	//cutfile19F = new TFile("cutjuan.root");		   // open file
	//CUTALFA19F = (TCutG *)cutfile19F->Get("CUTJUAN");	   // read CUTEC


//=============================


	file = new TFile(outfile.c_str(), "recreate");
	file->mkdir("raw/DSSD1");
  file->mkdir("raw/DSSD2");
  file->mkdir("raw/DSSD3");
  file->mkdir("raw/DSSD4");
  file->mkdir("cal/DSSD1_Cal");
  file->mkdir("cal/DSSD2_Cal");
  file->mkdir("cal/DSSD3_Cal");
  file->mkdir("cal/DSSD4_Cal");
  file->mkdir("raw/SiA");
  file->mkdir("raw/SiB");
	file->mkdir("cal/SSDs_Cal");
  file->mkdir("Tel1");
  file->mkdir("Tel2");
  file->mkdir("Tel3");
  file->mkdir("Tel4");
  file->mkdir("pos");
	file->mkdir("Ex");




  for(int i=0; i<32; i++){
	   hDSSD1X[i] = new TH1F(Form("DSSD1X_%d",i),Form("DSSD1X_%d",i),512,0,0);
     hDSSD1Y[i] = new TH1F(Form("DSSD1Y_%d",i),Form("DSSD1Y_%d",i),512,0,0);
     hDSSD2X[i] = new TH1F(Form("DSSD2X_%d",i),Form("DSSD2X_%d",i),512,0,0);
     hDSSD2Y[i] = new TH1F(Form("DSSD2Y_%d",i),Form("DSSD2Y_%d",i),512,0,0);
     hDSSD3X[i] = new TH1F(Form("DSSD3X_%d",i),Form("DSSD3X_%d",i),512,0,0);
     hDSSD3Y[i] = new TH1F(Form("DSSD3Y_%d",i),Form("DSSD3Y_%d",i),512,0,0);
     hDSSD4X[i] = new TH1F(Form("DSSD4X_%d",i),Form("DSSD4X_%d",i),512,0,0);
     hDSSD4Y[i] = new TH1F(Form("DSSD4Y_%d",i),Form("DSSD4Y_%d",i),512,0,0);

     hDSSD1X_Cal[i] = new TH1F(Form("DSSD1X_Cal_%d",i),Form("DSSD1X_Cal_%d",i),512,0,0);
     hDSSD1Y_Cal[i] = new TH1F(Form("DSSD1Y_Cal_%d",i),Form("DSSD1Y_Cal_%d",i),512,0,0);
     hDSSD2X_Cal[i] = new TH1F(Form("DSSD2X_Cal_%d",i),Form("DSSD2X_Cal_%d",i),512,0,0);
     hDSSD2Y_Cal[i] = new TH1F(Form("DSSD2Y_Cal_%d",i),Form("DSSD2Y_Cal_%d",i),512,0,0);
     hDSSD3X_Cal[i] = new TH1F(Form("DSSD3X_Cal_%d",i),Form("DSSD3X_Cal_%d",i),512,0,0);
     hDSSD3Y_Cal[i] = new TH1F(Form("DSSD3Y_Cal_%d",i),Form("DSSD3Y_Cal_%d",i),512,0,0);
     hDSSD4X_Cal[i] = new TH1F(Form("DSSD4X_Cal_%d",i),Form("DSSD4X_Cal_%d",i),512,0,0);
     hDSSD4Y_Cal[i] = new TH1F(Form("DSSD4Y_Cal_%d",i),Form("DSSD4Y_Cal_%d",i),512,0,0);
  }

  for(int i=0; i<16; i++){
     hSiAX[i] = new TH1F(Form("SiAX_%d",i),Form("SiAX_%d",i),512,0,0);
     hSiAY[i] = new TH1F(Form("SiAY_%d",i),Form("SiAY_%d",i),512,0,0);
     hSiBX[i] = new TH1F(Form("SiBX_%d",i),Form("SiBX_%d",i),512,0,0);
     hSiBY[i] = new TH1F(Form("SiBY_%d",i),Form("SiBY_%d",i),512,0,0);
   }

   hSSD1_Cal = new TH1F("SSD1_Cal","SSD1_Cal",512,0,0);
   hSSD2_Cal = new TH1F("SSD2_Cal","SSD2_Cal",512,0,0);
   hSSD3_Cal = new TH1F("SSD3_Cal","SSD3_Cal",512,0,0);
   hSSD4_Cal = new TH1F("SSD4_Cal","SSD4_Cal",512,0,0);
   htel1X_Cal = new TH1F("tel1X_Cal","tel1X_Cal",512,0,0);
   htel1Y_Cal = new TH1F("tel1Y_Cal","tel1Y_Cal",512,0,0);
   htel2X_Cal = new TH1F("tel2X_Cal","tel2X_Cal",512,0,0);
   htel2Y_Cal = new TH1F("tel2Y_Cal","tel2Y_Cal",512,0,0);
   htel3X_Cal = new TH1F("tel3X_Cal","tel3X_Cal",512,0,0);
   htel3Y_Cal = new TH1F("tel3Y_Cal","tel3Y_Cal",512,0,0);
   htel4X_Cal = new TH1F("tel4X_Cal","tel4X_Cal",512,0,0);
   htel4Y_Cal = new TH1F("tel4Y_Cal","tel4Y_Cal",512,0,0);
   hPID1 = new TH2F("dE-E_1","dE-E_1",512,0,0,512,0,0);
   hPID2 = new TH2F("dE-E_2","dE-E_2",512,0,0,512,0,0);
   hPID3 = new TH2F("dE-E_3","dE-E_3",512,0,0,512,0,0);
   hPID4 = new TH2F("dE-E_4","dE-E_4",512,0,0,512,0,0);
   angBeam = new TH1F("angBeam","angBeam",512,0,0);
   for(int i=0;i<4;i++){
      posDSSD[i] = new TH2F(Form("posDSSD_%d",i),Form("posDSSD_%d",i),512,0,0,512,0,0);
      angDSSD[i] = new TH1F(Form("angDSSD_%d",i),Form("angDSSD_%d",i),512,0,0);
      angCM[i] = new TH1F(Form("angCM_%d",i),Form("angCM_%d",i),512,0,0);
      ExRec[i] = new TH1F(Form("ExRec_%d",i),Form("ExRec_%d",i),512,0,0);
    }
   for(int i=0;i<2;i++){
      posSi[i] = new TH2F(Form("posSi_%d",i),Form("posSi_%d",i),512,0,0,512,0,0);
    }

	if (fChain == 0) return;

   	Long64_t nentries = fChain->GetEntriesFast();

	Long64_t nbytes = 0, nb = 0;
   	for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //for (Long64_t jentry=0; jentry<500;jentry++) {
      		Long64_t ientry = LoadTree(jentry);
      		if (ientry < 0) break;
      		nb = fChain->GetEntry(jentry);   nbytes += nb;




          //=======creating  few histograms
           for(int i=0; i<32; i++){
         	   if(DSSD1X_mul>0) hDSSD1X[i]->Fill(DSSD1X[i]);
             if(DSSD1Y_mul>0) hDSSD1Y[i]->Fill(DSSD1Y[i]);
             if(DSSD2X_mul>0) hDSSD2X[i]->Fill(DSSD2X[i]);
             if(DSSD2Y_mul>0) hDSSD2Y[i]->Fill(DSSD2Y[i]);
             if(DSSD3X_mul>0) hDSSD3X[i]->Fill(DSSD3X[i]);
             if(DSSD3Y_mul>0) hDSSD3Y[i]->Fill(DSSD3Y[i]);
             if(DSSD4X_mul>0) hDSSD4X[i]->Fill(DSSD4X[i]);
             if(DSSD4Y_mul>0) hDSSD4Y[i]->Fill(DSSD4Y[i]);
           }

           for(int i=0; i<16; i++){
         	   if(DSSD5X_mul>0) hSiAX[i]->Fill(DSSD5X[i]);
             if(DSSD5Y_mul>0) hSiAY[i]->Fill(DSSD5Y[i]);
             if(DSSD6X_mul>0) hSiBX[i]->Fill(DSSD6X[i]);
             if(DSSD6Y_mul>0) hSiBY[i]->Fill(DSSD6Y[i]);
           }

           int DSSD1i_int = (int)DSSD1i;
           if(DSSD1X_mul>0) hDSSD1X_Cal[DSSD1i_int]->Fill(DSSD1X_eng_MeV);
           int DSSD1j_int = (int)DSSD1j;
           if(DSSD1Y_mul>0) hDSSD1Y_Cal[DSSD1j_int]->Fill(DSSD1Y_eng_MeV);
           int DSSD2i_int = (int)DSSD2i;
           if(DSSD2X_mul>0) hDSSD2X_Cal[DSSD2i_int]->Fill(DSSD2X_eng_MeV);
           int DSSD2j_int = (int)DSSD2j;
           if(DSSD2Y_mul>0) hDSSD2Y_Cal[DSSD2j_int]->Fill(DSSD2Y_eng_MeV);
           int DSSD3i_int = (int)DSSD3i;
           if(DSSD3X_mul>0) hDSSD3X_Cal[DSSD3i_int]->Fill(DSSD3X_eng_MeV);
           int DSSD3j_int = (int)DSSD3j;
           if(DSSD3Y_mul>0) hDSSD3Y_Cal[DSSD3j_int]->Fill(DSSD3Y_eng_MeV);
           int DSSD4i_int = (int)DSSD4i;
           if(DSSD4X_mul>0) hDSSD4X_Cal[DSSD4i_int]->Fill(DSSD4X_eng_MeV);
           int DSSD4j_int = (int)DSSD4j;
           if(DSSD4Y_mul>0) hDSSD4Y_Cal[DSSD4j_int]->Fill(DSSD4Y_eng_MeV);

           if(DSSD1X_mul>0 && DSSD1Y_mul>0) hSSD1_Cal->Fill(SSD1_MeV);
           if(DSSD2X_mul>0 && DSSD2Y_mul>0) hSSD2_Cal->Fill(SSD2_MeV);
           if(DSSD3X_mul>0 && DSSD3Y_mul>0) hSSD3_Cal->Fill(SSD3_MeV);
           if(DSSD4X_mul>0 && DSSD4Y_mul>0) hSSD4_Cal->Fill(SSD4_MeV);

           if(DSSD1X_mul>0) htel1X_Cal->Fill(tel1X_MeV);
           if(DSSD1Y_mul>0) htel1Y_Cal->Fill(tel1Y_MeV);
           if(DSSD2X_mul>0) htel2X_Cal->Fill(tel2X_MeV);
           if(DSSD2Y_mul>0) htel2Y_Cal->Fill(tel2Y_MeV);
           if(DSSD3X_mul>0) htel3X_Cal->Fill(tel3X_MeV);
           if(DSSD3Y_mul>0) htel3Y_Cal->Fill(tel3Y_MeV);
           if(DSSD4X_mul>0) htel4X_Cal->Fill(tel4X_MeV);
           if(DSSD4Y_mul>0) htel4Y_Cal->Fill(tel4Y_MeV);

           if(DSSD1X_mul>0) hPID1->Fill(tel1X_MeV,DSSD1X_eng_MeV);
           if(DSSD2X_mul>0) hPID2->Fill(tel2X_MeV,DSSD2X_eng_MeV);
           if(DSSD3X_mul>0) hPID3->Fill(tel3X_MeV,DSSD3X_eng_MeV);
           if(DSSD4X_mul>0) hPID4->Fill(tel4X_MeV,DSSD4X_eng_MeV);


           //=======Filling  vectors
           vDSSD[0].SetXYZ(position_1x,position_1y,position_1z);
           vDSSD[1].SetXYZ(position_2x,position_2y,position_2z);
           vDSSD[2].SetXYZ(position_3x,position_3y,position_3z);
           vDSSD[3].SetXYZ(position_4x,position_4y,position_4z);
           vSi[0].SetXYZ(position_5x,position_5y,position_5z);
           vSi[1].SetXYZ(position_6x,position_6y,position_6z);
           vector<int> vMultX{DSSD1X_mul,DSSD2X_mul,DSSD3X_mul,DSSD4X_mul,DSSD5X_mul,DSSD6X_mul};
           vTarg.SetXYZ(target_x,target_y,target_z);
           vector<double> vK_ejec{tel1X_MeV,tel2X_MeV,tel3X_MeV,tel4X_MeV};

           for(int i=0;i<4;i++)
            if(vMultX[i]>0) posDSSD[i]->Fill(vDSSD[i].X(),vDSSD[i].Y());

           for(int i=0;i<2;i++)
             if(vMultX[i+4]>0) posSi[i]->Fill(vSi[i].X(),vSi[i].Y());

          TVector3 vBeam = vSi[0] - vSi[1];
          TVector3 vZaxis(0,0,1);
          if(vMultX[4]>0 && vMultX[5]>0) angBeam->Fill(vBeam.Angle(vZaxis)*180.0/3.1415);

          //=======Loop over DSSD hits
          //if(!CUTTEST->IsInside(tel3X_MeV,DSSD3X_eng_MeV)) continue;
          for(int i=0;i<4;i++){
           if(vMultX[i]>0){
              TVector3 vRdssd = vDSSD[i] - vTarg;
              //Lab angle of the reaction
              Float_t Rang = vRdssd.Angle(vBeam);
              angDSSD[i]->Fill(Rang*180.0/3.1415);

              //---Invariant mass kinematics
              double m1 = C11_mass;
            	double m2 = Pb208_mass;
            	double m3 = C11_mass;
            	double m4 = Pb208_mass;
            	double ebeam = 275.0 - 2.3;

              double *kinematics;
              //-----------------------------mass_proj, mass_target, mass_ejectile, mass_recoil, K_proj, thetalab, K_ejectile //-------All in MeV and theta in rad
              kinematics = kine_2b(C11_mass, Pb208_mass , C11_mass, Pb208_mass, ebeam, Rang, vK_ejec[i]);
              double theta_cm = *(kinematics+0);
              double Ex_recoil = *(kinematics+1);

              angCM[i]->Fill(theta_cm*180./3.1415);
              ExRec[i]->Fill(Ex_recoil);

              kinematics = NULL;
            }
         }


         /*Correction of 7Be punching through
            It is necessary to create a TCutG file for each telescope
            only selecting the region to be corrected*/
           //double corr7Be_tel1X_MeV, corr7Be_tel2X_MeV, corr7Be_tel3X_MeV, corr7Be_tel4X_MeV;
           //if(CUT7Be1->IsInside(SSD1_MeV,DSSD1X_eng_MeV)) corr7Be_tel1X_MeV = Correct7BeTotEner(tel1X_MeV, 0);
           //if(CUT7Be2->IsInside(SSD2_MeV,DSSD2X_eng_MeV)) corr7Be_tel2X_MeV = Correct7BeTotEner(tel2X_MeV, 1);
           //if(CUT7Be3->IsInside(SSD3_MeV,DSSD3X_eng_MeV)) corr7Be_tel3X_MeV = Correct7BeTotEner(tel3X_MeV, 2);
           //if(CUT7Be4->IsInside(SSD4_MeV,DSSD4X_eng_MeV)) corr7Be_tel4X_MeV = Correct7BeTotEner(tel4X_MeV, 3);






          //=======cleaning the  vectors
           vMultX.clear();
           vK_ejec.clear();
           for(int i=0;i<4;i++) vDSSD[i].Clear();
           for(int i=0;i<2;i++) vSi[i].Clear();

	}// events loop


//=======Saving the histograms
  for(int i=0; i<32; i++){
     file->cd("raw/DSSD1");
	   hDSSD1X[i]->Write();
     hDSSD1Y[i]->Write();
     file->cd("raw/DSSD2");
     hDSSD2X[i]->Write();
     hDSSD2Y[i]->Write();
     file->cd("raw/DSSD3");
     hDSSD3X[i]->Write();
     hDSSD3Y[i]->Write();
     file->cd("raw/DSSD4");
     hDSSD4X[i]->Write();
     hDSSD4Y[i]->Write();

     file->cd("cal/DSSD1_Cal");
     if(hDSSD1X_Cal[i]->GetEntries()>0) hDSSD1X_Cal[i]->Write();
     if(hDSSD1Y_Cal[i]->GetEntries()>0) hDSSD1Y_Cal[i]->Write();
     file->cd("cal/DSSD2_Cal");
     if(hDSSD2X_Cal[i]->GetEntries()>0) hDSSD2X_Cal[i]->Write();
     if(hDSSD2Y_Cal[i]->GetEntries()>0) hDSSD2Y_Cal[i]->Write();
     file->cd("cal/DSSD3_Cal");
     if(hDSSD3X_Cal[i]->GetEntries()>0) hDSSD3X_Cal[i]->Write();
     if(hDSSD3Y_Cal[i]->GetEntries()>0) hDSSD3Y_Cal[i]->Write();
     file->cd("cal/DSSD4_Cal");
     if(hDSSD4X_Cal[i]->GetEntries()>0) hDSSD4X_Cal[i]->Write();
     if(hDSSD4Y_Cal[i]->GetEntries()>0) hDSSD4Y_Cal[i]->Write();
  }

  for(int i=0; i<16; i++){
    file->cd("raw/SiA");
    hSiAX[i]->Write();
    hSiAY[i]->Write();
    file->cd("raw/SiB");
    hSiBX[i]->Write();
    hSiBY[i]->Write();
  }

  file->cd("cal/SSDs_Cal");
  hSSD1_Cal->Write();
  hSSD2_Cal->Write();
  hSSD3_Cal->Write();
  hSSD4_Cal->Write();
  file->cd("Tel1");
  htel1X_Cal->Write();
  htel1Y_Cal->Write();
  hPID1->Write();
  file->cd("Tel2");
  htel2X_Cal->Write();
  htel2Y_Cal->Write();
  hPID2->Write();
  file->cd("Tel3");
  htel3X_Cal->Write();
  htel3Y_Cal->Write();
  hPID3->Write();
  file->cd("Tel4");
  htel4X_Cal->Write();
  htel4Y_Cal->Write();
  hPID4->Write();
  for(int i=0;i<4;i++){
     file->cd("pos");
     posDSSD[i]->Write();
     angDSSD[i]->Write();
     file->cd("Ex");
     ExRec[i]->Write();
     angCM[i]->Write();
   }
  for(int i=0;i<2;i++) posSi[i]->Write();
  angBeam->Write();

	file->Close();



}
