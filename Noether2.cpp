#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TChain.h"

#include <math.h>
#include <string.h>
#include <iostream>
#include "particle.hpp"
using namespace std;

void Noether() {
  double a, c, A = 0, C = 0;
  int countC = 0;
  cout << "Enter value 1: ";
  cin >> a;
  cout << "Enter value 2: ";
  cin >> c;
  TFile* tfile = TFile::Open("BKstTauMu.root");
  TTree* ttree = (TTree*) tfile->Get("BKstTauMuTuple/MCDecayTree");

  // Defining 4-momenta of particles using hpp file.
  Particle< Double_t > Kplus_M(    "Kplus",    ttree );
  Particle< Double_t > piminus_M(  "piminus",  ttree );
  Particle< Double_t > muplus_M(   "muplus",   ttree );
  Particle< Double_t > piplus_M(   "piplus",   ttree );
  Particle< Double_t > piminus0_M( "piminus0", ttree );
  Particle< Double_t > piminus1_M( "piminus1", ttree );

  // hpp file does not currently work for end vertex positions, so old method being used.
  // Defining B origin vertex components
  double B0_TRUEORIGINVERTEX_X, B0_TRUEORIGINVERTEX_Y, B0_TRUEORIGINVERTEX_Z;
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_X",&B0_TRUEORIGINVERTEX_X);
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_Y",&B0_TRUEORIGINVERTEX_Y);
  ttree->SetBranchAddress("B0_TRUEORIGINVERTEX_Z",&B0_TRUEORIGINVERTEX_Z);
  // Defining tauminus origin vertex components
  double tauminus_TRUEORIGINVERTEX_X, tauminus_TRUEORIGINVERTEX_Y, tauminus_TRUEORIGINVERTEX_Z;
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_X", &tauminus_TRUEORIGINVERTEX_X);
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_Y", &tauminus_TRUEORIGINVERTEX_Y);
  ttree->SetBranchAddress("tauminus_TRUEORIGINVERTEX_Z", &tauminus_TRUEORIGINVERTEX_Z);
  // Defining tauminus end vertex components
  double tauminus_TRUEENDVERTEX_X, tauminus_TRUEENDVERTEX_Y, tauminus_TRUEENDVERTEX_Z;
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_X", &tauminus_TRUEENDVERTEX_X);
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_Y", &tauminus_TRUEENDVERTEX_Y);
  ttree->SetBranchAddress("tauminus_TRUEENDVERTEX_Z", &tauminus_TRUEENDVERTEX_Z);
  // x uncertainites widely scattered adjust hist_E0 for better estimate of standard deviation
  TH1D* hist_E0 = new TH1D("B Meson Mass", "Reconstruction of B Meson Mass", 70, 0, 20000);
  TH1D* StdUnc  = new TH1D("StdDev", "StdDev", 70, 0, 20000);
  TH1D* Again   = new TH1D("StdDev", "Again", 70, 0, 20000);

  double *AA = (double*) calloc(10000, sizeof(double));
  double *STD = (double*) calloc(10000, sizeof(double));
  double *CC = (double*) calloc(10000, sizeof(double));

  TRandom3* rand = new TRandom3();

  for(A = 0; A < a;) {
    for(C = 0; C < c;) {
      for(Long64_t i = 0; i < ttree->GetEntries(); i++ ) {
        ttree->GetEntry(i);
        TLorentzVector Kplus_Mvec    = Kplus_M.getVec();
        TLorentzVector piminus_Mvec  = piminus_M.getVec();
        TLorentzVector muplus_Mvec   = muplus_M.getVec();
        TLorentzVector piplus_Mvec   = piplus_M.getVec();
        TLorentzVector piminus0_Mvec = piminus0_M.getVec();
        TLorentzVector piminus1_Mvec = piminus1_M.getVec();

        // Particle 2-4 and 6-8 mass in tau or mode change mass!
        double s_Kplus    = rand->Gaus(493.677, 0.013);      // from B J et al (2012) particle listings
        double s_piminus  = rand->Gaus(139.570018, 0.00035); // from C Amsler et al (2008) particle listings
        double s_muplus   = 105.65837;                     // from Beringer J et al (particle data group) 2012 particle summary
        double s_piplus   = rand->Gaus(139.570018, 0.00035);
        double s_piminus0 = rand->Gaus(139.570018, 0.00035);
        double s_piminus1 = rand->Gaus(139.570018, 0.00035);

        // Adding errors to particles 2-4
        double Kplus_abs             = sqrt(pow(Kplus_Mvec.X(), 2) + pow(Kplus_Mvec.Y(), 2) + pow(Kplus_Mvec.Z(), 2));
        double Kplus_abs_SMEARED     = rand->Gaus(Kplus_abs, 0.005*A*Kplus_abs);
        double K_smear_factor        = Kplus_abs_SMEARED/Kplus_abs;
        double Kplus_X               = (Kplus_Mvec.X()*K_smear_factor);
        double Kplus_Y               = (Kplus_Mvec.Y()*K_smear_factor);
        double Kplus_Z               = (Kplus_Mvec.Z()*K_smear_factor);
        double Kplus_Mvec_E          = sqrt((s_Kplus)*(s_Kplus) + (Kplus_X)*(Kplus_X) + (Kplus_Y)*(Kplus_Y) + (Kplus_Z)*(Kplus_Z));
        
        /*if(i == 0) {
        std::cout << A << std:: endl;
				}*/
                

        double piminus_abs           = sqrt(pow(piminus_Mvec.X(), 2) + pow(piminus_Mvec.Y(), 2) + pow(piminus_Mvec.Z(), 2));
        double piminus_abs_SMEARED   = rand->Gaus(piminus_abs, 0.005*A*piminus_abs);
        double piminus_smear_factor  = piminus_abs_SMEARED/piminus_abs;
        double piminus_X             = (piminus_Mvec.X()*piminus_smear_factor);
        double piminus_Y             = (piminus_Mvec.Y()*piminus_smear_factor);
        double piminus_Z             = (piminus_Mvec.Z()*piminus_smear_factor);
        double piminus_Mvec_E        = sqrt((s_piminus)*(s_piminus) + (piminus_X)*(piminus_X) + (piminus_Y)*(piminus_Y) + (piminus_Z)*(piminus_Z));


        double muplus_abs            = sqrt(pow(muplus_Mvec.X(), 2) + pow(muplus_Mvec.Y(), 2) + pow(muplus_Mvec.Z(), 2));
        double muplus_abs_SMEARED    = rand->Gaus(muplus_abs, 0.005*A*muplus_abs);
        double muplus_smear_factor   = muplus_abs_SMEARED/muplus_abs;
        double muplus_X              = (muplus_Mvec.X()*muplus_smear_factor);
        double muplus_Y              = (muplus_Mvec.Y()*muplus_smear_factor);
        double muplus_Z              = (muplus_Mvec.Z()*muplus_smear_factor);
        double muplus_Mvec_E         = sqrt((s_muplus)*(s_muplus) + (muplus_X)*(muplus_X) + (muplus_Y)*(muplus_Y) + (muplus_Z)*(muplus_Z));


        // Adding errors to particles 6-8

        double piplus_abs            = sqrt(piplus_Mvec.X()*piplus_Mvec.X() + piplus_Mvec.Y()*piplus_Mvec.Y() + piplus_Mvec.Z()*piplus_Mvec.Z());
        double piplus_abs_smeared    = rand->Gaus(piplus_abs, 0.005*A*piplus_abs);
        double piplus_smear_factor   = piplus_abs_smeared/piplus_abs;
        double piplus_X              = (piplus_Mvec.X()*piplus_smear_factor);
        double piplus_Y              = (piplus_Mvec.Y()*piplus_smear_factor);
        double piplus_Z              = (piplus_Mvec.Z()*piplus_smear_factor);
        double piplus_Mvec_E         = sqrt((s_piplus)*(s_piplus) + (piplus_X)*(piplus_X) + (piplus_Y)*(piplus_Y) + (piplus_Z)*(piplus_Z));

                  

        double piminus0_abs          = sqrt(piminus0_Mvec.X()*piminus0_Mvec.X() + piminus0_Mvec.Y()*piminus0_Mvec.Y() + piminus0_Mvec.Z()*piminus0_Mvec.Z());
        double piminus0_abs_smeared  = rand->Gaus(piminus0_abs, 0.005*A*piminus0_abs);
        double piminus0_smear_factor = piminus0_abs_smeared/piminus0_abs;
        double piminus0_X            = (piminus0_Mvec.X()*piminus0_smear_factor);
        double piminus0_Y            = (piminus0_Mvec.Y()*piminus0_smear_factor);
        double piminus0_Z            = (piminus0_Mvec.Z()*piminus0_smear_factor);
        double piminus0_Mvec_E        = sqrt((s_piminus0)*(s_piminus0) + (piminus0_X)*(piminus0_X) + (piminus0_Y)*(piminus0_Y) + (piminus0_Z)*(piminus0_Z));


        double piminus1_abs          = sqrt(piminus1_Mvec.X()*piminus1_Mvec.X() + piminus1_Mvec.Y()*piminus1_Mvec.Y() + piminus1_Mvec.Z()*piminus1_Mvec.Z());
        double piminus1_abs_smeared  = rand->Gaus(piminus1_abs, 0.005*A*piminus0_abs);
        double piminus1_smear_factor = piminus1_abs_smeared/piminus1_abs;
        double piminus1_X            = (piminus1_Mvec.X()*piminus1_smear_factor);
        double piminus1_Y            = (piminus1_Mvec.Y()*piminus1_smear_factor);
        double piminus1_Z            = (piminus1_Mvec.Z()*piminus1_smear_factor);
        double piminus1_Mvec_E       = sqrt((s_piminus1)*(s_piminus1) + (piminus1_X)*(piminus1_X) + (piminus1_Y)*(piminus1_Y) + (piminus1_Z)*(piminus1_Z));


        // Errors on verticies

        double B0_VERTEX_X = rand->Gaus(B0_TRUEORIGINVERTEX_X, 0.04*C*sqrt(3/100));
        double B0_VERTEX_Y = rand->Gaus(B0_TRUEORIGINVERTEX_Y, 0.04*C*sqrt(3/100));
        double tauminus_VERTEX_X = rand->Gaus(tauminus_TRUEORIGINVERTEX_X, 0.04*C);
        double tauminus_VERTEX_Y = rand->Gaus(tauminus_TRUEORIGINVERTEX_Y, 0.04*C);
        double tauminus_VERTEX_Z = rand->Gaus(tauminus_TRUEORIGINVERTEX_Z, 0.2*C);
        double tauminus_ENDVERTEX_X = rand->Gaus(tauminus_TRUEENDVERTEX_X, 0.04*C);
        double tauminus_ENDVERTEX_Y = rand->Gaus(tauminus_TRUEENDVERTEX_Y, 0.04*C);
        double tauminus_ENDVERTEX_Z = rand->Gaus(tauminus_TRUEENDVERTEX_Z, 0.2*C);


        // Main Program.

        // Reconstruction calculations.

        double s5unit_X = tauminus_ENDVERTEX_X - tauminus_VERTEX_X;
        double s5unit_Y = tauminus_ENDVERTEX_Y - tauminus_VERTEX_Y;
        double s5unit_Z = tauminus_ENDVERTEX_Z - tauminus_VERTEX_Z;

        double s1unit_X = tauminus_VERTEX_X - B0_VERTEX_X;
        double s1unit_Y = tauminus_VERTEX_Y - B0_VERTEX_Y;

        double xi = ((Kplus_X + piminus_X + muplus_X)*s1unit_Y - (Kplus_Y + piminus_Y + muplus_Y)*s1unit_X)/((s5unit_Y*s1unit_X) - (s5unit_X*s1unit_Y));

                                

        double p9_X = xi*s5unit_X - (piplus_X + piminus0_X + piminus1_X);
        double p9_Y = xi*s5unit_Y - (piplus_Y + piminus0_Y + piminus1_Y);
        double p9_Z = xi*s5unit_Z - (piplus_Z + piminus0_Z + piminus1_Z);
        double p9_E = sqrt(pow(p9_X, 2) + pow(p9_Y, 2) + pow(p9_Z, 2));

        double s5_0 = pow(p9_E + piplus_Mvec_E + piminus0_Mvec_E + piminus1_Mvec_E, 2);
        double s5_1 = pow(p9_X + piplus_X + piminus0_X + piminus1_X, 2) + pow(p9_Y + piplus_Y + piminus0_Y + piminus1_Y, 2) + pow(p9_Z + piplus_Z + piminus0_Z + piminus1_Z, 2);

        double s5 = sqrt(s5_0-s5_1);

        double s4_0 = pow((sqrt(pow(s5, 2) + (pow(s5unit_X, 2) + pow(s5unit_Y, 2) + pow(s5unit_Z, 2))*pow(xi, 2)) + Kplus_Mvec_E + piminus_Mvec_E +muplus_Mvec_E), 2);
        double s4_1 = pow((s5unit_X*xi + Kplus_X + piminus_X + muplus_X), 2) + pow((s5unit_Y*xi + Kplus_Y + piminus_Y + muplus_Y), 2) + pow((s5unit_Z*xi + Kplus_Z + piminus_Z + muplus_Z), 2);

        double s1 = sqrt((s4_0-s4_1));

        hist_E0->Fill(s1);
      }

      int n = hist_E0->GetStdDev();
      std::cout<< n <<std::endl;
      STD[countC] = n;
      CC[countC] = C;
      AA[countC] = A;
      countC = countC+1;
      C = C + 0.1*c;
      hist_E0->Reset();
    }
    A = A + 0.1*a;
  }

  TCanvas *c1 = new TCanvas("c1","h2smooth",200,10,600,400);
  TGraph2D *g = new TGraph2D(10000,CC,AA,STD);
  gStyle->SetPalette(kBird);
  g->SetTitle("Standard Deviation;;;");
  g->Draw("surf1");
  
  c1->SaveAs("testRead.pdf");
  return;
}
