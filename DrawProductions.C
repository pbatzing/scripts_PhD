//File to collect the results from the different productions
#include "TString.h"
#include "TFile.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TCanvas.h"
#include "iostream"

using namespace std;

void Drawpp(TFile* file, TString dir){
  cout << dir<<endl;
  TDirectory * binstats = file->GetDirectory("bin_stats");
  TCanvas * statcanvas = dynamic_cast<TCanvas*>(binstats->Get("samestatscanvas"));
  statcanvas->Print(Form("%s/1stats.pdf",dir.Data()));
  TFile * stats = TFile::Open(Form("%s/1stats.root",dir.Data()),"RECREATE");
  stats->cd();
  statcanvas->Write("stats");
  stats->Close();
  delete stats; delete statcanvas;
  TDirectory * uncorrecteddir = file->GetDirectory("divided");
  TCanvas * samecanvasPhiPhi = dynamic_cast<TCanvas*>(uncorrecteddir->Get("DPHIDPHI"));	  
  samecanvasPhiPhi->Print(Form("%s/2phiphisame.pdf",dir.Data()));
  TFile * phiphisame = TFile::Open(Form("%s/2phiphisame.root",dir.Data()),"RECREATE");
  phiphisame->cd();
  samecanvasPhiPhi->Write("phiphisame");
  samecanvasPhiPhi->Close();  
  delete samecanvasPhiPhi;delete phiphisame;
  TCanvas * samecanvasPhiEta = dynamic_cast<TCanvas*>(uncorrecteddir->Get("DPhi_1_DEta_12_SameSideCanvas"));	  
  samecanvasPhiEta->Print(Form("%s/2phietasame.pdf",dir.Data()));
  TFile * phietasame = TFile::Open(Form("%s/2phietasame.root",dir.Data()),"RECREATE");
  phietasame->cd();
  samecanvasPhiEta->Write("phietasame");
  samecanvasPhiEta->Close();  
  delete samecanvasPhiEta;delete phietasame; 
  TDirectory * correcteddir = file->GetDirectory("iteration1");
  TCanvas * corcanvasPhiPhi = dynamic_cast<TCanvas*>(correcteddir->Get("DPHIDPHI"));	  
  corcanvasPhiPhi->Print(Form("%s/3phiphicor.pdf",dir.Data()));
  TFile * phiphicor = TFile::Open(Form("%s/3phiphicor.root",dir.Data()),"RECREATE");
  phiphicor->cd();
  corcanvasPhiPhi->Write("phiphicor");
  corcanvasPhiPhi->Close();  
  delete corcanvasPhiPhi;delete phiphicor;
  TCanvas * corcanvasPhiEta = dynamic_cast<TCanvas*>(correcteddir->Get("DPhi_1_DEta_12_SameSideCanvas"));	  
  corcanvasPhiEta->Print(Form("%s/3phietacor.pdf",dir.Data()));
  TFile * phietacor = TFile::Open(Form("%s/3phietacor.root",dir.Data()),"RECREATE");
  phietacor->cd();
  corcanvasPhiEta->Write("phietacor");
  corcanvasPhiEta->Close();  
  delete corcanvasPhiEta;delete phietacor;  
//   TDirectory * yield = file->GetDirectory("yield");  
  
  
}
void DrawPbPb(TFile* file, TString dir){
  cout << dir<<endl;  
}


void DrawProductions(){
  TStringToken prods("LHC10b LHC10c LHC10d LHC10e LHC10h LHC11a LHC11h"," ");
  while(prods.NextToken()){
    TFile * rfile =  TFile::Open(Form("%s/results.root",prods.Data()),"READ");
    if(rfile){
      if(!rfile->IsZombie()){
// 	rfile->ls();
	TString dir(Form("results/%s",prods.Data()));
	if(prods.CompareTo("LHC10h")==0||prods.CompareTo("LHC11h")==0){DrawPbPb(rfile,dir);}
	else{Drawpp(rfile,dir);}
      }
    }
  }
}