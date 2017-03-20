#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TMath.h"
#include <TRandom3.h>

using namespace std;

Double_t CGaus(Double_t * x, Double_t * par){
  //Parameterization of signal, par[0] gives the integral of the gaussian.
  Double_t c=par[1] ;
  Double_t s=par[2] ;
  Double_t dx=(x[0]-c)/s ;
  return par[0]*TMath::Exp(-dx*dx/2.)/TMath::Sqrt(TMath::TwoPi())/s ;
}

void TestGaussian(){
//   TRandom3 * rndm = new TRandom3();
//   rndm->SetSeed(0);
// //   gRandom->rndm;
  TF1 * gaus = new TF1("fgs",CGaus,-TMath::Pi(),TMath::Pi(),3);
  gaus->SetParameters(1.0,0.0,0.1);
  gaus->SetParLimits(0,0.0,1.0E10);
  gaus->SetParLimits(2,0.0,6.0);
  TF1 * gausd = new TF1("fgsd",CGaus,-TMath::Pi(),TMath::Pi(),3);
  gausd->SetParameters(1.0,0.0,0.1);
  gausd->SetParLimits(0,0.0,1.0E10);
  gausd->SetParLimits(2,0.0,3.0);
  TF1 * gauseta = new TF1("fgseta",CGaus,-0.8,0.8,3);
  gauseta->SetParameters(1.0,0.0,0.1);
  gauseta->SetParLimits(0,0.0,1.0E10);
  gauseta->SetParLimits(2,0.0,1.0);
  TF1 * gausdeta = new TF1("fgsdeta",CGaus,-1.5,1.5,3);
  gausdeta->SetParameters(1.0,0.0,0.5);
  gausdeta->SetParLimits(0,0.0,1.0E12);
  gausdeta->SetParLimits(2,0.0,1.0);
  
  TH1D* theta = new TH1D("hTheta","Distribution of #theta",1000,0,TMath::Pi());
  TH1D* thetat = new TH1D("hTheta_t","Distribution of #theta for the triggers",1000,0,TMath::Pi());

  TH1D* theta12 = new TH1D("hTheta12","Distribution of #theta_12",1000,-TMath::Pi()/2.0,TMath::Pi()/2.0);
  TH1D* eta = new TH1D("heta","Distribution of #eta",1000,-1.5,1.5);
  TH1D* etat = new TH1D("heta_t","Distribution of #eta_t",1000,-1.5,1.5);

  TH1D* eta12 = new TH1D("heta12","Distribution of #eta_12",1000,-1.6,1.6);
  for(int i = 0; i<1000000;i++){
    double eta_t = 3.0*(gRandom->Rndm()-0.5);
    double theta_t = 2.0*TMath::ATan(TMath::Exp(-1.0*eta_t));
    if(theta_t>TMath::Pi()) theta_t-=TMath::Pi();
    if(theta_t<0.0) theta_t += TMath::Pi();
    
    double dtheta1 = gaus->GetRandom();
    double theta1 = theta_t + dtheta1;
    double dtheta2 = gaus->GetRandom();
    double theta2 = theta_t + dtheta2;
    double dtheta3 = gaus->GetRandom();
    double theta3 = theta_t + dtheta3;
  

    double dtheta12;
    dtheta12 = theta2-theta1;
    if(dtheta12>=  TMath::Pi()) dtheta12 -=TMath::Pi();
    if(dtheta12<= -TMath::Pi()) dtheta12 +=TMath::Pi();
    double dtheta13;
    dtheta13 = theta3-theta1;
    if(dtheta13>=  TMath::Pi()) dtheta13 -=TMath::Pi();
    if(dtheta13<= -TMath::Pi()) dtheta13 +=TMath::Pi();
    double dtheta23;
    dtheta23 = theta3-theta2;
    if(dtheta23>=  TMath::Pi()) dtheta23 -=TMath::Pi();
    if(dtheta23<= -TMath::Pi()) dtheta23 +=TMath::Pi();
    
    double eta1 = - TMath::Log(TMath::Tan(0.5*theta1));
    double eta2 = - TMath::Log(TMath::Tan(0.5*theta2));
    double eta3 = - TMath::Log(TMath::Tan(0.5*theta3));
    double etad12 = eta1-eta2;
    double etad13 = eta1-eta3;
    double etad23 = eta2-eta3;
    
//     cout << eta1<<" "<<eta2<<" "<<eta_t<<endl;
    if(abs(eta_t)<0.8){
      etat->Fill(eta_t);
      thetat->Fill(theta_t);
      if(abs(eta1)<0.8&&abs(eta2)<0.8){
	theta->Fill(theta1);
	theta->Fill(theta2);
	theta12->Fill(dtheta12);
	eta->Fill(eta1);
	eta->Fill(eta2);
	eta12->Fill(etad12);
      }
      if(abs(eta1)<0.8&&abs(eta3)<0.8){
	theta->Fill(theta3);
	theta12->Fill(dtheta13);
	eta->Fill(eta3);
	eta12->Fill(etad13);
	if(!abs(eta2)<0.8){
	  theta->Fill(theta1);
	  eta->Fill(eta1);
	}
      }
      if(abs(eta2)<0.8&&abs(eta3)<0.8){
	theta12->Fill(dtheta13);
	eta12->Fill(etad13);
	if(!abs(eta1)<0.8){
	  theta->Fill(theta3);
	  eta->Fill(eta3);
	  theta->Fill(theta2);
	  eta->Fill(eta2);	  
	}
      }
    }
    if(abs(eta1)<0.8||abs(eta2)<0.8) continue;
    theta->Fill(theta1);
    theta->Fill(theta2);
    theta12->Fill(dtheta12);
    eta->Fill(eta1);
    eta->Fill(eta2);
    eta12->Fill(etad12);
  }
  theta12->Fit(gausd,"MSQ","",-TMath::Pi()/4,5.0/4.0*TMath::Pi());
//   eta->Fit(gauseta,"MSQ","",-0.8,0.8);
  eta12->Fit(gausdeta,"MSQ","",-1.5,1.5);

  TFile * file =  TFile::Open("thetaeta.root","RECREATE");
  file->cd();
  theta->Write();
  thetat->Write();
  theta12->Write();
  eta->Write();
  etat->Write();
  eta12->Write();
  file->Close();
  delete theta; delete theta12;delete gaus;
}