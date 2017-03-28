#include <iostream>
#include "TFile.h"
#include "TRegexp.h"
#include "TPRegexp.h"
#include "TStyle.h"
#include "TString.h"
#include "TObject.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TNamed.h"
#include "TList.h"
#include "TDirectory.h"
#include <X3DDefs.h>
#include "TPaveText.h"
#include "TMinuit.h"
#include "TLegend.h"
#include "TPaveLabel.h"
#include "TObjString.h"
#include "THn.h"
#include "TArrayD.h"
#include "TSystemDirectory.h"
#include "TParameter.h"
#include "THashList.h"

using namespace std;

void Mergehnhists(THnF* input1, THnF* input2 , THnF* input3, THnF* output){
  output->Reset();
  output->Sumw2();
  TAxis* multaxis = input1->GetAxis(0);
  TAxis* vzaxis = input1->GetAxis(1);
  TAxis* phiaxis = input1->GetAxis(2);
  TAxis* etaaxis = input1->GetAxis(3);
  TAxis* ptaxis = input1->GetAxis(4);
  int bin[5] = {0,0,0,0,0};
//   int outbin[5] = {0,0,0,0,0};
  double value[5] = {0,0,0,0,0};
  for(int mult = 1;mult<=multaxis->GetNbins();mult++){
    double multc = multaxis->GetBinCenter(mult);
    for(int vz = 1; vz<= vzaxis->GetNbins();vz++){
      double vzc = vzaxis->GetBinCenter(vz);
      for(int phi = 1; phi<= phiaxis->GetNbins();phi++){
      double phic = phiaxis->GetBinCenter(phi);
	for(int eta = 1; eta<= etaaxis->GetNbins();eta++){
	double etac = etaaxis->GetBinCenter(eta);      
	  for(int pt = 1; pt<= ptaxis->GetNbins();pt++){
	  double ptc = ptaxis->GetBinCenter(pt);      
	  bin[0] = mult;
	  bin[1] = vz;
	  bin[2] = phi;
	  bin[3] = eta;
	  bin[4] = pt;
// 	  outbin[0] = output->GetAxis(0)->FindBin(multaxis->GetBinCenter(mult));
// 	  outbin[1] = output->GetAxis(1)->FindBin(vzaxis->GetBinCenter(vz));
// 	  outbin[2] = output->GetAxis(2)->FindBin(phiaxis->GetBinCenter(phi));
// 	  outbin[3] = output->GetAxis(3)->FindBin(etaaxis->GetBinCenter(eta));
// 	  outbin[4] = output->GetAxis(4)->FindBin(ptaxis->GetBinCenter(pt));
	  value[0] = multc;
	  value[1] = vzc;
	  value[2] = phic;
	  value[3] = etac;
	  value[4] = ptc;
// 	  output->SetBinContent(outbin,input1->GetBinContent(bin));
	  for(int val = 1;val<=input1->GetBinContent(bin);val++){
	      output->Fill(value);
	   }    
	  }
	}
      }
    }
  }
  
  multaxis = input2->GetAxis(0);
  vzaxis = input2->GetAxis(1);
  phiaxis = input2->GetAxis(2);
  etaaxis = input2->GetAxis(3);
  ptaxis = input2->GetAxis(4);
  for(int mult = 1;mult<=multaxis->GetNbins();mult++){
    double multc = multaxis->GetBinCenter(mult);
    for(int vz = 1; vz<= vzaxis->GetNbins();vz++){
      double vzc = vzaxis->GetBinCenter(vz);
      for(int phi = 1; phi<= phiaxis->GetNbins();phi++){
      double phic = phiaxis->GetBinCenter(phi);
	for(int eta = 1; eta<= etaaxis->GetNbins();eta++){
	double etac = etaaxis->GetBinCenter(eta);      
	  for(int pt = 1; pt<= ptaxis->GetNbins();pt++){
	  double ptc = ptaxis->GetBinCenter(pt);      
	  bin[0] = mult;
	  bin[1] = vz;
	  bin[2] = phi;
	  bin[3] = eta;
	  bin[4] = pt;
// 	  outbin[0] = output->GetAxis(0)->FindBin(multaxis->GetBinCenter(mult));
// 	  outbin[1] = output->GetAxis(1)->FindBin(vzaxis->GetBinCenter(vz));
// 	  outbin[2] = output->GetAxis(2)->FindBin(phiaxis->GetBinCenter(phi));
// 	  outbin[3] = output->GetAxis(3)->FindBin(etaaxis->GetBinCenter(eta));
// 	  outbin[4] = output->GetAxis(4)->FindBin(ptaxis->GetBinCenter(pt));
	  value[0] = multc;
	  value[1] = vzc;
	  value[2] = phic;
	  value[3] = etac;
	  value[4] = ptc;
// 	  output->SetBinContent(outbin,input2->GetBinContent(bin));
	  for(int val = 1;val<=input2->GetBinContent(bin);val++){
	      output->Fill(value);
	   }
	  }
	}
      }
    }
  }

  multaxis = input3->GetAxis(0);
  vzaxis = input3->GetAxis(1);
  phiaxis = input3->GetAxis(2);
  etaaxis = input3->GetAxis(3);
  ptaxis = input3->GetAxis(4);
  for(int mult = 1;mult<=multaxis->GetNbins();mult++){
    double multc = multaxis->GetBinCenter(mult);
    for(int vz = 1; vz<= vzaxis->GetNbins();vz++){
      double vzc = vzaxis->GetBinCenter(vz);
      for(int phi = 1; phi<= phiaxis->GetNbins();phi++){
      double phic = phiaxis->GetBinCenter(phi);
	for(int eta = 1; eta<= etaaxis->GetNbins();eta++){
	double etac = etaaxis->GetBinCenter(eta);      
	  for(int pt = 1; pt<= ptaxis->GetNbins();pt++){
	  double ptc = ptaxis->GetBinCenter(pt);      
	  bin[0] = mult;
	  bin[1] = vz;
	  bin[2] = phi;
	  bin[3] = eta;
	  bin[4] = pt;
// 	  outbin[0] = output->GetAxis(0)->FindBin(multaxis->GetBinCenter(mult));
// 	  outbin[1] = output->GetAxis(1)->FindBin(vzaxis->GetBinCenter(vz));
// 	  outbin[2] = output->GetAxis(2)->FindBin(phiaxis->GetBinCenter(phi));
// 	  outbin[3] = output->GetAxis(3)->FindBin(etaaxis->GetBinCenter(eta));
// 	  outbin[4] = output->GetAxis(4)->FindBin(ptaxis->GetBinCenter(pt));
	  value[0] = multc;
	  value[1] = vzc;
	  value[2] = phic;
	  value[3] = etac;
	  value[4] = ptc;
// 	  output->SetBinContent(outbin,input3->GetBinContent(bin));
	  for(int val = 1;val<=input3->GetBinContent(bin);val++){
	      output->Fill(value);
	   } 
	  }
	}
      }
    }
  }

}

Double_t myfunc(Double_t* x,Double_t* par){
  double limit = 8;
  if(x[0]<limit){
    return par[0]+par[1]*x[0]+par[2]*x[0]*x[0]+par[3]*x[0]*x[0]*x[0]+par[4]*x[0]*x[0]*x[0]*x[0];
  }
  else return par[0]+limit*par[1]+limit*limit*par[2]+limit*limit*limit*par[3]+par[4]*limit*limit*limit*limit;
}

void CollectCentBins(const char* runnr,TString basedir){
  TString ghidir = basedir + TString("LHC12a17ghi/");
  cout << basedir.Data()<< " " << ghidir << endl;
  TSystemDirectory * runsdir = new TSystemDirectory("LHC12a17ghi",ghidir.Data());
  TList * runs = runsdir->GetListOfFiles();
  Int_t nphi = 1;
  Double_t phimin = 0.0;
  Double_t phimax = 2.0*TMath::Pi();
  Int_t nEta = 63;
  Double_t EtaMin = -0.9;
  Double_t EtaMax =  0.9;
  Double_t pTmin = 0.5;
  Double_t pTmax = 100.0;
  Double_t fMaxPt = pTmax;
  
  Int_t npT = (16.0-pTmin)/0.1 - 0.5;//rounding up
  Int_t nextra = 0;
  if(pTmax > 16.0){
    nextra = 2 + (fMaxPt-20.0)/10.0 ;//rounding up, one bin 16 - 20, rest 10 GeV/c bins
  }
  Double_t pTbins[npT+nextra+1];
  Double_t nowpT = pTmin;
  for (int i = 0;i<=npT+nextra;i++){
    pTbins[i] = nowpT;
    if(nowpT<16.0){
      nowpT+=0.1;
    }
    if(nowpT>=20.0){
      nowpT +=10.0;
    }
    if(nowpT >= 16.0 && nowpT<20.0){
      nowpT=20.0;
    }

  }  
  
  

  Int_t  bins[5]   = {24, 10,nphi,nEta,npT+nextra};
  Double_t xmin[5] = {0.0,-10.0,phimin,EtaMin,pTmin};
  Double_t xmax[5] = {90.0,10.0,phimax,EtaMax,pTmax};
  Double_t Mbins[25] = {0.0,1.25,2.5,3.75,5.0,6.25,7.5,8.75,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0};
  Double_t VZbins[11] = {-10.0,-7.5,-5.0,-3.5,-2.0,0.0,2.0,3.5,5.0,7.5,10};
  Double_t etamin = -0.9;
  Double_t etamax = 0.9;
  Double_t deta = (etamax - etamin)/63;
  TArrayD * etaaxisAD = new TArrayD(63+1);
  for(int i = 0; i<=63;i++){
    etaaxisAD->AddAt(etamin+i*deta,i);
  }
  
  THnF * savehist = new THnF("hnsavehist1","",5,bins,xmin,xmax);
  savehist->GetAxis(0)->Set(24,Mbins);
  savehist->GetAxis(1)->Set(10,VZbins);
  savehist->GetAxis(4)->Set(npT+nextra,pTbins);
  THnF * savehist2 = new THnF("hnsavehist2","",5,bins,xmin,xmax);
  savehist2->GetAxis(0)->Set(24,Mbins);
  savehist2->GetAxis(1)->Set(10,VZbins);
  savehist2->GetAxis(4)->Set(npT+nextra,pTbins);

  for(int i = 0; i<runs->GetEntries();i++){
    TSystemFile *file = dynamic_cast<TSystemFile*>(runs->At(i)); 
    if(file->IsDirectory()){
      if(!TString(runs->At(i)->GetName()).BeginsWith(".")&&((!TString(runs->At(i)->GetName()).CompareTo(runnr))||(!TString(runnr).CompareTo("")))){
	TString LocFileg = basedir + TString("LHC12a17g/") + TString(runs->At(i)->GetName()) + TString("mcg/eff.root");
	TString LocFileh = basedir + TString("LHC12a17h/") + TString(runs->At(i)->GetName()) + TString("mch/eff.root");
	TString LocFilei = basedir + TString("LHC12a17i/") + TString(runs->At(i)->GetName()) + TString("mci/eff.root");
	TFile * effg = TFile::Open(LocFileg.Data(),"READ");
	TFile * effh = TFile::Open(LocFileh.Data(),"READ");
	TFile * effi = TFile::Open(LocFilei.Data(),"READ");
	if(effg&&effh&&effi){
	  cout << LocFileg.Data()<<" + "<<  LocFileh.Data()<<" + "<<  LocFilei.Data()<<" in run "<<   runs->At(i)->GetName() << endl;
	  TString Locfile = ghidir + TString("/") + TString(runs->At(i)->GetName()) + TString("/eff.root");
	  TFile * effr = TFile::Open(Locfile.Data(),"RECREATE");
	  TDirectory * rhists = effr->mkdir("hists");

	  TH2D* vertexvscentg = (TH2D*)(effg->Get("hnCentralityvsVertex"));
	  TH2D* vertexvscenth = (TH2D*)(effh->Get("hnCentralityvsVertex"));
	  TH2D* vertexvscenti =	(TH2D*)(effi->Get("hnCentralityvsVertex"));
	  TH2D* vertexvscent  = (TH2D*)(vertexvscentg->Clone("CentralityvsVertexofEvents"));
	  vertexvscent->Reset();
	  vertexvscent->SetTitle("Vertex vs Centrality distribution of the Events.");
	  vertexvscent->Add(vertexvscentg);
	  vertexvscent->Add(vertexvscenth);
	  vertexvscent->Add(vertexvscenti);
	  vertexvscent->Write();
	  delete vertexvscent;
	  delete vertexvscentg;
	  delete vertexvscenth;
	  delete vertexvscenti;	  
  
	  THnF * RecTracksg = (THnF*)(effg->Get("hnTracksReconstruced"));
	  THnF * RecTracksh = (THnF*)(effh->Get("hnTracksReconstruced"));
	  THnF * RecTracksi = (THnF*)(effi->Get("hnTracksReconstruced"));
	  Mergehnhists(RecTracksg,RecTracksh,RecTracksi,savehist);
	  delete RecTracksg; 
	  delete RecTracksi;
	  delete RecTracksh;
  	  rhists->cd();
	  TH1D * multrec = savehist->Projection(0);
	  TH1D * vzrec = savehist->Projection(1);
	  TH1D * etarec = savehist->Projection(3);
	  TH1D * ptrec = savehist->Projection(4);
// 	  for(int i=1; i<= ptrec->Integral();i++) ntrackstotal->Fill(0.5);
// 	  ptRecTracksgl->Add(ptrec);
	  savehist->Write("hnTracksReconstruced");
	  THnF * MCTracksg = (THnF*)(effg->Get("hnTracksProduced"));
	  THnF * MCTracksh = (THnF*)(effh->Get("hnTracksProduced"));
	  THnF * MCTracksi = (THnF*)(effi->Get("hnTracksProduced"));
	  Mergehnhists(MCTracksg,MCTracksh,MCTracksi,savehist2);
	  delete MCTracksg; 
	  delete MCTracksh;
	  delete MCTracksi;
	  TH1D * multmc = savehist2->Projection(0);
	  TH1D * vzmc = savehist2->Projection(1);
	  TH1D * etamc = savehist2->Projection(3);
	  TH1D * ptmc = savehist2->Projection(4);
	  savehist2->Write("hnTracksProduced");
	  
	  
	  multrec->Write("hnTracksReconstructed_Cent");
	  multmc->Write("hnTracksProduced_Cent");
	  vzrec->Write("hnTracksReconstructed_VZ");
	  vzmc->Write("hnTracksProduced_VZ");
	  etarec->Write("hnTracksReconstructed_eta");
	  etamc->Write("hnTracksProduced_eta");
	  ptrec->Write("hnTracksReconstructed_pT");
	  ptmc->Write("hnTracksProduced_pT");
	  
	  multrec->Divide(multmc);
	  multrec->SetTitle("Projection of the efficiency along the centrality axis.");
	  multrec->GetXaxis()->SetTitle("Centrality (%)");
	  multrec->Write("hnEfficiency_Cent");
	  vzrec->Divide(vzmc);
	  vzrec->SetTitle("Projection of the efficiency along the VZ axis.");
	  vzrec->GetXaxis()->SetTitle("VZ (cm)");
	  vzrec->Write("hnEfficiency_VZ");
	  etarec->Divide(etamc);
	  etarec->SetTitle("Projection of the efficiency along the eta axis.");
	  etarec->GetXaxis()->SetTitle("eta ");
	  etarec->Write("hnEfficiency_eta");
	  ptrec->Divide(ptmc);
	  ptrec->SetTitle("Projection of the efficiency along the pT axis.");
	  ptrec->GetXaxis()->SetTitle("pT (GeV/c)");
	  ptrec->Write("hnEfficiency_pT");
	  
	  delete multrec;delete multmc;delete vzrec;delete vzmc;delete etarec;delete etamc;delete ptrec;delete ptmc;	  
	  effr->cd();
	  
	  TH3D* RecTrackslptg = (TH3D*)(effg->Get("hist3drec"));
	  TH3D* RecTrackslpth = (TH3D*)(effh->Get("hist3drec"));
	  TH3D* RecTrackslpti = (TH3D*)(effi->Get("hist3drec"));	  
	  TH3D * RecTrackslpt = (TH3D*)RecTrackslptg->Clone("hist3drec");
	  RecTrackslpt->Reset();
	  RecTrackslpt->Add(RecTrackslptg);
	  RecTrackslpt->Add(RecTrackslpth);
	  RecTrackslpt->Add(RecTrackslpti);
	  delete RecTrackslptg; delete RecTrackslpth; delete RecTrackslpti;
	  RecTrackslpt->Write();
	  TH2D* RecTrackshptg = (TH2D*)(effg->Get("hist3drechpT"));
	  TH2D* RecTrackshpth = (TH2D*)(effh->Get("hist3drechpT"));
	  TH2D* RecTrackshpti = (TH2D*)(effi->Get("hist3drechpT"));		  
	  TH2D * RecTrackshpt = (TH2D*)(RecTrackshptg->Clone("hist3drechpT"));
	  RecTrackshpt->Reset();
	  RecTrackshpt->Add(RecTrackshptg);
	  RecTrackshpt->Add(RecTrackshpth);
	  RecTrackshpt->Add(RecTrackshpti);
	  delete RecTrackshptg;delete RecTrackshpth;delete RecTrackshpti;
	  RecTrackshpt->Write("hist2drechpT");
	  
	  TH3D* MCTrackslptg = (TH3D*)(effg->Get("hist3dMC"));
	  TH3D* MCTrackslpth = (TH3D*)(effh->Get("hist3dMC"));
	  TH3D* MCTrackslpti = (TH3D*)(effi->Get("hist3dMC"));	
	  TH3D * MCTrackslpt = (TH3D*)(MCTrackslptg->Clone("hist3dMC"));
	  MCTrackslpt->Reset();
	  MCTrackslpt->Add(MCTrackslptg);
	  MCTrackslpt->Add(MCTrackslpth);
	  MCTrackslpt->Add(MCTrackslpti);
	  delete MCTrackslptg; delete MCTrackslpth;delete MCTrackslpti;
	  MCTrackslpt->Write();
	  
	  TH2D* MCTrackshptg = (TH2D*)(effg->Get("hist3dMChpT"));
	  TH2D* MCTrackshpth = (TH2D*)(effh->Get("hist3dMChpT"));
	  TH2D* MCTrackshpti = (TH2D*)(effi->Get("hist3dMChpT"));
	  TH2D* MCTrackshpt = (TH2D*)(MCTrackshptg->Clone("hist3dMChpT"));
	  MCTrackshpt->Reset();
	  MCTrackshpt->Add(MCTrackshptg);
	  MCTrackshpt->Add(MCTrackshpth);
	  MCTrackshpt->Add(MCTrackshpti);
	  delete MCTrackshptg; delete MCTrackshpth; delete MCTrackshpti;
	  MCTrackshpt->Write("hist2dMChpT");
	  
	  MCTrackslpt->Divide(RecTrackslpt);
	  MCTrackslpt->Write("Weight_lpT");
	  MCTrackshpt->Divide(RecTrackshpt);
	  MCTrackshpt->Write("Weight_hpT");	  
	  TH1D * err = new TH1D("err","Errors in the 3d histograms",1000,0.0,20.0);
	  for(int x = 1; x<=MCTrackslpt->GetNbinsX();x++){
	    for(int y = 1; y<=MCTrackslpt->GetNbinsY();y++){
	      for(int z =1; z<=MCTrackslpt->GetNbinsZ();z++){
		if(MCTrackslpt->GetBinContent(x,y,z)>1.0E-10){
		  double content = MCTrackslpt->GetBinContent(x,y,z);
		  double error =  MCTrackslpt->GetBinError(x,y,z);
		  error = error/content;
		  err->Fill(100*error);
		}
	      }
	    }
	  }
	  for(int x = 1; x<=MCTrackshpt->GetNbinsX();x++){
	    for(int y = 1; y<=MCTrackshpt->GetNbinsY();y++){
	      if(MCTrackshpt->GetBinContent(x,y)>1.0E-10){  
		double content = MCTrackshpt->GetBinContent(x,y);
		double error =  MCTrackshpt->GetBinError(x,y);
		error = error/content;
		err->Fill(100*error);
		}
	      }
	    }
	  delete RecTrackslpt; delete RecTrackshpt;delete MCTrackslpt; delete MCTrackshpt;
	  
	  err->Write();
	  delete err;
	  


	  effr->Close();
	  delete effr;

	}
	if(effg)effg->Close();
	if(effh)effh->Close();
	if(effi)effi->Close();
	delete effg; delete effh; delete effi;
      }
    }
  }
  delete runs;
  delete savehist;delete savehist2;
  delete runsdir;
}

void CollectRuns(TString basedir,Double_t ptborder=3.0,Double_t MaxWeight=5.0){
  TString globfile = basedir + TString("/eff.root");
  TFile * effglob = TFile::Open(globfile.Data(),"RECREATE");
  TDirectory* ghists = effglob->mkdir("hists");
  TDirectory* gpruns = effglob->mkdir("pruns");

  
  Double_t pTmin = 0.5;
  Double_t pTmax = 100.0;
  Double_t fMaxPt = pTmax;
  Int_t npT = (16.0-pTmin)/0.1 - 0.5;//rounding up
  Int_t nextra = 0;
  if(pTmax > 16.0){
    nextra = 2 + (fMaxPt-20.0)/10.0 ;//rounding up, one bin 16 - 20, rest 10 GeV/c bins
  }
  Double_t pTbins[npT+nextra+1];
  Double_t nowpT = pTmin;
  for (int i = 0;i<=npT+nextra;i++){
    pTbins[i] = nowpT;
    if(nowpT<16.0){
      nowpT+=0.1;
    }
    if(nowpT>=20.0){
      nowpT +=10.0;
    }
    if(nowpT >= 16.0 && nowpT<20.0){
      nowpT=20.0;
    }

  }  
  
  TH1D* ptRecTracksgl  = new TH1D("hnTracksRec_pT","Number of reconstruced tracks as a function of pT.",npT+nextra,pTbins);
  TH1D* ptMCTracksgl   = new TH1D("hnTracksMC_pT","Number of produced tracks as a function of pT.",npT+nextra,pTbins);
  TH1D* nTracksTot = new TH1D("totaltracks","Total number of tracks in both cases.",2,0,2);
  nTracksTot->GetXaxis()->SetBinLabel(1,"reconstructed");
  nTracksTot->GetXaxis()->SetBinLabel(2,"produced");
  nTracksTot->GetXaxis()->LabelsOption("v");

  Double_t Mbins[25] = {0.0,1.25,2.5,3.75,5.0,6.25,7.5,8.75,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0};
  TH1D* MRecTracksgl  = new TH1D("hnTracksRec_MB","Number of reconstruced tracks as a function of Centrality.",24,Mbins);
  TH1D* MMCTracksgl  = new TH1D("hnTracksMC_MB","Number of produced tracks as a function of Centrality.",24,Mbins);
  Double_t VZbins[11] = {-10.0,-7.5,-5.0,-3.5,-2.0,0.0,2.0,3.5,5.0,7.5,10};
  TH1D* VZRecTracksgl  = new TH1D("hnTracksRec_VZ","Number of reconstructed tracks as a function of VZ.",10,VZbins);
  TH1D* VZMCTracksgl  = new TH1D("hnTracksMC_VZ","Number of produced tracks as a function of VZ.",10,VZbins);
  TH1D* EtaRecTraclsgl  = new TH1D("hnTracksRec_eta","Number of reconstructed tracks as a function of eta.",63,-0.9,0.9);
  TH1D* EtaMCTraclsgl  = new TH1D("hnTracksMC_eta","Number of produced tracks as a function of eta.",63,-0.9,0.9);


  Double_t multaxisArray[9] = {0.0,5.0,10.0,15.0,20.0,30.0,40.0,50.0,95.0};
  Double_t  vzaxisArray[8] = {-10.0,-7.5,-5.0,-2.0,2.0,5.0,7.5,10.0};
  Double_t  pTaxisArray[16] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0};
  TArrayD * pTaxisA = new TArrayD(16,pTaxisArray);
//   Double_t  pTaxisArray[26] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0};
  Double_t etamin = -0.9;
  Double_t etamax = 0.9;
  Double_t deta = (etamax - etamin)/63;
  TArrayD * etaaxisAD = new TArrayD(63+1);
  for(int i = 0; i<=63;i++){
    etaaxisAD->AddAt(etamin+i*deta,i);
  }
  
  TH3D* histlptrec 	= 	new TH3D("h3dlptrec","Reconstructed tracks in low pT",63,etaaxisAD->GetArray(),7,vzaxisArray,15,pTaxisArray);
  TH3D* histlptmc 	= 	new TH3D("h3dlptmc","Produced tracks in low pT",63,etaaxisAD->GetArray(),7,vzaxisArray,15,pTaxisArray);
  TH2D* histhptrec 	= 	new TH2D("h2dhptrec","Reconstructed tracks in high pT",63,etaaxisAD->GetArray(),7,vzaxisArray);
  TH2D* histhptmc 	= 	new TH2D("h2dhptmc","Produced tracks in high pT",63,etaaxisAD->GetArray(),7,vzaxisArray);
  int nruns = 107;
  Int_t runnumbersP11h[107] = {170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167987, 167985, 167920, 167915};
  TH1D * TracksperRun = new TH1D("TracksperRun", "# reconstruced tracks per Run", nruns, 0,1);
  TH1D * EventsperRun = new TH1D("EventsperRun", "# Events per Run", nruns, 0,1);
  TH2D * TracksperRunandVertex = new TH2D("TracksperRunandVertex", "# reconstruced tracks per Run and per vertex", nruns, 0,1,100,-10.0,10.0);
  TH2D * EventsperRunandVertex = new TH2D("EventsperRunandVertex", "# Events per Run and per vertex", nruns, 0,1,100,-10.0,10.0);
  TH2D * EventsperRunandCent = new TH2D("EventsperRunandVertex", "# Events per Run and per vertex", nruns, 0,1,100,-10.0,10.0);
  //   TH1D * TracksperRunMC = new TH1D("MCperRun", "# MC tracks per Run", nruns, 0,1);
  TH1D * EffperRun = new TH1D("EffperRun", "Efficiency per Run", nruns, 0,1);
  EffperRun->SetStats(kFALSE);
  for(int j=0; j<nruns; j++){
    TString lable = Form("%i",runnumbersP11h[j]);
    TracksperRun->GetXaxis()->SetBinLabel(j+1, lable);
    EventsperRun->GetXaxis()->SetBinLabel(j+1, lable);
//     TracksperRunMC->GetXaxis()->SetBinLabel(j+1, lable);
    EffperRun->GetXaxis()->SetBinLabel(j+1, lable);
    TracksperRunandVertex->GetXaxis()->SetBinLabel(j+1, lable);
    EventsperRunandVertex->GetXaxis()->SetBinLabel(j+1, lable);
  }
  TracksperRun->GetXaxis()->LabelsOption("v");
//   TracksperRunMC->GetXaxis()->LabelsOption("v");
  EffperRun->GetXaxis()->LabelsOption("v");
  
  TH1D * errg = new TH1D("errg","Errors in the 3d histograms",1000,0.0,20.0);
  for(int i =0;i<nruns;i++){
    int runnr = runnumbersP11h[106-i];
    int binnumber = -1;
    for(int k =0;k<=nruns;k++){
      TString binlable = TString(EffperRun->GetXaxis()->GetBinLabel(k));
      if(binlable.CompareTo(Form("%i",runnr)) == 0)binnumber = k;
    }
    
    TString Locfile = basedir + TString("/") + TString(Form("%i",runnr)) + TString("/eff.root");
    cout << Locfile.Data()<< " "  << runnr <<endl;
    TFile * effr = TFile::Open(Locfile.Data(),"READ");
    TH1D  * efforrun = (TH1D*)(effr->Get("err"));
    errg->Add(efforrun);
    
    TH2D * eventsvertcent = (TH2D*)(effr->Get("CentralityvsVertexofEvents"));
//     TParameter<double>* effofrun = (TParameter<double>*)(effr->Get("EffRun"));
//     EffperRun->SetBinContent(binnumber,effofrun->GetVal());
    
    TH1D* ntracksrecpt = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksReconstructed_pT"));
    TH1D* ntracksmcpt  = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksProduced_pT"));
    int nTracksRec = ntracksrecpt->Integral();
    int nTracksMC = ntracksmcpt->Integral();
    ptRecTracksgl->Add(ntracksrecpt);
    ptMCTracksgl->Add(ntracksmcpt);
    nTracksTot->Fill(0.5,nTracksRec);
    nTracksTot->Fill(1.5,nTracksMC);
    
    double eff = ntracksrecpt->Integral()/ntracksmcpt->Integral();
    EffperRun->SetBinContent(binnumber,eff);
    EventsperRun->SetBinContent(binnumber,eventsvertcent->Integral());
    TracksperRun->SetBinContent(binnumber,ntracksrecpt->Integral());
    for(int i =1; i<=EventsperRunandVertex->GetNbinsY();i++){
      EventsperRunandVertex->SetBinContent(binnumber,i,eventsvertcent->Integral(1,eventsvertcent->GetNbinsX(),i,i));
    
    }
    
    TH1D* ntracksrecm = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksReconstructed_Cent"));
    TH1D* ntracksmcm  = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksProduced_Cent"));
    MRecTracksgl->Add(ntracksrecm);
    MMCTracksgl->Add(ntracksmcm);
    
    TH1D* ntracksrecvz = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksReconstructed_VZ"));
    TH1D* ntracksmcvz  = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksProduced_VZ"));
    VZRecTracksgl->Add(ntracksrecvz);
    VZMCTracksgl->Add(ntracksmcvz);   

    TH1D* ntracksreceta = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksReconstructed_eta"));
    TH1D* ntracksmceta  = (TH1D*)(effr->GetDirectory("hists")->Get("hnTracksProduced_eta"));
    EtaRecTraclsgl->Add(ntracksreceta);
    EtaMCTraclsgl->Add(ntracksmceta);   
    
    TH3D* histlptrecl = (TH3D*)(effr->Get("hist3drec"));
    histlptrec->Add(histlptrecl);
    TH3D* histlptmcl = (TH3D*)(effr->Get("hist3dMC"));
    histlptmc->Add(histlptmcl);
    TH2D* histhptrecl = (TH2D*)(effr->Get("hist2drechpT"));
    histhptrec->Add(histhptrecl);
    TH2D* histhptmcl = (TH2D*)(effr->Get("hist2dMChpT"));
    histhptmc->Add(histhptmcl);
    effr->Close();
    delete effr;
  }
  gpruns->cd();
  errg->Write();
  EventsperRun->Write();
  EventsperRunandVertex->Write();  
  EffperRun->Write();
  TracksperRun->Write();
  TracksperRunandVertex->Write();
  TracksperRun->Divide(EventsperRun);
  TracksperRun->SetTitle("Tracks per event as a function of run.");
  TracksperRun->Write("TracksperEventperRun");
  ghists->cd();
  TCanvas * canvas = new TCanvas("canvas");
  nTracksTot->Write();
  ptRecTracksgl->Write();
  ptMCTracksgl->Write();
  histlptrec->Write();
  histlptmc->Write();
  histhptrec->Write();
  histhptmc->Write();
  MRecTracksgl->Write();
  MMCTracksgl->Write();
  VZRecTracksgl->Write();
  VZMCTracksgl->Write();
  EtaRecTraclsgl->Write();
  EtaMCTraclsgl->Write();
  MRecTracksgl->Divide(MMCTracksgl);
  MRecTracksgl->SetTitle("N_reconstructed/N_MC");
  MRecTracksgl->Write("Efficiency_Cent");
  canvas->cd();
  MRecTracksgl->SetStats(false);
  MRecTracksgl->SetMinimum(0.0);
  MRecTracksgl->GetXaxis()->SetTitle("Centrality [cm]");
  MRecTracksgl->GetYaxis()->SetTitle("#varepsilon []");
  MRecTracksgl->Draw();
  canvas->SaveAs("Efficiency_Mult.eps");
  canvas->Clear();
  VZRecTracksgl->Divide(VZMCTracksgl);
  VZRecTracksgl->SetTitle("N_reconstructed/N_MC");
  VZRecTracksgl->Write("Efficiency_VZ");
  VZRecTracksgl->SetStats(false);
  VZRecTracksgl->SetMinimum(0.0);
  VZRecTracksgl->GetXaxis()->SetTitle("Z-Vertex [cm]");
  VZRecTracksgl->GetYaxis()->SetTitle("#varepsilon []");
  VZRecTracksgl->Draw();
  canvas->SaveAs("Efficiency_VZ.eps");
  canvas->Clear();
  EtaRecTraclsgl->Divide(EtaMCTraclsgl);
  EtaRecTraclsgl->SetTitle("N_reconstructed/N_MC");
  EtaRecTraclsgl->Write("Efficiency_eta");
  EtaRecTraclsgl->SetStats(false);
  EtaRecTraclsgl->SetMinimum(0.0);
  EtaRecTraclsgl->GetXaxis()->SetTitle("#eta []");
  EtaRecTraclsgl->GetYaxis()->SetTitle("#varepsilon []");
  EtaRecTraclsgl->Draw();
  canvas->SaveAs("Efficiency_eta.eps");
  canvas->Clear();  
  TH1D * effpt = (TH1D*)ptRecTracksgl->Clone("Efficiency_pT");
  effpt->Divide(ptMCTracksgl);
  effpt->SetTitle("N_reconstructed/N_MC");
  effpt->Write();
  effpt->GetXaxis()->SetRangeUser(0.0,50.0);
  effpt->SetStats(false);
  effpt->SetMinimum(0.0);
  effpt->GetYaxis()->SetTitle("pT [GeV/c]");
  effpt->GetXaxis()->SetTitle("#varepsilon []");
  effpt->Draw();
  canvas->SaveAs("Efficiency_pT.eps");
  canvas->Clear();    effglob->cd();
  Double_t nbins0 = 0.0;
  Double_t content1 = 0.0;
  Double_t content2 = 0.0;
  Double_t average = 0.0;
  
  
  for(int x = 1;x<=histhptmc->GetXaxis()->GetNbins();x++){
    for(int y = 1;y<=histhptmc->GetYaxis()->GetNbins();y++){
      if(histhptmc->GetBinContent(x,y)>1.0E-10){
	if(histhptrec->GetBinContent(x,y)>1.0E-10){
	  if((histhptmc->GetBinContent(x,y)/histhptrec->GetBinContent(x,y))>MaxWeight){
	    histhptrec->SetBinContent(x,y,0.0);
	    histhptrec->SetBinError(x,y,0.0);
	    histhptmc->SetBinContent(x,y,0.0);
	    histhptmc->SetBinError(x,y,0.0);
	    continue;
	  }
	  nbins0 +=1.0;
	  content1+=histhptmc->GetBinContent(x,y);
	  content2+=histhptrec->GetBinContent(x,y);
	  average += content1/content2;
	}
      }
    }
  }
  cout << content1/content2<< " " << average/nbins0<<endl;
  
  histlptmc->Divide(histlptrec);
  histhptmc->Divide(histhptrec);
  Double_t nbins = 0.0;
  Double_t content = 0.0;
  TH1D* errnow = new TH1D("errnow","",1000,0.0,100.0);
  //remove 
  for(int x = 1;x<=histhptmc->GetXaxis()->GetNbins();x++){
    for(int y = 1;y<=histhptmc->GetYaxis()->GetNbins();y++){
      if(histhptmc->GetBinContent(x,y)>1.0E-10){
	nbins +=1.0;
	content+=histhptmc->GetBinContent(x,y);
	double err = histhptmc->GetBinError(x,y);
	double con = histhptmc->GetBinContent(x,y);
	errnow->Fill(100*err/con);
      }
    }
  }
  for(int x = 1;x<=histlptmc->GetXaxis()->GetNbins();x++){
    for(int y = 1;y<=histlptmc->GetYaxis()->GetNbins();y++){
      for(int z = 1;z<=histlptmc->GetZaxis()->GetNbins();z++){
	if(histlptmc->GetBinContent(x,y,z)>1.0E-10){
	  if(histlptmc->GetBinContent(x,y,z)>MaxWeight){
	    histlptmc->SetBinContent(x,y,z,0.0);
	    histlptmc->SetBinError(x,y,z,0.0);
	  }
	  else{
	    nbins +=1.0;
	    content+=histlptmc->GetBinContent(x,y,z);
	    double err = histlptmc->GetBinError(x,y,z);
	    double con = histlptmc->GetBinContent(x,y,z);
    // 	  if(con>10)cout << con<< " " << err <<" " <<err/con <<endl;
	    errnow->Fill(100*err/con);
	  }
	}
      }
    }
  }  
  histlptmc->Write("Weight_lpT_pT_vertex_mult");
  histhptmc->Write("Weight_hpT_eta_vertex_mult");
  canvas->cd();
  histhptmc->SetTitle("Weight (N_{MC}/N_{reconstructed})");
  histhptmc->SetStats(kFALSE);
  histhptmc->GetZaxis()->SetLabelSize(0.02);
  histhptmc->Draw("colz");
  canvas->SaveAs("Weight_eta_vz.eps");
  double averagept = 0;
  double ntracks = 0;
  for(int i = ptRecTracksgl->FindBin(2.0);i<ptRecTracksgl->GetNbinsX();i++){
    averagept += ptRecTracksgl->GetBinCenter(i)*ptRecTracksgl->GetBinContent(i);
    ntracks += ptRecTracksgl->GetBinContent(i);
  }
  averagept = averagept/ntracks;
  
  
  ptMCTracksgl->Divide(ptRecTracksgl);
  double weightataverage = ptMCTracksgl->GetBinContent(ptMCTracksgl->FindBin(averagept));
  ptMCTracksgl->Scale(1.0/weightataverage);
  
  TF1* ptfunc = new TF1("pT_function",myfunc,2,50,5);
  ptfunc->SetNpx(10000);
  ptMCTracksgl->Fit(ptfunc,"","",2.0,50);
  ptMCTracksgl->Write("Weight_pT");
  canvas->cd();
  ptMCTracksgl->GetXaxis()->SetTitle("pT [GeV/c]");
  ptMCTracksgl->Draw();
  ptMCTracksgl->GetXaxis()->SetRangeUser(0.0,50.0);
//   ptMCTracksgl->SetMaximum(1.2);
//   ptMCTracksgl->SetMinimum(0.6);
  ptMCTracksgl->SetTitle("Weight in pT scaled to 1 at average pT.");
  canvas->Update();
  canvas->SaveAs("Weight_pT.eps");
  canvas->Clear();
  
  errnow->Write();
  TString  outfile = basedir + TString("LHC11hWeight.root");
  TFile* outfile2 = TFile::Open(outfile.Data(),"RECREATE");
  histlptmc->Write("hnWeight");
  histhptmc->Write("hnWeight_highpt");
  ptfunc->Write();
  outfile2->Close();
  delete outfile2;
    effglob->Close();

}

void RunStatsPlus(){
  TString basedir = TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/MCproductions/LHC12a17ghi");
  TString globfile = basedir + TString("/eff.root");
  TFile * effglob = TFile::Open(globfile.Data(),"READ");

  TString runsfiles = basedir + TString("/AnalysisResults.root");
  TFile * statsglob = TFile::Open(runsfiles.Data(),"READ");
  
  TString dir2 = TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/LHC11h/QA/AnalysisResults.root");
  TFile * RealQA = TFile::Open(dir2.Data(),"READ");
  
  TString runsfile = basedir + TString("/runs.root");
  TFile * runs = TFile::Open(runsfile.Data(),"RECREATE");
  
  Double_t pTmin = 0.5;
  Double_t pTmax = 16.0;
  Int_t npT = (pTmax-pTmin)/0.1 - 0.5;//rounding up
  TH2D* ptEff  = new TH2D("hnEff_pT","Efficiency as a function of pT in runs.",(pTmax-pTmin)/0.1 - 0.5,0.5,16.0,15,0.5,0.9);
  TH2D* ptRelEff = new TH2D("hnEff_Rel_pT","Relative difference in Efficiency between runs and average as a function of pT.",(pTmax-pTmin)/0.1 - 0.5,0.5,16.0,15,0.85,1.15);
  Double_t Mbins[25] = {0.0,1.25,2.5,3.75,5.0,6.25,7.5,8.75,10.0,15.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,60.0,65.0,70.0,75.0,80.0,85.0,90.0};
  TH2D* MEff  = new TH2D("hnEff_Cent","Efficiency as a function of centrality in runs.",24,Mbins,15,0.7,0.86);
  TH2D* MRelEff  = new TH2D("hnEff_Rel_Cent","Relative difference in Efficiency between runs and average as a function of centrality.",24,Mbins,15,0.85,1.15);
  Double_t VZbins[11] = {-10.0,-7.5,-5.0,-3.5,-2.0,0.0,2.0,3.5,5.0,7.5,10};
  TH2D* VZEff  = new TH2D("hnEff_VZ","Efficiency as a function of VZ in runs.",10,VZbins,15,0.6,0.9);
  TH2D* VZRelEff  = new TH2D("hnEff_Rel_VZ","Relative difference in Efficiency between runs and average as a function of VZ.",10,VZbins,15,0.85,1.15);
  TH2D* EtaEff  = new TH2D("hnEff_eta","Efficiency as a function of eta in runs.",63,-0.9,0.9,15,0.4,0.9);
  TH2D* EtaRelEff  = new TH2D("hnEff_Rel_eta","Relative difference in Efficiency between runs and average as a function of eta.",63,-0.9,0.9,15,0.85,1.15);
  
  TH1D* Eff_Cent_glob = (TH1D*)(effglob->GetDirectory("hists")->Get("Efficiency_Cent"));
  TH1D* Eff_VZ_glob = (TH1D*)(effglob->GetDirectory("hists")->Get("Efficiency_VZ"));
  TH1D* Eff_eta_glob = (TH1D*)(effglob->GetDirectory("hists")->Get("Efficiency_eta"));
  TH1D* Eff_pt_glob = (TH1D*)(effglob->GetDirectory("hists")->Get("Efficiency_pT"));
  
  //3d hists:
  Double_t multaxisArray[8] = {0.0,5.0,10.0,15.0,20.0,30.0,40.0,95.0};
  Double_t  vzaxisArray[8] = {-10.0,-7.5,-5.0,-2.0,2.0,5.0,7.5,10.0};
  Double_t  pTaxisArray[28] = {0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9,3.0,3.5,4.0};
  Double_t etamin = -0.9;
  Double_t etamax = 0.9;
  Double_t deta = (etamax - etamin)/63;
  TArrayD * etaaxisAD = new TArrayD(63/3.0+1);
  for(int i = 0; i<=63/3.0;i++){
    etaaxisAD->AddAt(etamin+i*3.0*deta,i);
  }
  TH2D* ptWeight  = new TH2D("hnWeight_pT","Efficiency as a function of pT in runs.",27,pTaxisArray,15,1.0,2.0);
  TH2D* ptRelWeight = new TH2D("hnWeight_Rel_pT","Relative difference in Efficiency between runs and average as a function of pT.",27,pTaxisArray,30,0.0,2.0);
  TH2D* MWeight  = new TH2D("hnWeight_Cent","Efficiency as a function of centrality in runs.",7,multaxisArray,15,1.0,2.0);
  TH2D* MRelWeight  = new TH2D("hnWeight_Rel_Cent","Relative difference in Efficiency between runs and average as a function of centrality.",7,multaxisArray,30,0.0,2.0);
  TH2D* VZWeight  = new TH2D("hnWeight_VZ","Efficiency as a function of VZ in runs.",7,vzaxisArray,15,1.0,2.0);
  TH2D* VZRelWeight  = new TH2D("hnWeight_Rel_VZ","Relative difference in Efficiency between runs and average as a function of VZ.",7,vzaxisArray,30,0.0,2.0);
  TH2D* EtaWeight  = new TH2D("hnWeight_eta","Efficiency as a function of eta in runs.",21,etaaxisAD->GetArray(),15,1.0,2.0);
  TH2D* EtaRelWeight  = new TH2D("hnWeight_Rel_eta","Relative difference in Efficiency between runs and average as a function of eta.",21,etaaxisAD->GetArray(),30,0.0,2.0);
  
  
  TH3D* Weight_lpT = (TH3D*)(effglob->Get("Weight_lpT_pT_vertex_mult"));
  TH3D* Weight_hpT = (TH3D*)(effglob->Get("Weight_hpT_eta_vertex_mult"));
  

  Int_t runnumbersP11h[107] = {170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167987, 167985, 167920, 167915};

  for(int j=0; j<107; j++){
    int runnr = runnumbersP11h[106-j];
    TString Locfile = basedir + TString("/") + TString(Form("%i",runnr)) + TString("/eff.root");
    TFile * effr = TFile::Open(Locfile.Data(),"READ");
    
    //Centrality:
    TH1D* Cent_loc = (TH1D*)(effr->GetDirectory("hists")->Get("hnEfficiency_Cent"));
    for(int i = 1;i<=Cent_loc->GetXaxis()->GetNbins();i++){
      double fillv = Cent_loc->GetBinCenter(i);
      MEff->Fill(fillv, Cent_loc->GetBinContent(i));
      MRelEff->Fill(fillv,Cent_loc->GetBinContent(i)/Eff_Cent_glob->GetBinContent(i));
    }
    //Vertex:
    TH1D* VZ_loc = (TH1D*)(effr->GetDirectory("hists")->Get("hnEfficiency_VZ"));
    for(int i = 1;i<=VZ_loc->GetXaxis()->GetNbins();i++){
      double fillv = VZ_loc->GetBinCenter(i);
      VZEff->Fill(fillv, VZ_loc->GetBinContent(i));
      VZRelEff->Fill(fillv,VZ_loc->GetBinContent(i)/Eff_VZ_glob->GetBinContent(i));
    }
    //eta:
    TH1D* eta_loc = (TH1D*)(effr->GetDirectory("hists")->Get("hnEfficiency_eta"));
    for(int i = 1;i<=eta_loc->GetXaxis()->GetNbins();i++){
      double fillv = eta_loc->GetBinCenter(i);
      EtaEff->Fill(fillv, eta_loc->GetBinContent(i));
      EtaRelEff->Fill(fillv,eta_loc->GetBinContent(i)/Eff_eta_glob->GetBinContent(i));
    }
    //pT:
    TH1D* pT_loc = (TH1D*)(effr->GetDirectory("hists")->Get("hnEfficiency_pT"));
    for(int i = 1;i<=pT_loc->GetXaxis()->GetNbins();i++){
      double fillv = pT_loc->GetBinCenter(i);
      ptEff->Fill(fillv, pT_loc->GetBinContent(i));
      ptRelEff->Fill(fillv,pT_loc->GetBinContent(i)/Eff_pt_glob->GetBinContent(i));
    }      
    //3d lpT:
    TH3D* Weight_lpt_loc = (TH3D*)(effr->Get("Weight_lpT"));
    for(int x = 1; x<=Weight_lpt_loc->GetNbinsX();x++){
      for(int y = 1;y<=Weight_lpt_loc->GetNbinsY();y++){
	for(int z = 1;z<=Weight_lpt_loc->GetNbinsZ();z++){
	  double centval = Weight_lpt_loc->GetXaxis()->GetBinCenter(x);
	  double vzval = Weight_lpt_loc->GetYaxis()->GetBinCenter(y);
	  double ptval = Weight_lpt_loc->GetZaxis()->GetBinCenter(z);
	  ptWeight->Fill(ptval,Weight_lpt_loc->GetBinContent(x,y,z));
	  ptRelWeight->Fill(ptval,Weight_lpt_loc->GetBinContent(x,y,z)/Weight_lpT->GetBinContent(x,y,z));
	  MWeight->Fill(centval,Weight_lpt_loc->GetBinContent(x,y,z));
	  MRelWeight->Fill(centval,Weight_lpt_loc->GetBinContent(x,y,z)/Weight_lpT->GetBinContent(x,y,z));
	  VZWeight->Fill(vzval,Weight_lpt_loc->GetBinContent(x,y,z)/Weight_lpT->GetBinContent(x,y,z));
	  VZRelWeight->Fill(vzval,Weight_lpt_loc->GetBinContent(x,y,z)/Weight_lpT->GetBinContent(x,y,z));
	}
      }
    }
    //3d hpT:
    TH2D* Weight_hpt_loc = (TH2D*)(effr->Get("Weight_hpT"));
    for(int x = 1; x<=Weight_hpt_loc->GetNbinsX();x++){
      for(int y = 1;y<=Weight_hpt_loc->GetNbinsY();y++){
// 	for(int z = 1;z<=Weight_hpt_loc->GetNbinsZ();z++){
	  double centval = Weight_hpt_loc->GetXaxis()->GetBinCenter(x);
	  double vzval = Weight_hpt_loc->GetYaxis()->GetBinCenter(y);
// 	  double eta = Weight_hpt_loc->GetZaxis()->GetBinCenter(z);
// 	  EtaWeight->Fill(eta,Weight_hpt_loc->GetBinContent(x,y,z));
// 	  EtaRelWeight->Fill(eta,Weight_hpt_loc->GetBinContent(x,y,z)/Weight_hpT->GetBinContent(x,y,z));
	  MWeight->Fill(centval,Weight_hpt_loc->GetBinContent(x,y));
	  MRelWeight->Fill(centval,Weight_hpt_loc->GetBinContent(x,y)/Weight_hpT->GetBinContent(x,y));
	  VZWeight->Fill(vzval,Weight_hpt_loc->GetBinContent(x,y)/Weight_hpT->GetBinContent(x,y));
	  VZRelWeight->Fill(vzval,Weight_hpt_loc->GetBinContent(x,y)/Weight_hpT->GetBinContent(x,y));
// 	}
      }
    }
        
    effr->Close();
    delete effr;
  }
  THashList* Runslist = (THashList*)statsglob->GetDirectory("ThreePartTrackEfficienciesPbPb_0_16")->Get("ThreePartTrackEfficienciesPbPb_0_16_0_0Coutput1");
  TH3D* tracksvertexeta12a = (TH3D*)Runslist->FindObject("NTracksVertexEta");
  TH2D* tracksvertex12a = (TH2D*)tracksvertexeta12a->Project3D("yx");

  TH2D* eventsvertex12a = (TH2D*)Runslist->FindObject("NEventsVertex");

  TH1D* tracks12a = tracksvertex12a->ProjectionX();
  TH1D* events12a = eventsvertex12a->ProjectionX();

    cout << "6"<<endl;

  THashList* QAlist = (THashList*)RealQA->GetDirectory("ThreePartTracksQA11h_0_8")->Get("ThreePartTracksQA11h_0_8_8_16Coutput1");
  cout << "1"<<endl;
  TH2D* tracksvertex11h = (TH2D*)QAlist->FindObject("NTracksVertex");
  cout << "2"<<endl;
  TH2D* eventsvertex11h = (TH2D*)QAlist->FindObject("NEventsVertex");
  cout << "3"<<endl;
  TH1D* tracks11h = tracksvertex11h->ProjectionX();
  cout << "4"<<endl;
  TH1D* events11h = eventsvertex11h->ProjectionX();
  cout << "5"<<endl;
  TH1D* MCtracksperrunperevent  = (TH1D*)effglob->GetDirectory("pruns")->Get("TracksperEventperRun");
  cout << "6"<<endl;

  runs->cd();
  tracksvertex11h->Write("tracksvertex11h");
  eventsvertex11h->Write("eventsvertex11h");
  tracksvertex11h->Divide(eventsvertex11h);
  tracksvertex11h->SetTitle("#Tracks per event vs Run and Vertex.");
  tracksvertex11h->Write("trackspereventvertex11h");

  tracksvertex12a->Write("tracksvertex12a");
  eventsvertex12a->Write("eventsvertex12a");
  tracksvertex12a->Divide(eventsvertex12a);
  tracksvertex12a->SetTitle("#Tracks per event vs Run and Vertex.");
  tracksvertex12a->Write("trackspereventvertex12a");
  
  tracksvertex11h->Divide(tracksvertex12a);
  tracksvertex11h->SetTitle("NTracks per Event vs run data divided by MC.");
  tracksvertex11h->Write("TracksVEta_Data_div_MC");
  
  tracks11h->Write("tracks11h");
  events11h->Write("events11h");
  tracks11h->Divide(events11h);
  tracksvertex11h->SetTitle("#Tracks per event vs Run.");
  tracks11h->Write("tracksperrun11h");
  tracks11h->Divide(MCtracksperrunperevent);
  tracks11h->SetTitle("#tracks per event data divided by MC.");
  tracks11h->Write("TrackspEvent_DdivMC");
  MEff->Write();
  MRelEff->Write();
  VZEff->Write();
  VZRelEff->Write();
  EtaRelEff->Write();
  EtaEff->Write();
  ptRelEff->Write();
  ptEff->Write();

  MWeight->Write();
  MRelWeight->Write();
  VZWeight->Write();
  VZRelWeight->Write();
  EtaWeight->Write();
  EtaRelWeight->Write();
  ptWeight->Write();
  ptRelWeight->Write();
  runs->Close();
  effglob->Close();
  delete runs;delete effglob;
}