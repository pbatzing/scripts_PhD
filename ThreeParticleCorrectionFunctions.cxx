#include "ThreeParticleorrectionFunctions.h"

TH2D* gSameEvent =0x0;
TH2D* gMETA = 0x0;
TH2D* gMETA2 = 0x0;
TH2D* gMETrigger = 0x0;

TH2D* gSameEventETA =0x0;
TH2D* gMETAETA = 0x0;
TH2D* gMETA2ETA = 0x0;
TH2D* gMETriggerETA = 0x0;

bool gislowpTbin = false;

double gkEtaFitRange = 1.7;
double gEtasigRange = 1.5;

Double_t CGausPol0(Double_t * x, Double_t * par){
    //Parameterization of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1];
    Double_t s=par[2];
    Double_t dx=(x[0]-c);
    return par[0]*exp(-dx*dx/(2.0*s*s))/TMath::Sqrt(2.0*TMath::Pi())/s+par[3];
  }
Double_t CGaus(Double_t * x, Double_t * par){
    //Parameterization of signal, par[0] gives the integral of the gaussian.
    Double_t c=par[1];
    Double_t s=par[2];
    Double_t dx=(x[0]-c);
    return par[0]*exp(-dx*dx/(2.0*s*s))/TMath::Sqrt(2.0*TMath::Pi())/s;
  }  
Double_t C2Gaus(Double_t * x, Double_t * par){
    //Parameterization of signal, par[0] gives the integral of the gaussian.
    Double_t c1=par[1];
    Double_t s1=par[2];
    Double_t dx1=(x[0]-c1);
    Double_t c2=par[4];
    Double_t s2=par[5];
    Double_t dx2=(x[0]-c2);
    return par[0]*exp(-dx1*dx1/(2.0*s1*s1))/TMath::Sqrt(2.0*TMath::Pi())/s1 + par[3]*exp(-dx2*dx2/(2.0*s2*s2))/TMath::Sqrt(2.0*TMath::Pi())/s2;
  }  
Double_t CBKGPol0(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    //reject if inside +- gEtasigRange:
    Double_t c=par[1];
    Double_t dx=(x[0]-c);
    
    if(dx<gEtasigRange&&dx>-gEtasigRange){TF1::RejectPoint();return 0;}
    return par[0] ;
  }
Double_t CBKGPol0p(Double_t * x, Double_t * par){
    //Parameterization of background. Provide par[1]=c  for reference.
    if(TMath::Abs(x[0])>2.0)return 0.0;
    return par[0] ;
  }  

BinDirs::BinDirs(TFile* file, const char* bin, bool empty){
    if(file){
      Bin = TString(bin);
      file->cd();
      TDirectory * METAparentdir = file->GetDirectory("META");
      if(!METAparentdir)METAparentdir = file->mkdir("META");
      METAparentdir->cd();
      METAdir = METAparentdir->GetDirectory(bin);
      if(!METAdir)METAdir= METAparentdir->mkdir(bin);
      METAdir->cd();
      if(!METAdir->GetDirectory("bin_stats")) METAdir->mkdir("bin_stats");
      if(!METAdir->GetDirectory("divided")) METAdir->mkdir("divided");
      if(empty) ::resultsdirectory(METAdir,"bin_stats");
      if(empty) ::resultsdirectory(METAdir,"divided");
      TDirectory * META2parentdir = file->GetDirectory("META2");
      if(!META2parentdir)META2parentdir = file->mkdir("META2");
      META2parentdir->cd();
      META2dir = META2parentdir->GetDirectory(bin);
      if(!META2dir)META2dir= METAparentdir->mkdir(bin);
      META2dir->cd();
      if(!META2dir->GetDirectory("bin_stats")) META2dir->mkdir("bin_stats");
      if(!META2dir->GetDirectory("divided")) META2dir->mkdir("divided");
      if(empty) ::resultsdirectory(META2dir,"bin_stats");
      if(empty) ::resultsdirectory(META2dir,"divided");
      TDirectory * METriggerparentdir = file->GetDirectory("METrigger");
      if(!METriggerparentdir)METriggerparentdir = file->mkdir("METrigger");
      METriggerparentdir->cd();
      METriggerdir = METriggerparentdir->GetDirectory(bin);
      if(!METriggerdir)METriggerdir= METriggerparentdir->mkdir(bin);
      METriggerdir->cd();
      if(!METriggerdir->GetDirectory("bin_stats")) METriggerdir->mkdir("bin_stats");
      if(!METriggerdir->GetDirectory("divided")) METriggerdir->mkdir("divided");
      if(empty) ::resultsdirectory(METriggerdir,"bin_stats");
      if(empty) ::resultsdirectory(METriggerdir,"divided");
      file->cd();
      Samedir = file->GetDirectory(bin);
      if(!Samedir)Samedir =  file->mkdir(bin);
      Samedir->cd();
      if(!Samedir->GetDirectory("bin_stats")) Samedir->mkdir("bin_stats");
      if(!Samedir->GetDirectory("divided")) Samedir->mkdir("divided");
      if(empty) ::resultsdirectory(Samedir,"bin_stats");
      if(empty) ::resultsdirectory(Samedir,"divided");
      ready = true;
    }
    else{
      Samedir = NULL;
      METAdir = NULL;
      META2dir = NULL;
      METriggerdir = NULL;
      Bin = TString("");
      ready = false;
    }
  }
BinDirs::BinDirs(TDirectory* dir, const char* bin, bool empty){
    if(dir){
      Bin = TString(bin);
      dir->cd();
      TDirectory * METAparentdir = dir->GetDirectory("META");
      if(!METAparentdir)METAparentdir = dir->mkdir("META");
      METAparentdir->cd();
      METAdir = METAparentdir->GetDirectory(bin);
      if(!METAdir)METAdir= METAparentdir->mkdir(bin);
      METAdir->cd();
      if(!METAdir->GetDirectory("bin_stats")) METAdir->mkdir("bin_stats");
      if(!METAdir->GetDirectory("divided")) METAdir->mkdir("divided");
      if(empty) ::resultsdirectory(METAdir,"bin_stats");
      if(empty) ::resultsdirectory(METAdir,"divided");
    TDirectory * META2parentdir = dir->GetDirectory("META2");
      if(!META2parentdir)META2parentdir = dir->mkdir("META2");
      META2parentdir->cd();
      META2dir = META2parentdir->GetDirectory(bin);
      if(!META2dir)META2dir= META2parentdir->mkdir(bin);
      META2dir->cd();
      if(!META2dir->GetDirectory("bin_stats")) META2dir->mkdir("bin_stats");
      if(!META2dir->GetDirectory("divided")) META2dir->mkdir("divided");
      if(empty) ::resultsdirectory(META2dir,"bin_stats");
      if(empty) ::resultsdirectory(META2dir,"divided");
      TDirectory * METriggerparentdir = dir->GetDirectory("METrigger");
      if(!METriggerparentdir)METriggerparentdir = dir->mkdir("METrigger");
      METriggerparentdir->cd();
      METriggerdir = METriggerparentdir->GetDirectory(bin);
      if(!METriggerdir)METriggerdir= METriggerparentdir->mkdir(bin);
      METriggerdir->cd();
      if(!METriggerdir->GetDirectory("bin_stats")) METriggerdir->mkdir("bin_stats");
      if(!METriggerdir->GetDirectory("divided")) METriggerdir->mkdir("divided");
      if(empty) ::resultsdirectory(METriggerdir,"bin_stats");
      if(empty) ::resultsdirectory(METriggerdir,"divided");
      dir->cd();
      Samedir = dir->GetDirectory(bin);
      if(!Samedir)Samedir =  dir->mkdir(bin);
      Samedir->cd();
      if(!Samedir->GetDirectory("bin_stats")) Samedir->mkdir("bin_stats");
      if(!Samedir->GetDirectory("divided")) Samedir->mkdir("divided");
      if(empty) ::resultsdirectory(Samedir,"bin_stats");
      if(empty) ::resultsdirectory(Samedir,"divided");
      ready = true;
    }
    else{
      Samedir = NULL;
      METAdir = NULL;
      META2dir = NULL;
      METriggerdir = NULL;
      Bin = TString("");
      ready = false;
    }
  }
BinDirs::BinDirs(TDirectory* samed, TDirectory * METAd, TDirectory* META2d, TDirectory * METriggerd, bool empty){
    if(samed&&METriggerd&&METAd){
      Bin = TString(samed->GetName());
      samed->cd();
      Samedir = samed;
      if(!Samedir->GetDirectory("bin_stats")) Samedir->mkdir("bin_stats");      
      else if(empty){::resultsdirectory(Samedir,"bin_stats");}
      if(!Samedir->GetDirectory("divided")){Samedir->mkdir("divided");}
      else if(empty){::resultsdirectory(Samedir,"divided");}
      METriggerd->cd();
      METriggerdir = METriggerd;
      if(!METriggerdir->GetDirectory("bin_stats")) METriggerdir->mkdir("bin_stats");
      else if(empty) ::resultsdirectory(METriggerdir,"bin_stats");
      if(!METriggerdir->GetDirectory("divided")) METriggerdir->mkdir("divided");
      else if(empty) ::resultsdirectory(METriggerdir,"divided");
      METAd->cd();
      METAdir = METAd;
      if(!METAdir->GetDirectory("bin_stats")) METAdir->mkdir("bin_stats");      
      else if(empty) ::resultsdirectory(METAdir,"bin_stats");
      if(!METAdir->GetDirectory("divided")) METAdir->mkdir("divided");
      else if(empty) ::resultsdirectory(METAdir,"divided");
      if(META2d){
	META2d->cd();
	META2dir = META2d;
	if(!META2dir->GetDirectory("bin_stats")) META2dir->mkdir("bin_stats");      
	else if(empty) ::resultsdirectory(META2dir,"bin_stats");
	if(!META2dir->GetDirectory("divided")) META2dir->mkdir("divided");
	else if(empty) ::resultsdirectory(META2dir,"divided");
      }
      ready = true;
    }
    else{
      Samedir = NULL;
      METAdir = NULL;
      META2dir = NULL;
      METriggerdir = NULL;
      Bin = TString("");
      ready = false;
    }
  }  
  
BinDirs BinDirs::resultsdirectory(const char* dir){
  TDirectory * samedirdir = ::resultsdirectory(Samedir,dir);
  TDirectory * METAdirdir = ::resultsdirectory(METAdir,dir); 
  TDirectory * META2dirdir = NULL;
  if(META2dir) META2dirdir = ::resultsdirectory(META2dir,dir);
  TDirectory * METriggerdirdir = ::resultsdirectory(METriggerdir,dir);
  return BinDirs(samedirdir,METAdirdir,META2dirdir,METriggerdirdir);
}


TList * GetMZDirectories(TDirectory* same){
  TList * keys = same->GetListOfKeys();
  TList * directories = new TList();
  for(int i=0; i<keys->GetEntries();i++){
    if(TString(keys->At(i)->GetName()).BeginsWith("BinM(")&&TString(keys->At(i)->GetName()).Contains("Z(")){
      BinDirs* tmp = NULL;
      if(same->GetDirectory("META2")) tmp = new BinDirs(same->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("META")->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("META2")->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("METrigger")->GetDirectory(keys->At(i)->GetName()));
      else tmp = new BinDirs(same->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("META")->GetDirectory(keys->At(i)->GetName()),NULL,same->GetDirectory("METrigger")->GetDirectory(keys->At(i)->GetName()));
      directories->Add(tmp);
    }
    
  }
  return directories;
}
TList * GetMZDirectories(TDirectory* same, float vertexcut,bool debug){
  TList * keys = same->GetListOfKeys();
  TList * directories = new TList();
  for(int i=0; i<keys->GetEntries();i++){
      TString totbin = TString(keys->At(i)->GetName());
      if(totbin.BeginsWith("BinM(")&&totbin.Contains("Z(")){
      TString mbin = TString(totbin.Tokenize("Z")->At(0)->GetName());
      TString zbin = TString(totbin.Tokenize("Z")->At(1)->GetName());
      float zedge1 = TMath::Abs(TString(TString(zbin.Tokenize("(")->At(0)->GetName()).Tokenize(")")->At(0)->GetName()).Atof());
      float zedge2 = TMath::Abs(TString(TString(zbin.Tokenize("(")->At(1)->GetName()).Tokenize(")")->At(0)->GetName()).Atof());
      bool collect = true;
      if(zedge1>vertexcut||zedge2>vertexcut)collect = false;
      if(debug&&collect) cout << "Bin "<<totbin.Data() <<" will be added to " << mbin.Data() <<endl;
      if(debug&&!collect) cout << "Bin "<<totbin.Data() <<" will not be added to " << mbin.Data() <<endl;      
      if(collect){
	BinDirs* tmp = NULL;
	if(same->GetDirectory("META2")) tmp = new BinDirs(same->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("META")->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("META2")->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("METrigger")->GetDirectory(keys->At(i)->GetName()));
	else tmp = new BinDirs(same->GetDirectory(keys->At(i)->GetName()),same->GetDirectory("META")->GetDirectory(keys->At(i)->GetName()),NULL,same->GetDirectory("METrigger")->GetDirectory(keys->At(i)->GetName()));
	directories->Add(tmp);
      }
    }
    
  }
  return directories;
}
TH2D * PrepareHist(TDirectory* dir, const char* subdir, const char* name,const char* title,const char* rmplmode){
  //gets and prepares a TH2D from a directory, removes plateau if wished:
  TH2D * hist;
  TString mode = TString(rmplmode);
  hist = dynamic_cast<TH2D*>(dir->GetDirectory(subdir)->Get(name));
  hist->GetXaxis()->SetTitleOffset(1.0);
  hist->GetXaxis()->SetTitleSize(0.04);
  hist->GetYaxis()->SetTitleOffset(1.0);
  hist->GetYaxis()->SetTitleSize(0.04);
  hist->SetTitle(title);
  if(mode.CompareTo("same") ==0){
    RemovePlateau(hist->Integral(hist->GetXaxis()->FindBin(-1.0),hist->GetXaxis()->FindBin(1.0),hist->GetYaxis()->FindBin(TMath::Pi()*0.5),hist->GetYaxis()->FindBin(TMath::Pi()*0.5))/(hist->GetXaxis()->FindBin(1.0)-hist->GetXaxis()->FindBin(-1.0)),hist);
  }
  if(mode.CompareTo("ZYAM") ==0){
    RemovePlateau(hist->GetBinContent(hist->GetMinimumBin()),hist);
  }
  return hist;
}

TDirectory* resultsdirectory(TDirectory* motherdir, const char* name){
  motherdir->cd();
  TDirectory * outdir = motherdir->GetDirectory(name);
  if(outdir) motherdir->rmdir(name);
  outdir = motherdir->mkdir(name);
  return outdir;
}
TCanvas* Makecanvas(TH2D* hist, const char* name, Bool_t Stats,Bool_t remedge)
{
  TCanvas * Canvas = new TCanvas(name);
  Canvas->Divide(2,2);
  if(!Stats)hist->SetStats(0);
  if(remedge)hist->SetAxisRange(-1.65,1.65);
  hist->GetZaxis()->SetTitleSize(0.03);
  hist->GetZaxis()->SetTitleOffset(1.4);
  hist->GetZaxis()->SetLabelSize(0.025);
  Canvas->cd(1);
  hist->Draw("surf2");
  Canvas->cd(2);
  hist->Draw("surf3");
  Canvas->cd(3);
  Int_t binpihn = hist->GetYaxis()->FindBin(0+0.2);
  Int_t binpiln = hist->GetYaxis()->FindBin(0-0.2);
  TH1D* projY = hist->ProjectionX(Form("%s%s",hist->GetName(),"_nearside"),binpiln,binpihn);
  if(!Stats) projY->SetStats(0);
//   projY->GetYaxis()->SetRangeUser(0., 1.1*projY->GetBinContent(projY->GetMaximumBin()));
  projY->SetTitle("Integral over the near side peak with #Delta#eta_{12} = 0#pm 0.2");
  projY->GetYaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  projY->SetTitleSize(0.04,"x");
  projY->SetTitleOffset(1.05,"x");
  projY->SetTitleSize(0.03,"y");
  projY->SetTitleOffset(1.4,"y");
  projY->GetXaxis()->SetLabelSize(0.025);
  projY->GetYaxis()->SetLabelSize(0.025);
  projY->Draw("E");
  Canvas->cd(4);
  Int_t binpih = hist->GetYaxis()->FindBin(TMath::Pi()+0.2);
  Int_t binpil = hist->GetYaxis()->FindBin(TMath::Pi()-0.2);
  TH1D* projX = hist->ProjectionX(Form("%s%s",hist->GetName(),"_px"),binpil,binpih);
  if(!Stats) projX->SetStats(0);
  projX->SetTitle("Integral over a slice of the away side peak around with #Delta#eta_{12} = #pi#pm 0.2");
  projX->GetYaxis()->SetTitle(hist->GetZaxis()->GetTitle());
  projX->SetTitleSize(0.04,"x");
  projX->SetTitleOffset(1.05,"x");
  projX->SetTitleSize(0.03,"y");
  projX->SetTitleOffset(1.4,"y");
  projX->GetXaxis()->SetLabelSize(0.025);
  projX->GetYaxis()->SetLabelSize(0.025);
  projX->Draw("E");
  return Canvas;
}
TCanvas* Makecanvas(TH2D* histtopl, TH2D* histtopr, TH2D* histbotl, TH2D* histbotr, const char* name, Bool_t Stats)
{
  TCanvas * Canvas = new TCanvas(name);
  Canvas->Divide(2,2);
  Canvas->cd(1);
  histtopl->GetZaxis()->SetTitleOffset(1.3);
  histtopl->GetZaxis()->SetTitleSize(0.03);
  histtopl->GetZaxis()->SetLabelSize(0.025);
  histtopr->GetZaxis()->SetTitleOffset(1.3);
  histtopr->GetZaxis()->SetTitleSize(0.03);
  histtopr->GetZaxis()->SetLabelSize(0.025);
  histbotl->GetZaxis()->SetTitleOffset(1.3);
  histbotl->GetZaxis()->SetTitleSize(0.03);
  histbotl->GetZaxis()->SetLabelSize(0.025);
  histbotr->GetZaxis()->SetTitleOffset(1.3);
  histbotr->GetZaxis()->SetTitleSize(0.03);
  histbotr->GetZaxis()->SetLabelSize(0.025);
  if(!Stats)histtopl->SetStats(0);
  histtopl->Draw("surf3");
  Canvas->cd(2);
  if(!Stats)histtopr->SetStats(0);
  histtopr->Draw("surf3");
  Canvas->cd(3);
  if(!Stats)histbotl->SetStats(0);
  histbotl->Draw("surf3");
  Canvas->cd(4);
  if(!Stats)histbotr->SetStats(0);
  histbotr->Draw("surf3");
  return Canvas;
}

void canvasmaker(const char* histname, TObjArray* multdirlist, bool isPbPb){
  BinDirs * All=NULL;BinDirs * Bin1=NULL;BinDirs * Bin2=NULL;BinDirs * Bin3=NULL;BinDirs * Bin4=NULL;BinDirs * Bin5=NULL;BinDirs * Bin6=NULL;BinDirs * Bin7=NULL;
  if(!isPbPb){
   All = dynamic_cast<BinDirs*>(multdirlist->At(0));
  }
  else{
    Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(0));
    if(multdirlist->GetEntries()>1) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(1));
    if(multdirlist->GetEntries()>2) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(2));
    if(multdirlist->GetEntries()>3) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(3));
    if(multdirlist->GetEntries()>4) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(4));
    if(multdirlist->GetEntries()>5) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(5));
    if(multdirlist->GetEntries()>6) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  }
  if(TString(histname).CompareTo("DPhi_1_DPHI_2_far")==0){
    TH2D* hist;TH2D* histnear;TH2D* histmid; TH2D* histfar;
    TCanvas* canvas;
    if(!isPbPb){
      hist = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(All->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);All->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
    }
    else{
      hist = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin1->Samediv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin1->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      if(Bin2){
	hist = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin2->Samediv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin2->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin3){
	hist = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin3->Samediv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin3->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin4){
	hist = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin4->Samediv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin4->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin5){
	hist = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin5->Samediv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin5->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin6){
	hist = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin6->Samediv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin6->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin7){
	hist = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin7->Samediv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin7->Samediv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      
      hist = dynamic_cast<TH2D*>(Bin1->METAdiv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin1->METAdiv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin1->METAdiv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin1->METAdiv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin1->METAdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      if(Bin2){
	hist = dynamic_cast<TH2D*>(Bin2->METAdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin2->METAdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin2->METAdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin2->METAdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin2->METAdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin3){
	hist = dynamic_cast<TH2D*>(Bin3->METAdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin3->METAdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin3->METAdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin3->METAdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin3->METAdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin4){
	hist = dynamic_cast<TH2D*>(Bin4->METAdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin4->METAdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin4->METAdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin4->METAdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin4->METAdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin5){
	hist = dynamic_cast<TH2D*>(Bin5->METAdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin5->METAdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin5->METAdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin5->METAdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin5->METAdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin6){
	hist = dynamic_cast<TH2D*>(Bin6->METAdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin6->METAdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin6->METAdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin6->METAdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin6->METAdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin7){
	hist = dynamic_cast<TH2D*>(Bin7->METAdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin7->METAdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin7->METAdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin7->METAdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin7->METAdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      
            hist = dynamic_cast<TH2D*>(Bin1->META2div()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin1->META2div()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin1->META2div()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin1->META2div()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin1->META2div()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      if(Bin2){
	hist = dynamic_cast<TH2D*>(Bin2->META2div()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin2->META2div()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin2->META2div()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin2->META2div()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin2->META2div()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin3){
	hist = dynamic_cast<TH2D*>(Bin3->META2div()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin3->META2div()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin3->META2div()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin3->META2div()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin3->META2div()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin4){
	hist = dynamic_cast<TH2D*>(Bin4->META2div()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin4->META2div()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin4->META2div()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin4->META2div()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin4->META2div()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin5){
	hist = dynamic_cast<TH2D*>(Bin5->META2div()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin5->META2div()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin5->META2div()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin5->META2div()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin5->META2div()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin6){
	hist = dynamic_cast<TH2D*>(Bin6->META2div()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin6->META2div()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin6->META2div()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin6->META2div()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin6->META2div()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin7){
	hist = dynamic_cast<TH2D*>(Bin7->META2div()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin7->META2div()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin7->META2div()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin7->META2div()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin7->META2div()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      
      hist = dynamic_cast<TH2D*>(Bin1->METriggerdiv()->Get("DPhi_1_DPHI_2"));
      histnear = dynamic_cast<TH2D*>(Bin1->METriggerdiv()->Get("DPhi_1_DPHI_2_near"));
      histmid = dynamic_cast<TH2D*>(Bin1->METriggerdiv()->Get("DPhi_1_DPHI_2_mid"));
      histfar = dynamic_cast<TH2D*>(Bin1->METriggerdiv()->Get("DPhi_1_DPHI_2_far"));
      if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin1->METriggerdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      if(Bin2){
	hist = dynamic_cast<TH2D*>(Bin2->METriggerdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin2->METriggerdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin2->METriggerdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin2->METriggerdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin2->METriggerdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin3){
	hist = dynamic_cast<TH2D*>(Bin3->METriggerdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin3->METriggerdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin3->METriggerdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin3->METriggerdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin3->METriggerdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin4){
	hist = dynamic_cast<TH2D*>(Bin4->METriggerdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin4->METriggerdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin4->METriggerdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin4->METriggerdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin4->METriggerdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin5){
	hist = dynamic_cast<TH2D*>(Bin5->METriggerdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin5->METriggerdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin5->METriggerdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin5->METriggerdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin5->METriggerdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin6){
	hist = dynamic_cast<TH2D*>(Bin6->METriggerdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin6->METriggerdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin6->METriggerdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin6->METriggerdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin6->METriggerdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
      if(Bin7){
	hist = dynamic_cast<TH2D*>(Bin7->METriggerdiv()->Get("DPhi_1_DPHI_2"));
	histnear = dynamic_cast<TH2D*>(Bin7->METriggerdiv()->Get("DPhi_1_DPHI_2_near"));
	histmid = dynamic_cast<TH2D*>(Bin7->METriggerdiv()->Get("DPhi_1_DPHI_2_mid"));
	histfar = dynamic_cast<TH2D*>(Bin7->METriggerdiv()->Get("DPhi_1_DPHI_2_far"));
	if(hist&&histnear&&histmid&&histfar){canvas = Makecanvas(hist,histnear,histmid,histfar,"DPHIDPHI",false);Bin7->METriggerdiv()->cd();canvas->Write();delete hist;delete histnear;delete histmid; delete histfar;delete canvas;}
      }
    }
  }
  else if(TString(histname).Contains("DPhi_1_DEta")||TString(histname).Contains("DPhi_12A_DEta")||TString(histname).Contains("DPhi_DEta")){
    TH2D* hist;
    TCanvas* canvas;
    if(!isPbPb){
      hist = dynamic_cast<TH2D*>(All->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      if(All->SameDir("mixed_event")){
	hist = dynamic_cast<TH2D*>(All->SameDir("mixed_event")->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
      }
      hist = dynamic_cast<TH2D*>(All->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      if(All->META2()){
	hist = dynamic_cast<TH2D*>(All->META2div()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->META2div()->cd();canvas->Write();delete hist;delete canvas;}
      }
      hist = dynamic_cast<TH2D*>(All->METriggerdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false);All->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
    }
    else{
      hist = dynamic_cast<TH2D*>(Bin1->Samediv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin1->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
      if(Bin1->SameDir("mixed_event")){
	hist = dynamic_cast<TH2D*>(Bin1->SameDir("mixed_event")->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin1->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
      }
      hist = dynamic_cast<TH2D*>(Bin1->METAdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin1->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
      if(Bin1->META2()){
	hist = dynamic_cast<TH2D*>(Bin1->META2div()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin1->META2div()->cd();canvas->Write();delete hist;delete canvas;}
      }
      hist = dynamic_cast<TH2D*>(Bin1->METriggerdiv()->Get(histname));
      if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin1->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      if(Bin2){
	hist = dynamic_cast<TH2D*>(Bin2->Samediv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin2->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin2->SameDir("mixed_event")){
	  hist = dynamic_cast<TH2D*>(Bin2->SameDir("mixed_event")->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin2->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
	}
	hist = dynamic_cast<TH2D*>(Bin2->METAdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin2->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin2->META2()){
	  hist = dynamic_cast<TH2D*>(Bin2->META2div()->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin2->META2div()->cd();canvas->Write();delete hist;delete canvas;}	
	}	
	hist = dynamic_cast<TH2D*>(Bin2->METriggerdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin2->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      }
      if(Bin3){
	hist = dynamic_cast<TH2D*>(Bin3->Samediv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin3->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin3->SameDir("mixed_event")){
	  hist = dynamic_cast<TH2D*>(Bin3->SameDir("mixed_event")->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin3->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
	}
	hist = dynamic_cast<TH2D*>(Bin3->METAdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin3->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin3->META2()){
	  hist = dynamic_cast<TH2D*>(Bin3->META2div()->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin3->META2div()->cd();canvas->Write();delete hist;delete canvas;}	  
	}
	hist = dynamic_cast<TH2D*>(Bin3->METriggerdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin3->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      }
      if(Bin4){
	hist = dynamic_cast<TH2D*>(Bin4->Samediv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin4->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin4->SameDir("mixed_event")){
	  hist = dynamic_cast<TH2D*>(Bin4->SameDir("mixed_event")->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin4->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
	}
	hist = dynamic_cast<TH2D*>(Bin4->METAdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin4->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin4->META2()){
	  hist = dynamic_cast<TH2D*>(Bin4->META2div()->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin4->META2div()->cd();canvas->Write();delete hist;delete canvas;}	  
	}
	hist = dynamic_cast<TH2D*>(Bin4->METriggerdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin4->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      }
      if(Bin5){
	hist = dynamic_cast<TH2D*>(Bin5->Samediv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin5->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin5->SameDir("mixed_event")){
	  hist = dynamic_cast<TH2D*>(Bin5->SameDir("mixed_event")->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin5->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
	}
	hist = dynamic_cast<TH2D*>(Bin5->METAdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin5->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin5->META2()){
	  hist = dynamic_cast<TH2D*>(Bin5->META2div()->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin5->META2div()->cd();canvas->Write();delete hist;delete canvas;}	  
	}
	hist = dynamic_cast<TH2D*>(Bin5->METriggerdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin5->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      }
      if(Bin6){
	hist = dynamic_cast<TH2D*>(Bin6->Samediv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin6->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin6->SameDir("mixed_event")){
	  hist = dynamic_cast<TH2D*>(Bin6->SameDir("mixed_event")->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin6->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
	}
	hist = dynamic_cast<TH2D*>(Bin6->METAdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin6->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin6->META2()){
	  hist = dynamic_cast<TH2D*>(Bin6->META2div()->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin6->META2div()->cd();canvas->Write();delete hist;delete canvas;}	  
	}
	hist = dynamic_cast<TH2D*>(Bin6->METriggerdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin6->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      }
      if(Bin7){
	hist = dynamic_cast<TH2D*>(Bin7->Samediv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin7->Samediv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin7->SameDir("mixed_event")){
	  hist = dynamic_cast<TH2D*>(Bin7->SameDir("mixed_event")->Get(histname));
	  if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin7->SameDir("mixed_event")->cd();canvas->Write();delete hist;delete canvas;}
	}
	hist = dynamic_cast<TH2D*>(Bin7->METAdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin7->METAdiv()->cd();canvas->Write();delete hist;delete canvas;}
	if(Bin7->META2()){
	hist = dynamic_cast<TH2D*>(Bin7->META2div()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin7->META2div()->cd();canvas->Write();delete hist;delete canvas;}
	}
	hist = dynamic_cast<TH2D*>(Bin7->METriggerdiv()->Get(histname));
	if(hist){canvas = Makecanvas(hist, Form("%sCanvas",histname),false,isPbPb);Bin7->METriggerdiv()->cd();canvas->Write();delete hist;delete canvas;}
      }
    }
  }
  return;
}
void canvasruns(const char* canname, const char* title, TH1* hist, TDirectory* dir){
  dir->cd();
  TCanvas * canvas = new TCanvas(canname);
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  TH2D* hist2d = dynamic_cast<TH2D*>(hist);
  if(hist1d){
    hist1d->SetStats(kFALSE);
    hist1d->Draw("E");
    canvas->Write();
  }
  if(hist2d){
    hist2d->SetStats(kFALSE);
    hist2d->Draw("colz");
    canvas->Write();
  }
  delete canvas;
  return;
}
TCanvas* canvaspTbins(const char* canname, const char* title, TList* subdirlist){
  TCanvas * canvas = new TCanvas(canname);
  canvas->SetTitle(title);
  
  int nbins = subdirlist->GetEntries();
  for(int i=0;i<subdirlist->GetEntries();i++){
    if(TString(subdirlist->At(i)->GetName()).Contains("pTBin")) nbins -= 1;
  }
  
  if(nbins<4){
    canvas->Divide(2,2);
    canvas->SetTitle(title);
    for(int i=0;i<subdirlist->GetEntries();i++){
      if(i>nbins-1) continue;
      canvas->cd(i+1);
      canvas->SetTitle(subdirlist->At(i)->GetName());
    }
    
  }
  else{
    canvas->Divide(3,2);
    canvas->SetTitle(title);
    for(int i=0;i<subdirlist->GetEntries();i++){
      if(i>nbins-1) continue;
      canvas->cd(i+1);
      canvas->SetTitle(subdirlist->At(i)->GetName());
    }    
  }
  return canvas;
}

TList * GetMZDirectories(BinDirs* same){
  TList * keys = same->Same()->GetListOfKeys();
  TList * directories = new TList();
  for(int i=0; i<keys->GetEntries();i++){if(TString(keys->At(i)->GetName()).BeginsWith("BinM(")&&TString(keys->At(i)->GetName()).Contains("Z(")) directories->Add(same->Same()->GetDirectory(keys->At(i)->GetName()) );}
  return directories;
}
TStringToken GetHistTokens(TDirectory * dir){
  TList * histlist = dir->GetListOfKeys();
  TString * histlistS = new TString("");
  for(int i=0; i<histlist->GetEntries();i++)
  {
    TKey* key = dynamic_cast<TKey*>(histlist->At(i));
    if((!dynamic_cast<TH1D*>(key->ReadObj()))&&(!dynamic_cast<TH2D*>(key->ReadObj()))&&(!dynamic_cast<TH3D*>(key->ReadObj())))continue;//exclude non histograms
    if(dynamic_cast<TH3D*>(key->ReadObj()))continue;//exclude 3d histograms
    if(TString(key->GetName()).Contains("scaled"))continue;
    histlistS->Append(Form("%s ",histlist->At(i)->GetName()));
    }
  TStringToken histtokens(histlistS->Data()," ");
  return histtokens;
}
void savedircontent(BinDirs* same, double METAScale, double METriggerScale, const char* iteration){
  TPaveText * whatinthisbin;
  TParameter<double>* metascale = new TParameter<double>("METASCALE",METAScale);
  TParameter<double>* metriggerscale = new TParameter<double>("METRIGGERSCALE",METriggerScale);
  TCanvas * whatinthisbincanvas = new TCanvas("Thisversion");
  whatinthisbincanvas->cd();
  whatinthisbin = new TPaveText(.05,.1,.95,.8);
  whatinthisbin->AddText("The corrected same event histogram");
  whatinthisbin->AddText(Form("- %4.2f times META     ",METAScale));
  whatinthisbin->AddText(Form("- %4.2f times METrigger",METriggerScale));
  whatinthisbin->Draw();
  same->SameDir(iteration)->cd();
  whatinthisbincanvas->Write();
  metascale->Write();
  metriggerscale->Write();
  delete metascale;
  delete metriggerscale;
  delete whatinthisbincanvas;
  delete whatinthisbin;
}
void savedircontent(BinDirs* same, double METAScale, double META2Scale, double METriggerScale, const char* iteration){
  TPaveText * whatinthisbin;
  TParameter<double>* metascale = new TParameter<double>("METASCALE",METAScale);
  TParameter<double>* metriggerscale = new TParameter<double>("METRIGGERSCALE",METriggerScale);
  TParameter<double>* meta2scale = new TParameter<double>("META2SCALE",META2Scale);
  TCanvas * whatinthisbincanvas = new TCanvas("Thisversion");
  whatinthisbincanvas->cd();
  whatinthisbin = new TPaveText(.05,.1,.95,.8);
  whatinthisbin->AddText("The corrected same event histogram");
  whatinthisbin->AddText(Form("- %4.2f times META     ",METAScale));
  whatinthisbin->AddText(Form("- %4.2f times META2",META2Scale));
  whatinthisbin->AddText(Form("- %4.2f times METrigger",METriggerScale));
  whatinthisbin->Draw();
  same->SameDir(iteration)->cd();
  whatinthisbincanvas->Write();
  metascale->Write();
  meta2scale->Write();
  metriggerscale->Write();
  delete metascale;
  delete meta2scale;
  delete metriggerscale;
  delete whatinthisbincanvas;
  delete whatinthisbin;
}

void RemovePlateau(Double_t plateauheight, TH2D * hist){
  Double_t loccontent=0.0;
  Double_t locrelerr = 0.0;
  for(int x = hist->GetXaxis()->FindBin(-1.8)+1;x<=hist->GetXaxis()->FindBin(1.8)-1;x++){
    for(int y = 0; y<hist->GetNbinsY()+1;y++){
      locrelerr = hist->GetBinError(x,y);
      if(hist->GetBinContent(x,y)>1.0E-12) locrelerr = locrelerr/hist->GetBinContent(x,y);
      else locrelerr = 0.0;
      loccontent = hist->GetBinContent(x,y)-plateauheight;
      hist->SetBinContent(x,y,loccontent);
      hist->SetBinError(x,y,loccontent*locrelerr);
    }
  }
}

void CollectHistbinstats(const char* histname,TList * directories, TObjArray* multdirlist, Bool_t isPbPb){
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("bin_stats")->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  if(hist1d){cout << "TH1D " <<histname <<endl; CollectHistbs(hist1d,directories,multdirlist,isPbPb);}
  canvasmaker(histname,multdirlist);
}
void CollectHist(const char* histname,TList * directories, TObjArray* multdirlist,Bool_t pearsonserrors ,Bool_t collecdivfirst, Bool_t isPbPb ){
//   if(TString(histname).CompareTo("DPhi_DEta")==0) return;
  TH1* hist = dynamic_cast<TH1*>(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("divided")->Get(histname));
  TH1D* hist1d = dynamic_cast<TH1D*>(hist);
  TH2D* hist2d = dynamic_cast<TH2D*>(hist);
  TH3D* hist3d = dynamic_cast<TH3D*>(hist);
  if(hist1d){cout << "TH1D " <<histname<<endl;CollectHist(hist1d,directories,multdirlist, isPbPb);}
  if(hist2d){cout << "TH2D " <<histname<<endl;CollectHist(hist2d,directories,multdirlist, pearsonserrors,collecdivfirst,isPbPb,true);}
  if(hist3d){cout << "TH3D " <<histname<<endl;}//CollectHist(hist3d,directories,multdirlist);}
  canvasmaker(histname,multdirlist);
}


void CollectHistbs(TH1D* histo, TList * directories, TObjArray* multdirlist, Bool_t isPbPb){
  BinDirs * All=NULL;BinDirs * Bin1=NULL;BinDirs * Bin2=NULL;BinDirs * Bin3=NULL;BinDirs * Bin4=NULL;BinDirs * Bin5=NULL;BinDirs * Bin6=NULL;BinDirs * Bin7=NULL;
  if(!isPbPb){
   All = dynamic_cast<BinDirs*>(multdirlist->At(0));
  }
  else{
    Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(0));
    if(multdirlist->GetEntries()>1) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(1));
    if(multdirlist->GetEntries()>2) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(2));
    if(multdirlist->GetEntries()>3) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(3));
    if(multdirlist->GetEntries()>4) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(4));
    if(multdirlist->GetEntries()>5) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(5));
    if(multdirlist->GetEntries()>6) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  }
  
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //Bin 1:
  TH1D* histbin1 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH1D* histMETAbin1 	  	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH1D* histMETA2bin1 	  	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1META2"	,histo->GetName())));
  TH1D* histMETriggerbin1 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  //Bin 2:
  TH1D* histbin2 		= new TH1D("emptybin2","title",1,0,1);
  TH1D* histMETAbin2 		= new TH1D("emptybinMETA2","title",1,0,1);
  TH1D* histMETA2bin2 		= new TH1D("emptybinMETA22","title",1,0,1);
  TH1D* histMETriggerbin2 	= new TH1D("emptybinMETRIGGER2","title",1,0,1);
  if(Bin2){
    histbin2->Delete();
    histMETAbin2->Delete();histMETA2bin2->Delete();histMETriggerbin2->Delete();
    histbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histMETAbin2	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETA2bin2	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META2"	,histo->GetName())));
    histMETriggerbin2 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));    
  }
  //Bin 3:
  TH1D* histbin3 		= new TH1D("emptybin3","title",1,0,1);
  TH1D* histMETAbin3 		= new TH1D("emptybinMETA3","title",1,0,1);
  TH1D* histMETA2bin3 		= new TH1D("emptybinMETA33","title",1,0,1);
  TH1D* histMETriggerbin3 	= new TH1D("emptybinMETRIGGER3","title",1,0,1);
  if(Bin3){
    histbin3->Delete();histMETAbin3->Delete();histMETA2bin3->Delete();histMETriggerbin3->Delete();
    histbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histMETAbin3	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETA2bin3	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META2"	,histo->GetName())));
    histMETriggerbin3 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  }
  //Bin 4:
  TH1D* histbin4 		= new TH1D("emptybin4","title",1,0,1);
  TH1D* histMETAbin4 		= new TH1D("emptybinMETA4","title",1,0,1);
  TH1D* histMETA2bin4 		= new TH1D("emptybinMETA44","title",1,0,1);
  TH1D* histMETriggerbin4 	= new TH1D("emptybinMETRIGGER4","title",1,0,1);
  if(Bin4){
    histbin4->Delete();histMETAbin4->Delete();histMETA2bin4->Delete();histMETriggerbin4->Delete();
    histbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histMETAbin4	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETA2bin4	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META2"	,histo->GetName())));
    histMETriggerbin4 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  }
  //Bin 5:
  TH1D* histbin5 		= new TH1D("emptybin5","title",1,0,1);
  TH1D* histMETAbin5 		= new TH1D("emptybinMETA5","title",1,0,1);
  TH1D* histMETA2bin5 		= new TH1D("emptybinMETA55","title",1,0,1);
  TH1D* histMETriggerbin5 	= new TH1D("emptybinMETRIGGER5","title",1,0,1);
  if(Bin5){
    histbin5->Delete();histMETAbin5->Delete();histMETA2bin5->Delete();histMETriggerbin5->Delete();
    histbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histMETAbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETA2bin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META2"	,histo->GetName())));
    histMETriggerbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  }
  //Bin 6:
  TH1D* histbin6 		= new TH1D("emptybin6","title",1,0,1);
  TH1D* histMETAbin6 		= new TH1D("emptybinMETA6","title",1,0,1);
  TH1D* histMETA2bin6 		= new TH1D("emptybinMETA66","title",1,0,1);
  TH1D* histMETriggerbin6 	= new TH1D("emptybinMETRIGGER6","title",1,0,1);
  if(Bin6){
    histbin6->Delete();histMETAbin6->Delete();histMETA2bin6->Delete();histMETriggerbin6->Delete();
    histbin6		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6"		,histo->GetName())));
    histMETAbin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6META"	,histo->GetName())));
    histMETA2bin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6META2"	,histo->GetName())));
    histMETriggerbin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6METrigger"	,histo->GetName())));
  }
  //Bin 7:
  TH1D* histbin7 		= new TH1D("emptybin7","title",1,0,1);
  TH1D* histMETAbin7 		= new TH1D("emptybinMETA7","title",1,0,1);
  TH1D* histMETA2bin7 		= new TH1D("emptybinMETA77","title",1,0,1);
  TH1D* histMETriggerbin7 	= new TH1D("emptybinMETRIGGER7","title",1,0,1);
  if(Bin7){
    histbin7->Delete();histMETAbin7->Delete();histMETA2bin7->Delete();histMETriggerbin7->Delete();
    histbin7		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7"		,histo->GetName())));
    histMETAbin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7META"	,histo->GetName())));
    histMETA2bin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7META2"	,histo->GetName())));
    histMETriggerbin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7METrigger"	,histo->GetName())));
  }  
  //Doubles to hold the values until they are put into the hists.
  TH1D* histtmp			= new TH1D("tmp","title",1,0,1);
  TH1D* histMETAtmp		= new TH1D("tmpMETA","title",1,0,1);
  TH1D* histMETA2tmp		= new TH1D("tmpMETA2","title",1,0,1);
  TH1D* histMETriggertmp	= new TH1D("tmpMETrigger","title",1,0,1);
  //loop over the bins:
  for(int i=0;i<directories->GetEntries();i++){
    int Mbin = 0;
    BinDirs* bin = dynamic_cast<BinDirs*>(directories->At(i));
    //find the Multiplicity bin we are in
    for(int j = 1;j<multdirlist->GetEntries();j++){
      if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
      if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo("BinM(0.00)->(10.00)")){
	if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(0.00)->(5.00)")==0) Mbin = j;
	if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(5.00)->(10.00)")==0)Mbin = j;
      }
    }
    //test if the relevant histogram exists in this bin:
    if(!dynamic_cast<TH1D*>(bin->SameDir("bin_stats")->Get(histo->GetName())))
    { 	cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
      continue;}
    //extract histogram in relevant bin:
    if(i==0){histtmp->Delete();    histMETAtmp->Delete();histMETA2tmp->Delete(); histMETriggertmp->Delete();}
    histtmp 		= dynamic_cast<TH1D*>(bin->SameDir("bin_stats")->Get(histo->GetName()));
    histMETAtmp 	= dynamic_cast<TH1D*>(bin->METADir("bin_stats")->Get(histo->GetName()));
    histMETriggertmp 	= dynamic_cast<TH1D*>(bin->METriggerDir("bin_stats")->Get(histo->GetName()));
    if(bin->META()) histMETA2tmp 	= dynamic_cast<TH1D*>(bin->META2Dir("bin_stats")->Get(histo->GetName()));
//     hist->Add(histtmp); histMETA->Add(histMETAtmp); histMETrigger->Add(histMETriggertmp);
    if(Mbin == 1||!isPbPb){	histbin1->Add(histtmp);histMETAbin1->Add(histMETAtmp); histMETA2bin1->Add(histMETA2tmp); histMETriggerbin1->Add(histMETriggertmp);}
    else if(Mbin == 2&&Bin2){ histbin2->Add(histtmp);histMETAbin2->Add(histMETAtmp); histMETA2bin2->Add(histMETA2tmp); histMETriggerbin2->Add(histMETriggertmp);}	      
    else if(Mbin == 3&&Bin3){ histbin3->Add(histtmp);histMETAbin3->Add(histMETAtmp); histMETA2bin3->Add(histMETA2tmp); histMETriggerbin3->Add(histMETriggertmp);}
    else if(Mbin == 4&&Bin4){ histbin4->Add(histtmp);histMETAbin4->Add(histMETAtmp); histMETA2bin4->Add(histMETA2tmp); histMETriggerbin4->Add(histMETriggertmp);}	      
    else if(Mbin == 5&&Bin5){ histbin5->Add(histtmp);histMETAbin5->Add(histMETAtmp); histMETA2bin5->Add(histMETA2tmp); histMETriggerbin5->Add(histMETriggertmp);}
    else if(Mbin == 6&&Bin6){ histbin6->Add(histtmp);histMETAbin6->Add(histMETAtmp); histMETA2bin6->Add(histMETA2tmp); histMETriggerbin6->Add(histMETriggertmp);}
    else if(Mbin == 7&&Bin7){ histbin7->Add(histtmp);histMETAbin7->Add(histMETAtmp); histMETA2bin7->Add(histMETA2tmp); histMETriggerbin7->Add(histMETriggertmp);}

    
  }
  //save the histograms in the relevant directories:
  if(!isPbPb){
    All->Same()->GetDirectory("bin_stats")->cd();
    histbin1->Write(histo->GetName());
    All->META()->GetDirectory("bin_stats")->cd();
    histMETAbin1->Write(histo->GetName());
    if(All->META2()){
      All->META2()->GetDirectory("bin_stats")->cd();
      histMETA2bin1->Write(histo->GetName());    
    }
    All->METrigger()->GetDirectory("bin_stats")->cd();
    histMETriggerbin1->Write(histo->GetName());
  }
  else{
    Bin1->Same()->GetDirectory("bin_stats")->cd();
    histbin1->Write(histo->GetName());
    Bin1->META()->GetDirectory("bin_stats")->cd();
    histMETAbin1->Write(histo->GetName());
    if(Bin1->META2()){
      Bin1->META2()->GetDirectory("bin_stats")->cd();
      histMETA2bin1->Write(histo->GetName());
    }
    Bin1->METrigger()->GetDirectory("bin_stats")->cd();
    histMETriggerbin1->Write(histo->GetName());
    if(Bin2){
      Bin2->Same()->GetDirectory("bin_stats")->cd();
      histbin2->Write(histo->GetName());
      Bin2->META()->GetDirectory("bin_stats")->cd();
      histMETAbin2->Write(histo->GetName());
      if(Bin2->META2()){
	Bin2->META2()->GetDirectory("bin_stats")->cd();
	histMETA2bin2->Write(histo->GetName());
      }
      Bin2->METrigger()->GetDirectory("bin_stats")->cd();
      histMETriggerbin2->Write(histo->GetName());
    }
    if(Bin3){
      Bin3->Same()->GetDirectory("bin_stats")->cd();
      histbin3->Write(histo->GetName());
      Bin3->META()->GetDirectory("bin_stats")->cd();
      histMETAbin3->Write(histo->GetName());
      if(Bin3->META2()){
	Bin3->META2()->GetDirectory("bin_stats")->cd();
	histMETA2bin3->Write(histo->GetName());
      }
      Bin3->METrigger()->GetDirectory("bin_stats")->cd();
      histMETriggerbin3->Write(histo->GetName());
    }
    if(Bin4){
      Bin4->Same()->GetDirectory("bin_stats")->cd();
      histbin4->Write(histo->GetName());
      Bin4->META()->GetDirectory("bin_stats")->cd();
      histMETAbin4->Write(histo->GetName());
      if(Bin4->META2()){
	Bin4->META2()->GetDirectory("bin_stats")->cd();
	histMETA2bin4->Write(histo->GetName());
      }
      Bin4->METrigger()->GetDirectory("bin_stats")->cd();
      histMETriggerbin4->Write(histo->GetName());
    }
    if(Bin5){
      Bin5->Same()->GetDirectory("bin_stats")->cd();
      histbin5->Write(histo->GetName());
      Bin5->META()->GetDirectory("bin_stats")->cd();
      histMETAbin5->Write(histo->GetName());
      if(Bin5->META2()){
	Bin5->META2()->GetDirectory("bin_stats")->cd();
	histMETA2bin5->Write(histo->GetName());
      }
      Bin5->METrigger()->GetDirectory("bin_stats")->cd();
      histMETriggerbin5->Write(histo->GetName());
    }
    if(Bin6){
      Bin6->Same()->GetDirectory("bin_stats")->cd();
      histbin6->Write(histo->GetName());
      Bin6->META()->GetDirectory("bin_stats")->cd();
      histMETAbin6->Write(histo->GetName());
      if(Bin6->META2()){
	Bin6->META2()->GetDirectory("bin_stats")->cd();
	histMETA2bin6->Write(histo->GetName());
      }
      Bin6->METrigger()->GetDirectory("bin_stats")->cd();
      histMETriggerbin6->Write(histo->GetName());
    }
    if(Bin7){
      Bin7->Same()->GetDirectory("bin_stats")->cd();
      histbin7->Write(histo->GetName());
      Bin7->META()->GetDirectory("bin_stats")->cd();
      histMETAbin7->Write(histo->GetName());
      if(Bin7->META2()){
	Bin7->META2()->GetDirectory("bin_stats")->cd();
	histMETA2bin7->Write(histo->GetName());
      }
      Bin7->METrigger()->GetDirectory("bin_stats")->cd();
      histMETriggerbin7->Write(histo->GetName());
    }
  }
  delete histbin1; delete histMETAbin1; delete histMETA2bin1; delete histMETriggerbin1;   
  delete histbin2; delete histMETAbin2; delete histMETA2bin2; delete histMETriggerbin2;   
  delete histbin3; delete histMETAbin3; delete histMETA2bin3; delete histMETriggerbin3;   
  delete histbin4; delete histMETAbin4; delete histMETA2bin4; delete histMETriggerbin4;   
  delete histbin5; delete histMETAbin5; delete histMETA2bin5; delete histMETriggerbin5;   
  delete histbin6; delete histMETAbin6; delete histMETA2bin6; delete histMETriggerbin6;   
  delete histbin7; delete histMETAbin7; delete histMETA2bin7; delete histMETriggerbin7;   
}
void CollectHist(TH1D* histo, TList * directories, TObjArray* multdirlist, Bool_t isPbPb){
  BinDirs * All=NULL;BinDirs * Bin1=NULL;BinDirs * Bin2=NULL;BinDirs * Bin3=NULL;BinDirs * Bin4=NULL;BinDirs * Bin5=NULL;BinDirs * Bin6=NULL;BinDirs * Bin7=NULL;
if(!isPbPb){
   Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(0));
  }
  else{
    Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(0));
    if(multdirlist->GetEntries()>1) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(1));
    if(multdirlist->GetEntries()>2) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(2));
    if(multdirlist->GetEntries()>3) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(3));
    if(multdirlist->GetEntries()>4) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(4));
    if(multdirlist->GetEntries()>5) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(5));
    if(multdirlist->GetEntries()>6) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  }
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //Bin 1:
  TH1D* histbin1 		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH1D* histMETAbin1 	  	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH1D* histMETA2bin1 	  	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1META2"	,histo->GetName())));
  TH1D* histMETriggerbin1 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
   //Bin 2:
  TH1D* histbin2 		= new TH1D("emptybin2","title",1,0,1);
  TH1D* histMETAbin2 		= new TH1D("emptybinMETA2","title",1,0,1);
  TH1D* histMETA2bin2 		= new TH1D("emptybinMETA22","title",1,0,1);
  TH1D* histMETriggerbin2 	= new TH1D("emptybinMETRIGGER2","title",1,0,1);
  if(Bin2){
    histbin2->Delete();histMETAbin2->Delete();histMETA2bin2->Delete();histMETriggerbin2->Delete();
    histbin2		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histMETAbin2	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETA2bin2	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2META2"	,histo->GetName())));
    histMETriggerbin2 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));    
  }
  //Bin 3:
  TH1D* histbin3 		= new TH1D("emptybin3","title",1,0,1);
  TH1D* histMETAbin3 		= new TH1D("emptybinMETA3","title",1,0,1);
  TH1D* histMETA2bin3 		= new TH1D("emptybinMETA33","title",1,0,1);
  TH1D* histMETriggerbin3 	= new TH1D("emptybinMETRIGGER3","title",1,0,1);
  if(Bin3){
    histbin3->Delete();histMETAbin3->Delete();histMETA2bin3->Delete();histMETriggerbin3->Delete();
    histbin3		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histMETAbin3	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETA2bin3	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3META2"	,histo->GetName())));
    histMETriggerbin3 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  }
  //Bin 4:
  TH1D* histbin4 		= new TH1D("emptybin4","title",1,0,1);
  TH1D* histMETAbin4 		= new TH1D("emptybinMETA4","title",1,0,1);
  TH1D* histMETA2bin4 		= new TH1D("emptybinMETA44","title",1,0,1);
  TH1D* histMETriggerbin4 	= new TH1D("emptybinMETRIGGER4","title",1,0,1);
  if(Bin4){
    histbin4->Delete();histMETAbin4->Delete();histMETA2bin4->Delete();histMETriggerbin4->Delete();
    histbin4		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histMETAbin4	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETA2bin4	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4META2"	,histo->GetName())));
    histMETriggerbin4 	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  }
  //Bin 5:
  TH1D* histbin5 		= new TH1D("emptybin5","title",1,0,1);
  TH1D* histMETAbin5 		= new TH1D("emptybinMETA5","title",1,0,1);
  TH1D* histMETA2bin5 		= new TH1D("emptybinMETA55","title",1,0,1);
  TH1D* histMETriggerbin5 	= new TH1D("emptybinMETRIGGER5","title",1,0,1);
  if(Bin5){
    histbin5->Delete();histMETAbin5->Delete();histMETA2bin5->Delete();histMETriggerbin5->Delete();
    histbin5		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histMETAbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETA2bin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5META2"	,histo->GetName())));
    histMETriggerbin5	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  }
  //Bin 6:
  TH1D* histbin6 		= new TH1D("emptybin6","title",1,0,1);
  TH1D* histMETAbin6 		= new TH1D("emptybinMETA6","title",1,0,1);
  TH1D* histMETA2bin6 		= new TH1D("emptybinMETA66","title",1,0,1);
  TH1D* histMETriggerbin6 	= new TH1D("emptybinMETRIGGER6","title",1,0,1);
  if(Bin6){
    histbin6->Delete();histMETAbin6->Delete();histMETA2bin6->Delete();histMETriggerbin6->Delete();
    histbin6		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6"		,histo->GetName())));
    histMETAbin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6META"	,histo->GetName())));
    histMETA2bin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6META2"	,histo->GetName())));
    histMETriggerbin6	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin6METrigger"	,histo->GetName())));
  }
  //Bin 7:
  TH1D* histbin7 		= new TH1D("emptybin7","title",1,0,1);
  TH1D* histMETAbin7 		= new TH1D("emptybinMETA7","title",1,0,1);
  TH1D* histMETA2bin7 		= new TH1D("emptybinMETA77","title",1,0,1);
  TH1D* histMETriggerbin7 	= new TH1D("emptybinMETRIGGER7","title",1,0,1);
  if(Bin7){
    histbin7->Delete();histMETAbin7->Delete();histMETA2bin7->Delete();histMETriggerbin7->Delete();
    histbin7		= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7"		,histo->GetName())));
    histMETAbin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7META"	,histo->GetName())));
    histMETA2bin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7META2"	,histo->GetName())));
    histMETriggerbin7	= dynamic_cast<TH1D*>( histo->Clone(Form("%sdivbin7METrigger"	,histo->GetName())));
  }  
  //Doubles to hold the values until they are put into the hists.
  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;Double_t bincontlMETA2   = 0.0;Double_t binerrorlMETA2  = 0.0;
  Double_t BinContent = 0.0;Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETA2bin1 = 0.0;Double_t BinErrorMETA2bin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETA2bin2 = 0.0;Double_t BinErrorMETA2bin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETA2bin3 = 0.0;Double_t BinErrorMETA2bin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETA2bin4 = 0.0;Double_t BinErrorMETA2bin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETA2bin5 = 0.0;Double_t BinErrorMETA2bin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  Double_t BinContentbin6 = 0.0;Double_t BinErrorbin6   = 0.0;Double_t BinContentMETAbin6 = 0.0;Double_t BinErrorMETAbin6   = 0.0;Double_t BinContentMETA2bin6 = 0.0;Double_t BinErrorMETA2bin6   = 0.0;Double_t BinContentMETriggerbin6 = 0.0;Double_t BinErrorMEtriggerbin6   = 0.0; 
  Double_t BinContentbin7 = 0.0;Double_t BinErrorbin7   = 0.0;Double_t BinContentMETAbin7 = 0.0;Double_t BinErrorMETAbin7   = 0.0;Double_t BinContentMETA2bin7 = 0.0;Double_t BinErrorMETA2bin7   = 0.0;Double_t BinContentMETriggerbin7 = 0.0;Double_t BinErrorMEtriggerbin7   = 0.0; 

  //loop over the bins:
  for(int x=0;x<=histo->GetNbinsX()+1;x++){
    //loop over all Multiplicity-Vertex bins:
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
      BinDirs* bin = dynamic_cast<BinDirs*>(directories->At(i));
      //find the Multiplicity bin we are in
      for(int j = 1;j<multdirlist->GetEntries();j++){
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo("BinM(0.00)->(10.00)")){
	  if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(0.00)->(5.00)")==0) Mbin = j;
	  if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(5.00)->(10.00)")==0)Mbin = j;
	}
      }
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH1D*>(bin->SameDir("divided")->Get(histo->GetName())))
      {cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
	continue;}
      //extract bin content and error in the relevant bin:
      bincontl 			= dynamic_cast<TH1D*>(bin->SameDir("divided")->Get(histo->GetName()))->GetBinContent(x);
      binerrorl 		= dynamic_cast<TH1D*>(bin->SameDir("divided")->Get(histo->GetName()))->GetBinError(x);	
      bincontlMETA 		= dynamic_cast<TH1D*>(bin->METADir("divided")->Get(histo->GetName()))->GetBinContent(x);
      binerrorlMETA 		= dynamic_cast<TH1D*>(bin->METADir("divided")->Get(histo->GetName()))->GetBinError(x);	
      if(bin->META2()){
	bincontlMETA2 		= dynamic_cast<TH1D*>(bin->META2Dir("divided")->Get(histo->GetName()))->GetBinContent(x);
	binerrorlMETA2 		= dynamic_cast<TH1D*>(bin->META2Dir("divided")->Get(histo->GetName()))->GetBinError(x);	
      }
      bincontlMEtrigger 	= dynamic_cast<TH1D*>(bin->METriggerDir("divided")->Get(histo->GetName()))->GetBinContent(x);
      binerrorlMETrigger 	= dynamic_cast<TH1D*>(bin->METriggerDir("divided")->Get(histo->GetName()))->GetBinError(x);  
      if(bincontl>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContent += bincontl/(binerrorl*binerrorl);
	if(Mbin == 1||!isPbPb){		BinContentbin1 += bincontl/(binerrorl*binerrorl); BinErrorbin1   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 2&&isPbPb){	BinContentbin2 += bincontl/(binerrorl*binerrorl); BinErrorbin2   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 3&&isPbPb){	BinContentbin3 += bincontl/(binerrorl*binerrorl); BinErrorbin3   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 4&&isPbPb){	BinContentbin4 += bincontl/(binerrorl*binerrorl); BinErrorbin4   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 5&&isPbPb){ 	BinContentbin5 += bincontl/(binerrorl*binerrorl); BinErrorbin5   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 6&&isPbPb){ 	BinContentbin6 += bincontl/(binerrorl*binerrorl); BinErrorbin6   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 7&&isPbPb){ 	BinContentbin7 += bincontl/(binerrorl*binerrorl); BinErrorbin7   += 1.0/(binerrorl*binerrorl);}      }
      else{
	if(Mbin == 1||!isPbPb)		BinErrorbin1   += 1.0;
	else if(Mbin == 2&&isPbPb)	BinErrorbin2   += 1.0;
	else if(Mbin == 3&&isPbPb)	BinErrorbin3   += 1.0;
	else if(Mbin == 4&&isPbPb)	BinErrorbin4   += 1.0;	      
	else if(Mbin == 5&&isPbPb)	BinErrorbin5   += 1.0;
	else if(Mbin == 6&&isPbPb)	BinErrorbin6   += 1.0;
	else if(Mbin == 7&&isPbPb)	BinErrorbin7   += 1.0;
      }
        if(bincontlMETA>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	if(Mbin == 1||!isPbPb){			BinContentMETAbin1 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin1   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 2&&isPbPb){		BinContentMETAbin2 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin2   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 3&&isPbPb){		BinContentMETAbin3 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin3   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 4&&isPbPb){		BinContentMETAbin4 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin4   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 5&&isPbPb){		BinContentMETAbin5 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin5   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 6&&isPbPb){		BinContentMETAbin6 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin6   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 7&&isPbPb){		BinContentMETAbin7 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin7   += 1.0/(binerrorlMETA*binerrorlMETA);}
	    }
      else{
	if(Mbin == 1||!isPbPb)			BinErrorMETAbin1   += 1.0;
	else if(Mbin == 2&&isPbPb)		BinErrorMETAbin2   += 1.0;	      
	else if(Mbin == 3&&isPbPb)		BinErrorMETAbin3   += 1.0;
	else if(Mbin == 4&&isPbPb)		BinErrorMETAbin4   += 1.0;
	else if(Mbin == 5&&isPbPb)		BinErrorMETAbin5   += 1.0;
	else if(Mbin == 6&&isPbPb)		BinErrorMETAbin6   += 1.0;
	else if(Mbin == 7&&isPbPb)		BinErrorMETAbin7   += 1.0;
      }
        if(bincontlMETA2>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	if(Mbin == 1||!isPbPb){			BinContentMETA2bin1 += bincontlMETA2/(binerrorlMETA2*binerrorlMETA2); BinErrorMETA2bin1   += 1.0/(binerrorlMETA2*binerrorlMETA2);}
	else if(Mbin == 2&&isPbPb){		BinContentMETA2bin2 += bincontlMETA2/(binerrorlMETA2*binerrorlMETA2); BinErrorMETA2bin2   += 1.0/(binerrorlMETA2*binerrorlMETA2);}	      
	else if(Mbin == 3&&isPbPb){		BinContentMETA2bin3 += bincontlMETA2/(binerrorlMETA2*binerrorlMETA2); BinErrorMETA2bin3   += 1.0/(binerrorlMETA2*binerrorlMETA2);}
	else if(Mbin == 4&&isPbPb){		BinContentMETA2bin4 += bincontlMETA2/(binerrorlMETA2*binerrorlMETA2); BinErrorMETA2bin4   += 1.0/(binerrorlMETA2*binerrorlMETA2);}	      
	else if(Mbin == 5&&isPbPb){		BinContentMETA2bin5 += bincontlMETA2/(binerrorlMETA2*binerrorlMETA2); BinErrorMETA2bin5   += 1.0/(binerrorlMETA2*binerrorlMETA2);}
	else if(Mbin == 6&&isPbPb){		BinContentMETA2bin6 += bincontlMETA2/(binerrorlMETA2*binerrorlMETA2); BinErrorMETA2bin6   += 1.0/(binerrorlMETA2*binerrorlMETA2);}
	else if(Mbin == 7&&isPbPb){		BinContentMETA2bin7 += bincontlMETA2/(binerrorlMETA2*binerrorlMETA2); BinErrorMETA2bin7   += 1.0/(binerrorlMETA2*binerrorlMETA2);}
	}
      else{
	if(Mbin == 1||!isPbPb)			BinErrorMETA2bin1   += 1.0;
	else if(Mbin == 2&&isPbPb)		BinErrorMETA2bin2   += 1.0;	      
	else if(Mbin == 3&&isPbPb)		BinErrorMETA2bin3   += 1.0;
	else if(Mbin == 4&&isPbPb)		BinErrorMETA2bin4   += 1.0;
	else if(Mbin == 5&&isPbPb)		BinErrorMETA2bin5   += 1.0;
	else if(Mbin == 6&&isPbPb)		BinErrorMETA2bin6   += 1.0;
	else if(Mbin == 7&&isPbPb)		BinErrorMETA2bin7   += 1.0;
      }
      if(bincontlMEtrigger>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	if(Mbin == 1||!isPbPb){			BinContentMETriggerbin1 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin1   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 2&&isPbPb){		BinContentMETriggerbin2 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin2   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 3&&isPbPb){		BinContentMETriggerbin3 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin3   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 4&&isPbPb){		BinContentMETriggerbin4 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin4   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 5&&isPbPb){		BinContentMETriggerbin5 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin5   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 6&&isPbPb){		BinContentMETriggerbin6 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin6   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 7&&isPbPb){		BinContentMETriggerbin7 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin7   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
      }
      else{
	if(Mbin == 1||!isPbPb)			BinErrorMEtriggerbin1   += 1.0;
	else if(Mbin == 2&&isPbPb)		BinErrorMEtriggerbin2   += 1.0;
	else if(Mbin == 3&&isPbPb)		BinErrorMEtriggerbin3   += 1.0;
	else if(Mbin == 4&&isPbPb)		BinErrorMEtriggerbin4   += 1.0;
	else if(Mbin == 5&&isPbPb)		BinErrorMEtriggerbin5   += 1.0;
	else if(Mbin == 6&&isPbPb)		BinErrorMEtriggerbin6   += 1.0;
	else if(Mbin == 7&&isPbPb)		BinErrorMEtriggerbin7   += 1.0;
      }
      //reset the local bins to zero just to be shure:
      bincontl  = 0.0;		binerrorl = 0.0;
      bincontlMETA =0.0;	binerrorlMETA=0.0;
      bincontlMETA2 =0.0;	binerrorlMETA2=0.0;
      bincontlMEtrigger=0.0;	binerrorlMETrigger=0.0;
    }//end loop over M-V bins
    //normalize the bin and the error for same:    
    if(BinErrorbin1>1.0e-10){		BinContentbin1 = BinContentbin1/BinErrorbin1;					BinErrorbin1 = 1.0/TMath::Sqrt(BinErrorbin1);			}
    else{				BinContentbin1=0.0;								BinErrorbin1=0.0;						}
    if(BinErrorbin2>1.0e-10){		BinContentbin2 = BinContentbin2/BinErrorbin2;					BinErrorbin2 = 1.0/TMath::Sqrt(BinErrorbin2);			}
    else{				BinContentbin2=0.0;								BinErrorbin2=0.0;						}
    if(BinErrorbin3>1.0e-10){		BinContentbin3 = BinContentbin3/BinErrorbin3;					BinErrorbin3 = 1.0/TMath::Sqrt(BinErrorbin3);			}
    else{				BinContentbin3=0.0;								BinErrorbin3=0.0;						}
    if(BinErrorbin4>1.0e-10){		BinContentbin4 = BinContentbin4/BinErrorbin4;					BinErrorbin4 = 1.0/TMath::Sqrt(BinErrorbin4);			}
    else{				BinContentbin4=0.0;								BinErrorbin4=0.0;						}
    if(BinErrorbin5>1.0e-10){		BinContentbin5 = BinContentbin5/BinErrorbin5;					BinErrorbin5 = 1.0/TMath::Sqrt(BinErrorbin5);			}
    else{				BinContentbin5=0.0;								BinErrorbin5=0.0;						}
    if(BinErrorbin6>1.0e-10){		BinContentbin6 = BinContentbin6/BinErrorbin6;					BinErrorbin6 = 1.0/TMath::Sqrt(BinErrorbin6);			}
    else{				BinContentbin6=0.0;								BinErrorbin6=0.0;						}
    if(BinErrorbin7>1.0e-10){		BinContentbin7 = BinContentbin7/BinErrorbin7;					BinErrorbin7 = 1.0/TMath::Sqrt(BinErrorbin7);			}
    else{				BinContentbin7=0.0;								BinErrorbin7=0.0;						}
    //normalize the bin and the error for META:    
    if(BinErrorMETAbin1>1.0e-10){	BinContentMETAbin1 = BinContentMETAbin1/BinErrorMETAbin1;			BinErrorMETAbin1 = 1.0/TMath::Sqrt(BinErrorMETAbin1);		}
    else{				BinContentMETAbin1=0.0;								BinErrorMETAbin1=0.0;						}
    if(BinErrorMETAbin2>1.0e-10){	BinContentMETAbin2 = BinContentMETAbin2/BinErrorMETAbin2;			BinErrorMETAbin2 = 1.0/TMath::Sqrt(BinErrorMETAbin2);		}
    else{				BinContentMETAbin2=0.0;								BinErrorMETAbin2=0.0;						}
    if(BinErrorMETAbin3>1.0e-10){	BinContentMETAbin3 = BinContentMETAbin3/BinErrorMETAbin3;			BinErrorMETAbin3 = 1.0/TMath::Sqrt(BinErrorMETAbin3);		}
    else{				BinContentMETAbin3=0.0;								BinErrorMETAbin3=0.0;						}
    if(BinErrorMETAbin4>1.0e-10){	BinContentMETAbin4 = BinContentMETAbin4/BinErrorMETAbin4;			BinErrorMETAbin4 = 1.0/TMath::Sqrt(BinErrorMETAbin4);		}
    else{				BinContentMETAbin4=0.0;								BinErrorMETAbin4=0.0;						}
    if(BinErrorMETAbin5>1.0e-10){	BinContentMETAbin5 = BinContentMETAbin5/BinErrorMETAbin5;			BinErrorMETAbin5 = 1.0/TMath::Sqrt(BinErrorMETAbin5);		}
    else{				BinContentMETAbin5=0.0;								BinErrorMETAbin5=0.0;						}
    if(BinErrorMETAbin6>1.0e-10){	BinContentMETAbin6 = BinContentMETAbin6/BinErrorMETAbin6;			BinErrorMETAbin6 = 1.0/TMath::Sqrt(BinErrorMETAbin6);		}
    else{				BinContentMETAbin6=0.0;								BinErrorMETAbin6=0.0;						}
    if(BinErrorMETAbin7>1.0e-10){	BinContentMETAbin7 = BinContentMETAbin7/BinErrorMETAbin7;			BinErrorMETAbin7 = 1.0/TMath::Sqrt(BinErrorMETAbin7);		}
    else{				BinContentMETAbin7=0.0;								BinErrorMETAbin7=0.0;						}
    //normalize the bin and the error for META2:    
    if(BinErrorMETA2bin1>1.0e-10){	BinContentMETA2bin1 = BinContentMETA2bin1/BinErrorMETA2bin1;			BinErrorMETA2bin1 = 1.0/TMath::Sqrt(BinErrorMETA2bin1);		}
    else{				BinContentMETA2bin1=0.0;								BinErrorMETA2bin1=0.0;					}
    if(BinErrorMETA2bin2>1.0e-10){	BinContentMETA2bin2 = BinContentMETA2bin2/BinErrorMETA2bin2;			BinErrorMETA2bin2 = 1.0/TMath::Sqrt(BinErrorMETA2bin2);		}
    else{				BinContentMETA2bin2=0.0;								BinErrorMETA2bin2=0.0;					}
    if(BinErrorMETA2bin3>1.0e-10){	BinContentMETA2bin3 = BinContentMETA2bin3/BinErrorMETA2bin3;			BinErrorMETA2bin3 = 1.0/TMath::Sqrt(BinErrorMETA2bin3);		}
    else{				BinContentMETA2bin3=0.0;								BinErrorMETA2bin3=0.0;					}
    if(BinErrorMETA2bin4>1.0e-10){	BinContentMETA2bin4 = BinContentMETA2bin4/BinErrorMETA2bin4;			BinErrorMETA2bin4 = 1.0/TMath::Sqrt(BinErrorMETA2bin4);		}
    else{				BinContentMETA2bin4=0.0;								BinErrorMETA2bin4=0.0;					}
    if(BinErrorMETA2bin5>1.0e-10){	BinContentMETA2bin5 = BinContentMETA2bin5/BinErrorMETA2bin5;			BinErrorMETA2bin5 = 1.0/TMath::Sqrt(BinErrorMETA2bin5);		}
    else{				BinContentMETA2bin5=0.0;								BinErrorMETA2bin5=0.0;					}
    if(BinErrorMETA2bin6>1.0e-10){	BinContentMETA2bin6 = BinContentMETA2bin6/BinErrorMETA2bin6;			BinErrorMETA2bin6 = 1.0/TMath::Sqrt(BinErrorMETA2bin6);		}
    else{				BinContentMETA2bin6=0.0;								BinErrorMETA2bin6=0.0;					}
    if(BinErrorMETA2bin7>1.0e-10){	BinContentMETA2bin7 = BinContentMETA2bin7/BinErrorMETA2bin7;			BinErrorMETA2bin7 = 1.0/TMath::Sqrt(BinErrorMETA2bin7);		}
    else{				BinContentMETA2bin7=0.0;								BinErrorMETA2bin7=0.0;					}
    //normalize the bin and the error for METrigger:    
    if(BinErrorMEtriggerbin1>1.0e-10){	BinContentMETriggerbin1 = BinContentMETriggerbin1/BinErrorMEtriggerbin1;	BinErrorMEtriggerbin1 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin1);	}
    else{				BinContentMETriggerbin1=0.0;							BinErrorMEtriggerbin1=0.0;					}
    if(BinErrorMEtriggerbin2>1.0e-10){	BinContentMETriggerbin2 = BinContentMETriggerbin2/BinErrorMEtriggerbin2;	BinErrorMEtriggerbin2 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin2);	}
    else{				BinContentMETriggerbin2=0.0;							BinErrorMEtriggerbin2=0.0;					}
    if(BinErrorMEtriggerbin3>1.0e-10){	BinContentMETriggerbin3 = BinContentMETriggerbin3/BinErrorMEtriggerbin3;	BinErrorMEtriggerbin3 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin3);	}
    else{				BinContentMETriggerbin3=0.0;							BinErrorMEtriggerbin3=0.0;					}
    if(BinErrorMEtriggerbin4>1.0e-10){	BinContentMETriggerbin4 = BinContentMETriggerbin4/BinErrorMEtriggerbin4;	BinErrorMEtriggerbin4 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin4);	}
    else{				BinContentMETriggerbin4=0.0;							BinErrorMEtriggerbin4=0.0;					}
    if(BinErrorMEtriggerbin5>1.0e-10){	BinContentMETriggerbin5 = BinContentMETriggerbin5/BinErrorMEtriggerbin5;	BinErrorMEtriggerbin5 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin5);	}
    else{				BinContentMETriggerbin5=0.0;							BinErrorMEtriggerbin5=0.0;					}
    if(BinErrorMEtriggerbin6>1.0e-10){	BinContentMETriggerbin6 = BinContentMETriggerbin6/BinErrorMEtriggerbin6;	BinErrorMEtriggerbin6 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin6);	}
    else{				BinContentMETriggerbin6=0.0;							BinErrorMEtriggerbin6=0.0;					}
    if(BinErrorMEtriggerbin7>1.0e-10){	BinContentMETriggerbin7 = BinContentMETriggerbin7/BinErrorMEtriggerbin7;	BinErrorMEtriggerbin7 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin7);	}
    else{				BinContentMETriggerbin7=0.0;							BinErrorMEtriggerbin7=0.0;					}
    //Set the bin content and error in every histogram.

  if(BinContent>1.0e-10){
      histbin1->SetBinContent(x,BinContentbin1);
      histbin1->SetBinError(x,BinErrorbin1);
      if(Bin2&&isPbPb){
	histbin2->SetBinContent(x,BinContentbin2);
	histbin2->SetBinError(x,BinErrorbin2);
      }
      if(Bin3&&isPbPb){
	histbin3->SetBinContent(x,BinContentbin3);
	histbin3->SetBinError(x,BinErrorbin3);
      }
      if(Bin4&&isPbPb){
	histbin4->SetBinContent(x,BinContentbin4);
	histbin4->SetBinError(x,BinErrorbin4);
      }
      if(Bin5&&isPbPb){
	histbin5->SetBinContent(x,BinContentbin5);
	histbin5->SetBinError(x,BinErrorbin5);
      }
      if(Bin6&&isPbPb){
	histbin6->SetBinContent(x,BinContentbin6);
	histbin6->SetBinError(x,BinErrorbin6);
      }
      if(Bin7&&isPbPb){
	histbin7->SetBinContent(x,BinContentbin7);
	histbin7->SetBinError(x,BinErrorbin7);
      }
      histMETAbin1->SetBinContent(x,BinContentMETAbin1);
      histMETAbin1->SetBinError(x,BinErrorMETAbin1);
      if(Bin2&&isPbPb){
	histMETAbin2->SetBinContent(x,BinContentMETAbin2);
	histMETAbin2->SetBinError(x,BinErrorMETAbin2);
      }
      if(Bin3&&isPbPb){
	histMETAbin3->SetBinContent(x,BinContentMETAbin3);
	histMETAbin3->SetBinError(x,BinErrorMETAbin3);
      }
      if(Bin4&&isPbPb){
	histMETAbin4->SetBinContent(x,BinContentMETAbin4);
	histMETAbin4->SetBinError(x,BinErrorMETAbin4);
      }
      if(Bin5&&isPbPb){
	histMETAbin5->SetBinContent(x,BinContentMETAbin5);
	histMETAbin5->SetBinError(x,BinErrorMETAbin5);
      }
      if(Bin6&&isPbPb){
	histMETAbin6->SetBinContent(x,BinContentMETAbin6);
	histMETAbin6->SetBinError(x,BinErrorMETAbin6);
      }
      if(Bin7&&isPbPb){
	histMETAbin7->SetBinContent(x,BinContentMETAbin7);
	histMETAbin7->SetBinError(x,BinErrorMETAbin7);
      }
      histMETA2bin1->SetBinContent(x,BinContentMETA2bin1);
      histMETA2bin1->SetBinError(x,BinErrorMETA2bin1);
      if(Bin2&&isPbPb){
	histMETA2bin2->SetBinContent(x,BinContentMETA2bin2);
	histMETA2bin2->SetBinError(x,BinErrorMETA2bin2);
      }
      if(Bin3&&isPbPb){
	histMETA2bin3->SetBinContent(x,BinContentMETA2bin3);
	histMETA2bin3->SetBinError(x,BinErrorMETA2bin3);
      }
      if(Bin4&&isPbPb){
	histMETA2bin4->SetBinContent(x,BinContentMETA2bin4);
	histMETA2bin4->SetBinError(x,BinErrorMETA2bin4);
      }
      if(Bin5&&isPbPb){
	histMETA2bin5->SetBinContent(x,BinContentMETA2bin5);
	histMETA2bin5->SetBinError(x,BinErrorMETA2bin5);
      }
      if(Bin6&&isPbPb){
	histMETA2bin6->SetBinContent(x,BinContentMETA2bin6);
	histMETA2bin6->SetBinError(x,BinErrorMETA2bin6);
      }
      if(Bin7&&isPbPb){
	histMETA2bin7->SetBinContent(x,BinContentMETA2bin7);
	histMETA2bin7->SetBinError(x,BinErrorMETA2bin7);
      }
      histMETriggerbin1->SetBinContent(x,BinContentMETriggerbin1);
      histMETriggerbin1->SetBinError(x,BinErrorMEtriggerbin1);
      if(Bin2&&isPbPb){
      histMETriggerbin2->SetBinContent(x,BinContentMETriggerbin2);
      histMETriggerbin2->SetBinError(x,BinErrorMEtriggerbin2);
      }
      if(Bin3&&isPbPb){
      histMETriggerbin3->SetBinContent(x,BinContentMETriggerbin3);
      histMETriggerbin3->SetBinError(x,BinErrorMEtriggerbin3);
      }
      if(Bin4&&isPbPb){
      histMETriggerbin4->SetBinContent(x,BinContentMETriggerbin4);
      histMETriggerbin4->SetBinError(x,BinErrorMEtriggerbin4);
      }
      if(Bin5&&isPbPb){
      histMETriggerbin5->SetBinContent(x,BinContentMETriggerbin5);
      histMETriggerbin5->SetBinError(x,BinErrorMEtriggerbin5);
      }
      if(Bin6&&isPbPb){      
	histMETriggerbin6->SetBinContent(x,BinContentMETriggerbin6);
	histMETriggerbin6->SetBinError(x,BinErrorMEtriggerbin6);
      }
      if(Bin7&&isPbPb){      
	histMETriggerbin7->SetBinContent(x,BinContentMETriggerbin7);
	histMETriggerbin7->SetBinError(x,BinErrorMEtriggerbin7);
      }      
    }
    //Reset all:
    BinContent = 0;
    BinContentbin1 = 0.0;BinErrorbin1 = 0.0;BinContentMETAbin1 = 0.0;BinErrorMETAbin1 = 0.0;BinContentMETA2bin1 = 0.0;BinErrorMETA2bin1 = 0.0;BinContentMETriggerbin1 = 0.0;BinErrorMEtriggerbin1 = 0.0;
    BinContentbin2 = 0.0;BinErrorbin2 = 0.0;BinContentMETAbin2 = 0.0;BinErrorMETAbin2 = 0.0;BinContentMETA2bin2 = 0.0;BinErrorMETA2bin2 = 0.0;BinContentMETriggerbin2 = 0.0;BinErrorMEtriggerbin2 = 0.0;
    BinContentbin3 = 0.0;BinErrorbin3 = 0.0;BinContentMETAbin3 = 0.0;BinErrorMETAbin3 = 0.0;BinContentMETA2bin3 = 0.0;BinErrorMETA2bin3 = 0.0;BinContentMETriggerbin3 = 0.0;BinErrorMEtriggerbin3 = 0.0;
    BinContentbin4 = 0.0;BinErrorbin4 = 0.0;BinContentMETAbin4 = 0.0;BinErrorMETAbin4 = 0.0;BinContentMETA2bin4 = 0.0;BinErrorMETA2bin4 = 0.0;BinContentMETriggerbin4 = 0.0;BinErrorMEtriggerbin4 = 0.0;
    BinContentbin5 = 0.0;BinErrorbin5 = 0.0;BinContentMETAbin5 = 0.0;BinErrorMETAbin5 = 0.0;BinContentMETA2bin5 = 0.0;BinErrorMETA2bin5 = 0.0;BinContentMETriggerbin5 = 0.0;BinErrorMEtriggerbin5 = 0.0;	
    BinContentbin6 = 0.0;BinErrorbin6 = 0.0;BinContentMETAbin6 = 0.0;BinErrorMETAbin6 = 0.0;BinContentMETA2bin6 = 0.0;BinErrorMETA2bin6 = 0.0;BinContentMETriggerbin6 = 0.0;BinErrorMEtriggerbin6 = 0.0;	
    BinContentbin7 = 0.0;BinErrorbin7 = 0.0;BinContentMETAbin7 = 0.0;BinErrorMETAbin7 = 0.0;BinContentMETA2bin7 = 0.0;BinErrorMETA2bin7 = 0.0;BinContentMETriggerbin7 = 0.0;BinErrorMEtriggerbin7 = 0.0;	
  }//end binloop
  //save the histograms in the relevant directories:
  if(!isPbPb){
    Bin1->Samediv()->cd();
    histbin1->Write(histo->GetName());
    Bin1->METAdiv()->cd();
    histMETAbin1->Write(histo->GetName());
    if(Bin1->META2()){
      Bin1->META2div()->cd();
      histMETA2bin1->Write(histo->GetName());
    }
    Bin1->METriggerdiv()->cd();
    histMETriggerbin1->Write(histo->GetName());
  }
  else{
    Bin1->Samediv()->cd();
    histbin1->Write(histo->GetName());
    Bin1->METAdiv()->cd();
    histMETAbin1->Write(histo->GetName());
    if(Bin1->META2()){
      Bin1->META2div()->cd();
      histMETA2bin1->Write(histo->GetName());
    }
    Bin1->METriggerdiv()->cd();
    histMETriggerbin1->Write(histo->GetName());
    if(Bin2){
      Bin2->Samediv()->cd();
      histbin2->Write(histo->GetName());
      Bin2->METAdiv()->cd();
      histMETAbin2->Write(histo->GetName());
      if(Bin2->META2()){
	Bin2->META2div()->cd();
	histMETA2bin2->Write(histo->GetName());
      }
      Bin2->METriggerdiv()->cd();
      histMETriggerbin2->Write(histo->GetName());
    }
    if(Bin3){
      Bin3->Samediv()->cd();
      histbin3->Write(histo->GetName());
      Bin3->METAdiv()->cd();
      histMETAbin3->Write(histo->GetName());
      if(Bin3->META2()){
	Bin3->META2div()->cd();
	histMETA2bin3->Write(histo->GetName());
      }
      Bin3->METriggerdiv()->cd();
      histMETriggerbin3->Write(histo->GetName());
    }
    if(Bin4){
      Bin4->Samediv()->cd();
      histbin4->Write(histo->GetName());
      Bin4->METAdiv()->cd();
      histMETAbin4->Write(histo->GetName());
      if(Bin4->META2()){
	Bin4->META2div()->cd();
	histMETA2bin4->Write(histo->GetName());
      }
      Bin4->METriggerdiv()->cd();
      histMETriggerbin4->Write(histo->GetName());
    }
    if(Bin5){
      Bin5->Samediv()->cd();
      histbin5->Write(histo->GetName());
      Bin5->METAdiv()->cd();
      histMETAbin5->Write(histo->GetName());
      if(Bin5->META2()){
	Bin5->META2div()->cd();
	histMETA2bin5->Write(histo->GetName());
    }
      Bin5->METriggerdiv()->cd();
      histMETriggerbin5->Write(histo->GetName());
    }
    if(Bin6){
      Bin6->Samediv()->cd();
      histbin6->Write(histo->GetName());
      Bin6->METAdiv()->cd();
      histMETAbin6->Write(histo->GetName());
      if(Bin6->META2()){
	Bin6->META2div()->cd();
	histMETA2bin6->Write(histo->GetName());
      }
      Bin6->METriggerdiv()->cd();
      histMETriggerbin6->Write(histo->GetName());
    }
    if(Bin7){
      Bin7->Samediv()->cd();
      histbin7->Write(histo->GetName());
      Bin7->METAdiv()->cd();
      histMETAbin7->Write(histo->GetName());
      if(Bin7->META2()){
	Bin7->META2div()->cd();
	histMETA2bin7->Write(histo->GetName());
      }
      Bin7->METriggerdiv()->cd();
      histMETriggerbin7->Write(histo->GetName());
    }
  }
  delete histbin1; delete histMETAbin1; delete histMETA2bin1; delete histMETriggerbin1;   
  delete histbin2; delete histMETAbin2; delete histMETA2bin2; delete histMETriggerbin2;   
  delete histbin3; delete histMETAbin3; delete histMETA2bin3; delete histMETriggerbin3;   
  delete histbin4; delete histMETAbin4; delete histMETA2bin4; delete histMETriggerbin4;   
  delete histbin5; delete histMETAbin5; delete histMETA2bin5; delete histMETriggerbin5;   
  delete histbin6; delete histMETAbin6; delete histMETA2bin6; delete histMETriggerbin6;   
  delete histbin7; delete histMETAbin7; delete histMETA2bin7; delete histMETriggerbin7;   

}
void CollectHist(TH2D* histo, TList * directories, TObjArray* multdirlist,Bool_t pearsonserrors ,Bool_t collectdivfirst , Bool_t isPbPb,Bool_t divbybinwidth){
  BinDirs * Bin1=NULL;BinDirs * Bin2=NULL;BinDirs * Bin3=NULL;BinDirs * Bin4=NULL;BinDirs * Bin5=NULL;BinDirs * Bin6=NULL;BinDirs * Bin7=NULL;
  if(!isPbPb){
   Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(0));
  }
  if(isPbPb){
    Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(0));
    if(multdirlist->GetEntries()>1) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(1));
    if(multdirlist->GetEntries()>2) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(2));
    if(multdirlist->GetEntries()>3) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(3));
    if(multdirlist->GetEntries()>4) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(4));
    if(multdirlist->GetEntries()>5) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(5));
    if(multdirlist->GetEntries()>6) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  }
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  double xaxisbinwidth = histo->GetXaxis()->GetBinCenter(5)-histo->GetXaxis()->GetBinCenter(4);
  double yaxisbinwidth = histo->GetYaxis()->GetBinCenter(5)-histo->GetYaxis()->GetBinCenter(4);
  
  TString namezaxis = TString("");
  if(divbybinwidth&&TString(histo->GetName()).Contains("DPhi_1_DPHI"))		namezaxis.Append("#frac{dN^{12}_{pairs}}{N_{trig}d#phi_{1} d#phi_{2} }(rad)^{-2}");
  if(divbybinwidth&&TString(histo->GetName()).Contains("DPhi_1_DEta"))		namezaxis.Append("#frac{dN^{12}_{pairs}}{N_{trig}d#phi d#eta }(rad)^{-1}");
  if(divbybinwidth&&TString(histo->GetName()).Contains("DPhi_1_DEta_12")){
    namezaxis.Clear();
    namezaxis.Append("#frac{dN^{12}_{pairs}}{N_{trig}d#phi_{1} d#eta_{12} }(rad)^{-1}");
  }
  if(!divbybinwidth) 								namezaxis.Append("#frac{N^{12}_{pairs}}{N_{trig}}");
  //Bin 1:
  TH2D* histbin1 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  histbin1->GetZaxis()->SetTitle(namezaxis.Data());
  TH2D* histMETAbin1 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  histMETAbin1->GetZaxis()->SetTitle(namezaxis.Data());
  TH2D* histMETA2bin1 		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1META2"	,histo->GetName())));
  histMETA2bin1->GetZaxis()->SetTitle(namezaxis.Data());
  TH2D* histMETriggerbin1 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  histMETriggerbin1->GetZaxis()->SetTitle(namezaxis.Data());

  //Bin 2:
  TH2D* histbin2 		= new TH2D("emptybin2","title",1,0,1,1,0,1);
  TH2D* histMETAbin2 		= new TH2D("emptybinMETA2","title",1,0,1,1,0,1);
  TH2D* histMETA2bin2 		= new TH2D("emptybinMETA22","title",1,0,1,1,0,1);
  TH2D* histMETriggerbin2 	= new TH2D("emptybinMETRIGGER2","title",1,0,1,1,0,1);
  if(Bin2&&isPbPb){
    histbin2->Delete();
    histMETAbin2->Delete();
    histMETA2bin2->Delete();
    histMETriggerbin2->Delete();
    histbin2		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histbin2->GetZaxis()->SetTitle(namezaxis.Data());
    histMETAbin2	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETAbin2->GetZaxis()->SetTitle(namezaxis.Data());
    histMETA2bin2	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2META2"	,histo->GetName())));
    histMETA2bin2->GetZaxis()->SetTitle(namezaxis.Data());
    histMETriggerbin2 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));    
    histMETriggerbin2->GetZaxis()->SetTitle(namezaxis.Data());
  }
  //Bin 3:
  TH2D* histbin3 		= new TH2D("emptybin3","title",1,0,1,1,0,1);
  TH2D* histMETAbin3 		= new TH2D("emptybinMETA3","title",1,0,1,1,0,1);
  TH2D* histMETA2bin3 		= new TH2D("emptybinMETA33","title",1,0,1,1,0,1);
  TH2D* histMETriggerbin3 	= new TH2D("emptybinMETRIGGER3","title",1,0,1,1,0,1);
  if(Bin3&&isPbPb){
    histbin3->Delete();
    histMETAbin3->Delete();
    histMETA2bin3->Delete();
    histMETriggerbin3->Delete();
    histbin3		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histbin3->GetZaxis()->SetTitle(namezaxis.Data());
    histMETAbin3	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETAbin3->GetZaxis()->SetTitle(namezaxis.Data());
    histMETA2bin3	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3META2"	,histo->GetName())));
    histMETA2bin3->GetZaxis()->SetTitle(namezaxis.Data());
    histMETriggerbin3 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
    histMETriggerbin3->GetZaxis()->SetTitle(namezaxis.Data());
  }
  //Bin 4:
  TH2D* histbin4 		= new TH2D("emptybin4","title",1,0,1,1,0,1);
  TH2D* histMETAbin4 		= new TH2D("emptybinMETA4","title",1,0,1,1,0,1);
  TH2D* histMETA2bin4 		= new TH2D("emptybinMETA44","title",1,0,1,1,0,1);
  TH2D* histMETriggerbin4 	= new TH2D("emptybinMETRIGGER4","title",1,0,1,1,0,1);
  if(Bin4&&isPbPb){
    histbin4->Delete();
    histMETAbin4->Delete();
    histMETA2bin4->Delete();
    histMETriggerbin4->Delete();
    histbin4		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histbin4->GetZaxis()->SetTitle(namezaxis.Data());
    histMETAbin4	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETAbin4->GetZaxis()->SetTitle(namezaxis.Data());
    histMETA2bin4	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4META2"	,histo->GetName())));
    histMETA2bin4->GetZaxis()->SetTitle(namezaxis.Data());
    histMETriggerbin4 	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
    histMETriggerbin4->GetZaxis()->SetTitle(namezaxis.Data());
  }
  //Bin 5:
  TH2D* histbin5 		= new TH2D("emptybin5","title",1,0,1,1,0,1);
  TH2D* histMETAbin5 		= new TH2D("emptybinMETA5","title",1,0,1,1,0,1);
  TH2D* histMETA2bin5 		= new TH2D("emptybinMETA55","title",1,0,1,1,0,1);
  TH2D* histMETriggerbin5 	= new TH2D("emptybinMETRIGGER5","title",1,0,1,1,0,1);
  if(Bin5&&isPbPb){
    histbin5->Delete();
    histMETAbin5->Delete();
    histMETA2bin5->Delete();
    histMETriggerbin5->Delete();
    histbin5		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histbin5->GetZaxis()->SetTitle(namezaxis.Data());
    histMETAbin5	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETAbin5->GetZaxis()->SetTitle(namezaxis.Data());
    histMETA2bin5	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5META2"	,histo->GetName())));
    histMETA2bin5->GetZaxis()->SetTitle(namezaxis.Data());
    histMETriggerbin5	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
    histMETriggerbin5->GetZaxis()->SetTitle(namezaxis.Data());
  }
  //Bin 6:
  TH2D* histbin6 		= new TH2D("emptybin6","title",1,0,1,1,0,1);
  TH2D* histMETAbin6 		= new TH2D("emptybinMETA6","title",1,0,1,1,0,1);
  TH2D* histMETA2bin6 		= new TH2D("emptybinMETA66","title",1,0,1,1,0,1);
  TH2D* histMETriggerbin6 	= new TH2D("emptybinMETRIGGER6","title",1,0,1,1,0,1);
  if(Bin6&&isPbPb){
    histbin6->Delete();
    histMETAbin6->Delete();
    histMETA2bin6->Delete();
    histMETriggerbin6->Delete();
    histbin6		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin6"		,histo->GetName())));
    histbin6->GetZaxis()->SetTitle(namezaxis.Data());
    histMETAbin6	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin6META"	,histo->GetName())));
    histMETAbin6->GetZaxis()->SetTitle(namezaxis.Data());
    histMETA2bin6	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin6META2"	,histo->GetName())));
    histMETA2bin6->GetZaxis()->SetTitle(namezaxis.Data());
    histMETriggerbin6	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin6METrigger"	,histo->GetName())));
    histMETriggerbin6->GetZaxis()->SetTitle(namezaxis.Data());
  }
  //Bin 7:
  TH2D* histbin7 		= new TH2D("emptybin7","title",1,0,1,1,0,1);
  TH2D* histMETAbin7 		= new TH2D("emptybinMETA7","title",1,0,1,1,0,1);
  TH2D* histMETA2bin7 		= new TH2D("emptybinMETA77","title",1,0,1,1,0,1);
  TH2D* histMETriggerbin7 	= new TH2D("emptybinMETRIGGER7","title",1,0,1,1,0,1);
  if(Bin7&&isPbPb){
    histbin7->Delete();
    histMETAbin7->Delete();
    histMETA2bin7->Delete();
    histMETriggerbin7->Delete();
    histbin7		= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin7"		,histo->GetName())));
    histbin7->GetZaxis()->SetTitle(namezaxis.Data());
    histMETAbin7	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin7META"	,histo->GetName())));
    histMETAbin7->GetZaxis()->SetTitle(namezaxis.Data());
    histMETA2bin7	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin7META2"	,histo->GetName())));
    histMETA2bin7->GetZaxis()->SetTitle(namezaxis.Data());
    histMETriggerbin7	= dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin7METrigger"	,histo->GetName())));
    histMETriggerbin7->GetZaxis()->SetTitle(namezaxis.Data());
  }  
  //For mixed event:
  TH2D* histbin1m 		= new TH2D("emptybin1m","title",1,0,1,1,0,1);
  TH2D* histbin2m 		= new TH2D("emptybin2m","title",1,0,1,1,0,1);
  TH2D* histbin3m 		= new TH2D("emptybin3m","title",1,0,1,1,0,1);
  TH2D* histbin4m 		= new TH2D("emptybin4m","title",1,0,1,1,0,1);
  TH2D* histbin5m 		= new TH2D("emptybin5m","title",1,0,1,1,0,1);
  TH2D* histbin6m 		= new TH2D("emptybin6m","title",1,0,1,1,0,1);
  TH2D* histbin7m 		= new TH2D("emptybin7m","title",1,0,1,1,0,1);
  
  if(!collectdivfirst){
    cout << "here"<<endl;
    histbin1m->Delete();
    histbin1m = dynamic_cast<TH2D*>( histo->Clone(Form("%sdivbin1m",histo->GetName())));
    histbin1m->Reset();
    if(Bin2&&isPbPb){histbin2m->Delete();histbin2m = dynamic_cast<TH2D*>(histbin1m->Clone(Form("%sdivbin2m",histo->GetName())));}
    if(Bin3&&isPbPb){histbin2m->Delete();histbin3m = dynamic_cast<TH2D*>(histbin1m->Clone(Form("%sdivbin3m",histo->GetName())));}
    if(Bin4&&isPbPb){histbin4m->Delete();histbin4m = dynamic_cast<TH2D*>(histbin1m->Clone(Form("%sdivbin4m",histo->GetName())));}
    if(Bin5&&isPbPb){histbin5m->Delete();histbin5m = dynamic_cast<TH2D*>(histbin1m->Clone(Form("%sdivbin5m",histo->GetName())));}
    if(Bin6&&isPbPb){histbin6m->Delete();histbin6m = dynamic_cast<TH2D*>(histbin1m->Clone(Form("%sdivbin6m",histo->GetName())));}
    if(Bin7&&isPbPb){histbin7m->Delete();histbin7m = dynamic_cast<TH2D*>(histbin1m->Clone(Form("%sdivbin7m",histo->GetName())));} 
  }
  //Doubles to hold the values until they are put into the hists.
  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;Double_t bincontlMETA2   = 0.0;Double_t binerrorlMETA2  = 0.0;
  Double_t BinContent = 0.0;Double_t BinError   = 0.0; Double_t BinContentMETA = 0.0;Double_t BinErrorMETA   = 0.0; Double_t BinContentMETA2 = 0.0;Double_t BinErrorMETA2   = 0.0;
  Double_t BinContentMETrigger = 0.0;Double_t BinErrorMEtrigger   = 0.0; Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETA2bin1 = 0.0;Double_t BinErrorMETA2bin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETA2bin2 = 0.0;Double_t BinErrorMETA2bin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETA2bin3 = 0.0;Double_t BinErrorMETA2bin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETA2bin4 = 0.0;Double_t BinErrorMETA2bin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETA2bin5 = 0.0;Double_t BinErrorMETA2bin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  Double_t BinContentbin6 = 0.0;Double_t BinErrorbin6   = 0.0;Double_t BinContentMETAbin6 = 0.0;Double_t BinErrorMETAbin6   = 0.0;Double_t BinContentMETA2bin6 = 0.0;Double_t BinErrorMETA2bin6   = 0.0;Double_t BinContentMETriggerbin6 = 0.0;Double_t BinErrorMEtriggerbin6   = 0.0; 
  Double_t BinContentbin7 = 0.0;Double_t BinErrorbin7   = 0.0;Double_t BinContentMETAbin7 = 0.0;Double_t BinErrorMETAbin7   = 0.0;Double_t BinContentMETA2bin7 = 0.0;Double_t BinErrorMETA2bin7   = 0.0;Double_t BinContentMETriggerbin7 = 0.0;Double_t BinErrorMEtriggerbin7   = 0.0; 
  if(collectdivfirst){
    //in this case, first add all same_event/mixed_event, then scale by #triggers.
    Double_t scalingfactorbin1=0.0;Double_t scalingfactorbin2=0.0;Double_t scalingfactorbin3=0.0;Double_t scalingfactorbin4=0.0;Double_t scalingfactorbin5=0.0;Double_t scalingfactorbin6=0.0;Double_t scalingfactorbin7=0.0;
    Double_t scalingfactorMETAbin1=0.0;Double_t scalingfactorMETAbin2=0.0;Double_t scalingfactorMETAbin3=0.0;Double_t scalingfactorMETAbin4=0.0;Double_t scalingfactorMETAbin5=0.0;Double_t scalingfactorMETAbin6=0.0;Double_t scalingfactorMETAbin7=0.0;
    Double_t scalingfactorMETA2bin1=0.0;Double_t scalingfactorMETA2bin2=0.0;Double_t scalingfactorMETA2bin3=0.0;Double_t scalingfactorMETA2bin4=0.0;Double_t scalingfactorMETA2bin5=0.0;Double_t scalingfactorMETA2bin6=0.0;Double_t scalingfactorMETA2bin7=0.0;
    Double_t scalingfactorMETriggerbin1=0.0;Double_t scalingfactorMETriggerbin2=0.0;Double_t scalingfactorMETriggerbin3=0.0;Double_t scalingfactorMETriggerbin4=0.0;Double_t scalingfactorMETriggerbin5=0.0;Double_t scalingfactorMETriggerbin6=0.0;Double_t scalingfactorMETriggerbin7=0.0;
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
      BinDirs * bin = dynamic_cast<BinDirs*>(directories->At(i));
      
      //find the Multiplicity bin we are in
      for(int j = 1;j<multdirlist->GetEntries();j++){
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo("BinM(0.00)->(10.00)")){
	  if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(0.00)->(5.00)")==0) Mbin = j;
	  if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(5.00)->(10.00)")==0)Mbin = j;
	}
      }
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH2D*>(bin->SameDir("same_event")->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(bin->SameDir("mixed_event")->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.
      if(!dynamic_cast<TH2D*>(bin->METADir("same_event")->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(bin->METADir("mixed_event")->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.
      if(!dynamic_cast<TH2D*>(bin->METriggerDir("same_event")->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(bin->METriggerDir("mixed_event")->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.      
      TH2D* histnow = dynamic_cast<TH2D*>(bin->SameDir("same_event")->Get(histo->GetName())->Clone("nowhist"));
      TH2D* histmixed = dynamic_cast<TH2D*>(bin->SameDir("mixed_event")->Get(histo->GetName()));
      Double_t scale = dynamic_cast<TParameter<double>*>(bin->SameDir("mixed_event")->Get(Form("%s_scale",histo->GetName())))->GetVal();
      histmixed->Scale(scale);
      scale = histmixed->GetMaximum();
      if(scale>1.0E-10)histmixed->Scale(1.0/scale);
      else histmixed->Scale(0.0);
      histnow->Divide(histmixed);
      if(Mbin == 0&&isPbPb){histbin1->Add(histnow);scalingfactorbin1+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(!isPbPb){histbin1->Add(histnow);scalingfactorbin1+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 1&&isPbPb){histbin2->Add(histnow);scalingfactorbin2+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 2&&isPbPb){histbin3->Add(histnow);scalingfactorbin3+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 3&&isPbPb){histbin4->Add(histnow);scalingfactorbin4+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 4&&isPbPb){histbin5->Add(histnow);scalingfactorbin5+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 5&&isPbPb){histbin6->Add(histnow);scalingfactorbin6+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 6&&isPbPb){histbin7->Add(histnow);scalingfactorbin7+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      delete histnow;
      TH2D* histnowMETA = dynamic_cast<TH2D*>(bin->METADir("same_event")->Get(histo->GetName())->Clone("nowhistMETA"));
      histnowMETA->Divide(histmixed);
      if(Mbin == 0){histMETAbin1->Add(histnowMETA);scalingfactorMETAbin1+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(!isPbPb){histMETAbin1->Add(histnowMETA);scalingfactorMETAbin1+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 1&&isPbPb){histMETAbin2->Add(histnowMETA);scalingfactorMETAbin2+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 2&&isPbPb){histMETAbin3->Add(histnowMETA);scalingfactorMETAbin3+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 3&&isPbPb){histMETAbin4->Add(histnowMETA);scalingfactorMETAbin4+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 4&&isPbPb){histMETAbin5->Add(histnowMETA);scalingfactorMETAbin5+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 5&&isPbPb){histMETAbin6->Add(histnowMETA);scalingfactorMETAbin6+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 6&&isPbPb){histMETAbin7->Add(histnowMETA);scalingfactorMETAbin7+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();} 
      delete histnowMETA;
      if(bin->META2()){
	TH2D* histnowMETA2 = dynamic_cast<TH2D*>(bin->META2Dir("same_event")->Get(histo->GetName())->Clone("nowhistMETA2"));
	histnowMETA2->Divide(histmixed);
	if(Mbin == 0){histMETA2bin1->Add(histnowMETA2);scalingfactorMETA2bin1+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(!isPbPb){histMETA2bin1->Add(histnowMETA2);scalingfactorMETA2bin1+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 1&&isPbPb){histMETA2bin2->Add(histnowMETA2);scalingfactorMETA2bin2+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 2&&isPbPb){histMETA2bin3->Add(histnowMETA2);scalingfactorMETA2bin3+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 3&&isPbPb){histMETA2bin4->Add(histnowMETA2);scalingfactorMETA2bin4+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 4&&isPbPb){histMETA2bin5->Add(histnowMETA2);scalingfactorMETA2bin5+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 5&&isPbPb){histMETA2bin6->Add(histnowMETA2);scalingfactorMETA2bin6+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 6&&isPbPb){histMETA2bin7->Add(histnowMETA2);scalingfactorMETA2bin7+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	delete histnowMETA2;
      }
      TH2D* histnowMETrigger = dynamic_cast<TH2D*>(bin->METriggerDir("same_event")->Get(histo->GetName())->Clone("nowhistMETrigger"));
      histnowMETrigger->Divide(histmixed);
      if(Mbin == 0){histMETriggerbin1->Add(histnowMETrigger);scalingfactorMETriggerbin1+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(!isPbPb){histMETriggerbin1->Add(histnowMETrigger);scalingfactorMETriggerbin1+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 1&&isPbPb){histMETriggerbin2->Add(histnowMETrigger);scalingfactorMETriggerbin2+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 2&&isPbPb){histMETriggerbin3->Add(histnowMETrigger);scalingfactorMETriggerbin3+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 3&&isPbPb){histMETriggerbin4->Add(histnowMETrigger);scalingfactorMETriggerbin4+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 4&&isPbPb){histMETriggerbin5->Add(histnowMETrigger);scalingfactorMETriggerbin5+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 5&&isPbPb){histMETriggerbin6->Add(histnowMETrigger);scalingfactorMETriggerbin6+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 6&&isPbPb){histMETriggerbin7->Add(histnowMETrigger);scalingfactorMETriggerbin7+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      delete histnowMETrigger;
    }

    if(!scalingfactorbin1 <1.0e-10)histbin1->Scale(1.0/scalingfactorbin1);
    if(divbybinwidth&&histbin1)histbin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin2 <1.0e-10)histbin2->Scale(1.0/scalingfactorbin2);
    if(divbybinwidth&&histbin1)histbin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin3 <1.0e-10)histbin3->Scale(1.0/scalingfactorbin3);
    if(divbybinwidth&&histbin3)histbin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin4 <1.0e-10)histbin4->Scale(1.0/scalingfactorbin4);
    if(divbybinwidth&&histbin4)histbin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin5 <1.0e-10)histbin5->Scale(1.0/scalingfactorbin5);
    if(divbybinwidth&&histbin5)histbin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin6 <1.0e-10)histbin6->Scale(1.0/scalingfactorbin6);
    if(divbybinwidth&&histbin6)histbin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin7 <1.0e-10)histbin7->Scale(1.0/scalingfactorbin7);
    if(divbybinwidth&&histbin7)histbin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));

    if(!scalingfactorMETAbin1 <1.0e-10)histMETAbin1->Scale(1.0/scalingfactorMETAbin1);
    if(divbybinwidth&&histbin1)histMETAbin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin2 <1.0e-10)histMETAbin2->Scale(1.0/scalingfactorMETAbin2);
    if(divbybinwidth&&histbin2)histMETAbin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin3 <1.0e-10)histMETAbin3->Scale(1.0/scalingfactorMETAbin3);
    if(divbybinwidth&&histbin3)histMETAbin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin4 <1.0e-10)histMETAbin4->Scale(1.0/scalingfactorMETAbin4);
    if(divbybinwidth&&histbin4)histMETAbin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin5 <1.0e-10)histMETAbin5->Scale(1.0/scalingfactorMETAbin5);
    if(divbybinwidth&&histbin5)histMETAbin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin6 <1.0e-10)histMETAbin6->Scale(1.0/scalingfactorMETAbin6);
    if(divbybinwidth&&histbin6)histMETAbin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin7 <1.0e-10)histMETAbin7->Scale(1.0/scalingfactorMETAbin7);
    if(divbybinwidth&&histbin7)histMETAbin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    
    if(!scalingfactorMETA2bin1 <1.0e-10)histMETA2bin1->Scale(1.0/scalingfactorMETA2bin1);
    if(divbybinwidth&&histbin1)histMETA2bin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin2 <1.0e-10)histMETA2bin2->Scale(1.0/scalingfactorMETA2bin2);
    if(divbybinwidth&&histbin2)histMETA2bin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin3 <1.0e-10)histMETA2bin3->Scale(1.0/scalingfactorMETA2bin3);
    if(divbybinwidth&&histbin3)histMETA2bin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin4 <1.0e-10)histMETA2bin4->Scale(1.0/scalingfactorMETA2bin4);
    if(divbybinwidth&&histbin4)histMETA2bin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin5 <1.0e-10)histMETA2bin5->Scale(1.0/scalingfactorMETA2bin5);
    if(divbybinwidth&&histbin5)histMETA2bin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin6 <1.0e-10)histMETA2bin6->Scale(1.0/scalingfactorMETA2bin6);
    if(divbybinwidth&&histbin6)histMETA2bin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin7 <1.0e-10)histMETA2bin7->Scale(1.0/scalingfactorMETA2bin7);
    if(divbybinwidth&&histbin7)histMETA2bin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    
    if(!scalingfactorMETriggerbin1 <1.0e-10)histMETriggerbin1->Scale(1.0/scalingfactorMETriggerbin1);
    if(divbybinwidth&&histbin1)histMETriggerbin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin2 <1.0e-10)histMETriggerbin2->Scale(1.0/scalingfactorMETriggerbin2);
    if(divbybinwidth&&histbin2)histMETriggerbin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin3 <1.0e-10)histMETriggerbin3->Scale(1.0/scalingfactorMETriggerbin3);
    if(divbybinwidth&&histbin3)histMETriggerbin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin4 <1.0e-10)histMETriggerbin4->Scale(1.0/scalingfactorMETriggerbin4);
    if(divbybinwidth&&histbin4)histMETriggerbin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin5 <1.0e-10)histMETriggerbin5->Scale(1.0/scalingfactorMETriggerbin5);
    if(divbybinwidth&&histbin5)histMETriggerbin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin6 <1.0e-10)histMETriggerbin6->Scale(1.0/scalingfactorMETriggerbin6);
    if(divbybinwidth&&histbin6)histMETriggerbin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin7 <1.0e-10)histMETriggerbin7->Scale(1.0/scalingfactorMETriggerbin7);
    if(divbybinwidth&&histbin7)histMETriggerbin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
  }
 if(!collectdivfirst){
   //sum same & mixed each for itself, then divide, then scale by triggers
    Double_t scalingfactorbin1=0.0;Double_t scalingfactorbin2=0.0;Double_t scalingfactorbin3=0.0;Double_t scalingfactorbin4=0.0;Double_t scalingfactorbin5=0.0;Double_t scalingfactorbin6=0.0;Double_t scalingfactorbin7=0.0;
    Double_t scalingfactormbin1=0.0;Double_t scalingfactormbin2=0.0;Double_t scalingfactormbin3=0.0;Double_t scalingfactormbin4=0.0;Double_t scalingfactormbin5=0.0;Double_t scalingfactormbin6=0.0;Double_t scalingfactormbin7=0.0;
    Double_t scalingfactorMETAbin1=0.0;Double_t scalingfactorMETAbin2=0.0;Double_t scalingfactorMETAbin3=0.0;Double_t scalingfactorMETAbin4=0.0;Double_t scalingfactorMETAbin5=0.0;Double_t scalingfactorMETAbin6=0.0;Double_t scalingfactorMETAbin7=0.0;
    Double_t scalingfactorMETA2bin1=0.0;Double_t scalingfactorMETA2bin2=0.0;Double_t scalingfactorMETA2bin3=0.0;Double_t scalingfactorMETA2bin4=0.0;Double_t scalingfactorMETA2bin5=0.0;Double_t scalingfactorMETA2bin6=0.0;Double_t scalingfactorMETA2bin7=0.0;
    Double_t scalingfactorMETriggerbin1=0.0;Double_t scalingfactorMETriggerbin2=0.0;Double_t scalingfactorMETriggerbin3=0.0;Double_t scalingfactorMETriggerbin4=0.0;Double_t scalingfactorMETriggerbin5=0.0;Double_t scalingfactorMETriggerbin6=0.0;Double_t scalingfactorMETriggerbin7=0.0;
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
      BinDirs * bin = dynamic_cast<BinDirs*>(directories->At(i));
      //find the Multiplicity bin we are in
      for(int j = 1;j<multdirlist->GetEntries();j++){
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo("BinM(0.00)->(10.00)")){
	  if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(0.00)->(5.00)")==0) Mbin = j;
	  if(TString(TString(bin->Same()->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(5.00)->(10.00)")==0)Mbin = j;
	}
      }
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH2D*>(bin->SameDir("same_event")->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(bin->SameDir("mixed_event")->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.
      if(!dynamic_cast<TH2D*>(bin->METADir("same_event")->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(bin->METADir("mixed_event")->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.
      if(!dynamic_cast<TH2D*>(bin->METriggerDir("same_event")->Get(histo->GetName()))){continue;}//check if the hist is in same_event
      if(!dynamic_cast<TH2D*>(bin->METriggerDir("mixed_event")->Get(histo->GetName()))){continue;}//check if the hist is in mixed_event
      if(!dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))){continue;}//check if the ntrigger hist is accessible.      
      TH2D* histnow = dynamic_cast<TH2D*>(bin->SameDir("same_event")->Get(histo->GetName())->Clone("nowhist"));

      if(Mbin == 0&&isPbPb){histbin1->Add(histnow);scalingfactorbin1+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(!isPbPb){histbin1->Add(histnow);scalingfactorbin1+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 1&&isPbPb){histbin2->Add(histnow);scalingfactorbin2+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 2&&isPbPb){histbin3->Add(histnow);scalingfactorbin3+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 3&&isPbPb){histbin4->Add(histnow);scalingfactorbin4+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 4&&isPbPb){histbin5->Add(histnow);scalingfactorbin5+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 5&&isPbPb){histbin6->Add(histnow);scalingfactorbin6+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 6&&isPbPb){histbin7->Add(histnow);scalingfactorbin7+= dynamic_cast<TH1D*>(bin->SameDir("same_event")->Get("number_of_triggers"))->Integral();}
      delete histnow;
      TH2D* histmixed = dynamic_cast<TH2D*>(bin->SameDir("mixed_event")->Get(histo->GetName()));
      Double_t scale = dynamic_cast<TParameter<double>*>(bin->SameDir("mixed_event")->Get(Form("%s_scale",histo->GetName())))->GetVal();
      histmixed->Scale(scale);

      if(Mbin == 0&&isPbPb){histbin1m->Add(histmixed);scalingfactormbin1+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      if(!isPbPb){histbin1m->Add(histmixed);scalingfactormbin1+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 1&&isPbPb){histbin2m->Add(histmixed);scalingfactormbin2+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 2&&isPbPb){histbin3m->Add(histmixed);scalingfactormbin3+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 3&&isPbPb){histbin4m->Add(histmixed);scalingfactormbin4+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 4&&isPbPb){histbin5m->Add(histmixed);scalingfactormbin5+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 5&&isPbPb){histbin6m->Add(histmixed);scalingfactormbin6+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 6&&isPbPb){histbin7m->Add(histmixed);scalingfactormbin7+= dynamic_cast<TH1D*>(bin->SameDir("mixed_event")->Get("number_of_triggers"))->Integral();}
      delete histmixed;
      TH2D* histnowMETA = dynamic_cast<TH2D*>(bin->METADir("same_event")->Get(histo->GetName())->Clone("nowhistMETA"));
      if(Mbin == 0){histMETAbin1->Add(histnowMETA);scalingfactorMETAbin1+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(!isPbPb){histMETAbin1->Add(histnowMETA);scalingfactorMETAbin1+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 1&&isPbPb){histMETAbin2->Add(histnowMETA);scalingfactorMETAbin2+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 2&&isPbPb){histMETAbin3->Add(histnowMETA);scalingfactorMETAbin3+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 3&&isPbPb){histMETAbin4->Add(histnowMETA);scalingfactorMETAbin4+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 4&&isPbPb){histMETAbin5->Add(histnowMETA);scalingfactorMETAbin5+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 5&&isPbPb){histMETAbin6->Add(histnowMETA);scalingfactorMETAbin6+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 6&&isPbPb){histMETAbin7->Add(histnowMETA);scalingfactorMETAbin7+= dynamic_cast<TH1D*>(bin->METADir("same_event")->Get("number_of_triggers"))->Integral();} 
      delete histnowMETA;
      if(bin->META2()){
	TH2D* histnowMETA2 = dynamic_cast<TH2D*>(bin->META2Dir("same_event")->Get(histo->GetName())->Clone("nowhistMETA2"));
	if(Mbin == 0){histMETA2bin1->Add(histnowMETA2);scalingfactorMETA2bin1+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(!isPbPb){histMETA2bin1->Add(histnowMETA2);scalingfactorMETA2bin1+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 1&&isPbPb){histMETA2bin2->Add(histnowMETA2);scalingfactorMETA2bin2+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 2&&isPbPb){histMETA2bin3->Add(histnowMETA2);scalingfactorMETA2bin3+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 3&&isPbPb){histMETA2bin4->Add(histnowMETA2);scalingfactorMETA2bin4+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 4&&isPbPb){histMETA2bin5->Add(histnowMETA2);scalingfactorMETA2bin5+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 5&&isPbPb){histMETA2bin6->Add(histnowMETA2);scalingfactorMETA2bin6+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	if(Mbin == 6&&isPbPb){histMETA2bin7->Add(histnowMETA2);scalingfactorMETA2bin7+= dynamic_cast<TH1D*>(bin->META2Dir("same_event")->Get("number_of_triggers"))->Integral();}
	delete histnowMETA2;
      }
      TH2D* histnowMETrigger = dynamic_cast<TH2D*>(bin->METriggerDir("same_event")->Get(histo->GetName())->Clone("nowhistMETrigger"));
      if(Mbin == 0){histMETriggerbin1->Add(histnowMETrigger);scalingfactorMETriggerbin1+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(!isPbPb){histMETriggerbin1->Add(histnowMETrigger);scalingfactorMETriggerbin1+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 1&&isPbPb){histMETriggerbin2->Add(histnowMETrigger);scalingfactorMETriggerbin2+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 2&&isPbPb){histMETriggerbin3->Add(histnowMETrigger);scalingfactorMETriggerbin3+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 3&&isPbPb){histMETriggerbin4->Add(histnowMETrigger);scalingfactorMETriggerbin4+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 4&&isPbPb){histMETriggerbin5->Add(histnowMETrigger);scalingfactorMETriggerbin5+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 5&&isPbPb){histMETriggerbin6->Add(histnowMETrigger);scalingfactorMETriggerbin6+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      if(Mbin == 6&&isPbPb){histMETriggerbin7->Add(histnowMETrigger);scalingfactorMETriggerbin7+= dynamic_cast<TH1D*>(bin->METriggerDir("same_event")->Get("number_of_triggers"))->Integral();}
      delete histnowMETrigger;

      
    }

    Double_t max = histbin1m->GetMaximum();
    histbin1->Divide(histbin1,histbin1m,max,1.0);
    histMETAbin1->Divide(histMETAbin1,histbin1m,max,1.0);
    histMETA2bin1->Divide(histMETA2bin1,histbin1m,max,1.0);
    histMETriggerbin1->Divide(histMETriggerbin1,histbin1m,max,1.0);

    if(Bin2){
      max = histbin2m->GetMaximum();
      histbin2->Divide(histbin2,histbin2m,max,1.0);
      histMETAbin2->Divide(histMETAbin2,histbin2m,max,1.0);
      histMETA2bin2->Divide(histMETA2bin2,histbin2m,max,1.0);
      histMETriggerbin2->Divide(histMETriggerbin2,histbin2m,max,1.0);
    }
    if(Bin3){
      max = histbin3m->GetMaximum();
      histbin3->Divide(histbin3,histbin3m,max,1.0);
      histMETAbin3->Divide(histMETAbin3,histbin3m,max,1.0);
      histMETA2bin3->Divide(histMETA2bin3,histbin3m,max,1.0);
      histMETriggerbin3->Divide(histMETriggerbin3,histbin3m,max,1.0);
    }    
    if(Bin4){
      max = histbin4m->GetMaximum();
      histbin4->Divide(histbin4,histbin4m,max,1.0);
      histMETAbin4->Divide(histMETAbin4,histbin4m,max,1.0);
      histMETA2bin4->Divide(histMETA2bin4,histbin4m,max,1.0);
      histMETriggerbin4->Divide(histMETriggerbin4,histbin4m,max,1.0);
    }
    if(Bin5){
      max = histbin5m->GetMaximum();
      histbin5->Divide(histbin5,histbin5m,max,1.0);
      histMETAbin5->Divide(histMETAbin5,histbin5m,max,1.0);
      histMETA2bin5->Divide(histMETA2bin5,histbin5m,max,1.0);
      histMETriggerbin5->Divide(histMETriggerbin5,histbin5m,max,1.0);
    }
    if(Bin6){
      max = histbin6m->GetMaximum();
      histbin6->Divide(histbin6,histbin6m,max,1.0);
      histMETAbin6->Divide(histMETAbin6,histbin6m,max,1.0);
      histMETA2bin6->Divide(histMETA2bin6,histbin6m,max,1.0);
      histMETriggerbin6->Divide(histMETriggerbin6,histbin6m,max,1.0);
    }    
    if(Bin7){
      max = histbin7m->GetMaximum();
      histbin7->Divide(histbin7,histbin7m,max,1.0);
      histMETAbin7->Divide(histMETAbin7,histbin7m,max,1.0);
      histMETA2bin7->Divide(histMETA2bin7,histbin7m,max,1.0);
      histMETriggerbin7->Divide(histMETriggerbin7,histbin7m,max,1.0);
    }    
    if(!scalingfactorbin1 <1.0e-10)histbin1->Scale(1.0/scalingfactorbin1);
    if(divbybinwidth&&histbin1)histbin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin2 <1.0e-10)histbin2->Scale(1.0/scalingfactorbin2);
    if(divbybinwidth&&histbin1)histbin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin3 <1.0e-10)histbin3->Scale(1.0/scalingfactorbin3);
    if(divbybinwidth&&histbin3)histbin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin4 <1.0e-10)histbin4->Scale(1.0/scalingfactorbin4);
    if(divbybinwidth&&histbin4)histbin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin5 <1.0e-10)histbin5->Scale(1.0/scalingfactorbin5);
    if(divbybinwidth&&histbin5)histbin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin6 <1.0e-10)histbin6->Scale(1.0/scalingfactorbin6);
    if(divbybinwidth&&histbin6)histbin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorbin7 <1.0e-10)histbin7->Scale(1.0/scalingfactorbin7);
    if(divbybinwidth&&histbin7)histbin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    
    if(!scalingfactormbin1 <1.0e-10)histbin1m->Scale(1.0/scalingfactormbin1);
    if(divbybinwidth&&histbin1m)histbin1m->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactormbin2 <1.0e-10)histbin2m->Scale(1.0/scalingfactormbin2);
    if(divbybinwidth&&histbin2m)histbin2m->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactormbin3 <1.0e-10)histbin3m->Scale(1.0/scalingfactormbin3);
    if(divbybinwidth&&histbin3m)histbin3m->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactormbin4 <1.0e-10)histbin4m->Scale(1.0/scalingfactormbin4);
    if(divbybinwidth&&histbin4m)histbin4m->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactormbin5 <1.0e-10)histbin5m->Scale(1.0/scalingfactormbin5);
    if(divbybinwidth&&histbin5m)histbin5m->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactormbin6 <1.0e-10)histbin6m->Scale(1.0/scalingfactormbin6);
    if(divbybinwidth&&histbin6m)histbin6m->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactormbin7 <1.0e-10)histbin7m->Scale(1.0/scalingfactormbin7);
    if(divbybinwidth&&histbin7m)histbin7m->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    
    if(!scalingfactorMETAbin1 <1.0e-10)histMETAbin1->Scale(1.0/scalingfactorMETAbin1);
    if(divbybinwidth&&histbin1)histMETAbin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin2 <1.0e-10)histMETAbin2->Scale(1.0/scalingfactorMETAbin2);
    if(divbybinwidth&&histbin2)histMETAbin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin3 <1.0e-10)histMETAbin3->Scale(1.0/scalingfactorMETAbin3);
    if(divbybinwidth&&histbin3)histMETAbin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin4 <1.0e-10)histMETAbin4->Scale(1.0/scalingfactorMETAbin4);
    if(divbybinwidth&&histbin4)histMETAbin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin5 <1.0e-10)histMETAbin5->Scale(1.0/scalingfactorMETAbin5);
    if(divbybinwidth&&histbin5)histMETAbin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin6 <1.0e-10)histMETAbin6->Scale(1.0/scalingfactorMETAbin6);
    if(divbybinwidth&&histbin6)histMETAbin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETAbin7 <1.0e-10)histMETAbin7->Scale(1.0/scalingfactorMETAbin7);
    if(divbybinwidth&&histbin7)histMETAbin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    
    if(!scalingfactorMETA2bin1 <1.0e-10)histMETA2bin1->Scale(1.0/scalingfactorMETA2bin1);
    if(divbybinwidth&&histbin1)histMETA2bin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin2 <1.0e-10)histMETA2bin2->Scale(1.0/scalingfactorMETA2bin2);
    if(divbybinwidth&&histbin2)histMETA2bin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin3 <1.0e-10)histMETA2bin3->Scale(1.0/scalingfactorMETA2bin3);
    if(divbybinwidth&&histbin3)histMETA2bin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin4 <1.0e-10)histMETA2bin4->Scale(1.0/scalingfactorMETA2bin4);
    if(divbybinwidth&&histbin4)histMETA2bin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin5 <1.0e-10)histMETA2bin5->Scale(1.0/scalingfactorMETA2bin5);
    if(divbybinwidth&&histbin5)histMETA2bin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin6 <1.0e-10)histMETA2bin6->Scale(1.0/scalingfactorMETA2bin6);
    if(divbybinwidth&&histbin6)histMETA2bin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETA2bin7 <1.0e-10)histMETA2bin7->Scale(1.0/scalingfactorMETA2bin7);
    if(divbybinwidth&&histbin7)histMETA2bin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    
    if(!scalingfactorMETriggerbin1 <1.0e-10)histMETriggerbin1->Scale(1.0/scalingfactorMETriggerbin1);
    if(divbybinwidth&&histbin1)histMETriggerbin1->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin2 <1.0e-10)histMETriggerbin2->Scale(1.0/scalingfactorMETriggerbin2);
    if(divbybinwidth&&histbin2)histMETriggerbin2->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin3 <1.0e-10)histMETriggerbin3->Scale(1.0/scalingfactorMETriggerbin3);
    if(divbybinwidth&&histbin3)histMETriggerbin3->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin4 <1.0e-10)histMETriggerbin4->Scale(1.0/scalingfactorMETriggerbin4);
    if(divbybinwidth&&histbin4)histMETriggerbin4->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin5 <1.0e-10)histMETriggerbin5->Scale(1.0/scalingfactorMETriggerbin5);
    if(divbybinwidth&&histbin5)histMETriggerbin5->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin6 <1.0e-10)histMETriggerbin6->Scale(1.0/scalingfactorMETriggerbin6);
    if(divbybinwidth&&histbin6)histMETriggerbin6->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
    if(!scalingfactorMETriggerbin7 <1.0e-10)histMETriggerbin7->Scale(1.0/scalingfactorMETriggerbin7);
    if(divbybinwidth&&histbin7)histMETriggerbin7->Scale(1.0/(xaxisbinwidth*yaxisbinwidth));
   
}
  if(!isPbPb){
    Bin1->Samediv()->cd();
    histbin1->Write(histo->GetName());
    Bin1->METAdiv()->cd();
    histMETAbin1->Write(histo->GetName());
    if(Bin1->META2()){
      Bin1->META2div()->cd();
      histMETA2bin1->Write(histo->GetName());
    }
    Bin1->METriggerdiv()->cd();
    histMETriggerbin1->Write(histo->GetName());
    if(!collectdivfirst){
      Bin1->Same()->mkdir("mixed_event");
      Bin1->SameDir("mixed_event")->cd();
      histbin1m->Write(histo->GetName());
    }
  }
  if(isPbPb){
    Bin1->Samediv()->cd();
    histbin1->Write(histo->GetName());
    Bin1->METAdiv()->cd();
    histMETAbin1->Write(histo->GetName());
    if(Bin1->META2()){
      Bin1->META2div()->cd();
      histMETA2bin1->Write(histo->GetName());
    }
    Bin1->METriggerdiv()->cd();
    histMETriggerbin1->Write(histo->GetName());
    if(!collectdivfirst){
      Bin1->Same()->mkdir("mixed_event");
      Bin1->SameDir("mixed_event")->cd();
      histbin1m->Write(histo->GetName());
    }
    if(Bin2){
      Bin2->Samediv()->cd();
      histbin2->Write(histo->GetName());
      Bin2->METAdiv()->cd();
      histMETAbin2->Write(histo->GetName());
      if(Bin2->META2()){
	Bin2->META2div()->cd();
	histMETA2bin2->Write(histo->GetName());
      }
      Bin2->METriggerdiv()->cd();
      histMETriggerbin2->Write(histo->GetName());
      if(!collectdivfirst){
	Bin2->Same()->mkdir("mixed_event");
	Bin2->SameDir("mixed_event")->cd();
	histbin2m->Write(histo->GetName());
      }
    }
    if(Bin3){
      Bin3->Samediv()->cd();
      histbin3->Write(histo->GetName());
      Bin3->METAdiv()->cd();
      histMETAbin3->Write(histo->GetName());
      if(Bin3->META2()){
	Bin3->META2div()->cd();
	histMETA2bin3->Write(histo->GetName());
      }
      Bin3->METriggerdiv()->cd();
      histMETriggerbin3->Write(histo->GetName());
      if(!collectdivfirst){
	Bin3->Same()->mkdir("mixed_event");
	Bin3->SameDir("mixed_event")->cd();
	histbin3m->Write(histo->GetName());
      }
    }
    if(Bin4){
      Bin4->Samediv()->cd();
      histbin4->Write(histo->GetName());
      Bin4->METAdiv()->cd();
      histMETAbin4->Write(histo->GetName());
      if(Bin4->META2()){
	Bin4->META2div()->cd();
	histMETA2bin4->Write(histo->GetName());
      }
      Bin4->METriggerdiv()->cd();
      histMETriggerbin4->Write(histo->GetName());
      if(!collectdivfirst){
	Bin4->Same()->mkdir("mixed_event");
	Bin4->SameDir("mixed_event")->cd();
	histbin4m->Write(histo->GetName());
      }
    }
    if(Bin5){
      Bin5->Samediv()->cd();
      histbin5->Write(histo->GetName());
      Bin5->METAdiv()->cd();
      histMETAbin5->Write(histo->GetName());
      if(Bin5->META2()){
        Bin5->META2div()->cd();
	histMETA2bin5->Write(histo->GetName());
    }
      Bin5->METriggerdiv()->cd();
      histMETriggerbin5->Write(histo->GetName());
      if(!collectdivfirst){
	Bin5->Same()->mkdir("mixed_event");
	Bin5->SameDir("mixed_event")->cd();
	histbin5m->Write(histo->GetName());
      }
    }
    if(Bin6){
      Bin6->Samediv()->cd();
      histbin6->Write(histo->GetName());
      Bin6->METAdiv()->cd();
      histMETAbin6->Write(histo->GetName());
      if(Bin6->META2()){
	Bin6->META2div()->cd();
	histMETA2bin6->Write(histo->GetName());
      }
      Bin6->METriggerdiv()->cd();
      histMETriggerbin6->Write(histo->GetName());
      if(!collectdivfirst){
	Bin6->Same()->mkdir("mixed_event");
	Bin6->SameDir("mixed_event")->cd();
	histbin6m->Write(histo->GetName());
      }
    }
    if(Bin7){
      Bin7->Samediv()->cd();
      histbin7->Write(histo->GetName());
      Bin7->METAdiv()->cd();
      histMETAbin7->Write(histo->GetName());
      if(Bin7->META2()){
	Bin7->META2div()->cd();
	histMETA2bin7->Write(histo->GetName());
      }
      Bin7->METriggerdiv()->cd();
      histMETriggerbin7->Write(histo->GetName());
      if(!collectdivfirst){
	Bin7->Same()->mkdir("mixed_event");
	Bin7->SameDir("mixed_event")->cd();
	histbin7m->Write(histo->GetName());
      }
    }
  }
  delete histbin1;delete histMETAbin1;delete histMETA2bin1;delete histMETriggerbin1;delete histbin1m;
  delete histbin2;delete histMETAbin2;delete histMETA2bin2;delete histMETriggerbin2;delete histbin2m;
  delete histbin3;delete histMETAbin3;delete histMETA2bin3;delete histMETriggerbin3;delete histbin3m;
  delete histbin4;delete histMETAbin4;delete histMETA2bin4;delete histMETriggerbin4;delete histbin4m;
  delete histbin5;delete histMETAbin5;delete histMETA2bin5;delete histMETriggerbin5;delete histbin5m;
  delete histbin6;delete histMETAbin6;delete histMETA2bin6;delete histMETriggerbin6;delete histbin6m;
  delete histbin7;delete histMETAbin7;delete histMETA2bin7;delete histMETriggerbin7;delete histbin7m;
}
void CollectHist(TH3D* histo, TList * directories, TObjArray* multdirlist, Bool_t isPbPb){
  BinDirs * All=NULL;BinDirs * Bin1=NULL;BinDirs * Bin2=NULL;BinDirs * Bin3=NULL;BinDirs * Bin4=NULL;BinDirs * Bin5=NULL;BinDirs * Bin6=NULL;BinDirs * Bin7=NULL;
  if(!isPbPb){
   All = dynamic_cast<BinDirs*>(multdirlist->At(0));
  }
  else{
    Bin1 = dynamic_cast<BinDirs*>(multdirlist->At(0));
    if(multdirlist->GetEntries()>1) Bin2 = dynamic_cast<BinDirs*>(multdirlist->At(1));
    if(multdirlist->GetEntries()>2) Bin3 = dynamic_cast<BinDirs*>(multdirlist->At(2));
    if(multdirlist->GetEntries()>3) Bin4 = dynamic_cast<BinDirs*>(multdirlist->At(3));
    if(multdirlist->GetEntries()>4) Bin5 = dynamic_cast<BinDirs*>(multdirlist->At(4));
    if(multdirlist->GetEntries()>5) Bin6 = dynamic_cast<BinDirs*>(multdirlist->At(5));
    if(multdirlist->GetEntries()>6) Bin7 = dynamic_cast<BinDirs*>(multdirlist->At(6));
  }
  //reset the histogram and create clones for all the types we want.
  histo->Reset();
  histo->ResetStats();
  //all bins:
  TH3D* hist 			= dynamic_cast<TH3D*>( histo->Clone(Form("%sdiv"		,histo->GetName())));
  TH3D* histMETA 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivMETA"		,histo->GetName())));
  TH3D* histMETrigger 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivMETrigger"	,histo->GetName())));
  //Bin 1:
  TH3D* histbin1 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin1"		,histo->GetName())));
  TH3D* histMETAbin1 		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin1META"	,histo->GetName())));
  TH3D* histMETriggerbin1 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin1METrigger"	,histo->GetName())));
  //Bin 2:
  TH3D* histbin2 		= new TH3D("emptybin2","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETAbin2 		= new TH3D("emptybinMETA2","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETriggerbin2 	= new TH3D("emptybinMETRIGGER2","title",1,0,1,1,0,1,1,0,1);
  if(Bin2){
    delete histbin2;delete histMETAbin2; delete histMETriggerbin2;
    histbin2		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2"		,histo->GetName())));
    histMETAbin2	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2META"	,histo->GetName())));
    histMETriggerbin2 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin2METrigger"	,histo->GetName())));
  }
  //Bin 3:
  TH3D* histbin3 		= new TH3D("emptybin3","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETAbin3 		= new TH3D("emptybinMETA3","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETriggerbin3 	= new TH3D("emptybinMETRIGGER3","title",1,0,1,1,0,1,1,0,1);
  if(Bin3){
    delete histbin3;delete histMETAbin3; delete histMETriggerbin3;
    histbin3		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3"		,histo->GetName())));
    histMETAbin3	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3META"	,histo->GetName())));
    histMETriggerbin3 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin3METrigger"	,histo->GetName())));
  }
  //Bin 4:
  TH3D* histbin4 		= new TH3D("emptybin4","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETAbin4 		= new TH3D("emptybinMETA4","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETriggerbin4 	= new TH3D("emptybinMETRIGGER4","title",1,0,1,1,0,1,1,0,1);
  if(Bin4){
    delete histbin4;delete histMETAbin4; delete histMETriggerbin4;
    histbin4		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4"		,histo->GetName())));
    histMETAbin4	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4META"	,histo->GetName())));
    histMETriggerbin4 	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin4METrigger"	,histo->GetName())));
  }
  //Bin 5:
  TH3D* histbin5 		= new TH3D("emptybin5","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETAbin5 		= new TH3D("emptybinMETA5","title",1,0,1,1,0,1,1,0,1);
  TH3D* histMETriggerbin5 	= new TH3D("emptybinMETRIGGER5","title",1,0,1,1,0,1,1,0,1);
  if(Bin5){
    delete histbin5;delete histMETAbin5; delete histMETriggerbin5;
    histbin5		= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5"		,histo->GetName())));
    histMETAbin5	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5META"	,histo->GetName())));
    histMETriggerbin5	= dynamic_cast<TH3D*>( histo->Clone(Form("%sdivbin5METrigger"	,histo->GetName())));
  }
  //Doubles to hold the values until they are put into the hists.
  Double_t bincontl   = 0.0;Double_t binerrorl  = 0.0;Double_t bincontlMETA   = 0.0;Double_t binerrorlMETA  = 0.0;
  Double_t BinContent = 0.0;Double_t BinError   = 0.0; Double_t BinContentMETA = 0.0;Double_t BinErrorMETA   = 0.0;
  Double_t BinContentMETrigger = 0.0;Double_t BinErrorMEtrigger   = 0.0; Double_t bincontlMEtrigger   = 0.0;Double_t binerrorlMETrigger  = 0.0;
  Double_t BinContentbin1 = 0.0;Double_t BinErrorbin1   = 0.0;Double_t BinContentMETAbin1 = 0.0;Double_t BinErrorMETAbin1   = 0.0;Double_t BinContentMETriggerbin1 = 0.0;Double_t BinErrorMEtriggerbin1   = 0.0;
  Double_t BinContentbin2 = 0.0;Double_t BinErrorbin2   = 0.0;Double_t BinContentMETAbin2 = 0.0;Double_t BinErrorMETAbin2   = 0.0;Double_t BinContentMETriggerbin2 = 0.0;Double_t BinErrorMEtriggerbin2   = 0.0;
  Double_t BinContentbin3 = 0.0;Double_t BinErrorbin3   = 0.0;Double_t BinContentMETAbin3 = 0.0;Double_t BinErrorMETAbin3   = 0.0;Double_t BinContentMETriggerbin3 = 0.0;Double_t BinErrorMEtriggerbin3   = 0.0;
  Double_t BinContentbin4 = 0.0;Double_t BinErrorbin4   = 0.0;Double_t BinContentMETAbin4 = 0.0;Double_t BinErrorMETAbin4   = 0.0;Double_t BinContentMETriggerbin4 = 0.0;Double_t BinErrorMEtriggerbin4   = 0.0;
  Double_t BinContentbin5 = 0.0;Double_t BinErrorbin5   = 0.0;Double_t BinContentMETAbin5 = 0.0;Double_t BinErrorMETAbin5   = 0.0;Double_t BinContentMETriggerbin5 = 0.0;Double_t BinErrorMEtriggerbin5   = 0.0; 
  //loop over the bins:'
  for(int x=0;x<=histo->GetNbinsX()+1;x++){for(int y=0;y<=histo->GetNbinsY();y++){for(int z=0;z<=histo->GetNbinsZ();z++){
    //loop over all Multiplicity-Vertex bins:
    for(int i=0;i<directories->GetEntries();i++){
      int Mbin = 0;
//       if((TString(directories->At(i)->GetName()).Contains("Z(-10.00)"))&&(TString(directories->At(i)->GetName()).Contains("Z(5.00)")))continue;      
      //find the Multiplicity bin we are in
      for(int j = 1;j<=5;j++){
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()))Mbin = j;
	if(dynamic_cast<BinDirs*>(multdirlist->At(j))->CompareTo("BinM(0.00)->(10.00)")){
	  if(TString(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(0.00)->(5.00)")==0) Mbin = j;
	  if(TString(TString(directories->At(i)->GetName()).Tokenize("Z")->At(0)->GetName()).CompareTo("BinM(5.00)->(10.00)")==0)Mbin = j;
	}
      }
      //test if the relevant histogram exists in this bin:
      if(!dynamic_cast<TH3D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName())))
      {
// 	cout << "No Histogram with name "<< histo->GetName() <<" exists in M-V bin "<< directories->At(i)->GetName() <<endl;
	continue;}
      //extract bin content and error in the relevant bin:
      bincontl 			= dynamic_cast<TH3D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y,z);
      binerrorl 		= dynamic_cast<TH3D*>(All->Same()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y,z);	
      bincontlMETA 		= dynamic_cast<TH3D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y,z);
      binerrorlMETA 		= dynamic_cast<TH3D*>(All->META()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y,z);	      
      bincontlMEtrigger 	= dynamic_cast<TH3D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinContent(x,y,z);
      binerrorlMETrigger 	= dynamic_cast<TH3D*>(All->METrigger()->GetDirectory(Form("%s/divided",directories->At(i)->GetName()))->Get(histo->GetName()))->GetBinError(x,y,z); 
      
      if(bincontl>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContent += bincontl/(binerrorl*binerrorl);
	BinError   += 1.0/(binerrorl*binerrorl);
	if(Mbin == 1){		BinContentbin1 += bincontl/(binerrorl*binerrorl); BinErrorbin1   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 2){	BinContentbin2 += bincontl/(binerrorl*binerrorl); BinErrorbin2   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 3){	BinContentbin3 += bincontl/(binerrorl*binerrorl); BinErrorbin3   += 1.0/(binerrorl*binerrorl);}
	else if(Mbin == 4){	BinContentbin4 += bincontl/(binerrorl*binerrorl); BinErrorbin4   += 1.0/(binerrorl*binerrorl);}	      
	else if(Mbin == 5){ 	BinContentbin5 += bincontl/(binerrorl*binerrorl); BinErrorbin5   += 1.0/(binerrorl*binerrorl);}
      }
      else{
	BinError += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorbin1   += 1.0;
	else if(Mbin == 2)	BinErrorbin2   += 1.0;
	else if(Mbin == 3)	BinErrorbin3   += 1.0;
	else if(Mbin == 4)	BinErrorbin4   += 1.0;	      
	else if(Mbin == 5)	BinErrorbin5   += 1.0;
      }
      if(bincontlMETA>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETA += bincontlMETA/(binerrorlMETA*binerrorlMETA);
	BinErrorMETA   += 1.0/(binerrorlMETA*binerrorlMETA);
	if(Mbin == 1){		BinContentMETAbin1 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin1   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 2){	BinContentMETAbin2 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin2   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 3){	BinContentMETAbin3 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin3   += 1.0/(binerrorlMETA*binerrorlMETA);}
	else if(Mbin == 4){	BinContentMETAbin4 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin4   += 1.0/(binerrorlMETA*binerrorlMETA);}	      
	else if(Mbin == 5){	BinContentMETAbin5 += bincontlMETA/(binerrorlMETA*binerrorlMETA); BinErrorMETAbin5   += 1.0/(binerrorlMETA*binerrorlMETA);}
      }
      else{
	BinErrorMETA += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMETAbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMETAbin2   += 1.0;	      
	else if(Mbin == 3)	BinErrorMETAbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMETAbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMETAbin5   += 1.0;
      }
      if(bincontlMEtrigger>1.0e-10){//if not, there were no fills in the bin and the error is ill defined
	BinContentMETrigger += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger);
	BinErrorMEtrigger   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);
	if(Mbin == 1){		BinContentMETriggerbin1 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin1   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 2){	BinContentMETriggerbin2 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin2   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 3){	BinContentMETriggerbin3 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin3   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
	else if(Mbin == 4){	BinContentMETriggerbin4 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin4   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}	      
	else if(Mbin == 5){	BinContentMETriggerbin5 += bincontlMEtrigger/(binerrorlMETrigger*binerrorlMETrigger); BinErrorMEtriggerbin5   += 1.0/(binerrorlMETrigger*binerrorlMETrigger);}
      }
      else{
	BinErrorMEtrigger += 1.0;//??? Still unshure about this.
	if(Mbin == 1)		BinErrorMEtriggerbin1   += 1.0;
	else if(Mbin == 2)	BinErrorMEtriggerbin2   += 1.0;
	else if(Mbin == 3)	BinErrorMEtriggerbin3   += 1.0;
	else if(Mbin == 4)	BinErrorMEtriggerbin4   += 1.0;
	else if(Mbin == 5)	BinErrorMEtriggerbin5   += 1.0;
      }
      //reset the local bins to zero just to be shure:
      bincontl  = 0.0;		binerrorl = 0.0;
      bincontlMETA =0.0;	binerrorlMETA=0.0;
      bincontlMEtrigger=0.0;	binerrorlMETrigger=0.0;
    }//end loop over M-V bins.
    //normalize the bin and the error for same:
    if(BinError>1.0e-10){		BinContent = BinContent/BinError;						BinError = 1.0/TMath::Sqrt(BinError);				}
    else{				BinContent=0.0;									BinError=0.0;							}
    if(BinErrorbin1>1.0e-10){		BinContentbin1 = BinContentbin1/BinErrorbin1;					BinErrorbin1 = 1.0/TMath::Sqrt(BinErrorbin1);			}
    else{				BinContentbin1=0.0;								BinErrorbin1=0.0;						}
    if(BinErrorbin2>1.0e-10){		BinContentbin2 = BinContentbin2/BinErrorbin2;					BinErrorbin2 = 1.0/TMath::Sqrt(BinErrorbin2);			}
    else{				BinContentbin2=0.0;								BinErrorbin2=0.0;						}
    if(BinErrorbin3>1.0e-10){		BinContentbin3 = BinContentbin3/BinErrorbin3;					BinErrorbin3 = 1.0/TMath::Sqrt(BinErrorbin3);			}
    else{				BinContentbin3=0.0;								BinErrorbin3=0.0;						}
    if(BinErrorbin4>1.0e-10){		BinContentbin4 = BinContentbin4/BinErrorbin4;					BinErrorbin4 = 1.0/TMath::Sqrt(BinErrorbin4);			}
    else{				BinContentbin4=0.0;								BinErrorbin4=0.0;						}
    if(BinErrorbin5>1.0e-10){		BinContentbin5 = BinContentbin5/BinErrorbin5;					BinErrorbin5 = 1.0/TMath::Sqrt(BinErrorbin5);			}
    else{				BinContentbin5=0.0;								BinErrorbin5=0.0;						}
    //normalize the bin and the error for META:    
    if(BinErrorMETA>1.0e-10){		BinContentMETA = BinContentMETA/BinErrorMETA;					BinErrorMETA = 1.0/TMath::Sqrt(BinErrorMETA);			}
    else{				BinContentMETA=0.0;								BinErrorMETA=0.0;						}
    if(BinErrorMETAbin1>1.0e-10){	BinContentMETAbin1 = BinContentMETAbin1/BinErrorMETAbin1;			BinErrorMETAbin1 = 1.0/TMath::Sqrt(BinErrorMETAbin1);		}
    else{				BinContentMETAbin1=0.0;								BinErrorMETAbin1=0.0;						}
    if(BinErrorMETAbin2>1.0e-10){	BinContentMETAbin2 = BinContentMETAbin2/BinErrorMETAbin2;			BinErrorMETAbin2 = 1.0/TMath::Sqrt(BinErrorMETAbin2);		}
    else{				BinContentMETAbin2=0.0;								BinErrorMETAbin2=0.0;						}
    if(BinErrorMETAbin3>1.0e-10){	BinContentMETAbin3 = BinContentMETAbin3/BinErrorMETAbin3;			BinErrorMETAbin3 = 1.0/TMath::Sqrt(BinErrorMETAbin3);		}
    else{				BinContentMETAbin3=0.0;								BinErrorMETAbin3=0.0;						}
    if(BinErrorMETAbin4>1.0e-10){	BinContentMETAbin4 = BinContentMETAbin4/BinErrorMETAbin4;			BinErrorMETAbin4 = 1.0/TMath::Sqrt(BinErrorMETAbin4);		}
    else{				BinContentMETAbin4=0.0;								BinErrorMETAbin4=0.0;						}
    if(BinErrorMETAbin5>1.0e-10){	BinContentMETAbin5 = BinContentMETAbin5/BinErrorMETAbin5;			BinErrorMETAbin5 = 1.0/TMath::Sqrt(BinErrorMETAbin5);		}
    else{				BinContentMETAbin5=0.0;								BinErrorMETAbin5=0.0;						}
    //normalize the bin and the error for METrigger:    
    if(BinErrorMEtrigger>1.0e-10){	BinContentMETrigger = BinContentMETrigger/BinErrorMEtrigger;			BinErrorMEtrigger = 1.0/TMath::Sqrt(BinErrorMEtrigger);		}
    else{				BinContentMETrigger=0.0;							BinErrorMEtrigger=0.0;						}
    if(BinErrorMEtriggerbin1>1.0e-10){	BinContentMETriggerbin1 = BinContentMETriggerbin1/BinErrorMEtriggerbin1;	BinErrorMEtriggerbin1 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin1);	}
    else{				BinContentMETriggerbin1=0.0;							BinErrorMEtriggerbin1=0.0;					}
    if(BinErrorMEtriggerbin2>1.0e-10){	BinContentMETriggerbin2 = BinContentMETriggerbin2/BinErrorMEtriggerbin2;	BinErrorMEtriggerbin2 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin2);	}
    else{				BinContentMETriggerbin2=0.0;							BinErrorMEtriggerbin2=0.0;					}
    if(BinErrorMEtriggerbin3>1.0e-10){	BinContentMETriggerbin3 = BinContentMETriggerbin3/BinErrorMEtriggerbin3;	BinErrorMEtriggerbin3 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin3);	}
    else{				BinContentMETriggerbin3=0.0;							BinErrorMEtriggerbin3=0.0;					}
    if(BinErrorMEtriggerbin4>1.0e-10){	BinContentMETriggerbin4 = BinContentMETriggerbin4/BinErrorMEtriggerbin4;	BinErrorMEtriggerbin4 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin4);	}
    else{				BinContentMETriggerbin4=0.0;							BinErrorMEtriggerbin4=0.0;					}
    if(!BinErrorMEtriggerbin5>1.0e-10){	BinContentMETriggerbin5 = BinContentMETriggerbin5/BinErrorMEtriggerbin5;	BinErrorMEtriggerbin5 = 1.0/TMath::Sqrt(BinErrorMEtriggerbin5);	}
    else{				BinContentMETriggerbin5=0.0;							BinErrorMEtriggerbin5=0.0;					}
    //Set the bin content and error in every histogram.
    if(BinContent>1.0e-10){
      hist->SetBinContent(x,y,z,BinContent);
      hist->SetBinError(x,y,z,BinError);
      histbin1->SetBinContent(x,y,z,BinContentbin1);
      histbin1->SetBinError(x,y,z,BinErrorbin1);
      if(Bin2){
	histbin2->SetBinContent(x,y,z,BinContentbin2);
	histbin2->SetBinError(x,y,z,BinErrorbin2);
	histbin3->SetBinContent(x,y,z,BinContentbin3);
	histbin3->SetBinError(x,y,z,BinErrorbin3);
	histbin4->SetBinContent(x,y,z,BinContentbin4);
	histbin4->SetBinError(x,y,z,BinErrorbin4);
	histbin5->SetBinContent(x,y,z,BinContentbin5);
	histbin5->SetBinError(x,y,z,BinErrorbin5);
      }
      histMETA->SetBinContent(x,y,z,BinContentMETA);
      histMETA->SetBinError(x,y,z,BinErrorMETA);
      histMETAbin1->SetBinContent(x,y,z,BinContentMETAbin1);
      histMETAbin1->SetBinError(x,y,z,BinErrorMETAbin1);
      if(Bin2){
	histMETAbin2->SetBinContent(x,y,z,BinContentMETAbin2);
	histMETAbin2->SetBinError(x,y,z,BinErrorMETAbin2);
	histMETAbin3->SetBinContent(x,y,z,BinContentMETAbin3);
	histMETAbin3->SetBinError(x,y,z,BinErrorMETAbin3);
	histMETAbin4->SetBinContent(x,y,z,BinContentMETAbin4);
	histMETAbin4->SetBinError(x,y,z,BinErrorMETAbin4);
	histMETAbin5->SetBinContent(x,y,z,BinContentMETAbin5);
	histMETAbin5->SetBinError(x,y,z,BinErrorMETAbin5);
      }
      histMETrigger->SetBinContent(x,y,z,BinContentMETrigger);
      histMETrigger->SetBinError(x,y,z,BinErrorMEtrigger);
      histMETriggerbin1->SetBinContent(x,y,z,BinContentMETriggerbin1);
      histMETriggerbin1->SetBinError(x,y,z,BinErrorMEtriggerbin1);
      if(Bin2){
	histMETriggerbin2->SetBinContent(x,y,z,BinContentMETriggerbin2);
	histMETriggerbin2->SetBinError(x,y,z,BinErrorMEtriggerbin2);
	histMETriggerbin3->SetBinContent(x,y,z,BinContentMETriggerbin3);
	histMETriggerbin3->SetBinError(x,y,z,BinErrorMEtriggerbin3);
	histMETriggerbin4->SetBinContent(x,y,z,BinContentMETriggerbin4);
	histMETriggerbin4->SetBinError(x,y,z,BinErrorMEtriggerbin4);
	histMETriggerbin5->SetBinContent(x,y,z,BinContentMETriggerbin5);
	histMETriggerbin5->SetBinError(x,y,z,BinErrorMEtriggerbin5);
      }
    }
    //Reset all:
    BinContent = 0.0;BinError = 0.0;BinContentMETA = 0.0;BinErrorMETA = 0.0;BinContentMETrigger = 0.0;BinErrorMEtrigger = 0.0;
    BinContentbin1 = 0.0;BinErrorbin1 = 0.0;BinContentMETAbin1 = 0.0;BinErrorMETAbin1 = 0.0;BinContentMETriggerbin1 = 0.0;BinErrorMEtriggerbin1 = 0.0;
    BinContentbin2 = 0.0;BinErrorbin2 = 0.0;BinContentMETAbin2 = 0.0;BinErrorMETAbin2 = 0.0;BinContentMETriggerbin2 = 0.0;BinErrorMEtriggerbin2 = 0.0;
    BinContentbin3 = 0.0;BinErrorbin3 = 0.0;BinContentMETAbin3 = 0.0;BinErrorMETAbin3 = 0.0;BinContentMETriggerbin3 = 0.0;BinErrorMEtriggerbin3 = 0.0;
    BinContentbin4 = 0.0;BinErrorbin4 = 0.0;BinContentMETAbin4 = 0.0;BinErrorMETAbin4 = 0.0;BinContentMETriggerbin4 = 0.0;BinErrorMEtriggerbin4 = 0.0;
    BinContentbin5 = 0.0;BinErrorbin5 = 0.0;BinContentMETAbin5 = 0.0;BinErrorMETAbin5 = 0.0;BinContentMETriggerbin5 = 0.0;BinErrorMEtriggerbin5 = 0.0;	
  }}}//end binloop
  //save the histograms in the relevant directories:
  //save the histograms in the relevant directories:
  All->Samediv()->cd();
  hist->Write(histo->GetName());
  All->METAdiv()->cd();
  histMETA->Write(histo->GetName());
  All->METriggerdiv()->cd();
  histMETrigger->Write(histo->GetName());
  Bin1->Samediv()->cd();
  histbin1->Write(histo->GetName());
  Bin1->METAdiv()->cd();
  histMETAbin1->Write(histo->GetName());
  Bin1->METriggerdiv()->cd();
  histMETriggerbin1->Write(histo->GetName());
  if(Bin2){
    Bin2->Samediv()->cd();
    histbin2->Write(histo->GetName());
    Bin2->METAdiv()->cd();
    histMETAbin2->Write(histo->GetName());
    Bin2->METriggerdiv()->cd();
    histMETriggerbin2->Write(histo->GetName());
    Bin3->Samediv()->cd();
    histbin3->Write(histo->GetName());
    Bin3->METAdiv()->cd();
    histMETAbin3->Write(histo->GetName());
    Bin3->METriggerdiv()->cd();
    histMETriggerbin3->Write(histo->GetName());
    Bin4->Samediv()->cd();
    histbin4->Write(histo->GetName());
    Bin4->METAdiv()->cd();
    histMETAbin4->Write(histo->GetName());
    Bin4->METriggerdiv()->cd();
    histMETriggerbin4->Write(histo->GetName());
    Bin5->Samediv()->cd();
    histbin5->Write(histo->GetName());
    Bin5->METAdiv()->cd();
    histMETAbin5->Write(histo->GetName());
    Bin5->METriggerdiv()->cd();
    histMETriggerbin5->Write(histo->GetName());
  }
  delete hist;
  delete histMETA;
  delete histMETrigger;
  delete histbin1;
  delete histMETAbin1;
  delete histMETriggerbin1;  
  delete histbin2;
  delete histMETAbin2;
  delete histMETriggerbin2;  
  delete histbin3;
  delete histMETAbin3;
  delete histMETriggerbin3;  
  delete histbin4;
  delete histMETAbin4;
  delete histMETriggerbin4;  
  delete histbin5;
  delete histMETAbin5;
  delete histMETriggerbin5;  

  
}


  
//function for minimizing the histogram edge:
Double_t Chi2(Double_t scaleMETA, Double_t scaleMETrigger, Double_t scaleMETA2){
  Double_t grad = 0.0;
  TH2D* testhist = dynamic_cast<TH2D*>(gSameEvent->Clone("testhist"));
  testhist->Add(gMETA,-1.0*scaleMETA);
  testhist->Add(gMETrigger,-1.0*scaleMETrigger); 
  testhist->Add(gMETA2,-1.0*scaleMETA2); 
  int pihalfbin = testhist->GetXaxis()->FindBin(TMath::Pi()/2.0);
  double gradloc = 0.0;
  double errgradloc = 0.0;
  for(int x=1;x<testhist->GetNbinsX();x++){
    gradloc = testhist->GetBinContent(x,1) - testhist->GetBinContent(x+1,1);
    errgradloc = testhist->GetBinError(x,1)*testhist->GetBinError(x,1) + testhist->GetBinError(x+1,1)* testhist->GetBinError(x+1,1);
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;
    gradloc = testhist->GetBinContent(x,testhist->GetNbinsY()) - testhist->GetBinContent(x+1,testhist->GetNbinsY());
    errgradloc = testhist->GetBinError(x,testhist->GetNbinsY())*testhist->GetBinError(x,testhist->GetNbinsY()) + testhist->GetBinError(x+1,testhist->GetNbinsY())* testhist->GetBinError(x+1,testhist->GetNbinsY());
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;
    gradloc =testhist->GetBinContent(1,x) - testhist->GetBinContent(1,x+1);
    errgradloc = testhist->GetBinError(1,x)*testhist->GetBinError(1,x) + testhist->GetBinError(1,x+1)* testhist->GetBinError(1,x+1);
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;
    gradloc = testhist->GetBinContent(testhist->GetNbinsY(),x) - testhist->GetBinContent(testhist->GetNbinsY(),x+1);
    errgradloc = testhist->GetBinError(testhist->GetNbinsY(),x)*testhist->GetBinError(testhist->GetNbinsY(),x) + testhist->GetBinError(testhist->GetNbinsY(),x+1)* testhist->GetBinError(testhist->GetNbinsY(),x+1);
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;
    gradloc =testhist->GetBinContent(x,pihalfbin) - testhist->GetBinContent(x+1,pihalfbin);
    errgradloc = testhist->GetBinError(x,pihalfbin)*testhist->GetBinError(x,pihalfbin) + testhist->GetBinError(x+1,pihalfbin)* testhist->GetBinError(x+1,pihalfbin);
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;
    gradloc = testhist->GetBinContent(x,pihalfbin) - testhist->GetBinContent(x+1,pihalfbin);
    errgradloc = testhist->GetBinError(x,pihalfbin)*testhist->GetBinError(x,pihalfbin) + testhist->GetBinError(x+1,pihalfbin)* testhist->GetBinError(x+1,pihalfbin);
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;
    gradloc =testhist->GetBinContent(pihalfbin,x) - testhist->GetBinContent(pihalfbin,x+1);
    errgradloc = testhist->GetBinError(pihalfbin,x)*testhist->GetBinError(pihalfbin,x) + testhist->GetBinError(pihalfbin,x+1)* testhist->GetBinError(pihalfbin,x+1);
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;
    gradloc = testhist->GetBinContent(pihalfbin,x) - testhist->GetBinContent(pihalfbin,x+1);
    errgradloc = testhist->GetBinError(pihalfbin,x)*testhist->GetBinError(pihalfbin,x) + testhist->GetBinError(pihalfbin,x+1)* testhist->GetBinError(pihalfbin,x+1);
    if(errgradloc>1.0E-25)grad += gradloc*gradloc;    
  }
  delete testhist;
  return 100000*grad;
}
Double_t Chi2(Double_t scaleMETA, Double_t scaleMETrigger){
  TH2D* testhist = dynamic_cast<TH2D*>(gSameEvent->Clone("testhist"));
  testhist->Add(gMETA,-1.0*scaleMETA);
  testhist->Add(gMETrigger,-1.0*scaleMETrigger); 
  testhist->Add(gMETA2,-1.0*scaleMETA); 
  int pihalfbin = testhist->GetXaxis()->FindBin(TMath::Pi()/2.0);
  double grad = 0.0;
  double gradloc = 0.0;
  double errgradloc = 0.0;
  for(int x=1;x<testhist->GetNbinsX();x++){
    gradloc = testhist->GetBinContent(x,1) - testhist->GetBinContent(x+1,1);
    errgradloc = testhist->GetBinError(x,1)*testhist->GetBinError(x,1) + testhist->GetBinError(x+1,1)*testhist->GetBinError(x+1,1);
    if(errgradloc>1.0E-25)grad += (gradloc*gradloc)/(errgradloc);
    gradloc = testhist->GetBinContent(x,testhist->GetNbinsY()) - testhist->GetBinContent(x+1,testhist->GetNbinsY());
    errgradloc = testhist->GetBinError(x,testhist->GetNbinsY())*testhist->GetBinError(x,testhist->GetNbinsY()) + testhist->GetBinError(x+1,testhist->GetNbinsY())*testhist->GetBinError(x+1,testhist->GetNbinsY());
    if(errgradloc>1.0E-25)grad += (gradloc*gradloc)/(errgradloc);
    gradloc = testhist->GetBinContent(x,pihalfbin) - testhist->GetBinContent(x+1,pihalfbin);
    errgradloc = testhist->GetBinError(x,pihalfbin)*testhist->GetBinError(x,pihalfbin) + testhist->GetBinError(x+1,pihalfbin)*testhist->GetBinError(x+1,pihalfbin);
    if(errgradloc>1.0E-25)grad += (gradloc*gradloc)/(errgradloc);
    gradloc = testhist->GetBinContent(pihalfbin,x) - testhist->GetBinContent(pihalfbin,x+1);
    errgradloc = testhist->GetBinError(pihalfbin,x)*testhist->GetBinError(pihalfbin,x) + testhist->GetBinError(pihalfbin,x+1)*testhist->GetBinError(pihalfbin,x+1);
    if(errgradloc>1.0E-25)grad += (gradloc*gradloc)/(errgradloc);  
  }
  delete testhist;
  return 100000*grad;
}


Double_t Chi2(Double_t scaleMETrigger){
  Double_t grad = 0.0;
  TH2D* testhist = dynamic_cast<TH2D*>(gSameEventETA->Clone("testhist"));
  testhist->Add(gMETriggerETA,-1.0*scaleMETrigger); 
  int pihalfbin = testhist->GetYaxis()->FindBin(TMath::Pi()/2.0);
  int firstbin = testhist->GetXaxis()->FindBin(-1.2);
  int lastbin = testhist->GetXaxis()->FindBin(1.2);
  Double_t gradloc = 0.0;
  Double_t errgradloc = 0.0;
  for(int x=1;x<testhist->GetNbinsX();x++){
    gradloc = testhist->GetBinContent(x,pihalfbin) - testhist->GetBinContent(x+1,pihalfbin);
    errgradloc = testhist->GetBinError(x,pihalfbin)*testhist->GetBinError(x,pihalfbin)+testhist->GetBinError(x+1,pihalfbin)*testhist->GetBinError(x+1,pihalfbin);
    if(errgradloc>1.0E-25){
      grad += (gradloc*gradloc)/(errgradloc);
    }
    gradloc = testhist->GetBinContent(x,pihalfbin-1) - testhist->GetBinContent(x+1,pihalfbin-1);
    errgradloc = testhist->GetBinError(x,pihalfbin-1)*testhist->GetBinError(x,pihalfbin-1)+testhist->GetBinError(x+1,pihalfbin-1)*testhist->GetBinError(x+1,pihalfbin-1);
    if(errgradloc>1.0E-25){
      grad += (gradloc*gradloc)/(errgradloc);
    }    
    gradloc = testhist->GetBinContent(x,pihalfbin+1) - testhist->GetBinContent(x+1,pihalfbin+1);
    errgradloc = testhist->GetBinError(x,pihalfbin+1)*testhist->GetBinError(x,pihalfbin+1)+testhist->GetBinError(x+1,pihalfbin+1)*testhist->GetBinError(x+1,pihalfbin+1);
    if(errgradloc>1.0E-25){
      grad += (gradloc*gradloc)/(errgradloc);
    }
    gradloc = testhist->GetBinContent(x,pihalfbin-2) - testhist->GetBinContent(x+1,pihalfbin-2);
    errgradloc = testhist->GetBinError(x,pihalfbin-2)*testhist->GetBinError(x,pihalfbin-2)+testhist->GetBinError(x+1,pihalfbin-2)*testhist->GetBinError(x+1,pihalfbin-2);
    if(errgradloc>1.0E-25){
      grad += (gradloc*gradloc)/(errgradloc);
    }    
    gradloc = testhist->GetBinContent(x,pihalfbin+2) - testhist->GetBinContent(x+1,pihalfbin+2);
    errgradloc = testhist->GetBinError(x,pihalfbin+2)*testhist->GetBinError(x,pihalfbin+2)+testhist->GetBinError(x+1,pihalfbin+2)*testhist->GetBinError(x+1,pihalfbin+2);
    if(errgradloc>1.0E-25){
      grad += (gradloc*gradloc)/(errgradloc);
    }  }
  delete testhist;
  return 100000*grad;
}


void FcnFitScale(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t) {
  //
  // Chi2 function
  // par[0] - scale factor for the META 	background histogram
  // par[1] - scale factor for the METrigger 	background histogram
  // par[2] - scale factor for the META2 	background histogram
  //
  Double_t chi2 = Chi2(par[0],par[1],par[2]);
  f = chi2;
} 
void FcnFitScale1(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t) {
  //
  // Chi2 function
  // par[0] - scale factor for the META 	background histogram
  // par[1] - scale factor for the METrigger 	background histogram
  //
  Double_t chi2 = Chi2(par[0],par[1]);
  f = chi2;
} 

void FcnFitScale2(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t) {
  //
  // Chi2 function
  // par[0] - scale factor for the METrigger 	background histogram
  //
  Double_t chi2 = Chi2(par[0]);
  f = chi2;
} 

void GetScalingFactors(BinDirs * BinDir,const char* in, const char* out,Double_t * METAscale, Double_t * METriggerscale,bool ispp, bool isPbPb, int iter)
{
  //Function to find the relevant scaling factors.
  //iter == 1  means finding the scaling factors in PhiPhi
  //iter == 2  means finding them in PhiEta12SS, but correcting only METrigger
  //Clean and create the output directory
  resultsdirectory(BinDir->Same(),out);
  if(!BinDir->SameDir(in)) return;

  
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * SameHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METAHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_META")));
  tmp = BinDir->METriggerdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METRiggerHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_METrigger")));
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DEta_12");if(!tmp){return;}
  TH2D * SameHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DEta_12");if(!tmp){return;}  
  TH2D * METAHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_META")));
   tmp = BinDir->METriggerdiv()->Get("DPhi_1_DEta_12");if(!tmp){return;}
  TH2D * METRiggerHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_METrigger")));  
  delete tmp;

  gSameEvent 	= SameHistPhiPhi;
  gMETA 	= METAHistPhiPhi;
  gMETrigger 	= METRiggerHistPhiPhi;
  gSameEventETA	= SameHistPhiEta;
  gMETAETA 	= METAHistPhiEta;
  gMETriggerETA	= METRiggerHistPhiEta;  
  

  TMinuit* minuit=new TMinuit(1);
  if(iter == 1){
    minuit->SetFCN(FcnFitScale1);
    Double_t arglist[2];
    Int_t ierflg=0;
    arglist[0] = 1;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg); 
    // Set starting values and step sizes for parameters
    Double_t vstart[2];
    Double_t step[2];
    vstart[0] = *METAscale;
    step[0] = 0.0001; 
    vstart[1] = *METriggerscale;
    step[1] = 0.0001;  
    
    minuit->mnparm(0, "METAscale", 	vstart[0], step[0], 0.0, 1000000.0, ierflg);
    minuit->mnparm(1, "METriggerscale", 	vstart[1], step[1], 0.0, 1000000.0, ierflg);
    
    arglist[0] = 500;
    arglist[1] = 1.;

    minuit->SetMaxIterations(10000);
    minuit->mnexcm("SIMPLEX", arglist ,2,ierflg);
    minuit->mnexcm("MIGRAD" , arglist ,2,ierflg); 
    Double_t temp =0.0;
    Double_t tempe = 0.0;
    minuit->GetParameter(0, temp, tempe); 
    *METAscale = temp;
    minuit->GetParameter(1, temp, tempe); 
    *METriggerscale = temp;
//     cout << *METAscale <<" " <<*METriggerscale<<endl;
  }
  if(iter == 2){
    minuit->SetFCN(FcnFitScale2);
  
    Double_t arglist[1];
    Int_t ierflg=0;
    arglist[0] = 1;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg); 
    // Set starting values and step sizes for parameters
    Double_t vstart[1];
    Double_t step[1];
    vstart[0] = *METriggerscale;
    step[0] = 0.0001; 
    
    minuit->mnparm(0, "METriggerscale", vstart[0], step[0], 0.0, 1000000.0, ierflg);
    
    arglist[0] = 500;
    arglist[1] = 1.;

    minuit->SetMaxIterations(10000);
    minuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
    minuit->mnexcm("MIGRAD" , arglist ,1,ierflg); 
    Double_t temp =0.0;
    Double_t tempe = 0.0;
    minuit->GetParameter(0, temp, tempe); 
    *METriggerscale = temp;
    *METAscale = 0.00;
//     cout << *METAscale <<" " <<*METriggerscale<<endl;    
    
  }
  
  
  
  
  delete SameHistPhiPhi;delete METAHistPhiPhi;delete METRiggerHistPhiPhi;delete SameHistPhiEta;delete METAHistPhiEta;delete METRiggerHistPhiEta;

  
}

void GetScalingFactors(BinDirs * BinDir,const char* in, const char* out,Double_t * METAscale, Double_t * META2scale, Double_t * METriggerscale,bool ispp, bool isPbPb, int iter)
{
  //Function to find the relevant scaling factors.
  //iter == 1  means finding the scaling factors in PhiPhi
  //iter == 2  means finding them in PhiEta12SS, but correcting only METrigger
  //Clean and create the output directory
  resultsdirectory(BinDir->Same(),out);
  if(!BinDir->SameDir(in)) return;

  
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * SameHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METAHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_META")));
  TH2D * META2HistPhiPhi = NULL;
  if(BinDir->META2()){
    tmp = BinDir->META2div()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
    META2HistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_META2")));
  }
  tmp = BinDir->METriggerdiv()->Get("DPhi_1_DPHI_2");if(!tmp){return;}
  TH2D * METRiggerHistPhiPhi = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DPHI_2_METrigger")));
  tmp = BinDir->SameDir(in)->Get("DPhi_1_DEta_12_SameSide");if(!tmp){return;}
  TH2D * SameHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_Same")));
  tmp = BinDir->METAdiv()->Get("DPhi_1_DEta_12_SameSide");if(!tmp){return;}  
  TH2D * METAHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_META")));
  TH2D * META2HistPhiEta = NULL;  
  if(BinDir->META2()){
    tmp = BinDir->METAdiv()->Get("DPhi_1_DEta_12_SameSide");if(!tmp){return;}  
    META2HistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_META2")));
  }
  tmp = BinDir->METriggerdiv()->Get("DPhi_1_DEta_12_SameSide");if(!tmp){return;}
  TH2D * METRiggerHistPhiEta = dynamic_cast<TH2D*>(tmp->Clone(Form("DPhi_1_DEta_12_METrigger")));  
  delete tmp;

  gSameEvent 	= SameHistPhiPhi;
  gMETA 	= METAHistPhiPhi;
  gMETA2 	= META2HistPhiPhi;
  gMETrigger 	= METRiggerHistPhiPhi;
  gSameEventETA	= SameHistPhiEta;
  gMETAETA 	= METAHistPhiEta;
  gMETA2ETA 	= META2HistPhiEta;
  gMETriggerETA	= METRiggerHistPhiEta;  
  

  TMinuit* minuit=new TMinuit(1);
  if(iter == 1){
    minuit->SetFCN(FcnFitScale);

    Double_t arglist[2];
    Int_t ierflg=0;
    arglist[0] = 1;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg); 
    // Set starting values and step sizes for parameters
    Double_t vstart[3];
    Double_t step[3];
    vstart[0] = *METAscale;
    step[0] = 0.0001; 
    vstart[1] = *METriggerscale;
    step[1] = 0.0001;  
    vstart[2] = *META2scale;
    step[2] = 0.0001;    
    minuit->mnparm(0, "METAscale", 	vstart[0], step[0], 0.0, 1000000.0, ierflg);
    minuit->mnparm(1, "METriggerscale", 	vstart[1], step[1], 0.0, 1000000.0, ierflg);
    minuit->mnparm(2, "META2scale", 	vstart[2], step[2], 0.0, 1000000.0, ierflg);

    
    arglist[0] = 500;
    arglist[1] = 1.;

    minuit->SetMaxIterations(10000);
    minuit->mnexcm("SIMPLEX", arglist ,3,ierflg);
    minuit->mnexcm("MIGRAD" , arglist ,3,ierflg); 
    Double_t temp =0.0;
    Double_t tempe = 0.0;
    minuit->GetParameter(0, temp, tempe); 
    *METAscale = temp;
    minuit->GetParameter(1, temp, tempe); 
    *METriggerscale = temp;
    minuit->GetParameter(2, temp, tempe); 
    *META2scale = temp;
//     cout << *METAscale <<" " << *META2scale <<" " <<*METriggerscale<<endl;


  }
  if(iter == 2){
    minuit->SetFCN(FcnFitScale2);
  
    Double_t arglist[1];
    Int_t ierflg=0;
    arglist[0] = 1;
    minuit->mnexcm("SET ERR", arglist, 1, ierflg); 
    // Set starting values and step sizes for parameters
    Double_t vstart[1];
    Double_t step[1];
    vstart[0] = *METriggerscale;
    step[0] = 0.0001; 
    
    minuit->mnparm(0, "METriggerscale", vstart[0], step[0], 0.0, 1000000.0, ierflg);
    
    arglist[0] = 500;
    arglist[1] = 1.;

    minuit->SetMaxIterations(10000);
    minuit->mnexcm("SIMPLEX", arglist ,1,ierflg);
    minuit->mnexcm("MIGRAD" , arglist ,1,ierflg); 
    Double_t temp =0.0;
    Double_t tempe = 0.0;
    minuit->GetParameter(0, temp, tempe); 
    *METriggerscale = temp;
    *METAscale = 0.00;
    *META2scale = 0.00;
//     cout << *METAscale <<" " << *META2scale <<" " <<*METriggerscale<<endl;    
    
  }

  delete SameHistPhiPhi;delete METAHistPhiPhi;delete META2HistPhiPhi;delete METRiggerHistPhiPhi;delete SameHistPhiEta;delete METAHistPhiEta;delete META2HistPhiEta;delete METRiggerHistPhiEta;

  
}





void AddSigBins(TH2D* hist, TH2D* METAhist, TH2D* METriggerhist){
  //Calculates hist -METAhist -METriggerhist for nonzero bins.
  if(!hist||!METAhist||!METriggerhist)return;
  for(int i = 0; i<=hist->GetNbinsX();i++){
    for(int j=0; j<=hist->GetNbinsY();j++){
      if(hist->GetBinContent(i,j)>1.0E-10){
	Double_t content = hist->GetBinContent(i,j);
	Double_t error = hist->GetBinError(i,j);
	Double_t METAcontent = METAhist->GetBinContent(i,j);
	Double_t METAerror = METAhist->GetBinError(i,j);
	Double_t METriggercontent = METriggerhist->GetBinContent(i,j);
	Double_t METriggererror = METriggerhist->GetBinError(i,j);
	content -= METAcontent + METriggercontent;
	hist->SetBinContent(i,j,content);
	hist->SetBinError(i,j,TMath::Sqrt(error*error + METAerror*METAerror+METriggererror*METriggererror));
      }
    }
  }
}

void AddSigBins(TH2D* hist, TH2D* METAhist, TH2D* META2hist, TH2D* METriggerhist){
  //Calculates hist -METAhist -METriggerhist for nonzero bins.
  if(!hist||!METAhist||!META2hist||!METriggerhist)return;
  for(int i = 0; i<=hist->GetNbinsX()+1;i++){
    for(int j=0; j<=hist->GetNbinsY()+1;j++){
      if(hist->GetBinContent(i,j)>1.0E-10){
	Double_t content = hist->GetBinContent(i,j);
	Double_t error = hist->GetBinError(i,j);
	Double_t METAcontent = METAhist->GetBinContent(i,j);
	Double_t METAerror = METAhist->GetBinError(i,j);
	Double_t META2content = META2hist->GetBinContent(i,j);
	Double_t META2error = META2hist->GetBinError(i,j);
	Double_t METriggercontent = METriggerhist->GetBinContent(i,j);
	Double_t METriggererror = METriggerhist->GetBinError(i,j);
	content -= METAcontent + METriggercontent + META2content;
	hist->SetBinContent(i,j,content);
	hist->SetBinError(i,j,TMath::Sqrt(error*error + METAerror*METAerror+METriggererror*METriggererror + META2error*META2error));
      }
    }
  }
}


void Correct(const char* histname, BinDirs * BinDir,const char* in,const char* out, Double_t  METAscale, Double_t  METriggerscale)
{
//   cout << histname << " "<<BinDir->Same()->GetName()<<"  "<< METAscale<<" "<< BinDir->META()->GetName()<< "  "<< METAscale<<" "<< BinDir->METrigger()->GetName()<<endl;
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get(histname);if(!tmp){return;}
  TH2D * SameHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname)));
  tmp = BinDir->METAdiv()->Get(histname);if(!tmp){return;}
  TH2D * METAHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname)));
  tmp = BinDir->METriggerdiv()->Get(histname);if(!tmp){return;}
  TH2D * METriggerHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname)));  
  delete tmp;
  //Scale the Mixed hists with the correct scales:
  METAHist->Scale(METAscale);
  METriggerHist->Scale(METriggerscale);
 //perform the removal:  
  AddSigBins(SameHist,METAHist,METriggerHist);
//   SameHist->Add(METAHist,-1.0*METAscale);SameHist->Add(METriggerHist,-1.0*METriggerscale);
//   RemovePlateau(SameHist->Integral(SameHist->GetXaxis()->FindBin(-1.6),SameHist->GetXaxis()->FindBin(1.6),SameHist->GetYaxis()->FindBin(TMath::Pi()*0.5),SameHist->GetYaxis()->FindBin(TMath::Pi()*0.5))/(SameHist->GetXaxis()->FindBin(1.6)-SameHist->GetXaxis()->FindBin(-1.6)),SameHist);
  TCanvas * canvas = Makecanvas(SameHist,Form("%sCanvas",histname),false);
  BinDir->SameDir(out)->cd();
  SameHist->Write(histname);
  canvas->Write();
  delete canvas;delete SameHist;delete METAHist; delete METriggerHist;
}

void Correct(const char* histname, BinDirs * BinDir,const char* in,const char* out, Double_t  METAscale,Double_t  META2scale, Double_t  METriggerscale)
{
//   cout << histname << " "<<BinDir->Same()->GetName()<<"  "<< METAscale<<" "<< BinDir->META()->GetName()<< "  "<< METAscale<<" "<< BinDir->METrigger()->GetName()<<endl;
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get(histname);if(!tmp){return;}
  TH2D * SameHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname)));
  tmp = BinDir->METAdiv()->Get(histname);if(!tmp){return;}
  TH2D * METAHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname)));
  tmp = BinDir->META2div()->Get(histname);if(!tmp){return;}
  TH2D * META2Hist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA2",histname)));  
  tmp = BinDir->METriggerdiv()->Get(histname);if(!tmp){return;}
  TH2D * METriggerHist = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname)));  
  delete tmp;
  //Scale the Mixed hists with the correct scales:
  METAHist->Scale(METAscale);
  META2Hist->Scale(META2scale);
  METriggerHist->Scale(METriggerscale);
 //perform the removal:  
  AddSigBins(SameHist,METAHist,META2Hist,METriggerHist);
//   SameHist->Add(METAHist,-1.0*METAscale);SameHist->Add(METriggerHist,-1.0*METriggerscale);
//   RemovePlateau(SameHist->Integral(SameHist->GetXaxis()->FindBin(-1.6),SameHist->GetXaxis()->FindBin(1.6),SameHist->GetYaxis()->FindBin(TMath::Pi()*0.5),SameHist->GetYaxis()->FindBin(TMath::Pi()*0.5))/(SameHist->GetXaxis()->FindBin(1.6)-SameHist->GetXaxis()->FindBin(-1.6)),SameHist);
  TCanvas * canvas = Makecanvas(SameHist,Form("%sCanvas",histname),false,true);
  BinDir->SameDir(out)->cd();
  SameHist->Write(histname);
  canvas->Write();
  delete canvas;delete SameHist;delete METAHist;delete META2Hist; delete METriggerHist;
}

void Correct(const char* histname1,const char* histname2,const char* histname3,const char* histname4, BinDirs * BinDir,const char* in,const char* out, Double_t  METAscale, Double_t  METriggerscale)
{
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get(histname1);if(!tmp){return;}
  TH2D * SameHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname1)));
  tmp = BinDir->METAdiv()->Get(histname1);if(!tmp){return;}
  TH2D * METAHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname1)));
  tmp = BinDir->METriggerdiv()->Get(histname1);if(!tmp){return;}
  TH2D * METriggerHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname1)));
  tmp = BinDir->SameDir(in)->Get(histname2);if(!tmp){return;}
  TH2D * SameHist2 = dynamic_cast<TH2D*>(BinDir->SameDir(in)->Get(histname2)->Clone(Form("%sSame",histname2)));
  tmp = BinDir->METAdiv()->Get(histname2);if(!tmp){return;}
  TH2D * METAHist2 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname2)));
  tmp = BinDir->METriggerdiv()->Get(histname2);if(!tmp){return;}
  TH2D * METriggerHist2 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname2)));
  tmp = BinDir->SameDir(in)->Get(histname3);if(!tmp){return;}
  TH2D * SameHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname3)));
  tmp = BinDir->METAdiv()->Get(histname3);if(!tmp){return;}
  TH2D * METAHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname3)));
  tmp = BinDir->METriggerdiv()->Get(histname3);if(!tmp){return;}
  TH2D * METriggerHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname3)));  
  tmp = BinDir->SameDir(in)->Get(histname4);if(!tmp){return;}
  TH2D * SameHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname4)));
  tmp = BinDir->METAdiv()->Get(histname4);if(!tmp){return;}
  TH2D * METAHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname4)));
  tmp = BinDir->METriggerdiv()->Get(histname4);if(!tmp){return;}
  TH2D * METriggerHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname4)));
  delete tmp;
  //Scale the Mixed hists with the correct scales:
  //perform the removal:  
  METAHist1->Scale(METAscale);
  METriggerHist1->Scale(METriggerscale);
  AddSigBins(SameHist1,METAHist1,METriggerHist1);
//   SameHist1->Add(METAHist1,-1.0*METAscale);SameHist1->Add(METriggerHist1,-1.0*METriggerscale);
  METAHist2->Scale(METAscale);
  METriggerHist2->Scale(METriggerscale);
  AddSigBins(SameHist2,METAHist2,METriggerHist2);
//   SameHist2->Add(METAHist2,-1.0*METAscale);SameHist2->Add(METriggerHist2,-1.0*METriggerscale);
  METAHist3->Scale(METAscale);
  METriggerHist3->Scale(METriggerscale);
  AddSigBins(SameHist3,METAHist3,METriggerHist3);
//   SameHist3->Add(METAHist3,-1.0*METAscale);SameHist3->Add(METriggerHist3,-1.0*METriggerscale);
  METAHist4->Scale(METAscale);
  METriggerHist4->Scale(METriggerscale);
  AddSigBins(SameHist4,METAHist4,METriggerHist4);
//   SameHist4->Add(METAHist4,-1.0*METAscale);SameHist4->Add(METriggerHist4,-1.0*METriggerscale);
  TCanvas * canvas = Makecanvas(SameHist1,SameHist2,SameHist3,SameHist4,"DPHIDPHI",false);
  BinDir->SameDir(out)->cd();
  SameHist1->Write(histname1);
  SameHist2->Write(histname2);
  SameHist3->Write(histname3);
  SameHist4->Write(histname4);
  canvas->Write();
  delete SameHist1;delete SameHist2;delete SameHist3;delete SameHist4;
  delete METAHist1;delete METAHist2;delete METAHist3;delete METAHist4;
  delete METriggerHist1;delete METriggerHist2;delete METriggerHist3;delete METriggerHist4;
  delete canvas;
}

void Correct(const char* histname1,const char* histname2,const char* histname3,const char* histname4, BinDirs * BinDir,const char* in,const char* out, Double_t  METAscale, Double_t  META2scale, Double_t  METriggerscale)
{
  TObject * tmp;
  tmp = BinDir->SameDir(in)->Get(histname1);if(!tmp){return;}
  TH2D * SameHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname1)));
  tmp = BinDir->METAdiv()->Get(histname1);if(!tmp){return;}
  TH2D * METAHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname1)));
  tmp = BinDir->META2div()->Get(histname1);if(!tmp){return;}
  TH2D * META2Hist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA2",histname1)));
  tmp = BinDir->METriggerdiv()->Get(histname1);if(!tmp){return;}
  TH2D * METriggerHist1 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname1)));
  tmp = BinDir->SameDir(in)->Get(histname2);if(!tmp){return;}
  TH2D * SameHist2 = dynamic_cast<TH2D*>(BinDir->SameDir(in)->Get(histname2)->Clone(Form("%sSame",histname2)));
  tmp = BinDir->METAdiv()->Get(histname2);if(!tmp){return;}
  TH2D * METAHist2 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname2)));
  tmp = BinDir->META2div()->Get(histname2);if(!tmp){return;}
  TH2D * META2Hist2 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA2",histname2)));
  tmp = BinDir->METriggerdiv()->Get(histname2);if(!tmp){return;}
  TH2D * METriggerHist2 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname2)));
  tmp = BinDir->SameDir(in)->Get(histname3);if(!tmp){return;}
  TH2D * SameHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname3)));
  tmp = BinDir->METAdiv()->Get(histname3);if(!tmp){return;}
  TH2D * METAHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname3)));
  tmp = BinDir->META2div()->Get(histname3);if(!tmp){return;}
  TH2D * META2Hist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA2",histname3)));
  tmp = BinDir->METriggerdiv()->Get(histname3);if(!tmp){return;}
  TH2D * METriggerHist3 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname3)));  
  tmp = BinDir->SameDir(in)->Get(histname4);if(!tmp){return;}
  TH2D * SameHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sSame",histname4)));
  tmp = BinDir->METAdiv()->Get(histname4);if(!tmp){return;}
  TH2D * METAHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA",histname4)));
  tmp = BinDir->META2div()->Get(histname4);if(!tmp){return;}
  TH2D * META2Hist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETA2",histname4)));
  tmp = BinDir->METriggerdiv()->Get(histname4);if(!tmp){return;}
  TH2D * METriggerHist4 = dynamic_cast<TH2D*>(tmp->Clone(Form("%sMETrigger",histname4)));
  delete tmp;
  //Scale the Mixed hists with the correct scales:
  //perform the removal:  
  METAHist1->Scale(METAscale);
  META2Hist1->Scale(META2scale);
  METriggerHist1->Scale(METriggerscale);
  AddSigBins(SameHist1,METAHist1,META2Hist1,METriggerHist1);
//   SameHist1->Add(METAHist1,-1.0*METAscale);SameHist1->Add(METriggerHist1,-1.0*METriggerscale);
  METAHist2->Scale(METAscale);
  META2Hist2->Scale(META2scale);
  METriggerHist2->Scale(METriggerscale);
  AddSigBins(SameHist2,METAHist2,META2Hist2,METriggerHist2);
//   SameHist2->Add(METAHist2,-1.0*METAscale);SameHist2->Add(METriggerHist2,-1.0*METriggerscale);
  METAHist3->Scale(METAscale);
  META2Hist3->Scale(META2scale);
  METriggerHist3->Scale(METriggerscale);
  AddSigBins(SameHist3,METAHist3,META2Hist3,METriggerHist3);
//   SameHist3->Add(METAHist3,-1.0*METAscale);SameHist3->Add(METriggerHist3,-1.0*METriggerscale);
  METAHist4->Scale(METAscale);
  META2Hist4->Scale(META2scale);
  METriggerHist4->Scale(METriggerscale);
  AddSigBins(SameHist4,METAHist4,META2Hist4,METriggerHist4);
//   SameHist4->Add(METAHist4,-1.0*METAscale);SameHist4->Add(METriggerHist4,-1.0*METriggerscale);
  TCanvas * canvas = Makecanvas(SameHist1,SameHist2,SameHist3,SameHist4,"DPHIDPHI",false);
  BinDir->SameDir(out)->cd();
  SameHist1->Write(histname1);
  SameHist2->Write(histname2);
  SameHist3->Write(histname3);
  SameHist4->Write(histname4);
  canvas->Write();
  delete SameHist1;delete SameHist2;delete SameHist3;delete SameHist4;
  delete METAHist1;delete METAHist2;delete METAHist3;delete METAHist4;
  delete METriggerHist1;delete METriggerHist2;delete METriggerHist3;delete METriggerHist4;
  delete canvas;
}

void Correctpp(BinDirs* BinDir,bool metaindependent)
{
  //If the bin is empty, skip the entire bin:
  if(!dynamic_cast<TH1D*>(BinDir->Samediv()->Get("number_of_triggers"))){cout << "skipping empty bin "<<BinDir->path()<<endl;return;}
  if(dynamic_cast<TH1D*>(BinDir->Samediv()->Get("number_of_triggers"))->GetEffectiveEntries()<1){cout << "skipping empty bin "<<BinDir->path()<<endl;return;}
  //create the scaling doubles.
  Double_t METAScale 		= 0.0;
  Double_t META2scale		= 0.0;
  Double_t METriggerScale 	= 0.0;
  //Normal estimation: 
  if(metaindependent)GetScalingFactors(BinDir,"divided","iteration1",&METAScale,&META2scale,&METriggerScale,true,false,1);
  else GetScalingFactors(BinDir,"divided","iteration1",&METAScale,&METriggerScale,true,false,1);
  if(!metaindependent)META2scale = METAScale;
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,META2scale,METriggerScale,"iteration1");
  Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12",BinDir,"divided","iteration1",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  
  //Alternative estimation:
  METAScale 		= 0.0;
  META2scale		= 0.0;
  METriggerScale 	= 0.0;
  GetScalingFactors(BinDir,"divided","Etasubstracted",&METAScale,&META2scale,&METriggerScale,true,false,2);
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,META2scale,METriggerScale,"Etasubstracted");
  Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12",BinDir,"divided","Etasubstracted",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);

}
void CorrectPbPb(BinDirs* BinDir,bool metaindependent)
{
  //If the bin is empty, skip the entire bin:
  if(!dynamic_cast<TH1D*>(BinDir->Samediv()->Get("number_of_triggers"))){cout << "skipping empty bin "<<BinDir->path()<<endl;return;}
  //create the scaling doubles.
  Double_t METAScale 		= 0.0;
  Double_t META2scale		= 0.0;
  Double_t METriggerScale 	= 0.0;
  //Normal estimation: 
  GetScalingFactors(BinDir,"divided","iteration1",&METAScale,&META2scale,&METriggerScale,false,true,1);
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,META2scale,METriggerScale,"iteration1");
  Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);  
  Correct("DPhi_12A_DEta_12",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);  
  Correct("DPhi_12A_DEta_12_SameSide",BinDir,"divided","iteration1",METAScale,META2scale,METriggerScale);  

  //Alternative estimation:
  METAScale 		= 0.0;
  META2scale		= 0.0;
  METriggerScale 	= 0.0;
  GetScalingFactors(BinDir,"divided","Etasubstracted",&METAScale,&META2scale,&METriggerScale,false,true,2);
  //Create a canvas that tells what is substracted.
  savedircontent(BinDir,METAScale,META2scale,METriggerScale,"Etasubstracted");
  Correct("DPhi_1_DPHI_2","DPhi_1_DPHI_2_near","DPhi_1_DPHI_2_mid","DPhi_1_DPHI_2_far",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12",BinDir,"divided","Etasubstracted",METAScale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_2PI",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_DPHI12_LESS_4PI",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);
  Correct("DPhi_1_DEta_12_SameSide",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);
  Correct("DPhi_12A_DEta_12",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);  
  Correct("DPhi_12A_DEta_12_SameSide",BinDir,"divided","Etasubstracted",METAScale,META2scale,METriggerScale);  
  
}   

void fitwith(TDirectory * dir, const char* type, TH2D * histo,Double_t etalimit){
  
  
  /*deactivated for now*/
//   //Fit the hist with the method specified in type.
//   //Types: fitfunc/rebin:
//   //Fitfunc: fgauspol1 , fgauspol2 , fgauspol0
//   //rebin: adds together n bins for /n.
//   TObjArray * types = TString(type).Tokenize("/");
//   int rebin = dynamic_cast<TObjString*>(types->At(1))->GetString().Atoi();
//   dir->cd();
//   TH2D* hist = dynamic_cast<TH2D*>(histo->Clone("DPhi_1_DEta_12_SameSidecc"));
//   if(rebin>1){hist->RebinX(3);hist->RebinY(2);}
//   hist->Write("DPhi_1_DEta_12_SameSide");
// 
//   TH1D* deta12ss = hist->ProjectionX("DEta12");
//   TH1D* deta12ssdraw = hist->ProjectionX("DEta12d");
//   TH1D* deta12ssbinc = hist->ProjectionX("DEta12dbc");
// 
//   deta12ss->Reset();deta12ssdraw->Reset();deta12ssbinc->Reset();
//   
//   TF1* fitsig;
//   TF1* fitbg;
//   TF1* rembg;
//   
//   TDirectory* typedir;
//   TDirectory* bindir;
//   
//   TString title = TString("#Delta#eta_{12} distribution in bin");
//   int color;
// //   const char* lable;
//   TH1D* yield = hist->ProjectionY("dphiyield");
//   yield->Reset();
//   yield->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
//   TH1D* yieldbc = hist->ProjectionY("dphiyieldbc");
//   yieldbc->Reset();
//   yield->SetTitle("Yield from bin counting vs #Delta#Phi_1");
//   yieldbc->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
//   
//   TH1D* background = dynamic_cast<TH1D*>(yield->Clone("dphibackground"));
//   background->SetTitle("Height of the background as a function of #Delta#Phi#.");
//   background->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
//   TH1D* width = dynamic_cast<TH1D*>(yield->Clone("dphiwidth")); 
//   width->SetTitle("Width of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
//   width->GetYaxis()->SetTitle("width of the peak (rad)");
//   TH1D* hRMS = dynamic_cast<TH1D*>(yield->Clone("dphiRMS")); 
//   width->SetTitle("RMS of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
//   width->GetYaxis()->SetTitle("RMS of the peak (rad)");
// 
//   TH1D* peakpos = dynamic_cast<TH1D*>(yield->Clone("dphipos")); 
//   peakpos->SetTitle("Position of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
//   peakpos->SetYTitle("Position of the peak (rad)");
//   TH1D* chisq = dynamic_cast<TH1D*>(yield->Clone("dphichisq"));
//   chisq->SetTitle("#Chi^2/NDF of the fit as a function of #Delta#Phi");
//   chisq->SetYTitle("#Chi^2/NDF ");
//   TH1D* prob = dynamic_cast<TH1D*>(yield->Clone("dphiprob"));
//   prob->SetTitle("Probability of the fit as a function of #Delta#Phi");
//   prob->SetYTitle("Probability ");
// 
//   TH1D * etasigrange = dynamic_cast<TH1D*>(yield->Clone("detasigrange"));
//   etasigrange->SetTitle("etarange bin by bin");
//   etasigrange->SetYTitle("eta range ");
//   
//   Double_t AvWidth = 1.0;
//   if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol0")==0){
//     typedir = resultsdirectory(dir,"GP0");    
//     bindir = resultsdirectory(typedir,"bins");
//     //flat background and a gaussian:
//     color = 2;
// //     lable = "flat bg + gaussian";
//     yield->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a flat background and a gaussian.");
//     yield->SetLineColor(color);
//     background->SetLineColor(color);
//     width->SetLineColor(color);
//     peakpos->SetLineColor(color);  
//     chisq->SetLineColor(color);  
//     prob->SetLineColor(color);    
//     //Initialize the fit function:
//     fitsig = new TF1("fgs",CGausPol0,-gkEtaFitRange,gkEtaFitRange,4) ;
//     fitsig->SetParNames("peakhight", "peakpos", "peakwidth", "B") ;
//     fitsig->SetLineColor(color);
//     fitsig->SetLineWidth(1);
//     fitsig->SetParameters(0.0,0.0,0.1,0.0);
//     fitsig->SetParLimits(0,0.0,0.3);
//     fitsig->SetParLimits(1,-0.1,0.1);
//     fitsig->SetParLimits(2,0.01,0.8);  
//     //Initialize the fit function for BKG:
//     fitbg = new TF1("fgs",CBKGPol0,-gkEtaFitRange,gkEtaFitRange,2) ;
//     fitbg->SetParNames("B", "peakpos" ) ;
//     fitbg->SetLineColor(color);
//     fitbg->SetLineWidth(1);
//     fitbg->SetParameters(0.0,0.0);
// //     fitbg->SetParLimits(0,0.0,1.0E10);
//     fitbg->SetParLimits(1,-0.1,0.1);    
//     //Initialize the fit function for BKG:
//     rembg = new TF1("fgs",CBKGPol0p,-gkEtaFitRange,gkEtaFitRange,1) ;
//     rembg->SetParNames("B", "peakpos" ) ;
//     rembg->SetLineColor(color);
//     rembg->SetLineWidth(1);
//     rembg->SetParameters(0.0,0.0);    
//     
//     //find the width at 0:
//     int bin1 = yield->GetXaxis()->FindBin(-0.001);
//     int bin2 = yield->GetXaxis()->FindBin(0.001);
//     
//     fillwithbinnr(hist,deta12ss,bin1);
//     fitsig->FixParameter(1,0.0);
//     TFitResultPtr fitresultwidth1 = deta12ss->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
//     Double_t width1 = fitsig->GetParameter(2);
//     fillwithbinnr(hist,deta12ss,bin2);
//     TFitResultPtr fitresultwidth2 = deta12ss->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
//     Double_t width2 = fitsig->GetParameter(2);  
//     deta12ss->Reset();
//     fitsig->ReleaseParameter(1);
//     
//     fitsig->FixParameter(2,width1+width2/2.0);
//     AvWidth = 3.0*fitsig->GetParameter(2);
// //     cout << AvWidth<<endl;
//     if(3.0*fitsig->GetParameter(2)<etalimit){AvWidth = 3.0*fitsig->GetParameter(2);}//set the range to 3*/sigma
//     else {AvWidth = etalimit;}
// //     cout << AvWidth<<endl;
//     fitsig->ReleaseParameter(2);
//     fitsig->SetParLimits(2,0.01,0.8);  
// 
//   }
//   /*
// //   if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol1")==0){
// //     typedir = resultsdirectory(dir,"GP1");
// //     bindir = resultsdirectory(typedir,"bins");
// //     //pol1 background and a gaussian:
// //     color = 3;
// // //     lable = "pol1 bg + gaussian";
// //     yield->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a pol1 background and a gaussian.");
// //     yield->SetLineColor(color);
// //     background->SetLineColor(color);
// //     width->SetLineColor(color);
// //     peakpos->SetLineColor(color);  
// //     chisq->SetLineColor(color);  
// //     prob->SetLineColor(color);
// //     //Initialize the fit function:
// //     fitsig = new TF1("p1gs",CGausAPol1,-gkEtaFitRange,gkEtaFitRange,5) ;
// //     fitsig->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C") ;
// //     fitsig->SetLineColor(color);
// //     fitsig->SetLineWidth(1);
// //     fitsig->SetParameters(0.0,0.0,0.1,0.0,0.0);
// //     fitsig->SetParLimits(0,0.0,0.3);
// //     fitsig->SetParLimits(1,-0.1,0.1);
// //     fitsig->SetParLimits(2,0.1,0.8);  
// //     fitsig->SetParLimits(4,0.0,1.0);
// //     //Initialize the fit function for BKG:
// //     fitbg = new TF1("p1gs",CBKGAPol1,-gkEtaFitRange,gkEtaFitRange,3) ;
// //     fitbg->SetParNames("B","peakpos", "C") ;
// //     fitbg->SetLineColor(color);
// //     fitbg->SetLineWidth(1);
// //     fitbg->SetParameters(0.0,0.0,0.0);
// // //     fitbg->SetParLimits(0,0.0,1.0E10);
// //     fitbg->SetParLimits(1,-0.1,0.1);
// //     fitbg->SetParLimits(2,0.0,1.0);
// //     //Initialize the fit function for BKG:
// //     rembg = new TF1("p1gs",CBKGAPol1p,-gkEtaFitRange,gkEtaFitRange,3) ;
// //     rembg->SetParNames("B","peakpos", "C") ;
// //     rembg->SetLineColor(color);
// //     rembg->SetLineWidth(1);
// //     rembg->SetParameters(0.0,0.0,0.0);
// // //     fitbg->SetParLimits(0,0.0,1.0E10);
// //     rembg->SetParLimits(1,-0.1,0.1);
// //     rembg->SetParLimits(2,0.0,1.0);
// //   }
// //   if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol2")==0){
// //     typedir = resultsdirectory(dir,"GP2");  
// //     bindir  = resultsdirectory(typedir,"bins");    
// //     //pol2 background and a gaussian:
// //     color = 4;
// // //     lable = "pol2 bg + gaussian";
// //     yield->SetTitle("Number of pairs per trigger as a function of #Delta#Phi extracted using a pol2 background and a gaussian.");
// //     yield->SetLineColor(color);
// //     background->SetLineColor(color);
// //     width->SetLineColor(color);
// //     peakpos->SetLineColor(color);  
// //     chisq->SetLineColor(color);  
// //     prob->SetLineColor(color);
// //     //Initialize the fit function:
// //     fitsig= new TF1("p2gs",CGausAPol2,-gkEtaFitRange,gkEtaFitRange,6) ;
// //     fitsig->SetParNames("peakhight", "peakpos", "peakwidth", "B", "C", "D") ;
// //     fitsig->SetLineColor(color);
// //     fitsig->SetLineWidth(1);
// //     fitsig->SetParameters(0.0,0.0,0.1,0.0,0.0,0.0);
// //     fitsig->SetParLimits(0,0.0,0.3);
// //     fitsig->SetParLimits(1,-0.1,0.1);
// //     fitsig->SetParLimits(2,0.1,0.8);  
// //     fitsig->SetParLimits(4,0.0,1.0);
// //     //Initialize the fit function:
// //     fitbg= new TF1("p2gs",CBKGAPol2,-gkEtaFitRange,gkEtaFitRange,4) ;
// //     fitbg->SetParNames("B", "peakpos",  "C", "D") ;
// //     fitbg->SetLineColor(color);
// //     fitbg->SetLineWidth(1);
// //     fitbg->SetParameters(0.0,0.0,0.0,0.0);
// // //     fitbg->SetParLimits(0,0.0,1.0E10);
// //     fitbg->SetParLimits(1,-0.1,0.1);
// //     fitbg->SetParLimits(2,0.0,1.0);
// //     //Initialize the fit function:
// //     rembg= new TF1("p2gs",CBKGAPol2p,-gkEtaFitRange,gkEtaFitRange,4) ;
// //     rembg->SetParNames("B", "peakpos",  "C", "D") ;
// //     rembg->SetLineColor(color);
// //     rembg->SetLineWidth(1);
// //     rembg->SetParameters(0.0,0.0,0.0,0.0);
// // //     fitbg->SetParLimits(0,0.0,1.0E10);
// //     rembg->SetParLimits(1,-0.1,0.1);
// //     rembg->SetParLimits(2,0.0,1.0);    
// //   }
// */
//   
//   double deltaphi = (yield->GetXaxis()->GetBinCenter(2)-yield->GetXaxis()->GetBinCenter(1));//Width of the bin in phi
//   bindir->cd();
//   for(int dphi=1;dphi<=hist->GetNbinsY();dphi++){
//     gEtasigRange = 1.0;//AvWidth;
// //     gEtasigRange = 0.6;
//     deta12ss->SetTitle(Form("%s %4.2f #pi < #Delta#Phi < %4.2f #pi",title.Data(),hist->GetYaxis()->GetBinLowEdge(dphi)/TMath::Pi(),hist->GetYaxis()->GetBinUpEdge(dphi)/TMath::Pi()));
//     fillwithbinnr(hist,deta12ss,dphi);
//     deta12ss->SetStats(false);
// //     fitsig->FixParameter(1,0.0);
//     TFitResultPtr fitresult = deta12ss->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
//     if(int(fitresult)!=4000)//if "error", try again with more:
//       fitresult = deta12ss->Fit(fitsig,"MSQ","",-gkEtaFitRange,gkEtaFitRange);
//     deta12ss->Write(Form("%sflatg_%i_first",deta12ss->GetName(),dphi));
//     fitbg->FixParameter(1,0.0);
// //     if(3.0*fitsig->GetParameter(2)<1.0) cout << fitsig->GetParameter(2);
//     
//     
// //     if(3.0*fitsig->GetParameter(2)<etalimit)gEtasigRange = 3.0*fitsig->GetParameter(2);//set the range to 3*/sigma
// //     else gEtasigRange = etalimit;
// //     if(fitsig->GetParameter(0)<1.0E-6) gEtasigRange = 0.3;//if there is no peak, most of it.
// 
//     TFitResultPtr bkgresult = deta12ss->Fit(fitbg,"SQ","",-gkEtaFitRange,gkEtaFitRange);
//     if(int(bkgresult)!=4000)//if "error", try again with more:
//       bkgresult = deta12ss->Fit(fitbg,"MSQ","",-gkEtaFitRange,gkEtaFitRange);    
//     deta12ss->Write(Form("%sflatg_%i_bkgonly",deta12ss->GetName(),dphi));
//     fitsig->FixParameter(3,0.0);//fitbg->GetParameter(0));
//     rembg->SetParameter(0,fitbg->GetParameter(0));
// /*
//     if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol1")==0){fitsig->FixParameter(4,fitbg->GetParameter(2));rembg->SetParameter(2,fitbg->GetParameter(2));}
//     if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol2")==0){fitsig->FixParameter(5,fitbg->GetParameter(3));rembg->SetParameter(2,fitbg->GetParameter(3));}*/
// //     
//     fillwithbinnr(hist,deta12ssbinc,dphi);
//     Double_t ErrBg = 0.0;
//     if(fitbg->GetNDF()!=0)ErrBg= fitbg->GetParError(0)*fitbg->GetChisquare()/fitbg->GetNDF();
//     removeconstant(deta12ssbinc,-1.0*fitbg->GetParameter(0),ErrBg);
//     
//     Double_t binerr;
//     etasigrange->SetBinContent(dphi,gEtasigRange);
//     Double_t binc = deta12ssbinc->IntegralAndError(deta12ssbinc->FindBin(-gEtasigRange),deta12ssbinc->FindBin(gEtasigRange),binerr);
//     Double_t rmsv = deta12ssbinc->GetRMS();
//     yieldbc->SetBinContent(dphi,binc/deltaphi);
//     yieldbc->SetBinError(dphi,binerr/deltaphi);
//     hRMS->SetBinContent(dphi,rmsv/deltaphi);
//     deta12ssbinc->Write(Form("%sminbgbin_%i",deta12ss->GetName(),dphi));
//     
//     TFitResultPtr fitresult2 = deta12ssbinc->Fit(fitsig,"SQ","",-gkEtaFitRange,gkEtaFitRange);
//     if(int(fitresult2)!=4000)//if "error", try again with more:
//       fitresult2 = deta12ssbinc->Fit(fitsig,"MSQ","",-gkEtaFitRange,gkEtaFitRange);    
//     deta12ssbinc->Write(Form("%sflatg_%i_full",deta12ss->GetName(),dphi));
//     fitsig->ReleaseParameter(1);fitsig->ReleaseParameter(3);
// //     if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol1")==0){fitsig->ReleaseParameter(4);}
// //     if(dynamic_cast<TObjString*>(types->At(0))->GetString().CompareTo("fgauspol2")==0){fitsig->ReleaseParameter(5);}    
// 
//     TCanvas * bincanvas = new TCanvas(Form("%sCanvas_%i",deta12ss->GetName(),dphi));
//     deta12ssdraw->SetTitle(Form("%s %4.2f #pi < #Delta#Phi < %4.2f #pi",title.Data(),hist->GetYaxis()->GetBinLowEdge(dphi)/TMath::Pi(),hist->GetYaxis()->GetBinUpEdge(dphi)/TMath::Pi()));
//     deta12ssdraw->Reset();
//     fillwithbinnr(hist,deta12ssdraw,dphi);
//     deta12ssdraw->SetStats(false);    
//     removeconstant(deta12ssdraw,-1.0*fitbg->GetParameter(0),ErrBg);
//     deta12ssdraw->Draw("ESAME");
//     fitsig->Draw("LSAME");
//     bincanvas->Update();
//     bincanvas->Write();
//     delete bincanvas;
// 
//     yield->SetBinContent(dphi,fitsig->GetParameter(0)/deltaphi);
//     yield->SetBinError(dphi,fitsig->GetParError(0)/deltaphi);
//     width->SetBinContent(dphi,fitsig->GetParameter(2)/deltaphi);
//     width->SetBinError(dphi,fitsig->GetParError(2)/deltaphi);
//     peakpos->SetBinContent(dphi,fitsig->GetParameter(1)/deltaphi);
//     peakpos->SetBinError(dphi,fitsig->GetParError(1)/deltaphi);
//     background->SetBinContent(dphi,fitbg->GetParameter(0)/deltaphi);
//     background->SetBinError(dphi,fitbg->GetParError(0)/deltaphi);
//     if(fitsig->GetNDF()>=1.0)chisq->SetBinContent(dphi,fitsig->GetChisquare()/fitsig->GetNDF());
//     prob->SetBinContent(dphi,fitsig->GetProb());
//     deta12ss->Reset();
//     
//     
//     fitresult.~TFitResultPtr();fitresult2.~TFitResultPtr();bkgresult.~TFitResultPtr();
//   }
//   typedir->cd();
//   yieldbc->Write();
//   hRMS->Write();
//   yield->Write();
//   background->Write();
//   width->Write();
//   peakpos->Write();
//   chisq->Write();
//   prob->Write();
//   etasigrange->Write();
//   
//   
//   delete types;delete yield; delete background;delete width; delete peakpos; delete chisq; delete prob;delete etasigrange;
}
void savedircontent(TDirectory* here, double yield_int_near,double yield_int_near_er,double yield_int_away,double yield_int_away_er,double width_near, double width_near_er,double width_away , double width_away_er,TF1 * fitfunc,TH1 * widthhist,TH1 * fittedwidthhist){
  here->cd();
  TPaveText * whatinthisbin;
  TParameter<double> widththetanear = TParameter<double>("width_thetanear",width_near);
  TParameter<double> widththetanear_error = TParameter<double>("width_thetanear_error",width_near_er);
  TParameter<double> widththetaaway = TParameter<double>("width_thetaaway",width_away);
  TParameter<double> widththetaaway_error = TParameter<double>("width_thetaaway_error",width_away_er);
  
  TParameter<double> widththetanear_f = TParameter<double>("width_thetanear_fitted",fitfunc->GetParameter(2));
  TParameter<double> widththetanear_f_error = TParameter<double>("width_thetanear_fitted_error",fitfunc->GetParError(2));
  TParameter<double> widththetaaway_f = TParameter<double>("width_thetaaway_fitted",fitfunc->GetParameter(5));
  TParameter<double> widththetaaway_f_error = TParameter<double>("width_thetaaway_fitted_error",fitfunc->GetParError(5));
  //use phi width to find the average eta width on near and away side:
  Double_t errorwidtheta = 0;
  double nearwidtheta = widthhist->IntegralAndError(widthhist->GetXaxis()->FindBin(0.0-width_near),widthhist->GetXaxis()->FindBin(0.0+width_near),errorwidtheta);
  nearwidtheta = nearwidtheta/(widthhist->GetXaxis()->FindBin(0.0+width_near)-widthhist->GetXaxis()->FindBin(0.0-width_near)+1);
  errorwidtheta = errorwidtheta/(widthhist->GetXaxis()->FindBin(0.0+width_near)-widthhist->GetXaxis()->FindBin(0.0-width_near)+1);
  TParameter<double> widthetanear = TParameter<double>("width_etanear",nearwidtheta);
  TParameter<double> widthetanear_error = TParameter<double>("width_etanear_error",errorwidtheta);
  
  double awaywidtheta = widthhist->IntegralAndError(widthhist->GetXaxis()->FindBin(TMath::Pi()-width_away),widthhist->GetXaxis()->FindBin(TMath::Pi()+width_away),errorwidtheta);  
  awaywidtheta = awaywidtheta/(widthhist->GetXaxis()->FindBin(TMath::Pi()+width_away)-widthhist->GetXaxis()->FindBin(TMath::Pi()-width_away)+1);
  errorwidtheta = errorwidtheta/(widthhist->GetXaxis()->FindBin(TMath::Pi()+width_away)-widthhist->GetXaxis()->FindBin(TMath::Pi()-width_away)+1);
  TParameter<double> widthetaaway = TParameter<double>("width_etaaway",awaywidtheta);
  TParameter<double> widthetaaway_error = TParameter<double>("width_etaaway_error",errorwidtheta);
  
  Double_t errorwidtheta_f = 0;  
  double nearwidtheta_f = fittedwidthhist->IntegralAndError(fittedwidthhist->GetXaxis()->FindBin(0.0-width_near),fittedwidthhist->GetXaxis()->FindBin(0.0+width_near),errorwidtheta_f);
  nearwidtheta_f = nearwidtheta_f/(fittedwidthhist->GetXaxis()->FindBin(0.0+width_near)-fittedwidthhist->GetXaxis()->FindBin(0.0-width_near)+1);
  errorwidtheta_f = errorwidtheta_f/(fittedwidthhist->GetXaxis()->FindBin(0.0+width_near)-fittedwidthhist->GetXaxis()->FindBin(0.0-width_near)+1);
  TParameter<double> widthetanear_f = TParameter<double>("width_etanear_fitted",nearwidtheta_f);
  TParameter<double> widthetanear_f_error = TParameter<double>("width_etanear_fitted_error",errorwidtheta_f);
  
  double awaywidtheta_f = fittedwidthhist->IntegralAndError(fittedwidthhist->GetXaxis()->FindBin(TMath::Pi()-width_away),fittedwidthhist->GetXaxis()->FindBin(TMath::Pi()+width_away),errorwidtheta_f);  
  awaywidtheta_f = awaywidtheta_f/(fittedwidthhist->GetXaxis()->FindBin(TMath::Pi()+width_away)-fittedwidthhist->GetXaxis()->FindBin(TMath::Pi()-width_away)+1);
  errorwidtheta_f = errorwidtheta_f/(fittedwidthhist->GetXaxis()->FindBin(TMath::Pi()+width_away)-fittedwidthhist->GetXaxis()->FindBin(TMath::Pi()-width_away)+1);
  TParameter<double> widthetaaway_f = TParameter<double>("width_etaaway_fitted",awaywidtheta_f);
  TParameter<double> widthetaaway_f_error = TParameter<double>("width_etaaway_fitted_error",errorwidtheta_f);
  
  
  
  
  //extract integrated and fitted yields:
  TParameter<double> yieldnearfit = TParameter<double>("yield_near_fitted",fitfunc->GetParameter(0));
  TParameter<double> yieldnearfit_error = TParameter<double>("yield_near_fitted_error",fitfunc->GetParError(0));
  TParameter<double> yieldawayfit = TParameter<double>("yield_away_fitted",fitfunc->GetParameter(3));
  TParameter<double> yieldawayfit_error = TParameter<double>("yield_away_fitted_error",fitfunc->GetParError(3));
  TParameter<double> yieldnear = TParameter<double>("yield_near",yield_int_near);
  TParameter<double> yieldnear_error = TParameter<double>("yield_near_error",yield_int_near_er);
  TParameter<double> yieldaway = TParameter<double>("yield_away",yield_int_away);
  TParameter<double> yieldaway_error = TParameter<double>("yield_away_error",yield_int_away_er);

  TCanvas * whatinthisbincanvas = new TCanvas("Peaks");
  whatinthisbincanvas->cd();
  whatinthisbin = new TPaveText(.05,.1,.95,.8);
  whatinthisbin->AddText("The nearside peak has:");
  whatinthisbin->AddText(Form("A width of %4.2f #pm %4.2f in #Delta#theta_{1}, a width of %4.2f #pm %4.2f in #Delta#eta_{12},",widththetanear.GetVal(),widththetanear_error.GetVal(),widthetanear.GetVal(),widthetanear_error.GetVal()));
  whatinthisbin->AddText(Form("and an integrated yield of %4.4f #pm %4.4f",yield_int_near, yield_int_near_er));
  whatinthisbin->AddText("The awayside peak has:");
  whatinthisbin->AddText(Form("A width of %4.2f #pm %4.2f in #Delta#theta_{1}, a width of %4.2f #pm %4.2f in #Delta#eta_{12},",widththetaaway.GetVal(),widththetaaway_error.GetVal(),widthetaaway.GetVal(),widthetaaway_error.GetVal()));
  whatinthisbin->AddText(Form("and an integrated yield of %4.4f #pm %4.4f",yield_int_away, yield_int_away_er));
  whatinthisbin->Draw();
  
  
  whatinthisbincanvas->Write();
  widththetanear.Write();
  widththetanear_error.Write();
  widththetaaway.Write();
  widththetaaway_error.Write();
  widththetanear_f.Write();
  widththetanear_f_error.Write();
  widththetaaway_f.Write();
  widththetaaway_f_error.Write();
  
  widthetanear.Write();
  widthetanear_error.Write();
  widthetaaway.Write();
  widthetaaway_error.Write();
  widthetanear_f.Write();
  widthetanear_f_error.Write();
  widthetaaway_f.Write();
  widthetaaway_f_error.Write();

  yieldnearfit.Write();
  yieldnearfit_error.Write();
  yieldawayfit.Write();
  yieldawayfit_error.Write();
  yieldnear.Write();
  yieldnear_error.Write();
  yieldaway.Write();
  yieldaway_error.Write();

  delete whatinthisbincanvas;
  delete whatinthisbin;
}

void savedircontent(TDirectory* here, 
		    double yield_int_near,double yield_int_near_er,double yield_int_away,double yield_int_away_er, 
		    double yield_f_near,double yield_f_near_er,double yield_f_away,double yield_f_away_er, 
		    double widtheta_near, double widtheta_near_err,double widtheta_away,double widtheta_away_err,
		    double widtheta_near_f, double widtheta_near_f_err,double widtheta_away_f,double widtheta_away_f_err){
  here->cd();
  TPaveText * whatinthisbin;
  TParameter<double> widthetanear = TParameter<double>("width_etanear",widtheta_near);
  TParameter<double> widthetanear_error = TParameter<double>("width_etanear_error",widtheta_near_err);
  TParameter<double> widthetaaway = TParameter<double>("width_etaaway",widtheta_away);
  TParameter<double> widthetaaway_error = TParameter<double>("width_etaaway_error",widtheta_away_err);
  
  TParameter<double> widthetanear_f = TParameter<double>("width_etanear_fitted",widtheta_near_f);
  TParameter<double> widthetanear_f_error = TParameter<double>("width_etanear_fitted_error",widtheta_near_f_err);
  TParameter<double> widthetaaway_f = TParameter<double>("width_etaaway_fitted",widtheta_away_f);
  TParameter<double> widthetaaway_f_error = TParameter<double>("width_etaaway_fitted_error",widtheta_away_f_err);
  //use phi width to find the average eta width on near and away side:
    
  //extract integrated and fitted yields:
  TParameter<double> yieldnearfit = TParameter<double>("yield_near_fitted",yield_f_near);
  TParameter<double> yieldnearfit_error = TParameter<double>("yield_near_fitted_error",yield_f_near_er);
  TParameter<double> yieldawayfit = TParameter<double>("yield_away_fitted",yield_f_away);
  TParameter<double> yieldawayfit_error = TParameter<double>("yield_away_fitted_error",yield_f_away_er);
  TParameter<double> yieldnear = TParameter<double>("yield_near",yield_int_near);
  TParameter<double> yieldnear_error = TParameter<double>("yield_near_error",yield_int_near_er);
  TParameter<double> yieldaway = TParameter<double>("yield_away",yield_int_away);
  TParameter<double> yieldaway_error = TParameter<double>("yield_away_error",yield_int_away_er);

  TCanvas * whatinthisbincanvas = new TCanvas("Peaks");
  whatinthisbincanvas->cd();
  whatinthisbin = new TPaveText(.05,.1,.95,.8);
  whatinthisbin->AddText("The nearside peak has:");
  whatinthisbin->AddText(Form("A width of %4.2f #pm %4.2f in #Delta#eta_{12}",widthetanear.GetVal(),widthetanear_error.GetVal()));
  whatinthisbin->AddText(Form("and an integrated yield of %4.4f #pm %4.4f",yield_int_near, yield_int_near_er));
  whatinthisbin->AddText("The awayside peak has:");
  whatinthisbin->AddText(Form("A width of %4.2f #pm %4.2f in #Delta#eta_{12}",widthetaaway.GetVal(),widthetaaway_error.GetVal()));
  whatinthisbin->AddText(Form("and an integrated yield of %4.4f #pm %4.4f",yield_int_away, yield_int_away_er));
  whatinthisbin->Draw();
  
  
  whatinthisbincanvas->Write();
  
  widthetanear.Write();
  widthetanear_error.Write();
  widthetaaway.Write();
  widthetaaway_error.Write();
  widthetanear_f.Write();
  widthetanear_f_error.Write();
  widthetaaway_f.Write();
  widthetaaway_f_error.Write();

  yieldnearfit.Write();
  yieldnearfit_error.Write();
  yieldawayfit.Write();
  yieldawayfit_error.Write();
  yieldnear.Write();
  yieldnear_error.Write();
  yieldaway.Write();
  yieldaway_error.Write();

  delete whatinthisbincanvas;
  delete whatinthisbin;
}


void simpleyield(TDirectory * dir, TH2D * histo,double sidebandlow , double sidebandhigh){
  // estimate background in sidebandlow<|deta|<sidebandhigh
  dir->cd();
  TH2D* hist = dynamic_cast<TH2D*>(histo->Clone("DPhi_1_DEta_12_SameSidecc"));
  hist->Write("DPhi_1_DEta_12_SameSide");

  TH1D* deta12ss = hist->ProjectionX("DEta12");
  TH1D* deta12ssdraw = hist->ProjectionX("DEta12d");
  TH1D* deta12ssbinc = hist->ProjectionX("DEta12dbc");

  deta12ss->Reset();deta12ssdraw->Reset();deta12ssbinc->Reset();
  
  
  TDirectory* bindir = resultsdirectory(dir,"bins");
  
  TString title = TString("#Delta#eta_{12} distribution in bin");
  int color=2;
  TH1D* yield = hist->ProjectionY("dphiyield");
  yield->Reset();
  yield->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  yield->SetTitle("Yield from fitting the peak vs #Delta#Phi_1");
  TH1D* yieldbc = hist->ProjectionY("dphiyieldbc");
  yieldbc->Reset();
  yieldbc->SetTitle("Yield from bin counting vs #Delta#Phi_1");
  yieldbc->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  
  TH1D* background = dynamic_cast<TH1D*>(yield->Clone("dphibackground"));
  background->SetTitle("Height of the background as a function of #Delta#Phi#.");
  background->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#Phi_{1}} (rad)^{-1}");
  TH1D* width = dynamic_cast<TH1D*>(yield->Clone("dphiwidth")); 
  width->SetTitle("Width of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  width->GetYaxis()->SetTitle("width of the peak [rad]");
  TH1D* dphiRMS = dynamic_cast<TH1D*>(yield->Clone("dphiRMS")); 
  dphiRMS->SetTitle("RMS of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  dphiRMS->GetYaxis()->SetTitle("RMS of the peak [rad]");
  
  TH1D* peakpos = dynamic_cast<TH1D*>(yield->Clone("dphipos")); 
  peakpos->SetTitle("Position of the peak in #Delta#eta_{12} as a function of #Delta#Phi");
  peakpos->SetYTitle("Position of the peak (rad)");

  TH1D* backgroundhypothesis = dynamic_cast<TH1D*>(yield->Clone("chi2_background")); 
  backgroundhypothesis->SetTitle("chi2 for the background region");
  backgroundhypothesis->SetYTitle("chi2");
  
  double deltaphi = (yield->GetXaxis()->GetBinCenter(2)-yield->GetXaxis()->GetBinCenter(1));//Width of the bin in phi
  double deltabins = 2.0*(deta12ss->FindBin(sidebandhigh)-deta12ss->FindBin(sidebandlow) + 1);//how many bins between 
  Double_t error1;  Double_t error2;
  Double_t background1;Double_t background2;
  Double_t erroryield;
  Double_t yieldc;
  
  //function to characterize the peak:
  TF1* fitsig = new TF1("fgs",CGaus,-sidebandlow,sidebandlow,3);
  fitsig->SetParNames("peakhight", "peakpos", "peakwidth") ;
  fitsig->SetLineColor(color);
  fitsig->SetLineWidth(1);
  fitsig->SetParameters(0.0,0.0,0.1);
  fitsig->SetParLimits(0,0.0,1000);
  fitsig->SetParLimits(1,-0.2,0.2);
  fitsig->SetParLimits(2,0.1,0.5);  
  fitsig->SetNpx(500);
  
  //function to characterize the peaks  in dphi:
  TF1* fitsig2 = new TF1("fgs",C2Gaus,-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0,6);
  fitsig2->SetParNames("peakhight_near", "peakpos_near", "peakwidth_near","peakhight_away", "peakpos_away", "peakwidth_away") ;
  fitsig2->SetLineColor(color);
  fitsig2->SetLineWidth(1);
  fitsig2->SetParameters(0.0,0.0,0.1,0.0,TMath::Pi(),0.2);
  fitsig2->SetParLimits(0,0.0,1000);
  fitsig2->SetParLimits(1,-0.2,0.2);
  fitsig2->SetParLimits(2,0.05,1.0);  
  fitsig2->SetParLimits(3,0.0,1000);
  fitsig2->SetParLimits(4,TMath::Pi()-0.2,TMath::Pi()+0.2);
  fitsig2->SetParLimits(5,0.1,1.5);    
  fitsig2->SetNpx(500);
  
  bindir->cd();
  for(int dphi=1;dphi<=hist->GetNbinsY();dphi++){
    
    deta12ss->SetTitle(Form("%s %4.2f #pi < #Delta#Phi < %4.2f #pi",title.Data(),hist->GetYaxis()->GetBinLowEdge(dphi)/TMath::Pi(),hist->GetYaxis()->GetBinUpEdge(dphi)/TMath::Pi()));
    fillwithbinnr(hist,deta12ss,dphi,true);
    deta12ss->SetStats(false);
    //save the original histogram
    deta12ss->Write(Form("%sflatg_%i_first",deta12ss->GetName(),dphi));
    //Estimate the backlground from the sidebands, weighting by the errors:
    double backgroundestimate = 0.0;
    double weight = 0.0;
    for(int sidebin = deta12ss->FindBin(-sidebandhigh);sidebin<=deta12ss->FindBin(-sidebandlow);sidebin++){
      if(deta12ss->GetBinError(sidebin)>1.0E-10){
	backgroundestimate += deta12ss->GetBinContent(sidebin)/(deta12ss->GetBinError(sidebin)*deta12ss->GetBinError(sidebin));
	weight += 1/(deta12ss->GetBinError(sidebin)*deta12ss->GetBinError(sidebin));
      }
    }
   for(int sidebin = deta12ss->FindBin(sidebandlow);sidebin<=deta12ss->FindBin(sidebandhigh);sidebin++){
      if(deta12ss->GetBinError(sidebin)>1.0E-10){
	backgroundestimate += deta12ss->GetBinContent(sidebin)/(deta12ss->GetBinError(sidebin)*deta12ss->GetBinError(sidebin));
	weight += 1/(deta12ss->GetBinError(sidebin)*deta12ss->GetBinError(sidebin));
      }
    }

    Double_t estimatedbackground = 0.0;
    Double_t errorestimatedbackground = 0.0;
    if(weight > 1.0E-10){
      estimatedbackground = backgroundestimate/weight;
      errorestimatedbackground = TMath::Sqrt(1.0/weight);
      background->SetBinContent(dphi,estimatedbackground);
      background->SetBinError(dphi,errorestimatedbackground);
    }
    //average deviation from flat:
    double nbins = 0.0;
    double chi = 0.0;
    for(int k = deta12ss->FindBin(-sidebandhigh); k<= deta12ss->FindBin(-sidebandlow);k++){
      if(TMath::Abs(deta12ss->GetBinError(k))>0.0000001)
	chi += (deta12ss->GetBinContent(k)-estimatedbackground)*(deta12ss->GetBinContent(k)-estimatedbackground)/(deta12ss->GetBinError(k)*deta12ss->GetBinError(k));
	nbins +=1.0;
    }
    for(int k = deta12ss->FindBin(sidebandlow); k<= deta12ss->FindBin(sidebandhigh);k++){
      if(TMath::Abs(deta12ss->GetBinError(k))>0.0000001)
	chi += (deta12ss->GetBinContent(k)-estimatedbackground)*(deta12ss->GetBinContent(k)-estimatedbackground)/(deta12ss->GetBinError(k)*deta12ss->GetBinError(k));
	nbins +=1.0;
    }
    if(nbins>1.5)backgroundhypothesis->SetBinContent(dphi,chi/(nbins-1));//ndof = nbins-1

    //Remove the background and write our result:
    removeconstant(deta12ss,-1.0*estimatedbackground,errorestimatedbackground);
    deta12ss->Write(Form("%sflatg_%i_backgroundsubstracted",deta12ss->GetName(),dphi));    
    
    deta12ss->GetXaxis()->SetRange(deta12ss->FindBin(-sidebandlow),deta12ss->FindBin(sidebandlow));
    dphiRMS->SetBinContent(dphi,deta12ss->GetRMS(1));
    dphiRMS->SetBinError(dphi,deta12ss->GetRMS(11));
    deta12ss->GetXaxis()->SetRange(1,-1);
    
    //Estimate yield by bincounting and write to histogram.
    yieldc = deta12ss->IntegralAndError(deta12ss->FindBin(-sidebandlow),deta12ss->FindBin(sidebandlow),erroryield);
    if(erroryield>1.0E-25){
      yieldbc->SetBinContent(dphi,yieldc);
      yieldbc->SetBinError(dphi,erroryield);
    }
    //Try to fit the background substracted histogram with a gaussian:
    deta12ss->Fit(fitsig,"SQ","",-sidebandlow,sidebandlow);
    deta12ss->Write(Form("%sflatg_%i_gaussianfit",deta12ss->GetName(),dphi));
//     cout << fitsig->GetChisquare()/fitsig->GetNDF()<< " " << fitsig->GetParameter(0)<< " " << fitsig->GetParError(0)<<endl;
    if(fitsig->GetParameter(0)>0.00001){
	if(fitsig->GetParError(0)/fitsig->GetParameter(0)<0.7&&fitsig->GetParError(2)/fitsig->GetParameter(2)<0.7){
	yield->SetBinContent(dphi,fitsig->GetParameter(0));
	yield->SetBinError(dphi,fitsig->GetParError(0));
	peakpos->SetBinContent(dphi,fitsig->GetParameter(1));
	peakpos->SetBinError(dphi,fitsig->GetParError(1));
	width->SetBinContent(dphi,fitsig->GetParameter(2));
	width->SetBinError(dphi,fitsig->GetParError(2));
      }
    }
    deta12ss->Reset();
    error1 = 0.0;
    error2 = 0.0;
    background1 = 0.0;
    background2 = 0.0;
    erroryield = 0.0;
    yieldc = 0.0;
  }
  //fit the dphi picture:
  yieldbc->Fit(fitsig2,"SQ","",-TMath::Pi()/2.0,3.0*TMath::Pi()/2.0);
  double yield_int_near_er = 0.0;
  double yield_int_near = yieldbc->IntegralAndError(yieldbc->FindBin(-TMath::Pi()/2.0),yieldbc->FindBin(TMath::Pi()/2.0),yield_int_near_er)*deltaphi;
  yield_int_near_er *= deltaphi;
  double yield_int_away_er = 0.0;
  double yield_int_away = yieldbc->IntegralAndError(yieldbc->FindBin(TMath::Pi()/2.0),yieldbc->FindBin(3.0*TMath::Pi()/2.0),yield_int_away_er)*deltaphi;
  yieldbc->GetXaxis()->SetRange(yieldbc->FindBin(-TMath::Pi()/2.0),yieldbc->FindBin(TMath::Pi()/2.0));
  double RMS_near = yieldbc->GetRMS(1);
  double RMS_near_er = yieldbc->GetRMS(11);
  yieldbc->GetXaxis()->SetRange(yieldbc->FindBin(TMath::Pi()/2.0),yieldbc->FindBin(3.0*TMath::Pi()/2.0));
  double RMS_away = yieldbc->GetRMS(1);
  double RMS_away_er = yieldbc->GetRMS(11);  
  yieldbc->GetXaxis()->SetRange(1,-1);
  yield_int_away_er *= deltaphi;  
  savedircontent(dir,  yield_int_near, yield_int_near_er, yield_int_away, yield_int_away_er,RMS_near,RMS_near_er,RMS_away,RMS_away_er,fitsig2, dphiRMS,width);
    
  dir->cd();
  yieldbc->Write();
  background->Write();
  yield->Write();
  peakpos->Write();
  width->Write();
  dphiRMS->Write();
  backgroundhypothesis->Write();

  
  delete yield; delete background;delete width; delete peakpos;
  delete deta12ss; delete deta12ssbinc; delete deta12ssdraw;
  delete backgroundhypothesis;
}

void singleyield(TDirectory * dir, TH2D * histo, double sidebandlow , double sidebandhigh){
  // estimate background in sidebandlow<|deta|<sidebandhigh in one slice per side.
  dir->cd();
  TH2D* hist = dynamic_cast<TH2D*>(histo->Clone("DPhi_1_DEta_12_SameSidecc"));
  hist->Write("DPhi_1_DEta_12_SameSide");

  //function to characterize the peak:
  TF1* fitsig1 = new TF1("fgs",CGaus,-sidebandlow,sidebandlow,3);
  fitsig1->SetParNames("peakhight", "peakpos", "peakwidth") ;
  fitsig1->SetLineColor(kRed);
  fitsig1->SetLineWidth(1);
  fitsig1->SetParameters(0.0,0.0,0.1);
  fitsig1->SetParLimits(0,0.0,1000);
  fitsig1->SetParLimits(1,-0.2,0.2);
  fitsig1->SetParLimits(2,0.001,1.0);  
  fitsig1->SetNpx(500);
  

  TH1D* deta12ns = hist->ProjectionX("DEta12_near");
  deta12ns->Reset();
  deta12ns->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#eta_{12}} ");
  TH1D* deta12as = hist->ProjectionX("DEta12_away");
  deta12as->Reset();
  deta12as->GetYaxis()->SetTitle("#frac{dN_{pairs}}{N_{trig}d#eta_{12}} ");
  deta12ns->Write("deta12_Near_empty");
  deta12as->Write("deta12_Away_empty");  
  double tmpnumber = 0.0;
  double tmperror = 0.0;
  double deltaphi = hist->GetYaxis()->GetBinCenter(5)-hist->GetYaxis()->GetBinCenter(4);
  //for each eta bin:
  for(int deta=1;deta<=deta12ns->GetNbinsX();deta++){
    //run over dphi bins and add all bins for near and away side;
    for(int dphi=1;dphi<=hist->GetNbinsY();dphi++){
      tmpnumber 	= hist->GetBinContent(deta,dphi)*deltaphi;
      tmperror 		= hist->GetBinError(deta,dphi)*deltaphi;
      if(hist->GetYaxis()->GetBinUpEdge(dphi)<=(TMath::Pi()/2.0)){
	//near side
	tmpnumber 	= tmpnumber + deta12ns->GetBinContent(deta);
	tmperror 	= tmperror*tmperror + deta12ns->GetBinError(deta)*deta12ns->GetBinError(deta);
	deta12ns->SetBinContent(deta,tmpnumber);
	deta12ns->SetBinError(deta,TMath::Sqrt(tmperror));
	
      }
      if(hist->GetYaxis()->GetBinUpEdge(dphi)>(TMath::Pi()/2.0)){
	//away side
	tmpnumber 	= tmpnumber + deta12as->GetBinContent(deta);
	tmperror 	= tmperror*tmperror + deta12as->GetBinError(deta)*deta12as->GetBinError(deta);
	deta12as->SetBinContent(deta,tmpnumber);
	deta12as->SetBinError(deta,TMath::Sqrt(tmperror));
      }
      tmperror = 0.0;
      tmpnumber = 0.0;
    }
  }
  
  deta12ns->Write("deta12_Near");
  deta12as->Write("deta12_Away");
  
  
  //integrate over sideband and substract:
  double bkgnear    	= 0.0;
  double errbkgnear 	= 0.0;
  
  double bkgaway    	= 0.0;
  double errbkgaway 	= 0.0;
  
  double nsigbinsn 	= 0.0;
  double nsigbinsa 	= 0.0;
  
  for(int dint=deta12ns->FindBin(-sidebandhigh);dint<=deta12ns->FindBin(-sidebandlow);dint++){
    bkgnear	+=deta12ns->GetBinContent(dint);
    errbkgnear 	= TMath::Sqrt(errbkgnear*errbkgnear + deta12ns->GetBinError(dint)*deta12ns->GetBinError(dint));
    bkgaway	+=deta12as->GetBinContent(dint);
    errbkgaway 	= TMath::Sqrt(errbkgaway*errbkgaway + deta12as->GetBinError(dint)*deta12as->GetBinError(dint));
    if(deta12ns->GetBinError(dint)>1.0E-25)nsigbinsn +=1;//only count
    if(deta12as->GetBinError(dint)>1.0E-25)nsigbinsa +=1;//significant bins
  }
  for(int dint = deta12ns->FindBin(sidebandlow);dint<=deta12ns->FindBin(sidebandhigh);dint++){
    bkgnear	+=deta12ns->GetBinContent(dint);
    errbkgnear 	= TMath::Sqrt(errbkgnear*errbkgnear + deta12ns->GetBinError(dint)*deta12ns->GetBinError(dint));
    bkgaway	+=deta12as->GetBinContent(dint);
    errbkgaway 	= TMath::Sqrt(errbkgaway*errbkgaway + deta12as->GetBinError(dint)*deta12as->GetBinError(dint));
    if(deta12ns->GetBinError(dint)>1.0E-25)nsigbinsn +=1;//only count
    if(deta12as->GetBinError(dint)>1.0E-25)nsigbinsa +=1;//significant bins
  }
  if(nsigbinsn>0.1){
    
    bkgnear = bkgnear/nsigbinsn;
    errbkgnear = errbkgnear/nsigbinsn;
    removeconstant(deta12ns,-bkgnear,errbkgnear);
  }


  if(nsigbinsa>0.1){
    bkgaway = bkgaway/nsigbinsa;
    errbkgaway = errbkgaway/nsigbinsa;
    removeconstant(deta12as,-bkgaway,errbkgaway);
  }


  
  deta12ns->Write("deta12_Near_bs");
  deta12as->Write("deta12_Away_bs");
  
  //extract width and hight by bincounting and fit:
  double nearsidemeanx 		= 0.0;
  double nearsideRMSx 		= 0.0;
  double nearsideRMSxerr	= 0.0;
  double nearsidewidth		= 0.0;
  double nearsidewidth_err	= 0.0;
  double nearsideyield 		= 0.0;
  double nearsideyielderr	= 0.0;
  double nearsideyield_f 	= 0.0;
  double nearsideyield_ferr	= 0.0;
  
  double awaysidemeanx 		= 0.0;
  double awaysideRMSx	 	= 0.0;
  double awaysideRMSxerr	= 0.0;
  double awaysidewidth		= 0.0;
  double awaysidewidth_err	= 0.0;
  double awaysideyield 		= 0.0;
  double awaysideyielderr	= 0.0;
  double awaysideyield_f 	= 0.0;
  double awaysideyield_ferr	= 0.0;
  
//   deta12ns->GetXaxis()->SetRange(deta12ns->FindBin(-sidebandlow),deta12ns->FindBin(sidebandlow));
  nearsidemeanx 	= deta12ns->GetMean(1);
  nearsideRMSx  	= deta12ns->GetRMS(1);
  nearsideRMSxerr 	= deta12ns->GetRMS(11);
  nearsideyield 	= deta12ns->IntegralAndError(deta12ns->FindBin(-sidebandlow),deta12ns->FindBin(sidebandlow),nearsideyielderr,"width");

//   deta12as->GetXaxis()->SetRange(deta12as->FindBin(-sidebandlow),deta12as->FindBin(sidebandlow));
  awaysidemeanx 	= deta12as->GetMean(1);
  awaysideRMSx  	= deta12as->GetRMS(1);
  awaysideRMSxerr 	= deta12as->GetRMS(11);
  awaysideyield 	= deta12as->IntegralAndError(deta12as->FindBin(-sidebandlow),deta12as->FindBin(sidebandlow),awaysideyielderr,"width");

  //fit nearside:
  deta12ns->Fit(fitsig1,"SQ","",-sidebandlow,sidebandlow);
  deta12ns->GetXaxis()->SetRange(1,-1);  
  nearsidewidth 	= fitsig1->GetParameter(2);
  nearsidewidth_err 	= fitsig1->GetParError(2);
  nearsideyield_f	= fitsig1->GetParameter(1);
  nearsideyield_ferr	= fitsig1->GetParError(1);

  //fit awayside:
  deta12as->Fit(fitsig1,"SQ","",-sidebandlow,sidebandlow);
  deta12as->GetXaxis()->SetRange(1,-1);
  awaysidewidth 	= fitsig1->GetParameter(2);
  awaysidewidth_err 	= fitsig1->GetParError(2);
  awaysideyield_f	= fitsig1->GetParameter(1);
  awaysideyield_ferr	= fitsig1->GetParError(1);
  
  savedircontent(dir,nearsideyield,nearsideyielderr,awaysideyield,awaysideyielderr
		    ,nearsideyield_f,nearsideyield_ferr,awaysideyield_f,awaysideyield_ferr
		    ,nearsideRMSx,nearsideRMSxerr,awaysideRMSx,awaysideRMSxerr
		    ,nearsidewidth,nearsidewidth_err,awaysidewidth,awaysidewidth_err);

  
  deta12ns->Write("deta12_Near_fitted");
  deta12as->Write("deta12_Away_fitted");

  
}

void fillwithbinnr(TH2D* fillfrom, TH1D* fillto, int binnr,bool width){
  double widthdb = fillfrom->GetXaxis()->GetBinCenter(5) - fillfrom->GetXaxis()->GetBinCenter(4);
  for(int x=1;x<=fillto->GetNbinsX();x++){
    if(!width){
      fillto->SetBinContent(x,fillfrom->GetBinContent(x,binnr));
      fillto->SetBinError(x,fillfrom->GetBinError(x,binnr));
    }
    else{
      fillto->SetBinContent(x,fillfrom->GetBinContent(x,binnr)*widthdb);
      fillto->SetBinError(x,fillfrom->GetBinError(x,binnr)*widthdb);
      
    }
  }
}
void removeconstant(TH1D * hist, Double_t plateau, Double_t erroronit){
  if(plateau == 0.0 && erroronit == 0.0 )return;
  for(int x =0;x<=hist->GetNbinsX();x++){
    if(TMath::Abs(hist->GetBinContent(x))>1.0E-10){
      hist->SetBinContent(x,hist->GetBinContent(x) + plateau);
      hist->SetBinError(x,TMath::Sqrt(hist->GetBinError(x)*hist->GetBinError(x) + erroronit*erroronit));
    }
  }
}

void sidebandcompare(TObjArray* sidebandarray,TDirectory*dir,TH2D*histo,Double_t * lowlimit, Double_t * uplimit){
  int nbins = sidebandarray->GetEntries();
  //make histogram to store variables:
  
  TH1D* chisq = new TH1D("chisq", "Average over #phi for #Chi^2 for flat sideband assumption in different sidebands ",nbins,0.0,nbins);
//   chisq->GetXaxis()->LabelsOption("v");
  chisq->SetXTitle("sideband (#Delta#eta_{12}");
  TH1D* bkg = new TH1D("bkge", "Average over #phi for the relative error on the background ",nbins,0.0,nbins);
//   bkg->GetXaxis()->LabelsOption("v");
  bkg->SetXTitle("sideband (#Delta#eta_{12}");
  TH1D* nearsideyield = new TH1D("nearyield", "nearside-yield for different sidebands ",nbins,0.0,nbins);
//   nearsideyield->GetXaxis()->LabelsOption("v");
  nearsideyield->SetXTitle("sideband (#Delta#eta_{12}");
  TH1D* awaysideyield = new TH1D("awayyield", "awayside-yield for different sidebands ",nbins,0.0,nbins);
//   awaysideyield->GetXaxis()->LabelsOption("v");
  awaysideyield->SetXTitle("sideband (|#Delta#eta_{12}|)");
  
  for(int i = 0; i<nbins;i++){
    TDirectory * nowdir = dynamic_cast<TDirectory*>(sidebandarray->At(i));
    TString name = TString(nowdir->GetName());
    TObjArray* stringar = name.Tokenize("_");
    double lowbinedgenow = dynamic_cast<TObjString*>(stringar->At(1))->GetString().Atof();
    double upbinedgenow = dynamic_cast<TObjString*>(stringar->At(2))->GetString().Atof();
    delete stringar;
    simpleyield(nowdir,histo,lowbinedgenow,upbinedgenow);
    TH1D* backgroundnow = dynamic_cast<TH1D*>(nowdir->Get("dphibackground"));
    TH1D* chisqnow = dynamic_cast<TH1D*>(nowdir->Get("chi2_background"));
    double chi2av=0.0;
    double bkgeav=0.0;
    double n = 0;
    for(int j = 1; j<=backgroundnow->GetNbinsX();j++){
      chi2av += chisqnow->GetBinContent(j);
      if(TMath::Abs(backgroundnow->GetBinError(j))>1.0E-10)bkgeav += TMath::Abs(backgroundnow->GetBinError(j)/backgroundnow->GetBinContent(j));
      n+=1.0;
    }
    TString binlable = TString(Form("%1.1f - %1.1f",lowbinedgenow, upbinedgenow));
    chi2av = chi2av/n;
    bkgeav = bkgeav/n;
    
    chisq->SetBinContent(i+1,chi2av);
    chisq->GetXaxis()->SetBinLabel(i+1,binlable.Data());
    bkg->SetBinContent(i+1,bkgeav);
    bkg->GetXaxis()->SetBinLabel(i+1,binlable.Data());
    nearsideyield->SetBinContent(i+1,dynamic_cast<TParameter<double>*>(nowdir->Get("yield_near"))->GetVal());
    nearsideyield->SetBinError(i+1,dynamic_cast<TParameter<double>*>(nowdir->Get("yield_near_error"))->GetVal());
    nearsideyield->GetXaxis()->SetBinLabel(i+1,binlable.Data());
    awaysideyield->SetBinContent(i+1,dynamic_cast<TParameter<double>*>(nowdir->Get("yield_away"))->GetVal());
    awaysideyield->SetBinError(i+1,dynamic_cast<TParameter<double>*>(nowdir->Get("yield_away_error"))->GetVal());
    awaysideyield->GetXaxis()->SetBinLabel(i+1,binlable.Data());
    
  }
  dir->cd();
  chisq->Write();
  bkg->Write();
  nearsideyield->Write();
  awaysideyield->Write();

}

void extractbinyield(TDirectory* dir, TDirectory* yielddir, Double_t etalimit,const char* options){
  if(!dir||!yielddir)return;
  gkEtaFitRange = 1.4;
//   TDirectory * bindir = resultsdirectory(yielddir,"original_binning");
//   TDirectory * binrebindir = resultsdirectory(yielddir,"Rebinned_3_to_1");
  TDirectory * simpleyielddir = resultsdirectory(yielddir,"simple");
  
  //Vary the sideband continiously for
  TDirectory* sidebands = resultsdirectory(yielddir,"sidebands");
  TObjArray * sidebandarray = new TObjArray();
  
  sidebandarray->Add(resultsdirectory(sidebands,"simple_1.1_1.5"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_1.0_1.5"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_1.0_1.4"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_1.0_1.3"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_1.0_1.2"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_1.0_1.1"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.9_1.5"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.9_1.4"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.9_1.3"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.9_1.2"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.9_1.1"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.8_1.5"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.8_1.4"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.8_1.3"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.8_1.2"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.8_1.1"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.7_1.5"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.7_1.4"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.7_1.3"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.7_1.2"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.7_1.1"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.6_1.5"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.6_1.4"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.6_1.3"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.6_1.2"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.6_1.1"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.5_1.5"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.5_1.4"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.5_1.3"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.5_1.2"));
  sidebandarray->Add(resultsdirectory(sidebands,"simple_0.5_1.1"));

  
  TDirectory * singleyielddir = resultsdirectory(yielddir,"singleyield");
  TH2D* dphideta12ss;
  if(dynamic_cast<TH2D*>(dir->Get("DPhi_1_DEta_12_SameSide"))) dphideta12ss = dynamic_cast<TH2D*>(dir->Get("DPhi_1_DEta_12_SameSide")->Clone("DPhi_1_DEta_12_SameSidec"));
  else return;
  dphideta12ss->GetXaxis()->SetRange(1,-1);//so that the ranges dont influence the procedure.
  yielddir->cd();
  //different sidebands:
  Double_t lowlimit = 0.7;
  Double_t uplimit  = 1.5;  
  if(TString(options).Contains("pp")){
    //the low pt stuff is wide in pp
    lowlimit = 1.1;
    uplimit = 1.5;
  }
  sidebandcompare(sidebandarray,sidebands,dphideta12ss,&lowlimit,&uplimit);
 
  simpleyield(simpleyielddir,dphideta12ss,lowlimit,uplimit);

  singleyield(singleyielddir,dphideta12ss,lowlimit,uplimit);

  
  
  delete dphideta12ss;
}

double GetVal(TDirectory * dir, TString name){
  if(dynamic_cast<TParameter<double>*>(dir->Get(name.Data())))
    return dynamic_cast<TParameter<double>*>(dir->Get(name.Data()))->GetVal();
  else return -100.0;  
}
TCanvas * Periodcanvas(TH1 * hist, int ncol){
  TCanvas * Canvas = new TCanvas(Form("%sCanvas",hist->GetName()));
//   TPaveLabel* title = new TPaveLabel(0.1,0.96,0.9,0.99,hist->GetTitle());
//   title->Draw();
//   TPad* graphPad = new TPad("Graphs","Graphs",0.01,0.05,0.95,0.95);
//   graphPad->Draw();
//   graphPad->cd();
  Canvas->Divide(ncol,ncol);
  return Canvas;
}

TObjArray * PadArray(TCanvas * canvas, int ntrig, int nass,TObjArray* trignames, TObjArray* assnames, const char* title){
  canvas->cd();
  if(ntrig!=trignames->GetEntries()||nass!=assnames->GetEntries()){ cout << "wrong dimensions on the names"<<endl; return 0x0;}
  TObjArray * arrayofpads = new TObjArray();
  
//   TPaveLabel* overtitle = new TPaveLabel(0.15,0.95,0.9,0.99,title,"nb");
//   overtitle->Draw();
  double textsize = 20;
  TPaveLabel* Trigtitle = new TPaveLabel(0.0,0.9,0.9,1.0,"p_{T}^{trigger}:","nb");
  Trigtitle->SetTextFont(63);
  Trigtitle->SetTextSize(textsize);
  Trigtitle->Draw();
  TPaveLabel* asstitle = new TPaveLabel(0.0,0.85,0.15,1.0,"p_{T}^{associated}:","nb");
  asstitle->SetTextFont(63);
  asstitle->SetTextSize(textsize);
  asstitle->Draw();
  
  Double_t xdiv = 0.75/ntrig;
  Double_t ydiv = 0.85/nass;
  
  
  //make trigger lables:
  for(int nx = 1;nx<=ntrig;nx++){
    canvas->cd();
    Double_t xmin = 0.15+(nx-1)*xdiv;
    Double_t xmax = xmin + xdiv;
    Double_t ymin = 0.85;
    Double_t ymax = 0.9;
    TPaveLabel* Triglable = new TPaveLabel(xmin,ymin,xmax,ymax,trignames->At(nx-1)->GetName(),"nb");
    Triglable->SetTextFont(63);
    Triglable->SetTextSize(textsize-6);
    Triglable->Draw();
  }
  //so the lables for associated:
  for(int ny = 1;ny<=nass;ny++){
    canvas->cd();
    Double_t xmin = 0.0;
    Double_t xmax = 0.15;
    Double_t ymax = 0.85-(ny-1)*ydiv;
    Double_t ymin = ymax-ydiv;
    TPaveText* Asslable = new TPaveText(xmin,ymin,xmax,ymax);
    Asslable->SetOption("nb");
    TString string = dynamic_cast<TObjString*>( assnames->At(ny-1))->GetString();
    Asslable->SetTextFont(63);
    Asslable->SetTextSize(textsize-6);
    Asslable->AddText(string.Data());
    Asslable->Draw();
  }
  
  for(int nx = 1;nx<=ntrig;nx++){for(int ny=1;ny<=nass;ny++){
    canvas->cd();
    Double_t xmin = 0.15+(nx-1)*xdiv;
    Double_t xmax = xmin + xdiv;
    Double_t ymax = 0.85-(ny-1)*ydiv;
    Double_t ymin = ymax-ydiv;
    TPad * thispad = new TPad(Form("padno_%i_%i",nx,ny),"",xmin,ymin,xmax,ymax);
    arrayofpads->Add(thispad);
    thispad->Draw();
  }}
  
  return arrayofpads;
  
  
}

TVirtualPad * cdppcanvas(TCanvas * canvas, int i){
  TPad * pad = dynamic_cast<TPad*>(canvas->cd(i));
  return pad;//->cd(i);
}
void PresentPbPb(TFile*infile,TDirectory*outpath,const char* type){
  //categorys that will be plotted:
  TObjArray * trigarray = new TObjArray();
  trigarray->Add(new TObjString("4.0-8.0 GeV/c"));
  trigarray->Add(new TObjString("8.0-16.0 GeV/c"));
  TObjArray * assarray  = new TObjArray();
  assarray->Add(new TObjString(" 0.5-1.0 GeV/c"));
  assarray->Add(new TObjString(" 1.0-2.0 GeV/c"));
  assarray->Add(new TObjString(" 2.0-3.0 GeV/c"));
  assarray->Add(new TObjString(" 3.0-4.0 GeV/c"));
  assarray->Add(new TObjString(" 4.0-6.0 GeV/c"));
  assarray->Add(new TObjString(" 6.0-8.0 GeV/c"));
  assarray->Add(new TObjString("8.0-16.0 GeV/c"));

  
  TString Ttype = TString(type);
  if(!Ttype.CompareTo("yield")==0){
    TString typetitle = TString("");
    TString histname1 = TString("DPhi_1_DEta_12_SameSide");
    TString histtitle1 = TString("#Delta#Phi_{1} vs #Delta#eta_{12}");
    TString histname2 = TString("DPhi_1_DPHI_2");
    TString histtitle2 = TString("#Delta#Phi_{1} vs #Delta#Phi_{2}");
    TString path = TString("");
    if(Ttype.CompareTo("divided")==0){path.Append("BinM#/divided/");typetitle.Append("Acceptance corrected");}
    if(Ttype.CompareTo("iteration1")==0){path.Append("BinM#/iteration1/");typetitle.Append("Fully corrected");}
    if(Ttype.CompareTo("META")==0){path.Append("META/BinM#/divided/");typetitle.Append("Acceptance corrected META");}
    if(Ttype.CompareTo("META2")==0){path.Append("META2/BinM#/divided/");typetitle.Append("Acceptance corrected META2");}
    if(Ttype.CompareTo("METrigger")==0){path.Append("METrigger/BinM#/divided/");typetitle.Append("Acceptance corrected METrigger");}
    TString path05 = TString(path.Data());
    path05.ReplaceAll("#","(0.00)->(5.00)");
    TString path510 = TString(path.Data());
    path510.ReplaceAll("#","(5.00)->(10.00)");
    TString path1020 = TString(path.Data());
    path1020.ReplaceAll("#","(10.00)->(20.00)");
    TString path2040 = TString(path.Data());
    path2040.ReplaceAll("#","(20.00)->(40.00)");
    TString path4060 = TString(path.Data());
    path4060.ReplaceAll("#","(40.00)->(60.00)");
    TString path6080 = TString(path.Data());
    path6080.ReplaceAll("#","(60.00)->(80.00)");
    //create the nessecary canvases:
    TCanvas * cent05_hist1 = new TCanvas(Form("%s_Cent05_%s",histname1.Data(),Ttype.Data()));
    TCanvas * cent510_hist1 = new TCanvas(Form("%sCent510_%s",histname1.Data(),Ttype.Data()));
    TCanvas * cent1020_hist1 = new TCanvas(Form("%sCent1020_%s",histname1.Data(),Ttype.Data()));
    TCanvas * cent2040_hist1 = new TCanvas(Form("%sCent2040_%s",histname1.Data(),Ttype.Data()));
    TCanvas * cent4060_hist1 = new TCanvas(Form("%sCent4060_%s",histname1.Data(),Ttype.Data()));
    TCanvas * cent6080_hist1 = new TCanvas(Form("%sCent6080_%s",histname1.Data(),Ttype.Data()));
    TCanvas * cent05_hist2 = new TCanvas(Form("%sCent05_%s",histname2.Data(),Ttype.Data()));
    TCanvas * cent510_hist2 = new TCanvas(Form("%sCent510_%s",histname2.Data(),Ttype.Data()));
    TCanvas * cent1020_hist2 = new TCanvas(Form("%sCent1020_%s",histname2.Data(),Ttype.Data()));
    TCanvas * cent2040_hist2 = new TCanvas(Form("%sCent2040_%s",histname2.Data(),Ttype.Data()));
    TCanvas * cent4060_hist2 = new TCanvas(Form("%sCent4060_%s",histname2.Data(),Ttype.Data()));
    TCanvas * cent6080_hist2 = new TCanvas(Form("%sCent6080_%s",histname2.Data(),Ttype.Data()));
    TObjArray* arrayofpadscor05hist1 = PadArray(cent05_hist1,2,7,trigarray,assarray,Form("%s correlation function %s in 0%%-5%% events in different pT ranges",typetitle.Data(),histtitle1.Data()));
    TObjArray* arrayofpadscor510hist1 = PadArray(cent510_hist1,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle1.Data()));
    TObjArray* arrayofpadscor1020hist1 = PadArray(cent1020_hist1,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle1.Data()));
    TObjArray* arrayofpadscor2040hist1 = PadArray(cent2040_hist1,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle1.Data()));
    TObjArray* arrayofpadscor4060hist1 = PadArray(cent4060_hist1,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle1.Data()));
    TObjArray* arrayofpadscor6080hist1 = PadArray(cent6080_hist1,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle1.Data()));
    TObjArray* arrayofpadscor05hist2 = PadArray(cent05_hist2,2,7,trigarray,assarray,Form("%s correlation function %s in 0%%-5%% events in different pT ranges",typetitle.Data(),histtitle2.Data()));
    TObjArray* arrayofpadscor510hist2 = PadArray(cent510_hist2,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle2.Data()));
    TObjArray* arrayofpadscor1020hist2 = PadArray(cent1020_hist2,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle2.Data()));
    TObjArray* arrayofpadscor2040hist2 = PadArray(cent2040_hist2,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle2.Data()));
    TObjArray* arrayofpadscor4060hist2 = PadArray(cent4060_hist2,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle2.Data()));
    TObjArray* arrayofpadscor6080hist2 = PadArray(cent6080_hist2,2,7,trigarray,assarray,Form("%s correlation function %s in 5%%-10%% events in different pT ranges",typetitle.Data(),histtitle2.Data()));
    
    
    TList * triggerlist = infile->GetListOfKeys();
    for(int i = 0;i<triggerlist->GetEntries();i++){
      int triggerindex = 0;
      if(TString(triggerlist->At(i)->GetName()).Contains("8_16")) triggerindex = 7;
      TDirectory * trigdir = dynamic_cast<TDirectory*>(infile->GetDirectory(triggerlist->At(i)->GetName()));
      TList* Asslist = trigdir->GetListOfKeys();
      for(int j =0;j<Asslist->GetEntries();j++){
	int index = triggerindex;
	if(TString(Asslist->At(j)->GetName()).Contains("1_2"))index+=1;
	if(TString(Asslist->At(j)->GetName()).Contains("2_3"))index+=2;
	if(TString(Asslist->At(j)->GetName()).Contains("3_4"))index+=3;
	if(TString(Asslist->At(j)->GetName()).Contains("4_6"))index+=4;
	if(TString(Asslist->At(j)->GetName()).Contains("6_8"))index+=5;
	if(TString(Asslist->At(j)->GetName()).Contains("8_1"))index+=6;
	TDirectory * assdir = trigdir->GetDirectory(Asslist->At(j)->GetName());
	//go through cents:
	drawinpad(assdir, path05,histname1,arrayofpadscor05hist1,index);	
	drawinpad(assdir, path510,histname1,arrayofpadscor510hist1,index);	
	drawinpad(assdir, path1020,histname1,arrayofpadscor1020hist1,index);	
	drawinpad(assdir, path2040,histname1,arrayofpadscor2040hist1,index);	
	drawinpad(assdir, path4060,histname1,arrayofpadscor4060hist1,index);	
	drawinpad(assdir, path6080,histname1,arrayofpadscor6080hist1,index);	
	drawinpad(assdir, path05,histname2,arrayofpadscor05hist2,index);	
	drawinpad(assdir, path510,histname2,arrayofpadscor510hist2,index);	
	drawinpad(assdir, path1020,histname2,arrayofpadscor1020hist2,index);	
	drawinpad(assdir, path2040,histname2,arrayofpadscor2040hist2,index);	
	drawinpad(assdir, path4060,histname2,arrayofpadscor4060hist2,index);	
	drawinpad(assdir, path6080,histname2,arrayofpadscor6080hist2,index);		
      }
    }
    outpath->cd();
    cent05_hist1->Update();
    cent05_hist1->Write();
    cent510_hist1->Update();
    cent510_hist1->Write();
    cent1020_hist1->Update();
    cent1020_hist1->Write();
    cent2040_hist1->Update();
    cent2040_hist1->Write();    
    cent4060_hist1->Update();
    cent4060_hist1->Write();  
    cent6080_hist1->Update();
    cent6080_hist1->Write();  
    cent05_hist2->Update();
    cent05_hist2->Write();
    cent510_hist2->Update();
    cent510_hist2->Write();
    cent1020_hist2->Update();
    cent1020_hist2->Write();
    cent2040_hist2->Update();
    cent2040_hist2->Write();    
    cent4060_hist2->Update();
    cent4060_hist2->Write();  
    cent6080_hist2->Update();
    cent6080_hist2->Write();  
  }
  
  
  
}
void Presentpp(TFile*infile,TDirectory*outpath,const char* type){
  //create the nessecary canvases:
  TCanvas * canvas = new TCanvas("canvas");
  
  
}


void drawinpad(TDirectory * dir, TString centpath,TString histname,TObjArray* padarray, int index,const char* drawop ){
  TH1*hist1 = dynamic_cast<TH1*>(dir->GetDirectory(centpath.Data())->Get(histname.Data())->Clone("clonehist"));
  double mineta= 0.0;
  double maxeta= 0.0;
  if(histname.Contains("DEta_12")){
    hist1->GetXaxis()->SetRangeUser(-1.0,1.0);
    mineta = hist1->GetMinimum();
    maxeta = hist1->GetMaximum();

    hist1->GetZaxis()->SetRangeUser(mineta,maxeta);
    
  }
  dynamic_cast<TPad*>(padarray->At(index))->cd();
  hist1->SetStats(false);
  hist1->GetXaxis()->SetLabelSize(0.08);
  hist1->GetXaxis()->SetTitleSize(0.09);
  hist1->GetXaxis()->SetTitleOffset(0.4);
  hist1->GetYaxis()->SetLabelSize(0.08);
  hist1->GetYaxis()->SetTitleSize(0.09);
  hist1->GetYaxis()->SetTitleOffset(0.4);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetTitle("");
  hist1->GetYaxis()->SetTitle("");
  if(hist1->GetMinimum()>0.0)hist1->GetYaxis()->SetRangeUser(0.0-hist1->GetMaximum()*0.1,hist1->GetMaximum()+hist1->GetBinError(hist1->GetMaximumBin())*1.1);
  hist1->GetXaxis()->SetRange(1,-1);
  hist1->GetYaxis()->SetRange(1,-1);
  hist1->GetZaxis()->SetTitle("");
  hist1->GetZaxis()->SetLabelSize(0.05);
  hist1->Draw(drawop);
  
  if(TString(drawop).CompareTo("colz")==0){
    dynamic_cast<TPad*>(padarray->At(index))->Update();
    TPaletteAxis * palette = (TPaletteAxis*)hist1->GetListOfFunctions()->FindObject("palette");
    if(palette){
      palette->SetX1NDC(0.92);
//       if(histname.Contains("DEta_12"))hist1->GetZaxis()->SetRangeUser(mineta,maxeta);
      hist1->GetZaxis()->SetLabelSize(0.04);
      hist1->GetZaxis()->SetLabelOffset(0.009);
      dynamic_cast<TPad*>(padarray->At(index))->Modified();
    }
  }
}
void writeinpad(TObjArray*padarray,int index, const char* text){
  dynamic_cast<TPad*>(padarray->At(index))->cd();
  TPaveText * textp = new TPaveText(0.1,0.1,0.9,0.9,"NB");
  textp->SetTextFont(63);
  textp->SetTextSize(15);
  textp->AddText(text);
  textp->Draw();
}
void getminmax(TDirectory * dir,TString centpath, TString histname, double  &min,double & max){
  TH1*hist1 = dynamic_cast<TH1*>(dir->GetDirectory(centpath.Data())->Get(histname.Data()));
  int minbin = hist1->GetMinimumBin();
  double min_tmp = hist1->GetBinContent(minbin) - hist1->GetBinError(minbin);
  int maxbin = hist1->GetMaximumBin();
  double max_tmp = hist1->GetBinContent(maxbin)+hist1->GetBinError(maxbin);
  if(min>min_tmp)min = min_tmp;
  if(max<max_tmp)max = max_tmp;
  
  
}
void getminmax(TH1*hist1, double  &min,double & max){
  int minbin = hist1->GetMinimumBin();
  double min_tmp = hist1->GetBinContent(minbin) - hist1->GetBinError(minbin);
  int maxbin = hist1->GetMaximumBin();
  double max_tmp = hist1->GetBinContent(maxbin)+hist1->GetBinError(maxbin);
  if(min>min_tmp)min = min_tmp;
  if(max<max_tmp)max = max_tmp;
  
  
}
void getminmax(TH1*hist1, double  &min,double & max, double firstbin, double lastbin){
  hist1->GetXaxis()->SetRangeUser(firstbin, lastbin);
  int minbin = hist1->GetMinimumBin();
  double min_tmp = hist1->GetBinContent(minbin) - hist1->GetBinError(minbin);
  int maxbin = hist1->GetMaximumBin();
  double max_tmp = hist1->GetBinContent(maxbin)+hist1->GetBinError(maxbin);
  if(min>min_tmp)min = min_tmp;
  if(max<max_tmp)max = max_tmp;
  
}
void getminmaxpos(TH1*hist,double  &min,double & max){
  int minbin = hist->GetMinimumBin();
  double min_tmp = hist->GetBinContent(minbin);
  if(min_tmp-hist->GetBinError(minbin)>0.0)min_tmp-=hist->GetBinError(minbin);
  int maxbin = hist->GetMaximumBin();
  double max_tmp = hist->GetBinContent(maxbin)+hist->GetBinError(maxbin);
  if(min>min_tmp&&min_tmp>0.0)min = min_tmp;
  if(max<max_tmp)max = max_tmp;  
  }

void drawcanvas(TDirectory * dir, TString centpath,TString histname,int histcolor,TCanvas * canvas, TLegend * legend, TString legendname,TString histtitle,const char* drawop, double minimum, double maximum){
  TH1*hist1 = dynamic_cast<TH1*>(dir->GetDirectory(centpath.Data())->Get(histname.Data())->Clone("clonehist"));
  hist1->SetLineColor(histcolor);
  canvas->cd();
  hist1->SetStats(false);
  if(hist1->GetFunction("fgs"))hist1->GetFunction("fgs")->SetBit(TF1::kNotDraw);
  hist1->SetTitle("");
  if(minimum>-100.0&&maximum>-100)hist1->GetYaxis()->SetRangeUser(minimum,maximum);
  hist1->SetTitle(histtitle.Data());
//   hist1->GetXaxis()->SetTitleSize(0.03);
  hist1->GetXaxis()->SetTitleOffset(0.9);
  hist1->GetYaxis()->SetTitleSize(0.035);
  hist1->Draw(drawop);
  legend->AddEntry(hist1,legendname.Data());
}

void draw(TDirectory * dir, TString centpath,TString histname,TString histtitle,const char* drawop,TString picdir){
  histname.ReplaceAll("_fitted","");
  TH1*hist1 = dynamic_cast<TH1*>(dir->GetDirectory(centpath.Data())->Get(histname.Data())->Clone("clonehist"));
  TCanvas * canvas = new TCanvas();
  hist1->SetStats(false);
  hist1->SetTitle("");
  hist1->GetXaxis()->SetTitleOffset(0.9);
  hist1->GetYaxis()->SetTitleSize(0.035);
  hist1->Draw(drawop);
  canvas->SaveAs(picdir.Data());
  delete canvas;delete hist1;
  
}


void drawptbins(TObjArray* indirs, TDirectory* outdir, const TString histname,TString type, TString Cent, TObjArray* trigarray, TObjArray* assarray){
  TString canvasname = TString(histname+TString("_")+type);
  if(Cent.CompareTo("")!=0) canvasname.Append(Form("_%s",Cent.Data()));
  TCanvas * canvas = new TCanvas(canvasname);
  TString path = TString("");
  TString singlename = TString("");
  bool drawsingle = false;
  TObjArray * dirarray = TString(gSystem->pwd()).Tokenize("/");
  TString typetitle = TString("");
  if(type.CompareTo("divided")==0){path.Append("BinM#/divided/");typetitle.Append("Acceptance corrected");}
  if(type.CompareTo("iteration1")==0){path.Append("BinM#/iteration1/");typetitle.Append("Fully corrected");}
  if(type.CompareTo("Etasubstracted")==0){path.Append("BinM#/Etasubstracted/");typetitle.Append("Flat in long range Eta12 corrected");}
  if(type.CompareTo("META")==0){path.Append("META/BinM#/divided/");typetitle.Append("Acceptance corrected META");}
  if(type.CompareTo("META2")==0){path.Append("META2/BinM#/divided/");typetitle.Append("Acceptance corrected META2");}
  if(type.CompareTo("METrigger")==0){path.Append("METrigger/BinM#/divided/");typetitle.Append("Acceptance corrected METrigger");}
  TString histtitle = TString("");http://www.sge4ever.de/cartoon-mit-andi-mike-folge-13/
  TString histobjname = TString("");
  TString drawop      = TString("");
  if(histname.CompareTo("DPHIDPHI")==0){
    histtitle.Append("correlation function d#Phi_{1} vs d#Phi_{2}");
    histobjname.Append("DPhi_1_DPHI_2");
    drawop.Append("colz");
  }
  if(histname.CompareTo("DETA_12")==0){
    histtitle.Append("correlation function d#Phi_{1} vs d#eta_{12} ");
    histobjname.Append("DPhi_1_DEta_12");
    drawop.Append("colz");
  }
  if(histname.CompareTo("DETA_12ss")==0){
    histtitle.Append("correlation function d#Phi_{1} vs d#eta_{12} with both associated on the same side");
    histobjname.Append("DPhi_1_DEta_12_SameSide");
    drawop.Append("colz");
  }
  if(histname.CompareTo("yield")==0){
    histtitle.Append("yield as a function of d#Phi_{1}");
    histobjname.Append("dphiyieldbc");
    path.ReplaceAll("divided/","yield/simple/");
    path.ReplaceAll("iteration1/","yield/simple/");
    path.ReplaceAll("Etasubstracted/","Etasubstracted/yield/simple/");
    drawop.Append("E");
  }
  if(histname.CompareTo("nearpeak")==0){
    histtitle.Append("Nearside yield as a function of d#eta_{12}");
    histobjname.Append("deta12_Near_fitted");
    path.ReplaceAll("divided/","yield/singleyield/");
    path.ReplaceAll("iteration1/","yield/singleyield/");
    path.ReplaceAll("Etasubstracted/","Etasubstracted/yield/singleyield/");
    drawop.Append("E");
    drawsingle = true;
  }
  if(histname.CompareTo("awaypeak")==0){
    histtitle.Append("Awayside yield as a function of d#eta_{12}");
    histobjname.Append("deta12_Away_fitted");
    path.ReplaceAll("divided/","yield/singleyield/");
    path.ReplaceAll("iteration1/","yield/singleyield/");
    path.ReplaceAll("Etasubstracted/","Etasubstracted/yield/singleyield/");
    drawop.Append("E");
    drawsingle = true;
  }
  TString canvastitle = TString(Form("%s %s in different pT ranges ",typetitle.Data(),histtitle.Data()));
  TString centtitle = TString("");
  TString centdir = TString("");
  if(Cent.CompareTo("")==0){centtitle.Append("in minimum bias events.");path.ReplaceAll("BinM#/","");}
  if(Cent.CompareTo("0_5")==0){centtitle.Append(Form("in Events that have C=0%%-5%% ."));path.ReplaceAll("#","(0.00)->(5.00)");}
  if(Cent.CompareTo("5_10")==0){centtitle.Append(Form("in Events that have C=5%%-10%% ."));path.ReplaceAll("#","(5.00)->(10.00)");}
  if(Cent.CompareTo("0_10")==0){centtitle.Append(Form("in Events that have C=0%%-10%% ."));path.ReplaceAll("#","(0.00)->(10.00)");}
  if(Cent.CompareTo("10_20")==0){centtitle.Append(Form("in Events that have C=10%%-20%% ."));path.ReplaceAll("#","(10.00)->(20.00)");}
  if(Cent.CompareTo("20_40")==0){centtitle.Append(Form("in Events that have C=20%%-40%% ."));path.ReplaceAll("#","(20.00)->(40.00)");}
  if(Cent.CompareTo("40_60")==0){centtitle.Append(Form("in Events that have C=40%%-60%% ."));path.ReplaceAll("#","(40.00)->(60.00)");}
  if(Cent.CompareTo("60_90")==0){centtitle.Append(Form("in Events that have C=60%%-90%% ."));path.ReplaceAll("#","(60.00)->(90.00)");}
  canvastitle.Append(centtitle);
  TString prodname = TString("");
  TObjArray* arrayofpads = PadArray(canvas,3,6,trigarray,assarray,canvastitle);
  if(!histobjname.CompareTo("")==0){
    for(int i = 0; i<indirs->GetEntries();i++){
      TDirectory * dir = dynamic_cast<TDirectory*>(indirs->At(i));
      TString whichbin =  TString(dir->GetTitle());
      TObjArray * tokens = whichbin.Tokenize("_");
      if(i == 0 )prodname.Append(dynamic_cast<TObjString*>(tokens->At(0))->GetString().Data());
      int indexonpad = 0;
      if(dynamic_cast<TObjString*>(tokens->At(2))->GetString().Atoi()==8)indexonpad+=6;
      if(dynamic_cast<TObjString*>(tokens->At(2))->GetString().Atoi()==16)indexonpad+=12;
      if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==2)indexonpad+=0;
      if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==3)indexonpad+=1;
      if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==4)indexonpad+=2;
      if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==6)indexonpad+=3;
      if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==8)indexonpad+=4;
      if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==16)indexonpad+=5;
      drawinpad(dir,path,histobjname,arrayofpads,indexonpad);
      if(drawsingle){
	TString ptname = TString();
	if(dynamic_cast<TObjString*>(tokens->At(2))->GetString().Atoi()==4)ptname.Append("_48");
	if(dynamic_cast<TObjString*>(tokens->At(2))->GetString().Atoi()==8)ptname.Append("_816");
	if(dynamic_cast<TObjString*>(tokens->At(2))->GetString().Atoi()==16)ptname.Append("_1650");
	if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==2)ptname.Append("_23");
	if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==3)ptname.Append("_34");
	if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==4)ptname.Append("_46");
	if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==6)ptname.Append("_68");
	if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==8)ptname.Append("_816");
	if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==16)ptname.Append("_1650");
	singlename.Append(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/%s/singles/%s_%s.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),prodname.Data(),canvasname.Data(),ptname.Data()));
	draw(dir,path,histobjname,"",drawop,singlename);
	singlename.Clear();
      }
      if(indexonpad == 0){
	TString xname = TString(dynamic_cast<TH1*>(dir->GetDirectory(path.Data())->Get(histobjname.Data()))->GetXaxis()->GetTitle());
	xname.ReplaceAll("[]","");
	writeinpad(arrayofpads,4,Form("x-axis: %s",xname.Data()));
	writeinpad(arrayofpads,5,Form("y-axis: %s",dynamic_cast<TH1*>(dir->GetDirectory(path.Data())->Get(histobjname.Data()))->GetYaxis()->GetTitle()));
	if(histobjname.Contains("DPhi_1"))writeinpad(arrayofpads,1,Form("z-axis: %s",dynamic_cast<TH1*>(dir->GetDirectory(path.Data())->Get(histobjname.Data()))->GetZaxis()->GetTitle()));
      }
      delete tokens;
    }
  }
  outdir->cd();
  canvas->Write();
  gSystem->mkdir(prodname.Data());
  gSystem->mkdir(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/%s",dirarray->At(dirarray->GetEntries()-1)->GetName(),prodname.Data()));
  canvas->SaveAs(Form("%s/%s.eps",prodname.Data(),canvasname.Data()));
  canvas->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/%s/%s.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),prodname.Data(),canvasname.Data()));
  delete dirarray;  
  delete arrayofpads;
  delete canvas;
}

void drawcentcom(TObjArray* indirs, TDirectory* outdir, TString histname, TString type){
  //Draw histograms for all centralities in one canvas for all pT bins  
  //Save under prodname/Centralities
  
  //define colors and names
  TString c010 = TString("C = 0% - 10%");
  int color010 = 1;
  TString c05 = TString("C = 0% - 5%");  
  int color05 = 1;
  TString c510 = TString("C = 5% - 10%");  
  int color510 = 7;
  TString c1020 = TString("C = 10% - 20%");
  int color1020 = 4;
  TString c2040 = TString("C = 20% - 40%");  
  int color2040 = 6;
  TString c4060 = TString("C = 40% - 60%");  
  int color4060 = 3;
  TString c6090 = TString("C = 60% - 90%");
  int color6090 = 28;

  TString canvasname = TString(histname+TString("_")+type);
  TCanvas * canvas = new TCanvas(canvasname);
  TLegend * legend = new TLegend(0.7,0.75,0.9,0.9);
  TString path = TString("");
  TString typetitle = TString("");
  if(type.CompareTo("iteration1")==0){path.Append("iteration1/");typetitle.Append("Fully corrected");}
  if(type.CompareTo("Etasubstracted")==0){path.Append("Etasubstracted/");typetitle.Append("Flat in long range Eta12 corrected");}

  TString histtitle = TString("");
  TString histobjname = TString("");
  TString drawop      = TString("");
  if(histname.CompareTo("yield")==0){
    histtitle.Append("yield as a function of d#Phi_{1}");
    histobjname.Append("dphiyieldbc");
    path.ReplaceAll("iteration1/","yield/simple/");
    path.ReplaceAll("Etasubstracted/","Etasubstracted/yield/simple/");
    drawop.Append("E");
  }
  if(histname.CompareTo("nearpeak")==0){
    histtitle.Append("Nearside yield as a function of d#eta_{12}");
    histobjname.Append("deta12_Near_fitted");
    path.ReplaceAll("iteration1/","yield/singleyield/");
    path.ReplaceAll("Etasubstracted/","Etasubstracted/yield/singleyield/");
    drawop.Append("E");
  }
  if(histname.CompareTo("awaypeak")==0){
    histtitle.Append("Awayside yield as a function of d#eta_{12}");
    histobjname.Append("deta12_Away_fitted");
    path.ReplaceAll("iteration1/","yield/singleyield/");
    path.ReplaceAll("Etasubstracted/","Etasubstracted/yield/singleyield/");
    drawop.Append("E");
  }
  TString canvastitle = TString(Form("%s %s in different centrality ranges ",typetitle.Data(),histtitle.Data()));


  TString prodname = TString("");
  bool hasbeendrawn = false;
  if(!histobjname.CompareTo("")==0){
    for(int i = 0; i<indirs->GetEntries();i++){
      TDirectory * dir = dynamic_cast<TDirectory*>(indirs->At(i));
      dir->pwd();
      //name of pT bin:
      TString triggername = TString("_");
      TString assname = TString("_");
      TString histtitenow = TString(histtitle.Data());
      
      TString whichbin =  TString(dir->GetTitle());
      TObjArray * tokens = whichbin.Tokenize("_");
      if(i == 0){
	prodname.Append(dynamic_cast<TObjString*>(tokens->At(0))->GetString().Data());
	gSystem->mkdir(Form("%s/Centralities",prodname.Data()));
      }
      triggername.Append(dynamic_cast<TObjString*>(tokens->At(2))->GetString().Data());
      triggername.Append("_");
      triggername.Append(dynamic_cast<TObjString*>(tokens->At(3))->GetString().Data());
      assname.Append(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Data());
      assname.Append("_");
      histtitenow.Append(Form(" in the bin %s GeV/c<p_{T-trigger}<%s GeV/c, %s GeV/c<p_{T-associated}<",dynamic_cast<TObjString*>(tokens->At(2))->GetString().Data(),dynamic_cast<TObjString*>(tokens->At(3))->GetString().Data(),dynamic_cast<TObjString*>(tokens->At(5))->GetString().Data()));
      if(dynamic_cast<TObjString*>(tokens->At(6))->GetString().CompareTo("51")){
	assname.Append(dynamic_cast<TObjString*>(tokens->At(6))->GetString().Data());
	histtitenow.Append(dynamic_cast<TObjString*>(tokens->At(6))->GetString().Data());
	histtitenow.Append(" GeV/c");
      }
      else{
	assname.Append("50");
	histtitenow.Append("50 GeV/c");
      }
      delete tokens;
      canvas->cd();
      

            
      double min = 0.0;
      double max = 0.0;
      
      TDirectory* centdir05 = dir->GetDirectory("BinM(0.00)->(5.00)");
      TDirectory* centdir010 = dir->GetDirectory("BinM(0.00)->(10.00)");
      TDirectory* centdir510 = dir->GetDirectory("BinM(5.00)->(10.00)");
      TDirectory* centdir1020 = dir->GetDirectory("BinM(10.00)->(20.00)");
      TDirectory* centdir2040 = dir->GetDirectory("BinM(20.00)->(40.00)");
      TDirectory* centdir4060 = dir->GetDirectory("BinM(40.00)->(60.00)");
      TDirectory* centdir6090 = dir->GetDirectory("BinM(60.00)->(90.00)");
      if(centdir05)getminmax(centdir05,path,histobjname,min,max);
      if(centdir010)getminmax(centdir010,path,histobjname,min,max);
      if(centdir510)getminmax(centdir510,path,histobjname,min, max);      
      if(centdir1020)getminmax(centdir1020,path,histobjname,min, max);
      if(centdir2040)getminmax(centdir2040,path,histobjname,min, max);
      if(centdir4060)getminmax(centdir4060,path,histobjname,min, max);
      if(centdir6090)getminmax(centdir6090,path,histobjname,min, max);
      if(min<0.0) min*=1.1; 
      if(min>=0.0&&max>0.0) min = 0.0 -0.2*max;
      if(min>=0.0&&max<=0.0) min = 0.0;
      if(max>0.0) max = max*1.1;
      else max = 0.0;
      
      if(centdir05)drawcanvas(centdir05,path,histobjname,color05,canvas,legend,c05,histtitenow,"E",min,max);
      if(centdir010)drawcanvas(centdir010,path,histobjname,color010,canvas,legend,c010,histtitenow,"E",min,max);
      if(centdir510)drawcanvas(centdir510,path,histobjname,color510,canvas,legend,c510,histtitenow);      
      if(centdir1020)drawcanvas(centdir1020,path,histobjname,color1020,canvas,legend,c1020,histtitenow);
      if(centdir2040)drawcanvas(centdir2040,path,histobjname,color2040,canvas,legend,c2040,histtitenow);
      if(centdir4060)drawcanvas(centdir4060,path,histobjname,color4060,canvas,legend,c4060,histtitenow);
      if(centdir6090)drawcanvas(centdir6090,path,histobjname,color6090,canvas,legend,c6090,histtitenow);
      outdir->cd();
      legend->Draw("same");
      canvas->SetTitle(histtitenow.Data());
      canvas->Write(Form("%s%s%s",canvas->GetName(),triggername.Data(),assname.Data()));
      TObjArray * dirarray = TString(gSystem->pwd()).Tokenize("/");
      gSystem->mkdir(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/%s/Centralities",dirarray->At(dirarray->GetEntries()-1)->GetName(),prodname.Data()));
      canvas->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/%s/Centralities/%s%s%s.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),prodname.Data(),canvasname.Data(),triggername.Data(),assname.Data()));
      delete dirarray;        
      canvas->SaveAs(Form("%s/Centralities/%s%s%s.eps",prodname.Data(),canvasname.Data(),triggername.Data(),assname.Data()));      
      legend->Clear();
      canvas->Clear();
    }
  }

  delete canvas;
  
  
    
}



void ComparePt(bool isPbPb,TObjArray* indirs, TDirectory * outdir, TString histname){
  //categorys that will be plotted:
  TObjArray * trigarray = new TObjArray();
  trigarray->Add(new TObjString(" 4.0- 8.0 GeV/c"));
  trigarray->Add(new TObjString(" 8.0-16.0 GeV/c"));
  trigarray->Add(new TObjString("16.0-50.0 GeV/c"));
  TObjArray * assarray  = new TObjArray();
  assarray->Add(new TObjString(" 2.0- 3.0 GeV/c"));
  assarray->Add(new TObjString(" 3.0- 4.0 GeV/c"));
  assarray->Add(new TObjString(" 4.0- 6.0 GeV/c"));
  assarray->Add(new TObjString(" 6.0- 8.0 GeV/c"));
  assarray->Add(new TObjString(" 8.0-16.0 GeV/c"));
  assarray->Add(new TObjString("16.0-50.0 GeV/c"));
  
  
  
  TObjArray * hists = histname.Tokenize("#");
  for(int k = 0; k< hists->GetEntries();k++){
    TString histnamenow =  dynamic_cast<TObjString*>(hists->At(k))->GetString();
    bool yield = false;
    if(histnamenow.CompareTo("yield")==0)	yield=true;
    if(histnamenow.CompareTo("nearpeak")==0)	yield=true;
    if(histnamenow.CompareTo("awaypeak")==0)	yield=true;
    if(!isPbPb){
      if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","",trigarray,assarray);
      drawptbins(indirs,outdir,histnamenow,"iteration1","",trigarray,assarray);
      drawptbins(indirs,outdir,histnamenow,"Etasubstracted","",trigarray,assarray);
      if(!yield) drawptbins(indirs,outdir,histnamenow,"META","",trigarray,assarray);
      if(!yield) drawptbins(indirs,outdir,histnamenow,"META2","",trigarray,assarray);
      if(!yield) drawptbins(indirs,outdir,histnamenow,"METrigger","",trigarray,assarray);
    }
    if(isPbPb){
	if(dynamic_cast<TDirectory*>(indirs->At(0))->GetDirectory("BinM(0.00)->(10.00)")){
	  if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","0_10",trigarray,assarray);
	  drawptbins(indirs,outdir,histnamenow,"iteration1","0_10",trigarray,assarray);
	  drawptbins(indirs,outdir,histnamenow,"Etasubstracted","0_10",trigarray,assarray);
	}
	else{
	  if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","0_5",trigarray,assarray);
	  drawptbins(indirs,outdir,histnamenow,"iteration1","0_5",trigarray,assarray);
	  drawptbins(indirs,outdir,histnamenow,"Etasubstracted","0_5",trigarray,assarray);
	  if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","5_10",trigarray,assarray);
	  drawptbins(indirs,outdir,histnamenow,"iteration1","5_10",trigarray,assarray);
	  drawptbins(indirs,outdir,histnamenow,"Etasubstracted","5_10",trigarray,assarray);
	}
        if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","10_20",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"iteration1","10_20",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"Etasubstracted","10_20",trigarray,assarray);
        if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","20_40",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"iteration1","20_40",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"Etasubstracted","20_40",trigarray,assarray);
        if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","40_60",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"iteration1","40_60",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"Etasubstracted","40_60",trigarray,assarray);
        if(!yield)drawptbins(indirs,outdir,histnamenow,"divided","60_90",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"iteration1","60_90",trigarray,assarray);
        drawptbins(indirs,outdir,histnamenow,"Etasubstracted","60_90",trigarray,assarray);
      
    }
  }
  delete hists;
  

  delete trigarray;
  delete assarray;
}

void CompareCent(TObjArray* indirs, TDirectory * outdir, TString histname){
  
  TObjArray * hists = histname.Tokenize("#");
  for(int k = 0; k< hists->GetEntries();k++){
    TString histnamenow =  dynamic_cast<TObjString*>(hists->At(k))->GetString();
    bool yield = false;
    if(histnamenow.CompareTo("yield")==0)	yield=true;
    if(histnamenow.CompareTo("nearpeak")==0)	yield=true;
    if(histnamenow.CompareTo("awaypeak")==0)	yield=true;
    if(yield){
      drawcentcom(indirs,outdir,histnamenow,"iteration1");
      drawcentcom(indirs,outdir,histnamenow,"Etasubstracted");
    }
  
  }
}

void Compare(TObjArray * indirs , TDirectory* outdir , TString histname , TString comparetype, TString titlecard){
  if(comparetype.Contains("pp")){
    if(comparetype.Contains("pTbins"))ComparePt(false,indirs, outdir,histname);
  }
  if(comparetype.Contains("PbPb")){
    if(comparetype.Contains("pTbins"))ComparePt(true,indirs,outdir,histname);    
    if(comparetype.Contains("pTbins"))CompareCent(indirs, outdir, histname);
  }
  
}

void fillwithvalues(TDirectory* dir,int bin,TH1D* nearyield, TH1D* awayhist,TH1D* nearyield_f, TH1D* awayhist_f, 
		    TH1D* width_etanear,   TH1D* width_etaaway, TH1D* width_thetanear, TH1D* width_thetaaway)
{
  double nearside 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near"))->GetVal();
  double nearside_e 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near_error"))->GetVal();
  double awayside 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away"))->GetVal();
  double awayside_e 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away_error"))->GetVal();

  double nearside_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near_fitted"))->GetVal();
  double nearside_e_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near_fitted_error"))->GetVal();
  double awayside_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away_fitted"))->GetVal();
  double awayside_e_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away_fitted_error"))->GetVal();

  double width_eta_near 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etanear"))->GetVal();
  double width_eta_near_e 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etanear_error"))->GetVal();
  double width_eta_away 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etaaway"))->GetVal();
  double width_eta_away_e 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etaaway_error"))->GetVal();

  double width_theta_near 		= dynamic_cast<TParameter<double>*>(dir->Get("width_thetanear"))->GetVal();
  double width_theta_near_e 		= dynamic_cast<TParameter<double>*>(dir->Get("width_thetanear_error"))->GetVal();
  double width_theta_away 		= dynamic_cast<TParameter<double>*>(dir->Get("width_thetaaway"))->GetVal();
  double width_theta_away_e 		= dynamic_cast<TParameter<double>*>(dir->Get("width_thetaaway_error"))->GetVal();
  
  double dpt = 1.0;
  if(bin==1) dpt = 1.0;
  if(bin==2) dpt = 1.0;
  if(bin==3) dpt = 2.0;
  if(bin==4) dpt = 2.0;
  if(bin==5) dpt = 8.0;
  if(bin==6) dpt = 34.0;
  
  nearyield->SetBinContent(bin,nearside/dpt);
  nearyield->SetBinError(bin,nearside_e/dpt);
  awayhist->SetBinContent(bin,awayside/dpt);
  awayhist->SetBinError(bin,awayside_e/dpt);

  nearyield_f->SetBinContent(bin,nearside_f/dpt);
  nearyield_f->SetBinError(bin,nearside_e_f/dpt);
  awayhist_f->SetBinContent(bin,awayside_f/dpt);
  awayhist_f->SetBinError(bin,awayside_e_f/dpt);

  width_etanear->SetBinContent(bin, width_eta_near);
  width_etanear->SetBinError(bin, width_eta_near_e);
  width_etaaway->SetBinContent(bin, width_eta_away);
  width_etaaway->SetBinError(bin, width_eta_away_e); 

  width_thetanear->SetBinContent(bin, width_theta_near);
  width_thetanear->SetBinError(bin, width_theta_near_e);
  width_thetaaway->SetBinContent(bin, width_theta_away);
  width_thetaaway->SetBinError(bin, width_theta_away_e);   
}

void fillwithvalues(TDirectory* dir,int bin,TH1D* nearyield, TH1D* awayhist,TH1D* nearyield_f, TH1D* awayhist_f, 
		    TH1D* width_etanear,   TH1D* width_etaaway)
{
  double nearside 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near"))->GetVal();
  double nearside_e 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near_error"))->GetVal();
  double awayside 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away"))->GetVal();
  double awayside_e 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away_error"))->GetVal();

  double nearside_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near_fitted"))->GetVal();
  double nearside_e_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_near_fitted_error"))->GetVal();
  double awayside_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away_fitted"))->GetVal();
  double awayside_e_f 			= dynamic_cast<TParameter<double>*>(dir->Get("yield_away_fitted_error"))->GetVal();

  double width_eta_near 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etanear"))->GetVal();
  double width_eta_near_e 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etanear_error"))->GetVal();
  double width_eta_away 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etaaway"))->GetVal();
  double width_eta_away_e 		= dynamic_cast<TParameter<double>*>(dir->Get("width_etaaway_error"))->GetVal();

  
  double dpt = 1.0;
  if(bin==1) dpt = 1.0;
  if(bin==2) dpt = 1.0;
  if(bin==3) dpt = 2.0;
  if(bin==4) dpt = 2.0;
  if(bin==5) dpt = 8.0;
  if(bin==6) dpt = 34.0;
  
  nearyield->SetBinContent(bin,nearside/dpt);
  nearyield->SetBinError(bin,nearside_e/dpt);
  awayhist->SetBinContent(bin,awayside/dpt);
  awayhist->SetBinError(bin,awayside_e/dpt);
  
  nearyield_f->SetBinContent(bin,nearside_f/dpt);
  nearyield_f->SetBinError(bin,nearside_e_f/dpt);
  awayhist_f->SetBinContent(bin,awayside_f/dpt);
  awayhist_f->SetBinError(bin,awayside_e_f/dpt);

  width_etanear->SetBinContent(bin, width_eta_near);
  width_etanear->SetBinError(bin, width_eta_near_e);
  width_etaaway->SetBinContent(bin, width_eta_away);
  width_etaaway->SetBinError(bin, width_eta_away_e); 

}

void extractvalues(TObjArray * indirs, TDirectory* outdir ,TString subdir,TH1D* hist,TString methodstring, TString ppperiod){
  //extracts nearside and awayside yields and widths
  TString yaxisnear = TString("#frac{N_{pairs near}}{N_{triggers}dp_{T}}");
  TString yaxisaway = TString("#frac{N_{pairs away}}{N_{triggers}dp_{T}}");
  TString yaxiswidtheta = TString("Width (#eta_{12}) []");
  TString yaxiswidththeta = TString("Width (#theta_{1}) [rad]");

  bool ispp = false;
  if(subdir.CompareTo("") == 0)ispp = true;
  
  TString titlestring = TString();
  if(subdir.CompareTo("BinM(0.00)->(5.00)")==0)titlestring.Append(Form("PbPb with C = 0%%-5%%") );
  if(subdir.CompareTo("BinM(0.00)->(10.00)")==0)titlestring.Append(Form("PbPb with C = 0%%-10%%") );
  if(subdir.CompareTo("BinM(5.00)->(10.00)")==0)titlestring.Append(Form("PbPb with C = 5%%-10%%") );
  if(subdir.CompareTo("BinM(10.00)->(20.00)")==0)titlestring.Append(Form("PbPb with C = 10%%-20%%") );
  if(subdir.CompareTo("BinM(20.00)->(40.00)")==0)titlestring.Append(Form("PbPb with C = 20%%-40%%") );
  if(subdir.CompareTo("BinM(40.00)->(60.00)")==0)titlestring.Append(Form("PbPb with C = 40%%-60%%") );
  if(subdir.CompareTo("BinM(60.00)->(90.00)")==0)titlestring.Append(Form("PbPb with C = 60%%-90%%") );
  TString extraname = TString("");
  if(ppperiod.Contains("LHC11a")&&ispp)	extraname.Append("_11a");
  if(methodstring.Contains("Eta"))extraname.Append("_etasubstracted");
  if(methodstring.Contains("singleyield"))extraname.Append("_single");
  
  TH1D* nearyieldhist48 = dynamic_cast<TH1D*>(hist->Clone(Form("NearYield48%s",extraname.Data())));
  nearyieldhist48->SetYTitle(yaxisnear.Data());
  nearyieldhist48->SetTitle(Form("Near-side peak yield in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* nearyieldhist816 = dynamic_cast<TH1D*>(hist->Clone(Form("NearYield816%s",extraname.Data())));
  nearyieldhist816->SetYTitle(yaxisnear.Data());
  nearyieldhist816->SetTitle(Form("Near-side peak yield in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data()));
  TH1D* nearyieldhist1650 = dynamic_cast<TH1D*>(hist->Clone(Form("NearYield1650%s",extraname.Data())));
  nearyieldhist1650->SetYTitle(yaxisnear.Data());
  nearyieldhist1650->SetTitle(Form("Near-side peak yield in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data()));
  TH1D* awayieldhist48 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayYield48%s",extraname.Data())));
  awayieldhist48->SetYTitle(yaxisaway.Data());
  awayieldhist48->SetTitle(Form("Away-side peak yield in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* awayieldhist816 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayYield816%s",extraname.Data())));
  awayieldhist816->SetYTitle(yaxisaway.Data());
  awayieldhist816->SetTitle(Form("Away-side peak yield in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data()));
  TH1D* awayieldhist1650 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayYield1650%s",extraname.Data())));
  awayieldhist1650->SetYTitle(yaxisaway.Data());
  awayieldhist1650->SetTitle(Form("Away-side peak yield in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data()));
  
  TH1D* nearyield_f_hist48 = dynamic_cast<TH1D*>(hist->Clone(Form("NearYield_f_48%s",extraname.Data())));
  nearyield_f_hist48->SetYTitle(yaxisnear.Data());
  nearyield_f_hist48->SetTitle(Form("Near-side peak yield in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* nearyield_f_hist816 = dynamic_cast<TH1D*>(hist->Clone(Form("NearYield_f_816%s",extraname.Data())));
  nearyield_f_hist816->SetYTitle(yaxisnear.Data());
  nearyield_f_hist816->SetTitle(Form("Near-side peak yield in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data()));
  TH1D* nearyield_f_hist1650 = dynamic_cast<TH1D*>(hist->Clone(Form("NearYield_f_1650%s",extraname.Data())));
  nearyield_f_hist1650->SetYTitle(yaxisnear.Data());
  nearyield_f_hist1650->SetTitle(Form("Near-side peak yield in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data()));
  TH1D* awayield_f_hist48 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayYield_f_48%s",extraname.Data())));
  awayield_f_hist48->SetYTitle(yaxisaway.Data());
  awayield_f_hist48->SetTitle(Form("Away-side peak yield in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* awayield_f_hist816 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayYield_f_816%s",extraname.Data())));
  awayield_f_hist816->SetYTitle(yaxisaway.Data());
  awayield_f_hist816->SetTitle(Form("Away-side peak yield in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data()));  
  TH1D* awayield_f_hist1650 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayYield_f_1650%s",extraname.Data())));
  awayield_f_hist1650->SetYTitle(yaxisaway.Data());
  awayield_f_hist1650->SetTitle(Form("Away-side peak yield in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data())); 
  
  TH1D* widthnear_eta_48 = dynamic_cast<TH1D*>(hist->Clone(Form("NearWidth_eta12_48%s",extraname.Data())));
  widthnear_eta_48->SetYTitle(yaxiswidtheta.Data());
  widthnear_eta_48->SetTitle(Form("Near-side peak width in #Delta#eta_{12} in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* widthnear_eta_816 = dynamic_cast<TH1D*>(hist->Clone(Form("NearWidth_eta12_816%s",extraname.Data())));
  widthnear_eta_816->SetYTitle(yaxiswidtheta.Data());
  widthnear_eta_816->SetTitle(Form("Near-side peak width in #Delta#eta_{12} in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data()));
  TH1D* widthnear_eta_1650 = dynamic_cast<TH1D*>(hist->Clone(Form("NearWidth_eta12_1650%s",extraname.Data())));
  widthnear_eta_1650->SetYTitle(yaxiswidtheta.Data());
  widthnear_eta_1650->SetTitle(Form("Near-side peak width in #Delta#eta_{12} in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data()));
  TH1D* widthaway_eta_48 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayWidth_eta12_48%s",extraname.Data())));
  widthaway_eta_48->SetYTitle(yaxiswidtheta.Data());
  widthaway_eta_48->SetTitle(Form("Away-side peak width in #Delta#eta_{12} in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* widthaway_eta_816 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayWidth_eta12_816%s",extraname.Data())));
  widthaway_eta_816->SetYTitle(yaxiswidtheta.Data());
  widthaway_eta_816->SetTitle(Form("Away-side peak width in #Delta#eta_{12} in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data())); 
  TH1D* widthaway_eta_1650 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayWidth_eta12_1650%s",extraname.Data())));
  widthaway_eta_1650->SetYTitle(yaxiswidtheta.Data());
  widthaway_eta_1650->SetTitle(Form("Away-side peak width in #Delta#eta_{12} in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data())); 
  
  TH1D* widthnear_theta_48 = dynamic_cast<TH1D*>(hist->Clone(Form("NearWidth_theta1_48%s",extraname.Data())));
  widthnear_theta_48->SetYTitle(yaxiswidththeta.Data());
  widthnear_theta_48->SetTitle(Form("Near-side peak width in #Delta#theta_{1} in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* widthnear_theta_816 = dynamic_cast<TH1D*>(hist->Clone(Form("NearWidth_theta1_816%s",extraname.Data())));
  widthnear_theta_816->SetYTitle(yaxiswidththeta.Data());
  widthnear_theta_816->SetTitle(Form("Near-side peak width in #Delta#theta_{1} in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data()));
  TH1D* widthnear_theta_1650 = dynamic_cast<TH1D*>(hist->Clone(Form("NearWidth_theta1_1650%s",extraname.Data())));
  widthnear_theta_1650->SetYTitle(yaxiswidththeta.Data());
  widthnear_theta_1650->SetTitle(Form("Near-side peak width in #Delta#theta_{1} in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data()));
  TH1D* widthaway_theta_48 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayWidth_theta1_48%s",extraname.Data())));
  widthaway_theta_48->SetYTitle(yaxiswidththeta.Data());
  widthaway_theta_48->SetTitle(Form("Away-side peak width in #Delta#theta_{1} in %s for 4GeV/c<=pT_{trigger}<=8GeV/c",titlestring.Data()));
  TH1D* widthaway_theta_816 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayWidth_theta1_816%s",extraname.Data())));
  widthaway_theta_816->SetYTitle(yaxiswidththeta.Data());
  widthaway_theta_816->SetTitle(Form("Away-side peak width in #Delta#theta_{1} in %s for 8GeV/c<=pT_{trigger}<=16GeV/c",titlestring.Data()));   
  TH1D* widthaway_theta_1650 = dynamic_cast<TH1D*>(hist->Clone(Form("AwayWidth_theta1_1650%s",extraname.Data())));
  widthaway_theta_1650->SetYTitle(yaxiswidththeta.Data());
  widthaway_theta_1650->SetTitle(Form("Away-side peak width in #Delta#theta_{1} in %s for 16GeV/c<=pT_{trigger}<=50GeV/c",titlestring.Data()));
  outdir->cd();
  for(int i = 0 ; i< indirs->GetEntries();i++){
    TDirectory * dir = dynamic_cast<TDirectory*>(indirs->At(i));
    TDirectory * dirPbPb = dynamic_cast<TDirectory*>(indirs->At(i))->GetDirectory(subdir.Data());
    TString bindis = TString(dir->GetTitle());
    TObjArray* binnra = bindis.Tokenize("_");
    int binedge = dynamic_cast<TObjString*>(binnra->At(5))->GetString().Atoi();
    delete binnra;
    int binnr = 0;
    double dptbin = 1 ;
    if(binedge ==2) {binnr = 1; dptbin = 1.0;}
    if(binedge ==3) {binnr = 2; dptbin = 1.0;}
    if(binedge ==4) {binnr = 3; dptbin = 2.0;}
    if(binedge ==6) {binnr = 4; dptbin = 2.0;}
    if(binedge ==8) {binnr = 5; dptbin = 8.0;}
    if(binedge ==16){binnr = 6; dptbin = 34.0;}

    if(bindis.Contains(ppperiod.Data())&&ispp){
      if(TString(dir->GetTitle()).Contains("er_4_8")){
	if(!methodstring.Contains("singleyield"))fillwithvalues(dir->GetDirectory(methodstring.Data()),binnr,nearyieldhist48,awayieldhist48,nearyield_f_hist48,awayield_f_hist48,widthnear_eta_48,widthaway_eta_48,widthnear_theta_48,widthaway_theta_48);
	else fillwithvalues(dir->GetDirectory(methodstring.Data()),binnr,nearyieldhist48,awayieldhist48,nearyield_f_hist48,awayield_f_hist48,widthnear_eta_48,widthaway_eta_48);
      }
      if(TString(dir->GetTitle()).Contains("er_8_16")){
	if(!methodstring.Contains("singleyield"))fillwithvalues(dir->GetDirectory(methodstring.Data()),binnr,nearyieldhist816,awayieldhist816,nearyield_f_hist816,awayield_f_hist816,widthnear_eta_816,widthaway_eta_816,widthnear_theta_816,widthaway_theta_816);
	else fillwithvalues(dir->GetDirectory(methodstring.Data()),binnr,nearyieldhist816,awayieldhist816,nearyield_f_hist816,awayield_f_hist816,widthnear_eta_816,widthaway_eta_816);
      }
      if(TString(dir->GetTitle()).Contains("er_16_")){
	if(!methodstring.Contains("singleyield"))fillwithvalues(dir->GetDirectory(methodstring.Data()),binnr,nearyieldhist1650,awayieldhist1650,nearyield_f_hist1650,awayield_f_hist1650,widthnear_eta_1650,widthaway_eta_1650,widthnear_theta_1650,widthaway_theta_1650);
	else fillwithvalues(dir->GetDirectory(methodstring.Data()),binnr,nearyieldhist1650,awayieldhist1650,nearyield_f_hist1650,awayield_f_hist1650,widthnear_eta_1650,widthaway_eta_1650);
      }
      
    }
    if(bindis.Contains("PbPb")&&!ispp){
      if(subdir.CompareTo("BinM(60.00)->(90.00)")==0)if(!TString(dir->GetTitle()).Contains("LHC10h")) continue;
      if(subdir.CompareTo("BinM(60.00)->(90.00)")!=0)if(!TString(dir->GetTitle()).Contains("LHC11h")) continue;
      if(TString(dir->GetTitle()).Contains("er_4_8")){
	if(!methodstring.Contains("singleyield"))fillwithvalues(dirPbPb->GetDirectory(methodstring.Data()),binnr,nearyieldhist48,awayieldhist48,nearyield_f_hist48,awayield_f_hist48,widthnear_eta_48,widthaway_eta_48,widthnear_theta_48,widthaway_theta_48);
	else fillwithvalues(dirPbPb->GetDirectory(methodstring.Data()),binnr,nearyieldhist48,awayieldhist48,nearyield_f_hist48,awayield_f_hist48,widthnear_eta_48,widthaway_eta_48);
      }
      if(TString(dir->GetTitle()).Contains("er_8_16")){
	if(!methodstring.Contains("singleyield"))fillwithvalues(dirPbPb->GetDirectory(methodstring.Data()),binnr,nearyieldhist816,awayieldhist816,nearyield_f_hist816,awayield_f_hist816,widthnear_eta_816,widthaway_eta_816,widthnear_theta_816,widthaway_theta_816);
	else fillwithvalues(dirPbPb->GetDirectory(methodstring.Data()),binnr,nearyieldhist816,awayieldhist816,nearyield_f_hist816,awayield_f_hist816,widthnear_eta_816,widthaway_eta_816);
      }
      if(TString(dir->GetTitle()).Contains("er_16_")){
	if(!methodstring.Contains("singleyield"))fillwithvalues(dirPbPb->GetDirectory(methodstring.Data()),binnr,nearyieldhist1650,awayieldhist1650,nearyield_f_hist1650,awayield_f_hist1650,widthnear_eta_1650,widthaway_eta_1650,widthnear_theta_1650,widthaway_theta_1650);
	else fillwithvalues(dirPbPb->GetDirectory(methodstring.Data()),binnr,nearyieldhist1650,awayieldhist1650,nearyield_f_hist1650,awayield_f_hist1650,widthnear_eta_1650,widthaway_eta_1650);
      }
      
    }    
  }
  TString suboutdir = TString("");
  if(subdir.CompareTo("")==0)suboutdir.Append("pp_min_bias");
  if(subdir.CompareTo("BinM(0.00)->(5.00)")==0)suboutdir.Append("PbPb_0-5");
  if(subdir.CompareTo("BinM(0.00)->(10.00)")==0)suboutdir.Append("PbPb_0-10");
  if(subdir.CompareTo("BinM(5.00)->(10.00)")==0)suboutdir.Append("PbPb_5-10");
  if(subdir.CompareTo("BinM(10.00)->(20.00)")==0)suboutdir.Append("PbPb_10-20");
  if(subdir.CompareTo("BinM(20.00)->(40.00)")==0)suboutdir.Append("PbPb_20-40");
  if(subdir.CompareTo("BinM(40.00)->(60.00)")==0)suboutdir.Append("PbPb_40-60");
  if(subdir.CompareTo("BinM(60.00)->(90.00)")==0)suboutdir.Append("PbPb_60-90");
  if(!outdir->GetDirectory(suboutdir.Data()))outdir->mkdir(suboutdir.Data())->cd();
  else outdir->GetDirectory(suboutdir.Data())->cd();
  nearyieldhist48->Write();
  awayieldhist48->Write();
  nearyieldhist816->Write();
  awayieldhist816->Write();
  nearyieldhist1650->Write();
  awayieldhist1650->Write();
  nearyield_f_hist48->Write();
  awayield_f_hist48->Write();
  nearyield_f_hist816->Write();
  awayield_f_hist816->Write();
  nearyield_f_hist1650->Write();
  awayield_f_hist1650->Write();
  widthnear_eta_48->Write();
  widthaway_eta_48->Write();
  widthnear_eta_816->Write();
  widthaway_eta_816->Write();
  widthnear_eta_1650->Write();
  widthaway_eta_1650->Write();
  if(!methodstring.Contains("singleyield")){
    widthnear_theta_48->Write();
    widthaway_theta_48->Write();
    widthnear_theta_816->Write();
    widthaway_theta_816->Write();
    widthnear_theta_1650->Write();
    widthaway_theta_1650->Write();
  }
  
}

void comparemethods(TObjArray * indirs, TDirectory* outdir, TString prods,TString typehist)
{
//compare methods, productions and periods
  TObjArray * arrayprods = prods.Tokenize("#");
  int nprods = arrayprods->GetEntries();
  bool ptcompare = false;
  bool centcompare = false;
  
  if(typehist.Contains("pT")||typehist.Contains("PT")||typehist.Contains("pt")){
    ptcompare = true;
  }
  if(typehist.Contains("cent")||typehist.Contains("Cent")){
    centcompare = true;
  }  
  
  int triggerpt = 0;
  if(typehist.Contains("t4"))triggerpt=4;
  if(typehist.Contains("t8"))triggerpt=8;
  if(typehist.Contains("t16"))triggerpt=16;
  
  TString parname = TString("");
  bool added = true;
  if(typehist.Contains("yield_near")){parname.Append("yield_near");added = false;}
  if(typehist.Contains("yield_away")&&added){parname.Append("yield_away");added = false;}
  if(typehist.Contains("widthThetanear")&&added){parname.Append("width_thetanear");added = false;}
  if(typehist.Contains("widthThetaaway")&&added){parname.Append("width_thetaaway");added = false;}
  if(typehist.Contains("widthetanear")&&added){parname.Append("width_etanear");added = false;}
  if(typehist.Contains("widthetaaway")&&added){parname.Append("width_etaaway");added = false;}
  if(added) return;

  int colors[6] = {1,2,3,4,6,7};
  //create pt histograms:
  
  TH1D * histpt  = new TH1D("histpt", "title",6,0.0,6.0);
  histpt->SetStats(false);
  histpt->GetXaxis()->SetBinLabel(1,"2.0-3.0");
  histpt->GetXaxis()->SetBinLabel(2,"3.0-4.0");
  histpt->GetXaxis()->SetBinLabel(3,"4.0-6.0");
  histpt->GetXaxis()->SetBinLabel(4,"6.0-8.0");
  histpt->GetXaxis()->SetBinLabel(5,"8.0-16.0");
  histpt->GetXaxis()->SetBinLabel(6,"16.0-50.0");
  histpt->GetXaxis()->SetTitle("Pt_{associated} [GeV/c]");
  TCanvas * ptcanvas = new TCanvas("ptCanvas");
  
  //create centrality histogram:
  TH1D * histcent  = new TH1D("histcent", "title",5,0.0,5.0);
  histcent->SetStats(false);
  histcent->GetXaxis()->SetBinLabel(1,"0.0-10.0");
  histcent->GetXaxis()->SetBinLabel(2,"10.0-20.0");
  histcent->GetXaxis()->SetBinLabel(3,"20.0-40.0");
  histcent->GetXaxis()->SetBinLabel(4,"40.0-60.0");
  histcent->GetXaxis()->SetBinLabel(5,"60.0-90.0");
  histcent->GetXaxis()->SetTitle("Centrality [%]");
  TCanvas * centcanvas = new TCanvas("centCanvas");

  TObjArray * pthists = new TObjArray();
  pthists->SetOwner(kTRUE);
  TObjArray * centhists = new TObjArray();
  centhists->SetOwner(kTRUE);
  
  for(int l = 1;l<=nprods;l++){
    pthists->Add(histpt->Clone(Form("histpt%i",nprods)));
    centhists->Add(histcent->Clone(Form("centhist%i",nprods)));
  }
  
  int NCentbins = 1;
  if(typehist.Contains("cloop")) NCentbins = 5;
  if(NCentbins>1&&typehist.Contains("Bin=")){cout << "Not possible to combine cloop and specific centrality bin."<<endl; return;}
  for(int centbin=1; centbin <=NCentbins; centbin++){
    for(int l = 1;l<=nprods;l++){
      dynamic_cast<TH1*>(pthists->At(l-1))->Reset();
      dynamic_cast<TH1*>(pthists->At(l-1))->SetTitle("");
      dynamic_cast<TH1*>(centhists->At(l-1))->Reset();
      dynamic_cast<TH1*>(centhists->At(l-1))->SetTitle("");
    }
    TString BinString = TString("");
    int nowcentbin = 0;
    if((typehist.Contains("Bin=1"))||((!typehist.Contains("Bin="))&&centbin==1)){BinString.Append("BinM(0.00)->(10.00)");nowcentbin=1;}
    if((typehist.Contains("Bin=2"))||((!typehist.Contains("Bin="))&&centbin==2)){BinString.Append("BinM(10.00)->(20.00)");nowcentbin=2;}
    if((typehist.Contains("Bin=3"))||((!typehist.Contains("Bin="))&&centbin==3)){BinString.Append("BinM(20.00)->(40.00)");nowcentbin=3;}
    if((typehist.Contains("Bin=4"))||((!typehist.Contains("Bin="))&&centbin==4)){BinString.Append("BinM(40.00)->(60.00)");nowcentbin=4;}
    if((typehist.Contains("Bin=5"))||((!typehist.Contains("Bin="))&&centbin==5)){BinString.Append("BinM(60.00)->(90.00)");nowcentbin=5;}
    
    
    double min=0.0;
    double max=0.0;
    for(int i = 1; i<=nprods;i++){
      TString prod = dynamic_cast<TObjString*>(arrayprods->At(i-1))->GetString();
      TString periodstr = TString("");
      if(prod.Contains("LHC10d"))periodstr.Append("ppLHC10d");
      if(prod.Contains("LHC10h"))periodstr.Append("PbPbLHC10h");
      if(prod.Contains("LHC11a"))periodstr.Append("ppLHC11a");
      if(prod.Contains("LHC11h"))periodstr.Append("PbPbLHC11h");
      
      TString subdirectory = TString("");
      TString title = TString("");
      if(prod.Contains("iteration")){
	subdirectory.Append("yield/");
	title.Append("Correction method 1, ");
      }
      if(prod.Contains("etasubstracted")){
	subdirectory.Append("Etasubstracted/yield/");	
	title.Append("Correction method 2, ");
      }
      
      if(prod.Contains("simple")){
	subdirectory.Append("simple/");
	title.Append("yield extracted in #phi slices.");
      }
      if(prod.Contains("single")){
	subdirectory.Append("singleyield/");
	title.Append("yield extracted in a single slice.");
      }
      
      dynamic_cast<TH1D*>(pthists->At(i-1))->SetTitle(title.Data());
      dynamic_cast<TH1D*>(centhists->At(i-1))->SetTitle(title.Data());

      TString BinStringloc = BinString;
      if(periodstr.Contains("pp"))BinStringloc.Clear();
      for(int k = 0;k<indirs->GetEntries();k++){
	TDirectory * dir = dynamic_cast<TDirectory*>(indirs->At(k));
	TString whichbin =  TString(dir->GetTitle());
	if(whichbin.Contains(periodstr.Data())){
	  TDirectory * locdir = dir->GetDirectory(BinStringloc.Data())->GetDirectory(subdirectory);
	  TObjArray * tokens = whichbin.Tokenize("_");
	  tokens->SetOwner(kTRUE);
  	  if(dynamic_cast<TObjString*>(tokens->At(2))->GetString().Atoi()==triggerpt)
	  {
	    int assbin = 0;	    
	    if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==2)assbin=1;
	    if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==3)assbin=2;
	    if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==4)assbin=3;
	    if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==6)assbin=4;
	    if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==8)assbin=5;
	    if(dynamic_cast<TObjString*>(tokens->At(5))->GetString().Atoi()==16)assbin=6;	
	    dynamic_cast<TH1D*>(pthists->At(i-1))->SetBinContent(assbin,GetVal(locdir,parname.Data()));
	    dynamic_cast<TH1D*>(pthists->At(i-1))->SetBinError(assbin,GetVal(locdir,Form("%s_error",parname.Data())));
	  }
	  delete tokens;
	}
      }
      getminmax(dynamic_cast<TH1D*>(pthists->At(i-1)),min,max);      
    }
    
    ptcanvas->Clear();
    ptcanvas->cd();
    TLegend * legend = new TLegend(0.7,0.75,0.9,0.9);
    for(int j = 1; j<=nprods;j++){
      TH1D* hist = dynamic_cast<TH1D*>(pthists->At(j-1));
      hist->SetLineColor(colors[j-1]);
      legend->AddEntry(hist,hist->GetTitle());
      if(j==1){
	hist->SetTitle("Comparison");
	if(min<0.0) min*=1.1; 
	if(min>=0.0&&max>0.0) min = 0.0 -0.2*max;
	if(min>=0.0&&max<=0.0) min = 0.0;
	if(max>0.0) max = max*1.1;	
	hist->GetYaxis()->SetRangeUser(min,max);
	hist->Draw("E");
      }
      else hist->Draw("same");
      
    }
    legend->Draw("same");
    
    TString picname = TString("Yields_Compared/");
    picname.Append(parname.Data());
    if(NCentbins!=1)picname.Append(BinString.Data());
    if(triggerpt==4)picname.Append("_pTtrigger_4_8.eps");
    if(triggerpt==8)picname.Append("_pTtrigger_8_16.eps");
    if(triggerpt==16)picname.Append("_pTtrigger_16_50.eps");
    ptcanvas->SaveAs(picname.Data());
  }
  delete pthists; delete centhists;
  delete histpt; delete histcent; delete ptcanvas; delete centcanvas;
  
}

  
void yieldcanvases(TObjArray * indirs, TDirectory* outdir,const char* histname,TString prodpp){
  //Compare different centralities to pp in one plot:
  //categorys that will be plotted:
  TObjArray * trigarray = new TObjArray();
  trigarray->Add(new TObjString(" 4.0- 8.0 GeV/c"));
  trigarray->Add(new TObjString(" 8.0-16.0 GeV/c"));
  trigarray->Add(new TObjString("16.0-50.0 GeV/c"));
  TObjArray * assarray  = new TObjArray();
  assarray->Add(new TObjString(" 2.0- 3.0 GeV/c"));
  assarray->Add(new TObjString(" 3.0- 4.0 GeV/c"));
  assarray->Add(new TObjString(" 4.0- 6.0 GeV/c"));
  assarray->Add(new TObjString(" 6.0- 8.0 GeV/c"));
  assarray->Add(new TObjString(" 8.0-16.0 GeV/c"));
  assarray->Add(new TObjString("16.0-50.0 GeV/c"));
  TString extraname = TString("");
  if(prodpp.CompareTo("ppLHC11a")==0) extraname.Append("11a");
  bool is010 = false;
  int i = 0;
  bool breakp = false;
  while(!breakp){
    if(dynamic_cast<TDirectory*>(indirs->At(i))->GetDirectory("BinM(0.00)->(10.00)"))is010 = true;
    if(TString(dynamic_cast<TDirectory*>(indirs->At(i))->GetTitle()).Contains("PbPb"))breakp = true;    
    i++;
  }
  TString objectfinder = TString("");
  TString titlehists   = TString("");
  if(TString(histname).CompareTo("dphiyieldbc")==0){
    objectfinder.Append("yield/simple/dphiyieldbc");
    titlehists.Append("Yields as a function of #Delta#phi_{1} in pp and");
  }
  if(TString(histname).CompareTo("near")==0) {
    objectfinder.Append("yield/singleyield/deta12_Near_bs");
    titlehists.Append("Near-side yield as a function of #Delta#eta_{12} in pp and");
  }
  if(TString(histname).CompareTo("away")==0){
    objectfinder.Append("yield/singleyield/deta12_Away_bs");    
    titlehists.Append("Away-side yield as a function of #Delta#eta_{12} in pp and");
  } 
  
  TCanvas * canvas05  = new TCanvas();
  TCanvas * canvas510 = new TCanvas();  
  TObjArray* arrayofpads05; TObjArray* arrayofpads510 ;    
    
  if(is010){
    arrayofpads05 = PadArray(canvas05,3,6,trigarray,assarray,"");//Form("%s 0%%-10%% PbPb",titlehists.Data()));
    arrayofpads510 = NULL;
    
  }
  else{
    arrayofpads05  = PadArray(canvas05,3,6,trigarray,assarray,"");//Form("%s 0%%-5%% PbPb",titlehists.Data()));
    arrayofpads510 = PadArray(canvas510,3,6,trigarray,assarray,"");//Form("%s 5%%-10%% PbPb",titlehists.Data()));
  }
  TCanvas * canvas1020 = new TCanvas();  
  TObjArray* arrayofpads1020 = PadArray(canvas1020,3,6,trigarray,assarray,"");//Form("%s  10%%-20%% PbPb",titlehists.Data()));
  TCanvas * canvas2040 = new TCanvas();
  TObjArray* arrayofpads2040 = PadArray(canvas2040,3,6,trigarray,assarray,"");//Form("%s  20%%-40%% PbPb",titlehists.Data()));
  TCanvas * canvas4060 = new TCanvas();
  TObjArray* arrayofpads4060 = PadArray(canvas4060,3,6,trigarray,assarray,"");//Form("%s  40%%-60%% PbPb",titlehists.Data()));
  TCanvas * canvas6090 = new TCanvas();
  TObjArray* arrayofpads6090 = PadArray(canvas6090,3,6,trigarray,assarray,"");//Form("%s  60%%-90%% PbPb",titlehists.Data()));
  
  TCanvas * canvas = new TCanvas();
  if((TString(histname).CompareTo("dphiyieldbc")==0)&&(prodpp.CompareTo("ppLHC10d")==0)) outdir->mkdir("yields_pp_PbPb")->cd();
  else outdir->GetDirectory("yields_pp_PbPb")->cd();
  TLegend * leg1 = new TLegend(0.1,0.1,0.9,0.9);
  TLegend * leg2 = new TLegend(0.1,0.1,0.9,0.9);
  TLegend * leg3 = new TLegend(0.1,0.1,0.9,0.9);
  TLegend * leg4 = new TLegend(0.1,0.1,0.9,0.9);
  TLegend * leg5 = new TLegend(0.1,0.1,0.9,0.9);
  TLegend * leg6 = new TLegend(0.1,0.1,0.9,0.9);
  for(int i = 0; i<indirs->GetEntries();i++){
    canvas->Clear();
    TDirectory * dir = dynamic_cast<TDirectory*>(indirs->At(i));
    TString bindis = TString(dir->GetTitle());
    TObjArray* binnra = bindis.Tokenize("_");
    int binedge = dynamic_cast<TObjString*>(binnra->At(5))->GetString().Atoi();
    delete binnra;
    TString binname = TString("");
    int bint = 0;
    if(TString(dir->GetTitle()).Contains("er_4_8")){binname.Append("Trigger_4_8");bint=4;}
    if(TString(dir->GetTitle()).Contains("er_8_16")){binname.Append("Trigger_8_16");bint=8;}
    if(TString(dir->GetTitle()).Contains("er_16_")){binname.Append("Trigger_16_50");bint=16;}
    if(binedge ==0) continue;
    if(binedge ==1) continue;
    if(binedge ==2) binname.Append("_Associated_2_3");
    if(binedge ==3) binname.Append("_Associated_3_4");
    if(binedge ==4) binname.Append("_Associated_4_6");
    if(binedge ==6) binname.Append("_Associated_6_8");
    if(binedge ==8) binname.Append("_Associated_8_16");
    if(binedge ==16) binname.Append("_Associated_16_50");
    if(bindis.Contains(prodpp.Data())){
      TH1D* yieldhistpp = dynamic_cast<TH1D*>(dir->Get(objectfinder.Data()));
      yieldhistpp->SetTitle("");
      yieldhistpp->SetStats(false);
      for(int j = 0; j<indirs->GetEntries();j++){
	TDirectory * dir1 = dynamic_cast<TDirectory*>(indirs->At(j));
	TString title2 = TString(dir1->GetTitle());
	if(!title2.Contains("PbPb"))continue;	
	TObjArray* binnra1 = title2.Tokenize("_");
	int binedge1 = dynamic_cast<TObjString*>(binnra1->At(5))->GetString().Atoi();
	delete binnra1;
	if(!(binedge==binedge1))continue;
	if(!((title2.Contains("er_4_8")&&bint==4)||(title2.Contains("er_8_16")&&bint==8)||(title2.Contains("er_16_50")&&bint==16)))continue;
	if(title2.Contains("LHC11h")){
	  TH1D* yieldhist05; 
	  TH1D* yieldhist510;
	  if(is010){
	      yieldhist05 	= dynamic_cast<TH1D*>(dir1->GetDirectory("BinM(0.00)->(10.00)")->Get(objectfinder.Data()));
	      yieldhist510 	= NULL;
	  }
	  else{
	    yieldhist05 	= dynamic_cast<TH1D*>(dir1->GetDirectory("BinM(0.00)->(5.00)")->Get(objectfinder.Data()));
	    yieldhist510 	= dynamic_cast<TH1D*>(dir1->GetDirectory("BinM(5.00)->(10.00)")->Get(objectfinder.Data()));
	  }

	  TH1D* yieldhist1020 	= dynamic_cast<TH1D*>(dir1->GetDirectory("BinM(10.00)->(20.00)")->Get(objectfinder.Data()));
	  TH1D* yieldhist2040 	= dynamic_cast<TH1D*>(dir1->GetDirectory("BinM(20.00)->(40.00)")->Get(objectfinder.Data()));
	  TH1D* yieldhist4060 	= dynamic_cast<TH1D*>(dir1->GetDirectory("BinM(40.00)->(60.00)")->Get(objectfinder.Data()));
	  yieldhist05->SetTitle(Form("0%%-5%% PbPb"));
	  TString titlestring = TString("Yield vs #phi with ");
	  if(bint==4)titlestring.Append(" 4 GeV/c#leq pT_{trigger}#leq 8 GeV/c and");
	  if(bint==8)titlestring.Append(" 8 GeV/c#leq pT_{trigger}#leq 16 GeV/c and");
	  if(bint==16)titlestring.Append(" 16 GeV/c#leq pT_{trigger}#leq 50 GeV/c and");
	  if(binedge==2)titlestring.Append(" 2 GeV/c#leq pT_{ass}#leq 3 GeV/c .");
	  if(binedge==3)titlestring.Append(" 3 GeV/c#leq pT_{ass}#leq 4 GeV/c .");
	  if(binedge==4)titlestring.Append(" 4 GeV/c#leq pT_{ass}#leq 6 GeV/c .");
	  if(binedge==6)titlestring.Append(" 6 GeV/c#leq pT_{ass}#leq 8 GeV/c .");
	  if(binedge==8)titlestring.Append(" 8 GeV/c#leq pT_{ass}#leq 16 GeV/c .");
	  if(binedge==8)titlestring.Append(" 16 GeV/c#leq pT_{ass}#leq 50 GeV/c .");
	  canvas->SetTitle(titlestring.Data());
	  double padarrayindex = 0;
	  if(bint==8)padarrayindex = 6;
	  if(bint==16)padarrayindex = 12;
	  if(binedge == 3) padarrayindex +=1;
	  if(binedge == 4) padarrayindex +=2;
	  if(binedge == 6) padarrayindex +=3;
	  if(binedge == 8) padarrayindex +=4;
	  if(binedge == 16) padarrayindex +=5;

	  if(padarrayindex==0){
	    TString xname = TString(Form("x-axis : %s",yieldhistpp->GetXaxis()->GetTitle()));
	    TString yname = TString(Form("y-axis : %s",yieldhistpp->GetYaxis()->GetTitle()));
	    xname.ReplaceAll("[]","");	    
	    if(is010){
	      dynamic_cast<TPad*>(arrayofpads05->At(4))->cd();
	      leg1->AddEntry(yieldhistpp,"pp");
	      leg1->AddEntry(yieldhist05,Form("0%%-10%% PbPb"));
	      leg1->Draw("same");	  
	      writeinpad(arrayofpads05,5,xname.Data());
	      writeinpad(arrayofpads05,11,yname.Data());
	    }
	    else{
	      dynamic_cast<TPad*>(arrayofpads05->At(4))->cd();
	      leg1->AddEntry(yieldhistpp,"pp");
	      leg1->AddEntry(yieldhist05,Form("0%%-5%% PbPb"));
	      leg1->Draw("same");
	      writeinpad(arrayofpads05,5,xname.Data());
	      writeinpad(arrayofpads05,11,yname.Data());
	   
	      dynamic_cast<TPad*>(arrayofpads510->At(4))->cd();
	      leg2->AddEntry(yieldhistpp,"pp");
	      leg2->AddEntry(yieldhist510,Form("5%%-10%% PbPb"));
	      leg2->Draw("same");	    
	      writeinpad(arrayofpads510,5,xname.Data());
	      writeinpad(arrayofpads510,11,yname.Data());
	      
	    }
	    dynamic_cast<TPad*>(arrayofpads1020->At(4))->cd();
	    leg3->Clear();  
	    leg3->AddEntry(yieldhistpp,"pp");
	    leg3->AddEntry(yieldhist1020,Form("10%%-20%% PbPb"));
	    leg3->Draw("same");	  
	    writeinpad(arrayofpads1020,5,xname.Data());
	    writeinpad(arrayofpads1020,11,yname.Data());
	    
	    dynamic_cast<TPad*>(arrayofpads2040->At(4))->cd();
	    leg4->Clear();  
	    leg4->AddEntry(yieldhistpp,"pp");
	    leg4->AddEntry(yieldhist2040,Form("20%%-40%% PbPb"));
	    leg4->Draw("same");
	    writeinpad(arrayofpads2040,5,xname.Data());
	    writeinpad(arrayofpads2040,11,yname.Data());
	    
	    dynamic_cast<TPad*>(arrayofpads4060->At(4))->cd();
	    leg5->Clear();  
	    leg5->AddEntry(yieldhistpp,"pp");
	    leg5->AddEntry(yieldhist4060,Form("40%%-60%% PbPb"));
	    leg5->Draw("same");	  
	    writeinpad(arrayofpads4060,5,xname.Data());
	    writeinpad(arrayofpads4060,11,yname.Data());	    
	    
	  }
	  double max = -1000;
	  double min = 1000;
	  yieldhistpp->SetLineColor(2);
	  if(yieldhistpp->GetFunction("fgs"))yieldhistpp->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	  yieldhistpp->GetYaxis()->SetLabelSize(0.04);
	  yieldhistpp->GetYaxis()->SetTitleSize(0.03);
	  yieldhistpp->GetYaxis()->SetTitleOffset(1.2);
	  yieldhistpp->GetXaxis()->SetLabelSize(0.04);
	  yieldhistpp->GetXaxis()->SetTitleSize(0.03);
	  yieldhistpp->GetXaxis()->SetTitleOffset(1.2);
	  getminmax(yieldhist05,min, max);
	  getminmax(yieldhistpp,min, max);
	  min = min - 0.1*max;
	  max = max*1.5;
	  canvas->cd();
	  yieldhistpp->GetYaxis()->SetRangeUser(min,max);
	  yieldhistpp->Draw("E");

	  yieldhist05->SetLineColor(4);
	  if(yieldhist05->GetFunction("fgs"))yieldhist05->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	  yieldhist05->GetYaxis()->SetRangeUser(min,max);
	  yieldhist05->Draw("Same");
	  TLegend * leg = new TLegend(0.7,0.7,0.9,0.9);
	  leg->AddEntry(yieldhistpp,"pp");
	  if(!is010)leg->AddEntry(yieldhist05,Form("0%%-5%% PbPb"));
	  if(is010)leg->AddEntry(yieldhist05,Form("0%%-10%% PbPb"));
	  leg->Draw("same");
	  canvas->SetTitle(titlestring);
	  if(!is010)canvas->Write(TString(histname)+TString("_")+binname+TString("_0_5")+extraname);
	  if(is010) canvas->Write(TString(histname)+TString("_")+binname+TString("_0_10")+extraname);
	  canvas->Clear();
	  dynamic_cast<TPad*>(arrayofpads05->At(padarrayindex))->cd();

	  yieldhistpp->Draw("E");
	  yieldhist05->Draw("Same");
	  
	  if(!is010){
	    canvas->cd();
	    min = 1000;
	    max =-1000;
	    getminmax(yieldhistpp,min, max);
	    getminmax(yieldhist510,min, max);
	    min = min - 0.1*max;
	    max = max*1.5;
	    yieldhistpp->GetYaxis()->SetRangeUser(min,max);
	    yieldhistpp->Draw("E");	  
	    yieldhist510->SetLineColor(4);
	     if(yieldhist510->GetFunction("fgs"))yieldhist510->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	    yieldhist510->Draw("Same");
	    leg->Clear();
	    leg->AddEntry(yieldhistpp,"pp");
	    leg->AddEntry(yieldhist510,Form("5%%-10%% PbPb"));
	    leg->Draw("same");
	    canvas->SetTitle(titlestring);
	    canvas->Write(TString(histname)+TString("_")+binname+TString("_5_10")+extraname);
	    canvas->Clear();
	    
	    dynamic_cast<TPad*>(arrayofpads510->At(padarrayindex))->cd();
	    yieldhistpp->Draw("E");
	    yieldhist510->Draw("Same");
	  }
	  
	  canvas->cd();
	  min = 1000;
	  max =-1000;
	  getminmax(yieldhistpp,min, max);
	  getminmax(yieldhist1020,min, max);
	  min = min - 0.1*max;
	  max = max*1.5;
	  yieldhistpp->GetYaxis()->SetRangeUser(min,max);
	  yieldhistpp->Draw("E");
	  yieldhist1020->SetLineColor(4);
	  if(yieldhist1020->GetFunction("fgs"))yieldhist1020->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	  yieldhist1020->Draw("Same");
	  leg->Clear();
	  leg->AddEntry(yieldhistpp,"pp");
	  leg->AddEntry(yieldhist1020,Form("10%%-20%% PbPb"));
	  leg->Draw("same");
	  canvas->SetTitle(titlestring);
	  canvas->Write(TString(histname)+TString("_")+binname+TString("_10_20")+extraname);
	  canvas->Clear();	
	  
	  dynamic_cast<TPad*>(arrayofpads1020->At(padarrayindex))->cd();
	  yieldhistpp->Draw("E");
	  yieldhist1020->Draw("same");

	  canvas->cd();
	  min = 1000;
	  max =-1000;
	  getminmax(yieldhistpp,min, max);
	  getminmax(yieldhist2040,min, max);
	  min = min - 0.1*max;
	  max = max*1.5;
	  yieldhistpp->GetYaxis()->SetRangeUser(min,max);
	  yieldhistpp->Draw("E");
	  yieldhist2040->SetLineColor(4);
	  if(yieldhist2040->GetFunction("fgs"))yieldhist2040->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	  yieldhist2040->Draw("Same");
	  leg->Clear();
	  leg->AddEntry(yieldhistpp,"pp");
	  leg->AddEntry(yieldhist2040,Form("20%%-40%% PbPb"));
	  leg->Draw("same");
	  canvas->SetTitle(titlestring);
	  canvas->Write(TString(histname)+TString("_")+binname+TString("_20_40")+extraname);
	  canvas->Clear();

	  dynamic_cast<TPad*>(arrayofpads2040->At(padarrayindex))->cd();
	  yieldhistpp->Draw("E");
	  yieldhist2040->Draw("Same");

	  canvas->cd();
	  min = 1000;
	  max =-1000;
	  getminmax(yieldhistpp,min, max);
	  getminmax(yieldhist2040,min, max);
	  min = min - 0.1*max;
	  max = max*1.5;
	  yieldhistpp->GetYaxis()->SetRangeUser(min,max);
	  yieldhistpp->Draw("E");
	  yieldhist4060->SetLineColor(4);
	  if(yieldhist4060->GetFunction("fgs"))yieldhist4060->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	  yieldhist4060->Draw("Same");
	  leg->Clear();
	  leg->AddEntry(yieldhistpp,"pp");
	  leg->AddEntry(yieldhist4060,Form("40%%-60%% PbPb"));
	  leg->Draw("same");
	  canvas->SetTitle(titlestring);
	  canvas->Write(TString(histname)+TString("_")+binname+TString("_40_60")+extraname);
	  canvas->Clear();

	  dynamic_cast<TPad*>(arrayofpads4060->At(padarrayindex))->cd();
	  yieldhistpp->Draw("E");
	  yieldhist4060->Draw("Same");
	  
	  delete leg;
	}
	if(title2.Contains("LHC10h")){

	  TH1D* yieldhist6090 	= dynamic_cast<TH1D*>( dir1->GetDirectory("BinM(60.00)->(90.00)/")->Get(objectfinder.Data())  );
	  TString titlestring = TString("Yield vs #phi with ");
	  if(bint==4)titlestring.Append(" 4 GeV/c#leq pT_{trigger}#leq 8 GeV/c and");
	  if(bint==8)titlestring.Append(" 8 GeV/c#leq pT_{trigger}#leq 16 GeV/c and");
	  if(bint==16)titlestring.Append(" 16 GeV/c#leq pT_{trigger}#leq 50 GeV/c and");
	  if(binedge==2)titlestring.Append(" 2 GeV/c#leq pT_{ass}#leq 3 GeV/c .");
	  if(binedge==3)titlestring.Append(" 3 GeV/c#leq pT_{ass}#leq 4 GeV/c .");
	  if(binedge==4)titlestring.Append(" 4 GeV/c#leq pT_{ass}#leq 6 GeV/c .");
	  if(binedge==6)titlestring.Append(" 6 GeV/c#leq pT_{ass}#leq 8 GeV/c .");
	  if(binedge==8)titlestring.Append(" 8 GeV/c#leq pT_{ass}#leq 16 GeV/c .");
	  if(binedge==8)titlestring.Append(" 16 GeV/c#leq pT_{ass}#leq 50 GeV/c .");
	  canvas->SetTitle(titlestring.Data());
	  
	  double padarrayindex = 0;
	  if(bint==8)padarrayindex = 6;
	  if(bint==16)padarrayindex = 12;
	  if(binedge == 3) padarrayindex +=1;
	  if(binedge == 4) padarrayindex +=2;
	  if(binedge == 6) padarrayindex +=3;
	  if(binedge == 8) padarrayindex +=4;
	  if(binedge == 16) padarrayindex +=5;

	  if(padarrayindex==0){
	    TPaveText * axislablesx = new TPaveText(0.1,0.1,0.9,0.9,"NB");
	    TPaveText * axislablesy = new TPaveText(0.1,0.1,0.9,0.9,"NB");
	    TString xname = TString(Form("x-axis : %s",yieldhistpp->GetXaxis()->GetTitle()));
	    TString yname = TString(Form("y-axis : %s",yieldhistpp->GetYaxis()->GetTitle()));
	    xname.ReplaceAll("[]","");
	    dynamic_cast<TPad*>(arrayofpads6090->At(4))->cd();
	    leg6->Clear();  
	    leg6->AddEntry(yieldhistpp,"pp");
	    leg6->AddEntry(yieldhist6090,Form("60%%-90%% PbPb"));
	    leg6->Draw("same");	 
	    writeinpad(arrayofpads6090,5,xname.Data());
	    writeinpad(arrayofpads6090,11,yname.Data());
	    
	  }
	  double max = -1000;
	  double min = 1000;
	  
	  canvas->cd();
	  getminmax(yieldhistpp,min,max);
	  getminmax(yieldhist6090,min,max);
	  
	  min = min - 0.1*max;
	  max = max*1.5;
  
	  yieldhistpp->GetYaxis()->SetRangeUser(min,max);
	  yieldhistpp->Draw("E");
	  yieldhist6090->SetLineColor(4);
	  if(yieldhist6090->GetFunction("fgs"))yieldhist6090->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	  yieldhist6090->Draw("Same");
	  TLegend * leg = new TLegend(0.7,0.7,0.9,0.9);
	  leg->Clear();
	  leg->AddEntry(yieldhistpp,"pp");
	  leg->AddEntry(yieldhist6090,Form("60%%-90%% PbPb"));
	  leg->Draw("same");
	  canvas->SetTitle(titlestring);
	  canvas->Write(TString(histname)+TString("_")+binname+TString("_60_90")+extraname);
	  canvas->Clear();
	  
	  dynamic_cast<TPad*>(arrayofpads6090->At(padarrayindex))->cd();
	  yieldhistpp->Draw("E");
	  yieldhist6090->Draw("Same");
	  
	  delete leg;
	  
	}
      }
    }
  }  
  TObjArray * dirarray = TString(gSystem->pwd()).Tokenize("/");

  if((TString(histname).CompareTo("dphiyieldbc")==0)&&(prodpp.CompareTo("ppLHC10d")==0)){
    outdir->mkdir("yields_pp_PbPb_OnePlot")->cd();
    gSystem->mkdir("Yields_Compared");
    
    gSystem->mkdir(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared",dirarray->At(dirarray->GetEntries()-1)->GetName()));

  }
  else outdir->GetDirectory("yields_pp_PbPb_OnePlot")->cd();

  
  if(!is010){
    canvas05->Write(Form("%s_pp_PbPb_0_5",histname));
    canvas05->SaveAs(Form("Yields_Compared/%s_pp%s_PbPb_0_5.eps",histname,extraname.Data()));
    canvas05->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared/%s_pp%s_PbPb_0_5.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),histname,extraname.Data()));
    canvas510->Write(Form("%s_pp%s_PbPb_5_10",histname,extraname.Data()));
    canvas510->SaveAs(Form("Yields_Compared/%s_pp%s_PbPb_5_10.eps",histname,extraname.Data()));
    canvas510->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared/%s_pp%s_PbPb_5_10.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),histname,extraname.Data()));
  }
  else {
    canvas05->Write(Form("%s_pp%s_PbPb_0_10",histname,extraname.Data()));
    canvas05->SaveAs(Form("Yields_Compared/%s_pp%s_PbPb_0_10.eps",histname,extraname.Data()));
    canvas05->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared/%s_pp%s_PbPb_0_10.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),histname,extraname.Data()));
  }
  canvas1020->Write(Form("%s_pp%s_PbPb_10_20",histname,extraname.Data()));
  canvas1020->SaveAs(Form("Yields_Compared/%s_pp%s_PbPb_10_20.eps",histname,extraname.Data()));
  canvas1020->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared/%s_pp%s_PbPb_10_20.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),histname,extraname.Data()));
  canvas2040->Write(Form("%s_pp%s_PbPb_20_40",histname,extraname.Data()));
  canvas2040->SaveAs(Form("Yields_Compared/%s_pp%s_PbPb_20_40.eps",histname,extraname.Data()));
  canvas2040->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared/%s_pp%s_PbPb_20_40.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),histname,extraname.Data()));
  canvas4060->Write(Form("%s_pp%s_PbPb_40_60",histname,extraname.Data()));
  canvas4060->SaveAs(Form("Yields_Compared/%s_pp%s_PbPb_40_60.eps",histname,extraname.Data()));
  canvas4060->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared/%s_pp%s_PbPb_40_60.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),histname,extraname.Data()));
  canvas6090->Write(Form("%s_pp%s_PbPb_60_90",histname,extraname.Data()));
  canvas6090->SaveAs(Form("Yields_Compared/%s_pp%s_PbPb_60_90.eps",histname,extraname.Data()));
  canvas6090->SaveAs(Form("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/%s/Yields_Compared/%s_pp%s_PbPb_60_90.eps",dirarray->At(dirarray->GetEntries()-1)->GetName(),histname,extraname.Data()));
  delete canvas;
  delete leg1; delete leg2; delete leg3; delete leg4; delete leg5; delete leg6; 
  delete canvas05;delete canvas510;delete canvas1020;delete canvas2040;delete canvas4060; delete canvas6090;
  delete arrayofpads05;
  if(!is010)delete arrayofpads510;
  delete arrayofpads1020;delete arrayofpads2040;delete arrayofpads4060;delete arrayofpads6090;
  delete dirarray;
}
  
void yieldcanvasespp(TObjArray * indirs, TDirectory* outdir, TString prod1 ,TString prod2 ){
  //Compare two different pps:
  
  //categories that will be plotted:
  TObjArray * trigarray = new TObjArray();
  trigarray->Add(new TObjString("4 GeV/c <= p_{T}^{trigger}<=8 GeV/c"));
  trigarray->Add(new TObjString("8 GeV/c <= p_{T}^{trigger}<=16 GeV/c"));
  TObjArray * assarray  = new TObjArray();
//   assarray->Add(new TObjString("0.5-1.0 GeV/c"));
//   assarray->Add(new TObjString("1.0-2.0 GeV/c"));
  assarray->Add(new TObjString("2.0-3.0 GeV/c"));
  assarray->Add(new TObjString("3.0-4.0 GeV/c"));
  assarray->Add(new TObjString("4.0-6.0 GeV/c"));
  assarray->Add(new TObjString("6.0-8.0 GeV/c"));
  assarray->Add(new TObjString("8.0-16.0 GeV/c"));
  
  TCanvas * canvascom  = new TCanvas();
  TObjArray* arrayofpads = PadArray(canvascom,2,5,trigarray,assarray,Form("Yields in %s and %s",prod1.Data(),prod2.Data()));

  TCanvas * canvas = new TCanvas();
  outdir->mkdir("yields_pp")->cd();
  
  TLegend * leg1 = new TLegend(0.1,0.1,0.9,0.9);


  for(int i = 0; i<indirs->GetEntries();i++){
    canvas->Clear();
    TDirectory * dir = dynamic_cast<TDirectory*>(indirs->At(i));
    TString bindis = TString(dir->GetTitle());
    TObjArray* binnra = bindis.Tokenize("_");
    int binedge = dynamic_cast<TObjString*>(binnra->At(5))->GetString().Atoi();
    delete binnra;
    TString binname = TString("");
    int bint = 0;
    if(TString(dir->GetTitle()).Contains("er_4_8")){binname.Append("Trigger_4_8");bint=4;}
    if(TString(dir->GetTitle()).Contains("er_8_16")){binname.Append("Trigger_8_16");bint=8;}
    if(binedge ==0) continue;//binname.Append("_Associated_0.5_1");
    if(binedge ==1) continue;//binname.Append("_Associated_1_2");
    if(binedge ==2) binname.Append("_Associated_2_3");
    if(binedge ==3) binname.Append("_Associated_3_4");
    if(binedge ==4) binname.Append("_Associated_4_6");
    if(binedge ==6) binname.Append("_Associated_6_8");
    if(binedge ==8) binname.Append("_Associated_8_16");
    
    if(bindis.Contains("pp")&&!bindis.Contains("pp2")){
      TH1D* yieldhistpp = dynamic_cast<TH1D*>(dir->GetDirectory("yield/simple/")->Get("dphiyieldbc"));
      yieldhistpp->SetTitle(prod1.Data());
      for(int j = 0; j<indirs->GetEntries();j++){
	TDirectory * dir1 = dynamic_cast<TDirectory*>(indirs->At(j));
	TString title2 = TString(dir1->GetTitle());
	if(!title2.Contains("pp2"))continue;	
	TObjArray* binnra1 = title2.Tokenize("_");
	int binedge1 = dynamic_cast<TObjString*>(binnra1->At(5))->GetString().Atoi();
	delete binnra1;
	if(!(binedge==binedge1))continue;
	if(!((title2.Contains("er_4_8")&&bint==4)||(title2.Contains("er_8_16")&&bint==8)))continue;
	TH1D* yieldhistpp2 	= dynamic_cast<TH1D*>(dir1->GetDirectory("yield/simple/")->Get("dphiyieldbc"));
	yieldhistpp2->SetTitle(prod2.Data());
	TString titlestring = TString("Yield vs #phi with ");
	if(bint==4)titlestring.Append(" 4 GeV/c#leq pT_{trigger}#leq 8 GeV/c and");
	if(bint==8)titlestring.Append(" 8 GeV/c#leq pT_{trigger}#leq 16 GeV/c and");
	if(binedge==2)titlestring.Append(" 2 GeV/c#leq pT_{ass}#leq 3 GeV/c .");
	if(binedge==3)titlestring.Append(" 3 GeV/c#leq pT_{ass}#leq 4 GeV/c .");
	if(binedge==4)titlestring.Append(" 4 GeV/c#leq pT_{ass}#leq 6 GeV/c .");
	if(binedge==6)titlestring.Append(" 6 GeV/c#leq pT_{ass}#leq 8 GeV/c .");
	if(binedge==8)titlestring.Append(" 8 GeV/c#leq pT_{ass}#leq 16 GeV/c .");
	canvas->SetTitle(titlestring.Data());
	
	double padarrayindex = 0;
	if(bint==8)padarrayindex = 5;
	if(binedge == 3) padarrayindex +=1;
	if(binedge == 4) padarrayindex +=2;
	if(binedge == 6) padarrayindex +=3;
	if(binedge == 8) padarrayindex +=4;

	if(padarrayindex==0){
	  dynamic_cast<TPad*>(arrayofpads->At(4))->cd();
	  leg1->AddEntry(yieldhistpp,prod1.Data());
	  leg1->AddEntry(yieldhistpp2,prod2.Data());
	  leg1->Draw("same");	  
	}
	double max = 0;
	double min = 0;
	
	yieldhistpp->SetLineColor(2);
	yieldhistpp->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	yieldhistpp->GetYaxis()->SetLabelSize(0.04);
	yieldhistpp->GetYaxis()->SetTitleSize(0.03);
	yieldhistpp->GetYaxis()->SetTitleOffset(1.2);
	yieldhistpp->GetXaxis()->SetLabelSize(0.04);
	yieldhistpp->GetXaxis()->SetTitleSize(0.03);
	yieldhistpp->GetXaxis()->SetTitleOffset(1.2);
	yieldhistpp->SetTitle(titlestring);
	max = yieldhistpp->GetBinContent(yieldhistpp->FindBin(0.0));
	if(yieldhistpp2->GetBinContent(yieldhistpp2->FindBin(0.0))>max)max = yieldhistpp2->GetBinContent(yieldhistpp2->FindBin(0.0));
	min = yieldhistpp->GetMinimum();
	if(yieldhistpp2->GetMinimum()>min)min = yieldhistpp2->GetMinimum();
	min = min - 0.1*max;
	max = max*1.5;
	canvas->cd();
	canvas->DrawFrame(-TMath::Pi()/2.0,min,3.0*TMath::Pi()/2.0, max);
	yieldhistpp->Draw("Same");
	
	yieldhistpp2->SetLineColor(4);
	yieldhistpp2->GetFunction("fgs")->SetBit(TF1::kNotDraw);
	yieldhistpp2->Draw("Same");
	TLegend * leg = new TLegend(0.7,0.7,0.9,0.9);
	leg->AddEntry(yieldhistpp,prod1.Data());
	leg->AddEntry(yieldhistpp2,prod2.Data());
	leg->Draw("same");
	canvas->SetTitle(titlestring);
	canvas->Write(binname);
	canvas->Clear();
	
	dynamic_cast<TPad*>(arrayofpads->At(padarrayindex))->cd();
	dynamic_cast<TPad*>(arrayofpads->At(padarrayindex))->DrawFrame(-TMath::Pi()/2.0,min,3.0*TMath::Pi()/2.0, max);
	yieldhistpp->SetAxisRange(min,max,"y");
	yieldhistpp->Draw("Same");
	yieldhistpp2->Draw("Same");
	
	delete leg;
      }
    }
  }  
  outdir->mkdir("yields_pp_OnePlot")->cd();

  canvascom->Write("yields_pp_compare");
  
  delete canvas;
  delete leg1; 
  delete canvascom;
  delete arrayofpads;
}
  
void CompareYields(TObjArray* indirs, TDirectory * outdir , TString options){
  //Create the relevant histograms and save them, collect the relevant values from the indirs.
  bool noPbPb = false;
  if(options.Contains("ppcompare"))noPbPb = true;
  
  //generally, we want a histogram with associated energy on the x axis:
  TH1D * hist  = new TH1D("hist", "title",6,0.0,6.0);
  hist->SetStats(false);
  hist->GetXaxis()->SetBinLabel(1,"2.0-3.0");
  hist->GetXaxis()->SetBinLabel(2,"3.0-4.0");
  hist->GetXaxis()->SetBinLabel(3,"4.0-6.0");
  hist->GetXaxis()->SetBinLabel(4,"6.0-8.0");
  hist->GetXaxis()->SetBinLabel(5,"8.0-16.0");
  hist->GetXaxis()->SetBinLabel(6,"16.0-50.0");
  hist->GetXaxis()->SetTitle("P_{T}^{associated} [GeV/c]");
  
  if(!noPbPb){
    yieldcanvases(indirs,outdir,"dphiyieldbc","ppLHC10d");
    yieldcanvases(indirs,outdir,"near","ppLHC10d");
    yieldcanvases(indirs,outdir,"away","ppLHC10d");
    yieldcanvases(indirs,outdir,"dphiyieldbc","ppLHC11a");
    yieldcanvases(indirs,outdir,"near","ppLHC11a");
    yieldcanvases(indirs,outdir,"away","ppLHC11a");
    
    bool is010 = false;
    int i = 0;
    bool breakp = false;
    while(!breakp){
      if(dynamic_cast<TDirectory*>(indirs->At(i))->GetDirectory("BinM(0.00)->(10.00)"))is010 = true;
      if(TString(dynamic_cast<TDirectory*>(indirs->At(i))->GetTitle()).Contains("PbPb"))breakp = true;    
      i++;
    }    
    
    
    extractvalues(indirs,outdir,"",hist,"yield/simple","ppLHC10d");
    extractvalues(indirs,outdir,"",hist,"yield/simple","ppLHC11a");
    extractvalues(indirs,outdir,"",hist,"yield/singleyield","ppLHC10d");
    extractvalues(indirs,outdir,"",hist,"yield/singleyield","ppLHC11a");
    extractvalues(indirs,outdir,"",hist,"Etasubstracted/yield/simple","ppLHC10d");
    extractvalues(indirs,outdir,"",hist,"Etasubstracted/yield/simple","ppLHC11a");
    extractvalues(indirs,outdir,"",hist,"Etasubstracted/yield/singleyield","ppLHC10d");
    extractvalues(indirs,outdir,"",hist,"Etasubstracted/yield/singleyield","ppLHC11a");    
    
    if(is010){
      extractvalues(indirs,outdir,"BinM(0.00)->(10.00)",hist,"yield/simple");
      extractvalues(indirs,outdir,"BinM(0.00)->(10.00)",hist,"yield/singleyield");
      extractvalues(indirs,outdir,"BinM(0.00)->(10.00)",hist,"Etasubstracted/yield/simple");
      extractvalues(indirs,outdir,"BinM(0.00)->(10.00)",hist,"Etasubstracted/yield/singleyield");
      
    }
    else{
      extractvalues(indirs,outdir,"BinM(0.00)->(5.00)",hist,"yield/simple");
      extractvalues(indirs,outdir,"BinM(0.00)->(5.00)",hist,"yield/singleyield");
      extractvalues(indirs,outdir,"BinM(0.00)->(5.00)",hist,"Etasubstracted/yield/simple");
      extractvalues(indirs,outdir,"BinM(0.00)->(5.00)",hist,"Etasubstracted/yield/singleyield");
      extractvalues(indirs,outdir,"BinM(5.00)->(10.00)",hist,"yield/simple");    
      extractvalues(indirs,outdir,"BinM(5.00)->(10.00)",hist,"yield/singleyield");
      extractvalues(indirs,outdir,"BinM(5.00)->(10.00)",hist,"Etasubstracted/yield/simple");
      extractvalues(indirs,outdir,"BinM(5.00)->(10.00)",hist,"Etasubstracted/yield/singleyield");
    }
    extractvalues(indirs,outdir,"BinM(10.00)->(20.00)",hist,"yield/simple");
    extractvalues(indirs,outdir,"BinM(10.00)->(20.00)",hist,"yield/singleyield");
    extractvalues(indirs,outdir,"BinM(10.00)->(20.00)",hist,"Etasubstracted/yield/simple");
    extractvalues(indirs,outdir,"BinM(10.00)->(20.00)",hist,"Etasubstracted/yield/singleyield");
    extractvalues(indirs,outdir,"BinM(20.00)->(40.00)",hist,"yield/simple");
    extractvalues(indirs,outdir,"BinM(20.00)->(40.00)",hist,"yield/singleyield");
    extractvalues(indirs,outdir,"BinM(20.00)->(40.00)",hist,"Etasubstracted/yield/simple");
    extractvalues(indirs,outdir,"BinM(20.00)->(40.00)",hist,"Etasubstracted/yield/singleyield");    
    extractvalues(indirs,outdir,"BinM(40.00)->(60.00)",hist,"yield/simple");
    extractvalues(indirs,outdir,"BinM(40.00)->(60.00)",hist,"yield/singleyield");
    extractvalues(indirs,outdir,"BinM(40.00)->(60.00)",hist,"Etasubstracted/yield/simple");
    extractvalues(indirs,outdir,"BinM(40.00)->(60.00)",hist,"Etasubstracted/yield/singleyield");
    extractvalues(indirs,outdir,"BinM(60.00)->(90.00)",hist,"yield/simple");
    extractvalues(indirs,outdir,"BinM(60.00)->(90.00)",hist,"yield/singleyield");
    extractvalues(indirs,outdir,"BinM(60.00)->(90.00)",hist,"Etasubstracted/yield/simple");
    extractvalues(indirs,outdir,"BinM(60.00)->(90.00)",hist,"Etasubstracted/yield/singleyield");
  
  }
  else{
    TObjArray* opar = options.Tokenize("_");
    TString ppprod = TString(dynamic_cast<TObjString*>(opar->At(1))->GetString().Data());
    TString pp2prod = TString(dynamic_cast<TObjString*>(opar->At(2))->GetString().Data());
    delete opar;
    yieldcanvasespp(indirs,outdir,ppprod,pp2prod);

  }
//   outdir->cd();
//   hist->Write();

}
void DrawHists(TH1* hist1 , TH1* hist2, TDirectory* outdir, TString name, TString name2, TString dropboxname ){
  TCanvas * canvas = new TCanvas("c1");
  double min = 1000;
  double max = -1000;
  TLegend * leg = new TLegend(0.8,0.8,0.9,0.9);
  leg->AddEntry(hist1,"PbPb");
  leg->AddEntry(hist2,"pp");
  getminmax(hist1,min,max);
  getminmax(hist2,min,max);
  if(min>0.0) min = -max*0.1;
  else min =-max*0.1;
  max *=1.1;
  hist1->SetAxisRange(min,max,"Y");
  hist1->GetXaxis()->SetTitleOffset(1.15);
  hist1->GetYaxis()->SetLabelSize(0.025);
  hist1->GetYaxis()->SetTitleOffset(1.2);
  hist1->GetYaxis()->SetTitleSize(0.035);
  canvas->cd();
  
  hist1->Draw("");
  hist2->Draw("same");
  leg->Draw("same");
  outdir->cd();
  canvas->Write(name.Data());
  canvas->SaveAs(dropboxname.Data());
  
  if(name2.CompareTo("")){
    canvas->Clear();
    TH1D* hist1c = dynamic_cast<TH1D*>(hist1->Clone("clone"));
    hist1c->Divide(hist2);
    hist1c->GetYaxis()->UnZoom();
    hist1c->Draw();
    outdir->cd();
    canvas->Write(name2.Data());
    delete hist1c;
  }
  delete canvas;
  delete leg;
}

void PlotHists(TDirectory * indir, TDirectory* ppdir, TDirectory* outdir,TObjArray* hists,TString methodstring, TString ppperiod,const char* options){
  cout << "";
  if(!indir)return;
  TString extranamePbPb = TString("");
  TString extraname = TString("");
  if(ppperiod.Contains("LHC11a"))		extraname.Append("_11a");
  if(methodstring.Contains("Eta"))		extraname.Append("_etasubstracted");
  if(methodstring.Contains("singleyield"))	extraname.Append("_single");  
  if(methodstring.Contains("Eta"))		extranamePbPb.Append("_etasubstracted");
  if(methodstring.Contains("singleyield"))	extranamePbPb.Append("_single");  
  
  TH1D* yieldn48 = dynamic_cast<TH1D*>(indir->Get(Form("NearYield48%s",extranamePbPb.Data())));
  yieldn48->SetLineColor(kBlue);
  yieldn48->SetTitle(Form("Near-Side peak yield for 4GeV/C #leq pT_{trigger}#leq8GeV/c"));
  TH1D* yieldn48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearYield48%s",extraname.Data())));
  yieldn48pp->SetLineColor(kRed);
  TH1D* yielda48 = dynamic_cast<TH1D*>(indir->Get(Form("AwayYield48%s",extranamePbPb.Data())));
  yielda48->SetLineColor(kBlue);
  yielda48->SetTitle(Form("Away-Side peak yield for 4GeV/C #leq pT_{trigger}#leq8GeV/c"));
  TH1D* yielda48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayYield48%s",extraname.Data())));
  yielda48pp->SetLineColor(kRed);  

  TH1D* yieldn816 = dynamic_cast<TH1D*>(indir->Get(Form("NearYield816%s",extranamePbPb.Data())));
  yieldn816->SetLineColor(kBlue);
  yieldn816->SetTitle(Form("Near-Side peak yield for 8GeV/C #leq pT_{trigger}#leq16GeV/c"));
  TH1D* yieldn816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearYield816%s",extraname.Data())));
  yieldn816pp->SetLineColor(kRed);
  TH1D* yielda816 = dynamic_cast<TH1D*>(indir->Get(Form("AwayYield816%s",extranamePbPb.Data())));
  yielda816->SetLineColor(kBlue);
  yielda816->SetTitle(Form("Away-Side peak yield for 8GeV/C #leq pT_{trigger}#leq16GeV/c"));
  TH1D* yielda816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayYield816%s",extraname.Data())));
  yielda816pp->SetLineColor(kRed);  

  TH1D* yieldn16 = dynamic_cast<TH1D*>(indir->Get(Form("NearYield1650%s",extranamePbPb.Data())));
  yieldn16->SetLineColor(kBlue);
  yieldn16->SetTitle(Form("Near-Side peak yield for 16GeV/C #leq pT_{trigger}#leq50GeV/c"));
  TH1D* yieldn16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearYield1650%s",extraname.Data())));
  yieldn16pp->SetLineColor(kRed);
  TH1D* yielda16 = dynamic_cast<TH1D*>(indir->Get(Form("AwayYield1650%s",extranamePbPb.Data())));
  yielda16->SetLineColor(kBlue);
  yielda16->SetTitle(Form("Away-Side peak yield for 16GeV/C #leq pT_{trigger}#leq50GeV/c"));
  TH1D* yielda16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayYield1650%s",extraname.Data())));
  yielda16pp->SetLineColor(kRed);  

  
  TH1D* widthen48 = dynamic_cast<TH1D*>(indir->Get(Form("NearWidth_eta12_48%s",extranamePbPb.Data())));
  widthen48->SetLineColor(kBlue);
  widthen48->SetTitle(Form("Near-Side peak width in #eta_{12} for 4GeV/C #leq pT_{trigger}#leq8GeV/c"));
  TH1D* widthen48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearWidth_eta12_48%s",extraname.Data())));
  widthen48pp->SetLineColor(kRed);
  TH1D* widthea48 = dynamic_cast<TH1D*>(indir->Get(Form("AwayWidth_eta12_48%s",extranamePbPb.Data())));
  widthea48->SetLineColor(kBlue);
  widthea48->SetTitle(Form("Away-Side peak width in #eta_{12} for 4GeV/C #leq pT_{trigger}#leq8GeV/c"));
  TH1D* widthea48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayWidth_eta12_48%s",extraname.Data())));
  widthea48pp->SetLineColor(kRed);  

  TH1D* widthen816 = dynamic_cast<TH1D*>(indir->Get(Form("NearWidth_eta12_816%s",extranamePbPb.Data())));
  widthen816->SetLineColor(kBlue);
  widthen816->SetTitle(Form("Near-Side peak width in #eta_{12} for 8GeV/C #leq pT_{trigger}#leq16GeV/c"));
  TH1D* widthen816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearWidth_eta12_816%s",extraname.Data())));
  widthen816pp->SetLineColor(kRed);
  TH1D* widthea816 = dynamic_cast<TH1D*>(indir->Get(Form("AwayWidth_eta12_816%s",extranamePbPb.Data())));
  widthea816->SetLineColor(kBlue);
  widthea816->SetTitle(Form("Away-Side peak width  in #eta_{12} for 8GeV/C #leq pT_{trigger}#leq16GeV/c"));
  TH1D* widthea816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayWidth_eta12_816%s",extraname.Data())));
  widthea816pp->SetLineColor(kRed);  

  TH1D* widthen16 = dynamic_cast<TH1D*>(indir->Get(Form("NearWidth_eta12_1650%s",extranamePbPb.Data())));
  widthen16->SetLineColor(kBlue);
  widthen16->SetTitle(Form("Near-Side peak width in #eta_{12} for 16GeV/C #leq pT_{trigger}#leq50GeV/c"));
  TH1D* widthen16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearWidth_eta12_1650%s",extraname.Data())));
  widthen16pp->SetLineColor(kRed);
  TH1D* widthea16 = dynamic_cast<TH1D*>(indir->Get(Form("AwayWidth_eta12_1650%s",extranamePbPb.Data())));
  widthea16->SetLineColor(kBlue);
  widthea16->SetTitle(Form("Away-Side peak width in #eta_{12} for 16GeV/C #leq pT_{trigger}#leq50GeV/c"));
  TH1D* widthea16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayWidth_eta12_1650%s",extraname.Data())));
  widthea16pp->SetLineColor(kRed); 
  
  
  
  
  TH1D* widthpn48 = dynamic_cast<TH1D*>(indir->Get("NearWidth_theta1_48"));
  if(!methodstring.Contains("single")) widthpn48 = dynamic_cast<TH1D*>(indir->Get(Form("NearWidth_theta1_48%s",extranamePbPb.Data())));
  widthpn48->SetLineColor(kBlue);
  widthpn48->SetTitle(Form("Near-Side peak width in #Phi_{1} for 4GeV/C #leq pT_{trigger}#leq8GeV/c"));
  TH1D* widthpn48pp = dynamic_cast<TH1D*>(ppdir->Get("NearWidth_theta1_48"));
  if(!methodstring.Contains("single")) widthpn48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearWidth_theta1_48%s",extraname.Data())));
  widthpn48pp->SetLineColor(kRed);
  TH1D* widthpa48 = dynamic_cast<TH1D*>(indir->Get("AwayWidth_theta1_48"));
  if(!methodstring.Contains("single")) widthpa48 = dynamic_cast<TH1D*>(indir->Get(Form("AwayWidth_theta1_48%s",extranamePbPb.Data())));
  widthpa48->SetLineColor(kBlue);
  widthpa48->SetTitle(Form("Away-Side peak width in #Phi_{1} for 4GeV/C #leq pT_{trigger}#leq8GeV/c"));
  TH1D* widthpa48pp = dynamic_cast<TH1D*>(ppdir->Get("AwayWidth_theta1_48"));
  if(!methodstring.Contains("single")) widthpa48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayWidth_theta1_48%s",extraname.Data())));
  widthpa48pp->SetLineColor(kRed);  

  TH1D* widthpn816 = dynamic_cast<TH1D*>(indir->Get("NearWidth_theta1_816"));
  if(!methodstring.Contains("single")) widthpn816 = dynamic_cast<TH1D*>(indir->Get(Form("NearWidth_theta1_816%s",extranamePbPb.Data())));
  widthpn816->SetLineColor(kBlue);
  widthpn816->SetTitle(Form("Near-Side peak width in #Phi_{1} for 8GeV/C #leq pT_{trigger}#leq16GeV/c"));
  TH1D* widthpn816pp = dynamic_cast<TH1D*>(ppdir->Get("NearWidth_theta1_816"));
  if(!methodstring.Contains("single")) widthpn816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearWidth_theta1_816%s",extraname.Data())));
  widthpn816pp->SetLineColor(kRed);
  TH1D* widthpa816 = dynamic_cast<TH1D*>(indir->Get("AwayWidth_theta1_816"));
  if(!methodstring.Contains("single")) widthpa816 = dynamic_cast<TH1D*>(indir->Get(Form("AwayWidth_theta1_816%s",extranamePbPb.Data())));
  widthpa816->SetLineColor(kBlue);
  widthpa816->SetTitle(Form("Away-Side peak width  in #Phi_{1} for 8GeV/C #leq pT_{trigger}#leq16GeV/c"));
  TH1D* widthpa816pp = dynamic_cast<TH1D*>(ppdir->Get("AwayWidth_theta1_816"));
  if(!methodstring.Contains("single")) widthpa816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayWidth_theta1_816%s",extraname.Data())));
  widthpa816pp->SetLineColor(kRed);  

  TH1D* widthpn16 = dynamic_cast<TH1D*>(indir->Get("NearWidth_theta1_1650"));
  if(!methodstring.Contains("single")) widthpn16 = dynamic_cast<TH1D*>(indir->Get(Form("NearWidth_theta1_1650%s",extranamePbPb.Data())));
  widthpn16->SetLineColor(kBlue);
  widthpn16->SetTitle(Form("Near-Side peak width in #Phi_{1} for 16GeV/C #leq pT_{trigger}#leq50GeV/c"));
  TH1D* widthpn16pp = dynamic_cast<TH1D*>(ppdir->Get("NearWidth_theta1_1650"));
  if(!methodstring.Contains("single")) widthpn16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearWidth_theta1_1650%s",extraname.Data())));
  widthpn16pp->SetLineColor(kRed);
  TH1D* widthpa16 = dynamic_cast<TH1D*>(indir->Get("AwayWidth_theta1_1650"));
  if(!methodstring.Contains("single")) widthpa16 = dynamic_cast<TH1D*>(indir->Get(Form("AwayWidth_theta1_1650%s",extranamePbPb.Data())));
  widthpa16->SetLineColor(kBlue);
  widthpa16->SetTitle(Form("Away-Side peak width in #Phi_{1} for 16GeV/C #leq pT_{trigger}#leq50GeV/c"));
  TH1D* widthpa16pp = dynamic_cast<TH1D*>(ppdir->Get("AwayWidth_theta1_1650"));
  if(!methodstring.Contains("single")) widthpa16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayWidth_theta1_1650%s",extraname.Data())));
  widthpa16pp->SetLineColor(kRed); 


  TString string = TString();
  TString name = TString();
  int bincent = 0;
  if(TString(indir->GetName()).Contains("0-10")){string.Append(Form("PbPb events with C = 0%%-10%%")); name.Append("0-10"); bincent = 1;}
  if(TString(indir->GetName()).Contains("0-5"))string.Append(Form("PbPb events with C = 0%%-5%%"));
  if(TString(indir->GetName()).Contains("5-10"))string.Append(Form("PbPb events with C = 5%%-10%%"));
  if(TString(indir->GetName()).Contains("10-20")){string.Append(Form("PbPb events with C = 10%%-20%%"));name.Append("10-20");bincent = 2;}
  if(TString(indir->GetName()).Contains("20-40")){string.Append(Form("PbPb events with C = 20%%-40%%"));name.Append("20-40");bincent = 3;}
  if(TString(indir->GetName()).Contains("40-60")){string.Append(Form("PbPb events with C = 40%%-60%%"));name.Append("40-60");bincent = 4;}
  if(TString(indir->GetName()).Contains("60-90")){string.Append(Form("PbPb events with C = 60%%-90%%"));name.Append("60-90");bincent = 5;}

  name.Append(extraname.Data());


  TString outfolder = "~/Dropbox/uio/PhD/analysisnote/figures/Yields_vs_pT/";
  if(TString(options).CompareTo("")){
    outfolder.Clear();
    outfolder.Append("~/Dropbox/uio/PhD/analysisnote/figures/");
    outfolder.Append(options);
    outfolder.Append("/Yields_vs_pT/");
  }
  
  DrawHists(yieldn48,yieldn48pp,outdir,Form("Near_yield_48_%s",indir->GetName()),Form("Near_PbPb/pp_48_%s",indir->GetName()),Form("%sNSY_48_pp_and_PbPb_C_%s_OnePlot.eps",outfolder.Data(),name.Data()));
  DrawHists(yielda48,yielda48pp,outdir,Form("Away_yield_48_%s",indir->GetName()),Form("Away_PbPb/pp_48_%s",indir->GetName()),Form("%sASY_48_pp_and_PbPb_C_%s_OnePlot.eps",outfolder.Data(),name.Data()));
  DrawHists(yieldn816,yieldn816pp,outdir,Form("Near_yield_816_%s",indir->GetName()),Form("Near_PbPb/pp_816_%s",indir->GetName()),Form("%sNSY_816_pp_and_PbPb_C_%s_OnePlot.eps",outfolder.Data(),name.Data()));
  DrawHists(yielda816,yielda816pp,outdir,Form("Away_yield_816_%s",indir->GetName()),Form("Away_PbPb/pp_816_%s",indir->GetName()),Form("%sASY_816_pp_and_PbPb_C_%s_OnePlot.eps",outfolder.Data(),name.Data()));
  DrawHists(yieldn16,yieldn16pp,outdir,Form("Near_yield_16_%s",indir->GetName()),Form("Near_PbPb/pp_16_%s",indir->GetName()),Form("%sNSY_1650_pp_and_PbPb_C_%s_OnePlot.eps",outfolder.Data(),name.Data()));
  DrawHists(yielda16,yielda16pp,outdir,Form("Away_yield_16_%s",indir->GetName()),Form("Away_PbPb/pp_16_%s",indir->GetName()),Form("%sASY_1650_pp_and_PbPb_C_%s_OnePlot.eps",outfolder.Data(),name.Data()));
  
  TString outfolder2 = "~/Dropbox/uio/PhD/analysisnote/figures/Widths_eta_vs_pT/";
  if(TString(options).CompareTo("")){
    outfolder2.Clear();
    outfolder2.Append("~/Dropbox/uio/PhD/analysisnote/figures/");
    outfolder2.Append(options);
    outfolder2.Append("/Widths_eta_vs_pT/");
  }

  DrawHists(widthen48,widthen48pp,outdir,Form("Near_widthe_48_%s",indir->GetName()),"",Form("%sNSWE_48_pp_and_PbPb_C_%s_OnePlot.eps",outfolder2.Data(),name.Data()));
  DrawHists(widthea48,widthea48pp,outdir,Form("Away_widthe_48_%s",indir->GetName()),"",Form("%sASWE_48_pp_and_PbPb_C_%s_OnePlot.eps",outfolder2.Data(),name.Data()));
  DrawHists(widthen816,widthen816pp,outdir,Form("Near_widthe_816_%s",indir->GetName()),"",Form("%sNSWE_816_pp_and_PbPb_C_%s_OnePlot.eps",outfolder2.Data(),name.Data()));
  DrawHists(widthea816,widthea816pp,outdir,Form("Away_widthe_816_%s",indir->GetName()),"",Form("%sASWE_816_pp_and_PbPb_C_%s_OnePlot.eps",outfolder2.Data(),name.Data()));
  DrawHists(widthen16,widthen16pp,outdir,Form("Near_widthe_16_%s",indir->GetName()),"",Form("%sNSWE_1650_pp_and_PbPb_C_%s_OnePlot.eps",outfolder2.Data(),name.Data()));
  DrawHists(widthea16,widthea16pp,outdir,Form("Away_widthe_16_%s",indir->GetName()),"",Form("%sASWE_1650_pp_and_PbPb_C_%s_OnePlot.eps",outfolder2.Data(),name.Data()));

  
  TString outfolder3 = "~/Dropbox/uio/PhD/analysisnote/figures/Widths_theta_vs_pT/";
  if(TString(options).CompareTo("")){
    outfolder3.Clear();
    outfolder3.Append("~/Dropbox/uio/PhD/analysisnote/figures/");
    outfolder3.Append(options);
    outfolder3.Append("/Widths_theta_vs_pT/");
  }
  if(!methodstring.Contains("single")){
    DrawHists(widthpn48,widthpn48pp,outdir,Form("Near_widtht_48_%s",indir->GetName()),"",Form("%sNSWT_48_pp_and_PbPb_C_%s_OnePlot.eps",outfolder3.Data(),name.Data()));
    DrawHists(widthpa48,widthpa48pp,outdir,Form("Away_widtht_48_%s",indir->GetName()),"",Form("%sASWT_48_pp_and_PbPb_C_%s_OnePlot.eps",outfolder3.Data(),name.Data()));
    DrawHists(widthpn816,widthpn816pp,outdir,Form("Near_widtht_816_%s",indir->GetName()),"",Form("%sNSWT_816_pp_and_PbPb_C_%s_OnePlot.eps",outfolder3.Data(),name.Data()));
    DrawHists(widthpa816,widthpa816pp,outdir,Form("Away_widtht_816_%s",indir->GetName()),"",Form("%sASWT_816_pp_and_PbPb_C_%s_OnePlot.eps",outfolder3.Data(),name.Data()));
    DrawHists(widthpn16,widthpn16pp,outdir,Form("Near_widtht_16_%s",indir->GetName()),"",Form("%sNSWT_1650_pp_and_PbPb_C_%s_OnePlot.eps",outfolder3.Data(),name.Data()));
    DrawHists(widthpa16,widthpa16pp,outdir,Form("Away_widtht_16_%s",indir->GetName()),"",Form("%sASWT_1650_pp_and_PbPb_C_%s_OnePlot.eps",outfolder3.Data(),name.Data()));
  }
  
  for(int i =0; i<hists->GetEntries();i++){
    TH1D* hist = dynamic_cast<TH1D*>(hists->At(i));
    if(TString(hist->GetName()).CompareTo("t816a816")==0){
      hist->SetBinContent(bincent,yielda816->GetBinContent(5));
      hist->SetBinError(bincent,yielda816->GetBinError(5));
    }
    if(TString(hist->GetName()).CompareTo("t816a68")==0){
      hist->SetBinContent(bincent,yielda816->GetBinContent(4));
      hist->SetBinError(bincent,yielda816->GetBinError(4));
    }    
    if(TString(hist->GetName()).CompareTo("t816a46")==0){
      hist->SetBinContent(bincent,yielda816->GetBinContent(3));
      hist->SetBinError(bincent,yielda816->GetBinError(3));
    }    
    if(TString(hist->GetName()).CompareTo("t816a34")==0){
      hist->SetBinContent(bincent,yielda816->GetBinContent(2));
      hist->SetBinError(bincent,yielda816->GetBinError(2));
    }        
    if(TString(hist->GetName()).CompareTo("t48a68")==0){
      hist->SetBinContent(bincent,yielda48->GetBinContent(4));
      hist->SetBinError(bincent,yielda48->GetBinError(4));
    }    
    if(TString(hist->GetName()).CompareTo("t48a46")==0){
      hist->SetBinContent(bincent,yielda48->GetBinContent(3));
      hist->SetBinError(bincent,yielda48->GetBinError(3));
    }    
    if(TString(hist->GetName()).CompareTo("t48a34")==0){
      hist->SetBinContent(bincent,yielda48->GetBinContent(2));
      hist->SetBinError(bincent,yielda48->GetBinError(2));
    }    
    
  }
  
  
  
}


void removenegativepoints(TGraphErrors* graph){
  double x = 0;
  double y = 0;
  for(int i = graph->GetN()-1;i>=0;i--){
    graph->GetPoint(i,x,y);
    if(y<0.0){graph->RemovePoint(i);}
  }
}

void DrawHists(TH1* histppo , TH1* hist010o, TH1* hist1020o, TH1* hist2040o, TH1* hist4060o, TH1* hist6090o, TDirectory* outdir, TString name, TString dropboxname , TString zoomname){
  TH1* histpp = dynamic_cast<TH1*>(histppo->Clone("pphistfordrawing"));  
  TH1* hist010 = dynamic_cast<TH1*>(hist010o->Clone("pphistfordrawing"));
  TH1* hist1020 = dynamic_cast<TH1*>(hist1020o->Clone("pphistfordrawing"));
  TH1* hist2040 = dynamic_cast<TH1*>(hist2040o->Clone("pphistfordrawing"));
  TH1* hist4060 = dynamic_cast<TH1*>(hist4060o->Clone("pphistfordrawing"));
  TH1* hist6090 = dynamic_cast<TH1*>(hist6090o->Clone("pphistfordrawing"));
  
  
  TCanvas * canvas = new TCanvas("c1");
  double min = 1000;
  double max = -1000;
  TH1D * range = dynamic_cast<TH1D*>(histpp->Clone("range"));
  range->Reset();
  range->SetStats(kFALSE);
  TString title = TString("Preliminary ");
  if(TString(range->GetTitle()).Contains("Away"))title.Append("away-side per-trigger yields for p_{T}^{trigger} = ");
  if(TString(range->GetTitle()).Contains("Near"))title.Append("near-side per-trigger yields for p_{T}^{trigger} = ");
  if(TString(range->GetTitle()).Contains("4GeV/c<="))title.Append("4-8 GeV/c");
  if(TString(range->GetTitle()).Contains("8GeV/c<="))title.Append("8-16 GeV/c");
  if(TString(range->GetTitle()).Contains("16GeV/c<="))title.Append("16-50 GeV/c");
  range->SetTitle(title.Data());
  histpp->SetLineColor(kBlack);
  histpp->SetLineWidth(2);
  histpp->SetMarkerStyle(20);
  histpp->SetMarkerSize(1.3);
  histpp->SetMarkerColor(kBlack);
  hist010->SetLineColor(kRed);
  hist010->SetLineWidth(2);
  hist010->SetMarkerStyle(24);
  hist010->SetMarkerSize(1.3);
  hist010->SetMarkerColor(kRed);
  hist1020->SetLineColor(kBlue);
  hist1020->SetLineWidth(2);
  hist1020->SetMarkerStyle(24);
  hist1020->SetMarkerSize(1.3);
  hist1020->SetMarkerColor(kBlue);
  hist2040->SetLineColor(kGreen);
  hist2040->SetLineWidth(2);
  hist2040->SetMarkerStyle(24);
  hist2040->SetMarkerSize(1.3);
  hist2040->SetMarkerColor(kGreen);
  hist4060->SetLineColor(kTeal);
  hist4060->SetLineWidth(2);
  hist4060->SetMarkerStyle(24);
  hist4060->SetMarkerSize(1.3);
  hist4060->SetMarkerColor(kTeal);
  hist6090->SetLineColor(kViolet);
  hist6090->SetLineWidth(2);
  hist6090->SetMarkerStyle(24);
  hist6090->SetMarkerSize(1.3);
  hist6090->SetMarkerColor(kViolet);
  TLegend * leg = new TLegend(0.7,0.6,0.9,0.9);
  leg->AddEntry(histpp,"pp");
  leg->AddEntry(hist010,"0%-10%");
  leg->AddEntry(hist1020,"10%-20%");
  leg->AddEntry(hist2040,"20%-40%");
  leg->AddEntry(hist4060,"40%-60%");
  leg->AddEntry(hist6090,"60%-90%");
  
  getminmax(histpp,min,max);
  getminmax(hist010,min,max);
  getminmax(hist1020,min,max);
  getminmax(hist2040,min,max);
  getminmax(hist4060,min,max);
  getminmax(hist6090,min,max);

  if(min>0.0) min = -max*0.1;
  else min =-max*0.1;
  max *=1.1;
  range->SetAxisRange(min,max,"Y");
  range->GetXaxis()->SetTitleOffset(1.15);
  range->GetYaxis()->SetLabelSize(0.025);
  range->GetYaxis()->SetTitleOffset(1.2);
  range->GetYaxis()->SetTitleSize(0.035);
  
  TGraphErrors * graphpp = new TGraphErrors(histpp);
  TGraphErrors * graph010 = new TGraphErrors(hist010);
  TGraphErrors * graph1020 = new TGraphErrors(hist1020);
  TGraphErrors * graph2040 = new TGraphErrors(hist2040);
  TGraphErrors * graph4060 = new TGraphErrors(hist4060);
  TGraphErrors * graph6090 = new TGraphErrors(hist6090);
  
  for(int i = histpp->GetNbinsX() -1; i>=0;i--){
    Double_t x = 0;
    Double_t y = 0;
    graphpp->GetPoint(i,x,y);
    graphpp->SetPointError(i,0.0,graphpp->GetErrorY(i));
    graphpp->SetPoint(i,x-0.3,y);
    if(histpp->GetBinError(i+1)<1.0E-30)graphpp->RemovePoint(i);
    graph010->GetPoint(i,x,y);
    graph010->SetPointError(i,0.0,graph010->GetErrorY(i));
    graph010->SetPoint(i,x-0.2,y);
    if(hist010->GetBinError(i+1)<1.0E-30)graph010->RemovePoint(i);
    graph1020->GetPoint(i,x,y);
    graph1020->SetPointError(i,0.0,graph1020->GetErrorY(i));    
    graph1020->SetPoint(i,x-0.1,y);
    if(hist1020->GetBinError(i+1)<1.0E-30)graph1020->RemovePoint(i);
    graph2040->GetPoint(i,x,y);
    graph2040->SetPointError(i,0.0,graph2040->GetErrorY(i));
    graph2040->SetPoint(i,x,y);
    if(hist2040->GetBinError(i+1)<1.0E-30)graph2040->RemovePoint(i);
    graph4060->GetPoint(i,x,y);
    graph4060->SetPointError(i,0.0,graph4060->GetErrorY(i));
    graph4060->SetPoint(i,x+0.1,y);
    if(hist4060->GetBinError(i+1)<1.0E-30)graph4060->RemovePoint(i);
    graph6090->GetPoint(i,x,y);
    graph6090->SetPointError(i,0.0,graph6090->GetErrorY(i));
    graph6090->SetPoint(i,x+0.2,y);
    if(hist6090->GetBinError(i+1)<1.0E-30)graph6090->RemovePoint(i);
    
  }

  canvas->cd();
  range->Draw("");
  graphpp->Draw("psame");
  graph010->Draw("psame");
  graph1020->Draw("psame");
  graph2040->Draw("psame");
  graph4060->Draw("psame");
  graph6090->Draw("psame");
  leg->Draw("same");
  outdir->cd();
  canvas->Write(name.Data());
  canvas->SaveAs(dropboxname.Data());

  canvas->Clear();
  min = 1000;
  max = -1000;
  getminmax(histpp,min,max,1.0,6.0);
  getminmax(hist010,min,max,1,6.0);
  getminmax(hist1020,min,max,1,6.0);
  getminmax(hist2040,min,max,1,6.0);
  getminmax(hist4060,min,max,1,6.0);
  getminmax(hist6090,min,max,1,6.0);
  if(min>0.0) min = -max*0.1;
  else min =-max*0.1;
  max *=1.1;
  range->SetAxisRange(min,max,"Y");
  range->SetAxisRange(1.0,6.0,"X");
  canvas->cd();
  range->Draw("");
  graphpp->Draw("psame");
  graph010->Draw("psame");
  graph1020->Draw("psame");
  graph2040->Draw("psame");
  graph4060->Draw("psame");
  graph6090->Draw("psame");
  leg->Draw("same");
  canvas->SaveAs(zoomname.Data());
  zoomname.ReplaceAll("zoom","log");
  removenegativepoints(graphpp);
  removenegativepoints(graph010);
  removenegativepoints(graph1020);
  removenegativepoints(graph2040);
  removenegativepoints(graph4060);
  removenegativepoints(graph6090);
  
  canvas->Clear();
  min = 1000;
  max = -1000;
  getminmaxpos(histpp,min,max);
  getminmaxpos(hist010,min,max);
  getminmaxpos(hist1020,min,max);
  getminmaxpos(hist2040,min,max);
  getminmaxpos(hist4060,min,max);
  getminmaxpos(hist6090,min,max);
  max *=1.1;
  
  range->SetAxisRange(min,max,"Y");
  range->Draw("");
  canvas->SetLogy();
  graphpp->Draw("psame");
  graph010->Draw("psame");
  graph1020->Draw("psame");
  graph2040->Draw("psame");
  graph4060->Draw("psame");
  graph6090->Draw("psame");
  leg->Draw("same");
  canvas->SaveAs(zoomname.Data());

  
  delete canvas;
  delete leg;
  delete graphpp;delete graph010; delete graph1020; delete graph2040; delete graph4060; delete graph6090;
  delete histpp;delete hist010; delete hist1020;delete hist2040;delete hist4060; delete hist6090;
}

void PlotHists(TDirectory * dir010 , TDirectory* dir1020 , TDirectory* dir2040 , TDirectory * dir4060 , TDirectory * dir6090, TDirectory * ppdir,TDirectory* outdir,TString methodstring, TString ppperiod,const char* options){
  cout << "";
  if(!dir010)return;
  if(!dir1020)return;
  if(!dir2040)return;
  if(!dir4060)return;
  if(!dir6090)return;
  if(!ppdir)return;
  TString extranamePbPb = TString("");
  TString extraname = TString("");
  if(ppperiod.Contains("LHC11a"))		extraname.Append("_11a");
  if(methodstring.Contains("Eta"))		extraname.Append("_etasubstracted");
  if(methodstring.Contains("singleyield"))	extraname.Append("_single");  
  if(methodstring.Contains("Eta"))		extranamePbPb.Append("_etasubstracted");
  if(methodstring.Contains("singleyield"))	extranamePbPb.Append("_single");  
  
  //pp hists:
  TH1D* yieldn48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearYield48%s",extraname.Data())));
  TH1D* yielda48pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayYield48%s",extraname.Data())));
  TH1D* yieldn816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearYield816%s",extraname.Data())));
  TH1D* yielda816pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayYield816%s",extraname.Data())));
  TH1D* yielda16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("AwayYield1650%s",extraname.Data())));
  TH1D* yieldn16pp = dynamic_cast<TH1D*>(ppdir->Get(Form("NearYield1650%s",extraname.Data())));
  
  TH1D* yieldn48_010 = dynamic_cast<TH1D*>(dir010->Get(Form("NearYield48%s",extranamePbPb.Data())));
  TH1D* yielda48_010 = dynamic_cast<TH1D*>(dir010->Get(Form("AwayYield48%s",extranamePbPb.Data())));
  TH1D* yieldn816_010 = dynamic_cast<TH1D*>(dir010->Get(Form("NearYield816%s",extranamePbPb.Data())));
  TH1D* yielda816_010 = dynamic_cast<TH1D*>(dir010->Get(Form("AwayYield816%s",extranamePbPb.Data())));
  TH1D* yieldn16_010 = dynamic_cast<TH1D*>(dir010->Get(Form("NearYield1650%s",extranamePbPb.Data())));
  TH1D* yielda16_010 = dynamic_cast<TH1D*>(dir010->Get(Form("AwayYield1650%s",extranamePbPb.Data())));

  TH1D* yieldn48_1020 = dynamic_cast<TH1D*>(dir1020->Get(Form("NearYield48%s",extranamePbPb.Data())));
  TH1D* yielda48_1020 = dynamic_cast<TH1D*>(dir1020->Get(Form("AwayYield48%s",extranamePbPb.Data())));
  TH1D* yieldn816_1020 = dynamic_cast<TH1D*>(dir1020->Get(Form("NearYield816%s",extranamePbPb.Data())));
  TH1D* yielda816_1020 = dynamic_cast<TH1D*>(dir1020->Get(Form("AwayYield816%s",extranamePbPb.Data())));
  TH1D* yieldn16_1020 = dynamic_cast<TH1D*>(dir1020->Get(Form("NearYield1650%s",extranamePbPb.Data())));
  TH1D* yielda16_1020 = dynamic_cast<TH1D*>(dir1020->Get(Form("AwayYield1650%s",extranamePbPb.Data())));
  
  TH1D* yieldn48_2040 = dynamic_cast<TH1D*>(dir2040->Get(Form("NearYield48%s",extranamePbPb.Data())));
  TH1D* yielda48_2040 = dynamic_cast<TH1D*>(dir2040->Get(Form("AwayYield48%s",extranamePbPb.Data())));
  TH1D* yieldn816_2040 = dynamic_cast<TH1D*>(dir2040->Get(Form("NearYield816%s",extranamePbPb.Data())));
  TH1D* yielda816_2040 = dynamic_cast<TH1D*>(dir2040->Get(Form("AwayYield816%s",extranamePbPb.Data())));
  TH1D* yieldn16_2040 = dynamic_cast<TH1D*>(dir2040->Get(Form("NearYield1650%s",extranamePbPb.Data())));
  TH1D* yielda16_2040 = dynamic_cast<TH1D*>(dir2040->Get(Form("AwayYield1650%s",extranamePbPb.Data())));
  
  TH1D* yieldn48_4060 = dynamic_cast<TH1D*>(dir4060->Get(Form("NearYield48%s",extranamePbPb.Data())));
  TH1D* yielda48_4060 = dynamic_cast<TH1D*>(dir4060->Get(Form("AwayYield48%s",extranamePbPb.Data())));
  TH1D* yieldn816_4060 = dynamic_cast<TH1D*>(dir4060->Get(Form("NearYield816%s",extranamePbPb.Data())));
  TH1D* yielda816_4060 = dynamic_cast<TH1D*>(dir4060->Get(Form("AwayYield816%s",extranamePbPb.Data())));
  TH1D* yieldn16_4060 = dynamic_cast<TH1D*>(dir4060->Get(Form("NearYield1650%s",extranamePbPb.Data())));
  TH1D* yielda16_4060 = dynamic_cast<TH1D*>(dir4060->Get(Form("AwayYield1650%s",extranamePbPb.Data())));  
  
  TH1D* yieldn48_6090 = dynamic_cast<TH1D*>(dir6090->Get(Form("NearYield48%s",extranamePbPb.Data())));
  TH1D* yielda48_6090 = dynamic_cast<TH1D*>(dir6090->Get(Form("AwayYield48%s",extranamePbPb.Data())));
  TH1D* yieldn816_6090 = dynamic_cast<TH1D*>(dir6090->Get(Form("NearYield816%s",extranamePbPb.Data())));
  TH1D* yielda816_6090 = dynamic_cast<TH1D*>(dir6090->Get(Form("AwayYield816%s",extranamePbPb.Data())));
  TH1D* yieldn16_6090 = dynamic_cast<TH1D*>(dir6090->Get(Form("NearYield1650%s",extranamePbPb.Data())));
  TH1D* yielda16_6090 = dynamic_cast<TH1D*>(dir6090->Get(Form("AwayYield1650%s",extranamePbPb.Data())));
  
  
  TString outfolder = "~/Dropbox/uio/PhD/analysisnote/figures/";
  if(TString(options).CompareTo("")){
    outfolder.Clear();
    outfolder.Append("~/Dropbox/uio/PhD/analysisnote/figures/");
    outfolder.Append(options);
    outfolder.Append("/");
    
  }
  TString name = TString(extraname);
  DrawHists(yieldn48pp,yieldn48_010,yieldn48_1020,yieldn48_2040,yieldn48_4060,yieldn48_6090,outdir,Form("Near_yield_48"),Form("%sNSY_48_Cents_%s.eps",outfolder.Data(),name.Data()),Form("%sNSY_48_Cents_%s_zoom.eps",outfolder.Data(),name.Data()));
  DrawHists(yielda48pp,yielda48_010,yielda48_1020,yielda48_2040,yielda48_4060,yielda48_6090,outdir,Form("Away_yield_48"),Form("%sASY_48_Cents_%s.eps",outfolder.Data(),name.Data()),Form("%sASY_48_Cents_%s_zoom.eps",outfolder.Data(),name.Data()));
  DrawHists(yieldn816pp,yieldn816_010,yieldn816_1020,yieldn816_2040,yieldn816_4060,yieldn816_6090,outdir,Form("Near_yield_816"),Form("%sNSY_816_Cents_%s.eps",outfolder.Data(),name.Data()),Form("%sNSY_816_Cents_%s_zoom.eps",outfolder.Data(),name.Data()));
  DrawHists(yielda816pp,yielda816_010,yielda816_1020,yielda816_2040,yielda816_4060,yielda816_6090,outdir,Form("Away_yield_816_"),Form("%sASY_816_Cents_%s.eps",outfolder.Data(),name.Data()),Form("%sASY_816_Cents_%s_zoom.eps",outfolder.Data(),name.Data()));
  DrawHists(yieldn16pp,yieldn16_010,yieldn16_1020,yieldn16_2040,yieldn16_4060,yieldn16_6090,outdir,Form("Near_yield_16"),Form("%sNSY_1650_Cents_%s.eps",outfolder.Data(),name.Data()),Form("%sNSY_1650_Cents_%s_zoom.eps",outfolder.Data(),name.Data()));
  DrawHists(yielda16pp,yielda16_010,yielda16_1020,yielda16_2040,yielda16_4060,yielda16_6090,outdir,Form("Away_yield_16_"),Form("%sASY_1650_Cents_%s.eps",outfolder.Data(),name.Data()),Form("%sASY_1650_Cents_%s_zoom.eps",outfolder.Data(),name.Data()));  
}

void PlotCollectedValues(TDirectory* directory,TDirectory* outdir,const char* options){
  TDirectory * ppdir = directory->GetDirectory("pp_min_bias");
  TDirectory * PbPb010 = directory->GetDirectory("PbPb_0-10");
  TDirectory * PbPb05 = directory->GetDirectory("PbPb_0-5");
  TDirectory * PbPb510 = directory->GetDirectory("PbPb_5-10");
  TDirectory * PbPb1020 = directory->GetDirectory("PbPb_10-20");
  TDirectory * PbPb2040 = directory->GetDirectory("PbPb_20-40");
  TDirectory * PbPb4060 = directory->GetDirectory("PbPb_40-60");
  TDirectory * PbPb6090 = directory->GetDirectory("PbPb_60-90");
  if(!ppdir){
    cout << "ERROR - no pp."<<endl;
    return;
  }
  //now we want a histogram with centrality on the x axis:
  TH1D * hist  = new TH1D("hist", "title",5,0.0,5.0);
  hist->SetStats(false);
  hist->GetXaxis()->SetBinLabel(1,Form("0.0%%-10.0%%"));
  hist->GetXaxis()->SetBinLabel(2,Form("10.0%%-20.0%%"));
  hist->GetXaxis()->SetBinLabel(3,Form("20.0%%-40.0%%"));
  hist->GetXaxis()->SetBinLabel(4,Form("40.0%%-60.0%%"));
  hist->GetXaxis()->SetBinLabel(5,Form("60.0%%-90.0%%"));
//   hist->GetXaxis()->LabelsOption("v");
  hist->GetXaxis()->SetTitle(Form("Centrality [%%]"));
  TObjArray* hists = new TObjArray();
  hists->Add(hist->Clone("t816a816"));
  hists->Add(hist->Clone("t816a68"));
  hists->Add(hist->Clone("t816a46"));
  hists->Add(hist->Clone("t816a34"));
  hists->Add(hist->Clone("t48a68"));
  hists->Add(hist->Clone("t48a46"));
  hists->Add(hist->Clone("t48a34"));

  
  PlotHists(PbPb010,PbPb1020,PbPb2040,PbPb4060,PbPb6090,ppdir,outdir,TString("yield/simple"),("ppLHC10d"),options);
  PlotHists(PbPb010,PbPb1020,PbPb2040,PbPb4060,PbPb6090,ppdir,outdir,TString("yield/simple"),("ppLHC11a"),options);
  PlotHists(PbPb010,PbPb1020,PbPb2040,PbPb4060,PbPb6090,ppdir,outdir,TString("yield/singleyield"),("ppLHC10d"),options);
  PlotHists(PbPb010,PbPb1020,PbPb2040,PbPb4060,PbPb6090,ppdir,outdir,TString("yield/singleyield"),("ppLHC11a"),options);
  
  PlotHists(PbPb010,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb010,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb010,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb010,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb010,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb010,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb010,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb010,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb05,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb510,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb1020,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb2040,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC11a"),options);  
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb4060,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC11a"),options);  
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("yield/singleyield"),TString("ppLHC11a"),options);
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC10d"),options);
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("Etasubstracted/yield/simple"),TString("ppLHC11a"),options);
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC10d"),options);
  PlotHists(PbPb6090,ppdir,outdir,hists,TString("Etasubstracted/yield/singleyield"),TString("ppLHC11a"),options);  
  outdir->mkdir("centhists")->cd();
  for(int i =0; i<hists->GetEntries();i++){
    TH1D* hist = dynamic_cast<TH1D*>(hists->At(i));
    if(TString(hist->GetName()).CompareTo("t816a816")==0){
      hist->SetTitle("Ratio of awayside yield in PbPb/pp for 8GeV/c<p_{T trigger}<16GeV/c and  8GeV/c<p_{T A}<16GeV/c");
    }
    if(TString(hist->GetName()).CompareTo("t816a68")==0){
      hist->SetTitle("Ratio of awayside yield in PbPb/pp for 8GeV/c<p_{T trigger}<16GeV/c and  6GeV/c<p_{T A}<8GeV/c");
    }    
    if(TString(hist->GetName()).CompareTo("t816a46")==0){
      hist->SetTitle("Ratio of awayside yield in PbPb/pp for 8GeV/c<p_{T trigger}<16GeV/c and  4GeV/c<p_{T A}<6GeV/c");
    }    
    if(TString(hist->GetName()).CompareTo("t816a34")==0){
      hist->SetTitle("Ratio of awayside yield in PbPb/pp for 8GeV/c<p_{T trigger}<16GeV/c and  3GeV/c<p_{T A}<4GeV/c");
    }        
    if(TString(hist->GetName()).CompareTo("t48a68")==0){
      hist->SetTitle("Ratio of awayside yield in PbPb/pp for 4GeV/c<p_{T trigger}<8GeV/c and  6GeV/c<p_{T A}<8GeV/c");
    }    
    if(TString(hist->GetName()).CompareTo("t48a46")==0){
      hist->SetTitle("Ratio of awayside yield in PbPb/pp for 4GeV/c<p_{T trigger}<8GeV/c and  4GeV/c<p_{T A}<6GeV/c");
    }    
    if(TString(hist->GetName()).CompareTo("t48a34")==0){
      hist->SetTitle("Ratio of awayside yield in PbPb/pp for 4GeV/c<p_{T trigger}<8GeV/c and  3GeV/c<p_{T A}<4GeV/c");
    }      
    hist->Write();
  }
}

void PlotMethods(TObjArray* indirs,TDirectory* outdir , TString options){
  //compare methods for same production at all pT bins:
  comparemethods(indirs,outdir->mkdir("LHC10d"),"LHC10d_iteration_simple#LHC10d_iteration_single#LHC10d_etasubstracted_simple#LHC10d_etasubstracted_single","pt_yield_away_t4");
  comparemethods(indirs,outdir->GetDirectory("LHC10d"),"LHC10d_iteration_simple#LHC10d_iteration_single#LHC10d_etasubstracted_simple#LHC10d_etasubstracted_single","pt_widthetaaway_t4");
  comparemethods(indirs,outdir->GetDirectory("LHC10d"),"LHC10d_iteration_simple#LHC10d_etasubstracted_simple","pt_widthThetaaway_t4");

  comparemethods(indirs,outdir->mkdir("LHC10h"),"LHC10h_iteration_simple#LHC10h_iteration_single#LHC10h_etasubstracted_simple#LHC10h_etasubstracted_single","pt_yield_away_cloop_t4");
  comparemethods(indirs,outdir->GetDirectory("LHC10h"),"LHC10h_iteration_simple#LHC10h_iteration_single#LHC10h_etasubstracted_simple#LHC10h_etasubstracted_single","pt_widthetaaway_cloop_t4");
  comparemethods(indirs,outdir->GetDirectory("LHC10h"),"LHC10h_iteration_simple#LHC10h_etasubstracted_simple","pt_widthThetaaway_cloop_t4");
  
}

void Plot2dHist(Double_t &METAscale, Double_t &META2scale, Double_t &METriggerscale,TString infile, TString histname , TString folder, TString outdir, TString plotconf, TString Xaxisname, TString Yaxisname , TString title){
  TFile * file = TFile::Open(infile.Data());
  TDirectory * dir = file->GetDirectory(folder.Data());
  TH2D * hist = dynamic_cast<TH2D*>(dir->Get(histname.Data()));
  TCanvas * canv = new TCanvas("name");
  hist->SetTitle(title.Data());
  if(title.CompareTo("")==0)hist->SetBit(TH1::kNoTitle);
  hist->SetStats(false);
  if(Xaxisname.CompareTo(""))hist->GetXaxis()->SetTitle(Xaxisname.Data());
  hist->GetXaxis()->SetTitleOffset(1.2);
  if(Yaxisname.CompareTo(""))hist->GetYaxis()->SetTitle(Yaxisname.Data());
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetZaxis()->SetLabelSize(0.02);
  hist->GetZaxis()->SetTitleSize(0.02);
  hist->GetZaxis()->SetTitleOffset(2.0);
  if(plotconf.CompareTo("colz")==0){
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(0.7);
    hist->GetZaxis()->SetLabelSize(0.02);
    hist->GetZaxis()->SetTitleSize(0.01);
    hist->GetZaxis()->SetTitleOffset(1.0);
    hist->GetZaxis()->SetTitle("");
  }
  hist->Draw(plotconf.Data());
  METAscale = dynamic_cast<TParameter<double>*>(dir->Get("METASCALE"))->GetVal();
  META2scale = dynamic_cast<TParameter<double>*>(dir->Get("META2SCALE"))->GetVal();
  METriggerscale = dynamic_cast<TParameter<double>*>(dir->Get("METRIGGERSCALE"))->GetVal();
  canv->SaveAs(outdir.Data());
  
  delete canv;
  delete hist;
  file->Close();
}

void Plot2dHist(Double_t &zmax, Double_t &zmin, TString infile, TString histname , TString folder, TString outdir, TString plotconf, TString Xaxisname, TString Yaxisname , TString title){
  TFile * file = TFile::Open(infile.Data());
  TDirectory * dir = file->GetDirectory(folder.Data());
  TH2D * hist = dynamic_cast<TH2D*>(dir->Get(histname.Data()));
  TCanvas * canv = new TCanvas("name");
  hist->SetTitle(title.Data());
  hist->SetStats(false);
  if(Xaxisname.CompareTo(""))hist->GetXaxis()->SetTitle(Xaxisname.Data());
  hist->GetXaxis()->SetTitleOffset(1.2);
  if(Yaxisname.CompareTo(""))hist->GetYaxis()->SetTitle(Yaxisname.Data());
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetZaxis()->SetLabelSize(0.02);
  zmax = 1.1*hist->GetMaximum();
  zmin = 0.8*hist->GetMinimum();
  hist->SetAxisRange(zmin,zmax,"Z");
  hist->GetZaxis()->SetTitleSize(0.02);
  hist->GetZaxis()->SetTitleOffset(2.0);
  if(plotconf.CompareTo("colz")==0){
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(0.7);
    hist->GetZaxis()->SetLabelSize(0.02);
    hist->GetZaxis()->SetTitleSize(0.01);
    hist->GetZaxis()->SetTitleOffset(1.0);
    hist->GetZaxis()->SetTitle("");
  }
  hist->Draw(plotconf.Data());
  canv->SaveAs(outdir.Data());
  delete canv;
  delete hist;
  file->Close();
}

void Plot2dHist(TString infile,Double_t zmax, Double_t zmin, TString histname , TString folder, TString outdir, TString plotconf, TString Xaxisname, TString Yaxisname , TString title,Double_t scale){
  TFile * file = TFile::Open(infile.Data());
  TDirectory * dir = file->GetDirectory(folder.Data());
  TH2D * hist = dynamic_cast<TH2D*>(dir->Get(histname.Data()));
  TCanvas * canv = new TCanvas("name");
  hist->SetTitle(title.Data());
  hist->SetStats(false);
  if(Xaxisname.CompareTo(""))hist->GetXaxis()->SetTitle(Xaxisname.Data());
  hist->GetXaxis()->SetTitleOffset(1.2);
  if(Yaxisname.CompareTo(""))hist->GetYaxis()->SetTitle(Yaxisname.Data());
  hist->Scale(scale);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetZaxis()->SetLabelSize(0.02);
  hist->SetAxisRange(zmin,zmax,"Z");
  hist->GetZaxis()->SetTitleSize(0.02);
  hist->GetZaxis()->SetTitleOffset(2.0);
  if(plotconf.CompareTo("colz")==0){
    hist->GetXaxis()->SetTitleOffset(0.9);
    hist->GetYaxis()->SetTitleOffset(0.7);
    hist->GetZaxis()->SetLabelSize(0.02);
    hist->GetZaxis()->SetTitleSize(0.01);
    hist->GetZaxis()->SetTitleOffset(1.0);
    hist->GetZaxis()->SetTitle("");
  }
  hist->Draw(plotconf.Data());
  canv->SaveAs(outdir.Data());
  
  delete canv;
  delete hist;
  file->Close();
}

void ETStats(TString infile1, TString infile2,TString infile3,bool isPb,TString outdir){
  TFile * file1 = TFile::Open(infile1.Data());
  TFile * file2 = TFile::Open(infile2.Data());
  TFile * file3 = TFile::Open(infile3.Data());

  
  TList* dir1 = file1->GetListOfKeys();
  TDirectory* Dir1 = file1->GetDirectory(dir1->At(0)->GetName());
  TList* dir11 = Dir1->GetListOfKeys();
  THashList* List1 = dynamic_cast<THashList*>(Dir1->Get(dir11->At(0)->GetName()));
  
  TList* dir2 = file2->GetListOfKeys();
  TDirectory* Dir2 = file2->GetDirectory(dir2->At(0)->GetName());
  TList* dir21 = Dir2->GetListOfKeys();
  THashList* List2 = dynamic_cast<THashList*>(Dir2->Get(dir21->At(0)->GetName()));

  TList* dir3 = file3->GetListOfKeys();
  TDirectory* Dir3 = file3->GetDirectory(dir3->At(0)->GetName());
  TList* dir31 = Dir3->GetListOfKeys();
  THashList* List3 = dynamic_cast<THashList*>(Dir3->Get(dir31->At(0)->GetName()));
  
  
  TCanvas * canvas = new TCanvas("c1");
  
  //Track distributions:
  TH1D * trackpt1  = dynamic_cast<TH1D*>(List1->FindObject("trackPt"));
  TH1D * trackpt2  = dynamic_cast<TH1D*>(List2->FindObject("trackPt"));
  TH1D * trackpt3  = dynamic_cast<TH1D*>(List3->FindObject("trackPt"));
  trackpt1->Add(trackpt2);
  trackpt1->Add(trackpt3);
  trackpt1->GetXaxis()->SetRangeUser(0.0,50.0);
  trackpt1->SetXTitle("pT [GeV/c]");
  trackpt1->GetYaxis()->SetTitle("Number of Tracks");
  canvas->SetLogy();
  trackpt1->SetBit(TH1::kNoTitle,true);
  trackpt1->Draw("");
  canvas->SaveAs(Form("%s/trackpt.eps",outdir.Data()));
  
  canvas->Clear();
  canvas->SetLogy(false);
  
  TH1D * trackphi1  = dynamic_cast<TH1D*>(List1->FindObject("trackPhi"));
  TH1D * trackphi2  = dynamic_cast<TH1D*>(List2->FindObject("trackPhi"));
  TH1D * trackphi3  = dynamic_cast<TH1D*>(List3->FindObject("trackPhi"));
  trackphi1->Add(trackphi2);
  trackphi1->Add(trackphi3);
  trackphi1->SetXTitle("#phi [rad]");
  trackphi1->GetYaxis()->SetTitle("Number of Tracks");
  trackphi1->SetStats(false);
  double max = trackphi1->GetMaximum();
  trackphi1->SetAxisRange(0.0,max*1.1,"y");
  trackphi1->SetBit(TH1::kNoTitle,true);
  trackphi1->Draw("");
  canvas->SaveAs(Form("%s/trackphi.eps",outdir.Data()));
  
  canvas->Clear();
  
  //Event characateristics
  TH2D * vertexvscent = dynamic_cast<TH2D*>(List1->FindObject("centVsZVertex"));
  vertexvscent->Draw("colz");
  vertexvscent->SetBit(TH1::kNoTitle,true);
  vertexvscent->SetStats(false);
  if(isPb){
    vertexvscent->SetYTitle("V_{Z} [cm]");
    vertexvscent->SetXTitle("Centrality [%]");
    vertexvscent->SetTitle("Number of Events in V_{Z} and Centrality bins. ");
  }
  if(!isPb){
    vertexvscent->SetYTitle("V_{Z} [cm]");
    vertexvscent->SetXTitle("Number of tracks per Event");
    vertexvscent->SetTitle("Number of Events in V_{Z} and number of tracks bins. ");
  }
  canvas->SaveAs(Form("%s/vz_mult.eps",outdir.Data()));
  canvas->Clear();
  
  TH1D * vertex = vertexvscent->ProjectionY("vertex");
  vertex->SetYTitle("Number of Events");
  TH1D * cent =   vertexvscent->ProjectionX("# events");
  cent->SetYTitle("Number of Events");  
  
  
  vertex->SetStats(false);
  vertex->SetAxisRange(0.0,vertex->GetMaximum()*1.1,"y");
  vertex->SetBit(TH1::kNoTitle,true);
  vertex->GetYaxis()->SetLabelSize(0.03);
  vertex->Draw();
  canvas->SaveAs(Form("%s/vz.eps",outdir.Data()));
  
  if(isPb)cent->SetAxisRange(0.0,cent->GetMaximum()*1.1,"y");
  if(!isPb)canvas->SetLogy();
  cent->SetBit(TH1::kNoTitle,true);
  cent->GetYaxis()->SetLabelSize(0.03);
  
  cent->Draw();
  canvas->SaveAs(Form("%s/cent.eps",outdir.Data()));
  if(!isPb)canvas->SetLogy(false);
  
  file1->Close();
  file2->Close();
  file3->Close();
//   delete dir1; delete dir2;
  delete canvas;
  
}

