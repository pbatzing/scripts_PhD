#ifndef THREEPARTICLECORRECTIONFUNCTIONS_H
#define THREEPARTICLECORRECTIONFUNCTIONS_H
#include <iostream>
#include <stdio.h>
#include "TMath.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TString.h"
#include "TList.h"
#include "TPRegexp.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH3D.h"
#include "TKey.h"
#include "TObjArray.h"
#include "TCanvas.h"
#include "TMinuit.h"
#include "TPaveText.h"
#include "TParameter.h"
#include "TPaveLabel.h"
#include "TObjString.h"
#include "TAxis.h"
#include "TF1.h"
#include "TLegend.h"
#include "THashList.h" 
#include "TSystem.h"
#include "TPaletteAxis.h"
#include "TGraphErrors.h"
TDirectory* resultsdirectory(TDirectory* motherdir, const char* name);



//class for containing bin directories
class BinDirs : public TObject{
public:
  BinDirs(TFile* file, const char* bin, bool empty = false);
  BinDirs(TDirectory* dir, const char* bin, bool empty = false );
  BinDirs(TDirectory* samed, TDirectory * METAd,TDirectory * META2d, TDirectory * METriggerd, bool empty = false );

  bool CompareTo(const char* name){return (Bin.CompareTo(name)==0);}
  TDirectory * Same(){return Samedir;}
  TDirectory * Samediv(){return Samedir->GetDirectory("divided");}
  TDirectory * SameDir(const char* name){return dynamic_cast<TDirectory*>(Samedir->Get(name));}
  TDirectory * META(){return METAdir;}
  TDirectory * METAdiv(){return METAdir->GetDirectory("divided");}
  TDirectory * METADir(const char* name){return dynamic_cast<TDirectory*>(METAdir->Get(name));}
  TDirectory * META2(){return META2dir;}
  TDirectory * META2div(){return META2dir->GetDirectory("divided");}
  TDirectory * META2Dir(const char* name){return dynamic_cast<TDirectory*>(META2dir->Get(name));}
  TDirectory * METrigger(){return METriggerdir;}
  TDirectory * METriggerdiv(){return METriggerdir->GetDirectory("divided");}
  TDirectory * METriggerDir(const char* name){return dynamic_cast<TDirectory*>(METriggerdir->Get(name));}
  BinDirs resultsdirectory(const char* dir);
  const char* path(){return Samedir->GetPath();}
  const char* pathMETA(){return METAdir->GetPath();}
  const char* pathMETA2(){return META2dir->GetPath();}
  const char* pathMETrigger(){return METriggerdir->GetPath();}
private:
  TDirectory * Samedir;
  TDirectory * METAdir;
  TDirectory * META2dir;
  TDirectory * METriggerdir;
  TString    Bin;
  bool ready;

};

TList * 	GetMZDirectories(TDirectory* same);
TList * 	GetMZDirectories(TDirectory* same,float vertexcut,bool debug);

TStringToken 	GetHistTokens(TDirectory* dir);
TH2D * PrepareHist(TDirectory* dir, const char* subdir, const char* name,const char* title,const char* rmplmode);
void canvasmaker(const char* histname, TObjArray* multdirlist,bool isPbPb =  true);
void canvasruns(const char* canname, const char* title, TH1* hist, TDirectory* dir);
TCanvas* canvaspTbins(const char* canname, const char* title, TList* subdirlist);
TCanvas* Makecanvas(TH2D* hist, const char* name, Bool_t Stats,Bool_t remedge = false);
TCanvas* Makecanvas(TH2D* histtopl, TH2D* histtopr, TH2D* histbotl, TH2D* histbotr, const char* name, Bool_t Stats);
void savedircontent(BinDirs* same, double METAScale, double METriggerScale, const char* iteration);
void savedircontent(BinDirs* same, double METAScale,double META2Scale, double METriggerScale, const char* iteration);
void RemovePlateau(Double_t plateauheight, TH2D * hist);
void CollectHistbinstats(const char* histname,TList * directories, TObjArray* multdirlist,Bool_t isPbPb = true);
void CollectHist(const char* histname,TList * directories, TObjArray* multdirlist,Bool_t pearsonserrors = false,Bool_t collecdivfirst = false, Bool_t isPbPb = true);
void CollectHistbs(TH1D* histo, TList * directories, TObjArray* multdirlist,Bool_t isPbPb = true);
void CollectHist(TH1D* histo, TList * directories, TObjArray* multdirlist, Bool_t isPbPb=true);
void CollectHist(TH2D* histo, TList * directories, TObjArray* multdirlist,Bool_t pearsonserrors = false,Bool_t collectdivfirst = false, Bool_t isPbPb=true,Bool_t divbybinwidth = false);
Double_t Chi2(Double_t scaleMETA, Double_t scaleMETrigger);
void GetScalingFactors(BinDirs * BinDir,const char* in, const char* out,Double_t * METAscale, Double_t * METriggerscale,bool ispp, bool isPbPb, int iter);

void GetScalingFactors(BinDirs * BinDir,const char* in, const char* out,Double_t * METAscale, Double_t * META2scale, Double_t * METriggerscale,bool ispp, bool isPbPb, int iter);
void AddSigBins(TH2D* hist, TH2D* METAhist, TH2D* METriggerhist);
void Correctpp(BinDirs* BinDir,bool metaindependent);
void CorrectPbPb(BinDirs* BinDir,bool metaindependent);
void fitwith(TDirectory * dir, const char* type, TH2D * histo,Double_t etalimit);
void simpleyield(TDirectory * dir, TH2D* histo,double sidebandlow , double sidebandhigh);
void singleyield(TDirectory * dir, TH2D* histo,double sidebandlow , double sidebandhigh);
void fillwithbinnr(TH2D* fillfrom, TH1D* fillto, int binnr,bool width=false);
void removeconstant(TH1D * hist, Double_t plateau, Double_t erroronit);
void extractbinyield(TDirectory* dir, TDirectory* yielddir, Double_t etalimit,const char* options);
double GetValue(TDirectory * dir, TString name);

//functions to create canvases comparing things:
void PresentPbPb(TFile*infile,TDirectory*outpath,const char* type);
void Presentpp(TFile*infile,TDirectory*outpath,const char* type);
void drawinpad(TDirectory * dir, TString centpath,TString histname,TObjArray* padarray, int index,const char* drawop = "colz" );
void writeinpad(TObjArray*padarray,int index, const char* text);
void getminmax(TDirectory * dir,TString centpath,TString histname, double  &min,double & max);
void getminmax(TH1*hist, double  &min,double & max);
void drawcanvas(TDirectory * dir, TString centpath,TString histname,int histcolor,TCanvas * canvas, TLegend * legend, TString legendname,TString histtitleconst,const char* drawop = "same", double minimum = -100, double maximum = -100 );

//Comparing histograms:
void drawptbins(TObjArray* indirs, TDirectory* outdir, TString histname, TString type, TString Cent, TObjArray* trigarray, TObjArray* assarray);
void drawcentcom(TObjArray* indirs, TDirectory* outdir, TString histname, TString type, int ptT, int ptA);
void ComparePt(bool isPbPb,TObjArray* indirs, TDirectory * outdir, TString histname);
void CompareCent(TObjArray* indirs, TDirectory * outdir, TString histname);
void Compare(TObjArray * indirs , TDirectory* outdir , TString histname , TString comparetype, TString titlecard);

void extractvalues(TObjArray * indirs, TDirectory* outdir ,TString subdir, TH1D* hist,TString methodstring, TString ppperiod=TString("ppLHC10d"));
void comparemethods(TObjArray * indirs, TDirectory* outdir, TString prods,TString typehist);
void yieldcanvases(TObjArray * indirs, TDirectory* outdir ,const char* namehist, TString prodpp);
void yieldcanvasespp(TObjArray * indirs, TDirectory* outdir ,TString prod1 , TString prod2);
void PlotHists(TDirectory * indir, TDirectory* ppdir, TDirectory* outdir,TObjArray* hists,TString methodstring, TString ppperiod=TString("ppLHC10d"),const char* options = "");
void CompareYields(TObjArray * indirs,TDirectory* outdir ,TString options);
void DrawHists(TH1* hist1 , TH1* hist2, TDirectory* outdir, TString name,TString name2, TString dropboxname );
void PlotCollectedValues(TDirectory* directory,TDirectory* outdir,const char* options = "");
void PlotMethods(TObjArray* indirs,TDirectory* outdir , TString options);
void Plot2dHist(Double_t &METAscale,Double_t &META2scale,Double_t &METriggerscale,TString infile, TString histname , TString folder, TString outdir, TString plotconf, TString Xaxisname, TString Yaxisname , TString title);
void Plot2dHist(Double_t &zmax, Double_t &zmin, TString infile, TString histname , TString folder, TString outdir, TString plotconf, TString Xaxisname, TString Yaxisname , TString title);
void Plot2dHist(TString infile,Double_t zmax, Double_t zmin, TString histname , TString folder, TString outdir, TString plotconf, TString Xaxisname, TString Yaxisname , TString title,Double_t scale=1.0);
void ETStats(TString infile1, TString infile2,TString infile3,bool isPbPb,TString outdir);

#endif