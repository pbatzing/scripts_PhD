#include <iostream>
#include <stdio.h>
#include "TString.h"
using namespace std;
void runThreePartDC(const char* mode = "Compare" ,const char* options = ""){
  //run macro to compile and run correction code
  //mode decides the correct running mode, options is passed on to the function.
  TString workingdir = TString(gSystem->WorkingDirectory());
  if(workingdir.Contains("paulbatzing")){workingdir.Clear();workingdir.Append("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");}
  else workingdir.Clear();
  gROOT->LoadMacro(Form("%sThreeParticleCorrectionFunctions.cxx+g",workingdir.Data()));
  gROOT->LoadMacro(Form("%sThreePartDC.cxx+g",workingdir.Data()));
  gPrintViaErrorHandler = kTRUE;  //makes minuit quiet
  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  TGaxis::SetMaxDigits(4);
  TString TSmode = TString(mode);
  if(TSmode.CompareTo("Draw")==0){Draw(options);}
  if(TSmode.CompareTo("Collect")==0){Collect(options);}
  if(TSmode.CompareTo("Correct")==0){Correct(options);}
  if(TSmode.CompareTo("RunsLHC11h")==0){Runs11h(options,workingdir +TString("LHC11h/"));}
  if(TSmode.CompareTo("pTBins")==0){pTBins(options);}
  if(TSmode.CompareTo("Yield")==0){yield(options);}
  if(TSmode.CompareTo("Present")==0){Present(options);}
  if(TSmode.CompareTo("BinsPeriods")==0){pTBinsPer(options);}
  if(TSmode.CompareTo("Compare")==0){ComparisonFile(options);}
  if(TSmode.CompareTo("Sign")==0){Sign(options);}
  if(TSmode.CompareTo("ExamplePlots")==0){ExamplePlots(options);}
  if(TSmode.CompareTo("Compile") == 0){cout <<"compiled"<<endl;}
}
