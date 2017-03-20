#include "TString.h"
#include "TSystem.h"
#include "THashList.h"
#include <iostream>
#include <fstream>
void Stats(){
  TString workingdir = TString(gSystem->WorkingDirectory());
  if(workingdir.Contains("paulbatzing")){workingdir.Clear();workingdir.Append("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");}
  else workingdir.Clear();
  gROOT->LoadMacro(Form("%sStatsfunc.cxx+g",workingdir.Data()));
  gPrintViaErrorHandler = kTRUE;  //makes minuit quiet
  
  TString dirs = TString(gSystem->GetFromPipe("ls -d */"));
  TString subdirs = TString();
  TString subsubdirs = TString();
  TString tok;
  TString tok1;
  TString tok2;
  Ssiz_t from  = 0;
  Ssiz_t from1 = 0;
  Ssiz_t from2 = 0;
  TFile * file;
  TFile * collectedfile;
  TList * folderlist;
  TString whattoprint = TString();  
  while (dirs.Tokenize(tok, from, "\n")) {
    if(gSystem->cd(tok.Strip(TString::kTrailing,'/').Data())){
      subdirs.Append(gSystem->GetFromPipe("ls -d */"));
      tok1.Clear();
      from1 = 0;
      while(subdirs.Tokenize(tok1,from1,"\n")){
	if(gSystem->cd(tok1.Strip(TString::kTrailing,'/').Data()))){
	    whattoprint.Append(gSystem->GetFromPipe("pwd"));
	    givemestats(whattoprint);
	    whattoprint.Clear();
	    if(gSystem->cd("runs")){
	    subsubdirs.Append(gSystem->GetFromPipe("ls -d */"));
	    tok2.Clear();
	    from2 = 0;
	    while(subsubdirs.Tokenize(tok2,from2,"\n")){
	      if(gSystem->cd(tok2.Strip(TString::kTrailing,'/').Data())){
		whattoprint.Append(gSystem->GetFromPipe("pwd"));
		givemestats(whattoprint);
		whattoprint.Clear();
		gSystem->cd("../");
	      }
	    }
	    subsubdirs.Clear();
	    gSystem->cd("../");
	  }
	  gSystem->cd("../");
	}
      }
      subdirs.Clear();
      gSystem->cd("../");
    }
  }

}