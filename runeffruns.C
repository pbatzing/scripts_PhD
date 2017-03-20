#include <iostream>
#include <stdio.h>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TSystemDirectory.h"


using namespace std;
void runeffruns(const char* options = "", const char* o2 = "")
{
  
  Double_t MaxlowpT = 2.0;
  Double_t Cutoffpoint = 5.0;
  TString basedir=TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  //   MakeTestHists();
  if(TString(options).Contains("effruns ")||TString(options).Contains("effrun=")||TString(options).Contains("effrunsg")||TString(options).Contains("effrunsh")||TString(options).Contains("effrunsi")){
    gROOT->LoadMacro(Form("%seffmaker.cxx+g",basedir.Data()));
    TObjArray * dir = new TObjArray();
    dir->SetOwner(kTRUE);
    if((TString(options).Contains("effruns ")||TString(options).Contains("effrun=")||TString(options).Contains("effrunsg"))&&TString(options).Contains("LHC11h")){
      dir->Add(new TSystemDirectory("LHC12a17g","/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/MCproductions/LHC12a17g"));
    }
    if((TString(options).Contains("effruns ")||TString(options).Contains("effrun=")||TString(options).Contains("effrunsh"))&&TString(options).Contains("LHC11h")){
      dir->Add(new TSystemDirectory("LHC12a17h","/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/MCproductions/LHC12a17h"));
    }
    if((TString(options).Contains("effruns ")||TString(options).Contains("effrun=")||TString(options).Contains("effrunsi"))&&TString(options).Contains("LHC11h")){
      dir->Add(new TSystemDirectory("LHC12a17i","/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/MCproductions/LHC12a17i"));
    }
    if(TString(options).Contains("LHC10h")){
      dir->Add(new TSystemDirectory("LHC11a10a","/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/MCproductions/LHC11a10a"));      
    }
    TSystemDirectory * sysdir;
    TString * filename;
    TString * nowfilename;
    for(int i = 0; i<dir->GetEntries();i++){

      sysdir = dynamic_cast<TSystemDirectory*>(dir->At(i));
      TList * directs = sysdir->GetListOfFiles();
      filename = new TString(basedir.Data());
      filename->Append("MCproductions/");
      filename->Append(dir->At(i)->GetName());
      filename->Append("/");
      for(int j =0; j<directs->GetEntries();j++){
	nowfilename = new TString(filename->Data());
	TSystemFile *file = dynamic_cast<TSystemFile*>(directs->At(j)); 
	if(file->IsDirectory()){
	  if(!TString(file->GetName()).BeginsWith(".")){
	    if(TString(options).Contains("effrun=")){
	      bool isit = false;
	      TObjArray * optionsArray = TString(options).Tokenize(" ");
	      for(int k = 0; k<optionsArray->GetEntries();k++){
		if(TString(optionsArray->At(k)->GetName()).BeginsWith("effrun=")){
		  TObjArray * subArray = TString(optionsArray->At(k)->GetName()).Tokenize("=");
		 if(TString(TString( subArray->At(1)->GetName())).CompareTo(file->GetName()) == 0) isit = true;
		 delete subArray;
		}
	      }
	      delete optionsArray;
	      if(!isit)continue;
	    }
	    cout << file->GetName()<<endl;	    
	    nowfilename->Append(file->GetName());
	    nowfilename->Append("/AnalysisResults.root");
	    cout << nowfilename->Data()<<endl;
	    MakeEffHistsPbPb(nowfilename->Data());
	  }
	}
      delete nowfilename;
      }
    
    delete filename;
    }
    delete dir;
  }
  if(TString(options).Contains("CollectCentBins")){
    gROOT->LoadMacro(Form("%seffruns.cxx+g",basedir.Data()));
    CollectCentBins(o2);
  }
  if(TString(options).Contains("CollectRuns")){
    gROOT->LoadMacro(Form("%seffruns.cxx+g",basedir.Data()));    
    CollectRuns(MaxlowpT,Cutoffpoint);
  }
  if(TString(options).Contains("RunStatsPlus")){
    gROOT->LoadMacro(Form("%seffruns.cxx+g",basedir.Data()));    
    RunStatsPlus();
  }
}