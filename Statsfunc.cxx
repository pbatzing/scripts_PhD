#include "TString.h"
#include "TArrayD.h"
#include "TFile.h"
#include "THn.h"
#include "TH3D.h"
#include "TF1.h"
#include <iosfwd>
void givemestats(TString locfile){
  TFile * file = TFile::Open(Form("%s/AnalysisResults.root",locfile.Data()), "READ");
  TString whattoprint = TString();
  char* s = new char[1];
  if(file){

    TList* folderlist = file->GetListOfKeys();
    for(int i = 0;i<folderlist->GetEntries();i++){
      if(!TString(folderlist->At(i)->GetName()).BeginsWith("ThreePart"))continue;	
      const char* foldername = folderlist->At(i)->GetName();
      TDirectory * folderdir = file->GetDirectory(foldername);
      TList * Containerlist = folderdir->GetListOfKeys();		  
      for(int j=0;j<Containerlist->GetEntries();j++){
	 const char* containername = Containerlist->At(j)->GetName();
	if(TString(containername).BeginsWith("ThreePart")){	
	  TList * containterdir = dynamic_cast<TList*>(folderdir->Get(containername));
	  TH1D* hist = dynamic_cast<TH1D*>(containterdir->FindObject("trackCount"));
	  whattoprint.Append("Trigger pT ");
	  if(locfile.Contains("t16"))whattoprint.Append("16 - 50 GeV/c");
	  if(locfile.Contains("t816"))whattoprint.Append(" 8 - 16 GeV/c");
	  if(locfile.Contains("t48"))whattoprint.Append(" 4 -  8 GeV/c");
	  
	  if(locfile.Contains("a23"))whattoprint.Append(" Associated pT  2 -  3 GeV/c");
	  if(locfile.Contains("a34"))whattoprint.Append(" Associated pT  3 -  4 GeV/c");
	  if(locfile.Contains("a46"))whattoprint.Append(" Associated pT  4 -  6 GeV/c");
	  if(locfile.Contains("a68"))whattoprint.Append(" Associated pT  6 -  8 GeV/c");
	  if(locfile.Contains("a816"))whattoprint.Append(" Associated pT  8 - 16 GeV/c");
	  if(locfile.Contains("a16"))whattoprint.Append(" Associated pT 16 - 50 GeV/c");

	  if(locfile.Contains("runs/01")) whattoprint.Append(" division 1:");
	  if(locfile.Contains("runs/02")) whattoprint.Append(" division 2:");
	  if(locfile.Contains("runs/03")) whattoprint.Append(" division 3:");
	  if(locfile.Contains("runs/04")) whattoprint.Append(" division 4:");
	  if(locfile.Contains("runs/05")) whattoprint.Append(" division 5:");
	  if(locfile.Contains("runs/06")) whattoprint.Append(" division 6:");
	  cout <<whattoprint.Data()<< " " << hist->GetEntries()<<endl;
	  delete containterdir;
	  whattoprint.Clear();
	}
      }
    }
//     delete folderlist;
    file->Close();
  }
}