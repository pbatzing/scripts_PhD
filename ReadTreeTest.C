#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "AliFilteredEvent.h"
#include "AliFilteredTrack.h"

using namespace std;
void ReadTreeTest(){
//   gROOT->LoadMacro("AliFilteredEvent.cxx+"); 
//   gROOT->ProcessLine(".include $ALICE_PHYSICS/include $ALICE_PHYSICS/lib");
//   gSystem->Load("libPWGCFCorrelationsThreePart.so")
//   gROOT->ProcessLine(".L AliFilteredTrack.cxx");
  TFile *f = TFile::Open("dstTree.root","READ");
  if(!f)return;
  TDirectory* dir = f->GetDirectory("");
  TTree* T;
  dir->GetObject("DstTree",T);
  TBranch *events = T->GetBranch("Event");
  cout << T->GetEntries()<<endl;
  AliFilteredEvent* event = new AliFilteredEvent();
  events->SetAddress(&event);

  
  for(Int_t i=0;i<T->GetEntries();i++){
    T->GetEvent(i);
    TClonesArray*arr = event->GetTracks();
    int globt = 0;
    int bit4 = 0;
    int bit5 = 0;
    int bit6 = 0;
    for(int j=0;j<arr->GetEntries();j++){
      if(dynamic_cast<AliFilteredTrack*>(arr->At(j))->IsGlobalHybrid())globt+=1;
      if(dynamic_cast<AliFilteredTrack*>(arr->At(j))->IsBIT4())bit4+=1;
      if(dynamic_cast<AliFilteredTrack*>(arr->At(j))->IsBIT5())bit5+=1;
      if(dynamic_cast<AliFilteredTrack*>(arr->At(j))->IsBIT6())bit6+=1;
    }
    cout << event->GetNtrks()<<" "<< globt<< " " << bit4 << " " << bit5 << " " << bit6 <<endl;
  }
  f->Close();
  
  
}