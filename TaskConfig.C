#include <AliFilteredEvent.h>
#include <AliFilteredTrack.h>
#include <AliCorrelation3p.h>
#include <AliCorrelation3p_noQA.h>
#include <AliThreeParticleCorrelator.h>
#include <AliLog.h>

#include <TFile.h>
#include <TArrayD.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
using namespace std;
class ThreePartTreeCorrelations;

void TaskConfig(int connr, ThreePartTreeCorrelations * task){
  //different configurations for the task on the tree:
  const char* name = "ThreePartTracksPbPb";
  Double_t MinTriggerPt = 4.0;
  Double_t MaxTriggerPt = 8.0;
  Double_t MinAssociatedPt = 0.5;
  Double_t MaxAssociatedPt = 1.0;
  Double_t Acceptancecut = 0.9;
  const char* period = "11h";
  Int_t MaxNEventsMix = 10;
  Int_t MinNTracksMix = 2000;
  Int_t NMBins = 6;
  Double_t Mbin0 = 0.;
  Double_t Mbin1 = 5.;
  Double_t Mbin2 = 10.;
  Double_t Mbin3 = 20.;
  Double_t Mbin4 = 40.;
  Double_t Mbin5 = 60.;
  Double_t Mbin6 = 90.;
  Double_t Mbin7 = 90.;
  const char * file = "LHC11hWeight.root";
  const char * cutmask = "GlobalHybrid";
  Int_t binver = 1;
  Int_t Zbinver =1;
  if(connr == 1){
   cout << connr<<endl;
  }
  
  
  
  //Actually apply all of this:
  task->SetMinTriggerPt(MinTriggerPt);
  task->SetMaxTriggerPt(MaxTriggerPt);
  task->SetMinAssociatedPt(MinAssociatedPt);
  task->SetMaxAssociatedPt(MaxAssociatedPt);
  task->SetAcceptanceCut(Acceptancecut);
  task->SetTrackCut(cutmask);
  task->SetBinVer(binver);
  task->SetQAtask(true);
   if(TString(file).CompareTo("")!=0)   task->SetWeights(Form("%s",file));

  //Mixing scheme:
  Double_t *Mbin = new Double_t[NMBins+1];
  Mbin[0] = Mbin0;
  Mbin[1] = Mbin1;
  if(NMBins>1) Mbin[2] = Mbin2;
  if(NMBins>2) Mbin[3] = Mbin3;
  if(NMBins>3) Mbin[4] = Mbin4;
  if(NMBins>4) Mbin[5] = Mbin5;
  if(NMBins>5) Mbin[6] = Mbin6;
  if(NMBins>6) Mbin[7] = Mbin7;
  TArrayD tMbin(NMBins+1, Mbin);
  //z vertex binning.
  Int_t NZBins   =   19;
  Double_t Zbin0 = -10.;
  Double_t Zbin1 = -8.5;
  Double_t Zbin2 = -7.5;
  Double_t Zbin3 = -6.5;
  Double_t Zbin4 = -5.5;
  Double_t Zbin5 = -4.5;
  Double_t Zbin6 = -3.5;
  Double_t Zbin7 = -2.5;
  Double_t Zbin8 = -1.5;
  Double_t Zbin9 = -0.5;
  Double_t Zbin10 = 0.5;
  Double_t Zbin11 = 1.5;						      
  Double_t Zbin12 = 2.5;						      
  Double_t Zbin13 = 3.5;						      
  Double_t Zbin14 = 4.5;						      
  Double_t Zbin15 = 5.5;						      
  Double_t Zbin16 = 6.5;						      
  Double_t Zbin17 = 7.5;						      
  Double_t Zbin18 = 8.5;						      
  Double_t Zbin19 = 10.;  
  if(Zbinver==2){
    NZBins=    9;
    Zbin0 = -10.;
    Zbin1 = -7.5;
    Zbin2 = -5.5;
    Zbin3 = -3.5;
    Zbin4 = -1.5;
    Zbin5 =  1.5;
    Zbin6 =  3.5;
    Zbin7 =  5.5;
    Zbin8 =  7.5;
    Zbin9 =  10.;
  }
  Double_t *Zbin = new Double_t[NZBins+1];
  Zbin[0] = Zbin0;
  Zbin[1] = Zbin1;
  if(NZBins>1) Zbin[2] = Zbin2;
  if(NZBins>2) Zbin[3] = Zbin3;
  if(NZBins>3) Zbin[4] = Zbin4;
  if(NZBins>4) Zbin[5] = Zbin5;
  if(NZBins>5) Zbin[6] = Zbin6;
  if(NZBins>6) Zbin[7] = Zbin7;
  if(NZBins>7) Zbin[8] = Zbin8;
  if(NZBins>8) Zbin[9] = Zbin9;
  if(NZBins>9) Zbin[10] = Zbin10;
  if(NZBins>10) Zbin[11] = Zbin11;
  if(NZBins>11) Zbin[12] = Zbin12;
  if(NZBins>12) Zbin[13] = Zbin13;
  if(NZBins>13) Zbin[14] = Zbin14;
  if(NZBins>14) Zbin[15] = Zbin15;
  if(NZBins>15) Zbin[16] = Zbin16;
  if(NZBins>16) Zbin[17] = Zbin17;
  if(NZBins>17) Zbin[18] = Zbin18;
  if(NZBins>18) Zbin[19] = Zbin19;
  TArrayD tZbin(NZBins+1, Zbin);  
  task->SetMixingScheme(MaxNEventsMix,MinNTracksMix,tMbin,tZbin);
  
  if( TString(period).Contains("10h") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P10h);
  if( TString(period).Contains("11h") )
    task->SetPeriod(AliAnalysisTaskCorrelation3p::P11h);
  
}