#include "TreeCorrelations.h"

//Macro to run over tree of AliFilteredEvents filled with AliFilteredTracks and build 3 particle correlations.

ThreePartTreeCorrelations::ThreePartTreeCorrelations()
  : TObject()
  , fname(TString("AliAnalysisTaskCorrelation3p"))
  , fOutfile(NULL)
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption("")
  , fRun(130848)
  , fVertex()
  , fperiod(ThreePartTreeCorrelations::P10h)
  , fCollisionType(ThreePartTreeCorrelations::PbPb)
  , fQA(kFALSE)
  , fqatask(kFALSE)
  , fWeights(NULL)
  , fWeightshpt(NULL)
  , fpTfunction(NULL)
  , fCorrelator(NULL)
  , fRunNumberList(NULL)
  , fNruns(1000)
  , fRunFillValue(0.0)
  , fMBinEdges(TArrayD())  
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fBinVer(0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)  
  , fAcceptancecut(0.8)
  , fMinTriggerPt(3.0)
  , fMaxTriggerPt(8.0)
  , fMinAssociatedPt(1.0)
  , fMaxAssociatedPt(3.0)
  , fMinNClustersTPC(70)
  , fCutMask(0)
{
  // default constructor
  // 
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-5,-2.5,2.5,5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
}

ThreePartTreeCorrelations::ThreePartTreeCorrelations(const char *name,const char* opt)
  : TObject()
  , fname(TString(name))
  , fOutfile(NULL)
  , fOutput(NULL)
  , fTextBox(NULL)
  , fOption(opt)
  , fRun(130848)
  , fVertex()
  , fperiod(ThreePartTreeCorrelations::P10h)
  , fCollisionType(ThreePartTreeCorrelations::PbPb)
  , fQA(kFALSE)
  , fqatask(kFALSE)
  , fWeights(NULL)
  , fWeightshpt(NULL)
  , fpTfunction(NULL)
  , fCorrelator(NULL)
  , fRunNumberList(NULL)
  , fNruns(1000)
  , fRunFillValue(0.0)
  , fMBinEdges(TArrayD())  
  , fZBinEdges(TArrayD())  
  , fMaxNEventMix(100)
  , fMinNofTracksMix(10)
  , fCentralityPercentile(0)
  , fMultiplicity(0)
  , fBinVer(0)
  , fMaxVz(10.0)
  , fMaxMult(0)
  , fMaxNumberOfTracksInPPConsidered(200)
  , fNTriggers(0)
  , fNAssociated(0)  
  , fAcceptancecut(0.8)
  , fMinTriggerPt(3.0)
  , fMaxTriggerPt(8.0)
  , fMinAssociatedPt(1.0)
  , fMaxAssociatedPt(3.0)
  , fMinNClustersTPC(70)
  , fCutMask(0)
{
  // default constructor
  // 
  Double_t Medges[6] = {0,20,40,60,80,90};
  TArrayD MBEdges(6, Medges);
  Double_t Zedges[6] = {-10,-5,-2.5,2.5,5,10};
  TArrayD ZBedges(6, Zedges);
  fMBinEdges = MBEdges;
  fZBinEdges = ZBedges;
}

int ThreePartTreeCorrelations::GetInput()
{
  fInfile = TFile::Open("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/LHC11h/167915/dstTree.root","READ");
  if(!fInfile)return 0;
  fInfile->GetDirectory("")->GetObject("DstTree",fTree);
  TBranch *events = fTree->GetBranch("Event");
  cout << fTree->GetEntries()<<endl;
  fEvent = new AliFilteredEvent();
  events->SetAddress(&fEvent);
  return 1;
}


void ThreePartTreeCorrelations::CreateOutputObjects()
{
  fOutfile = TFile::Open("AnalysisResults.root","RECREATE");
  // create result objects and add to output file
  MakeRunNumbers();//Needs to be done once in the beginning
  TH1::SetDefaultSumw2(kTRUE);//want the collection of weights on all histograms.
  TString collisiontype;
  if(fCollisionType==pp) collisiontype.Append("pp");
  if(fCollisionType==PbPb) collisiontype.Append("PbPb");
  fOutput = new THashList;
  fOutput->SetOwner();
  if(!fQA&&!fqatask){
    //Create the appropriate ThreeParticleCorrelators and add the used one to be fCorrelator.
    AliThreeParticleCorrelator<AliCorrelation3p>* correlator=new AliThreeParticleCorrelator<AliCorrelation3p>;
    fCorrelator=correlator;

    //Initialize QA histograms and add them to fOutput
    InitializeQAhistograms();
    //Intitialize the Multiplicity and ZVertex bins.
    const Int_t    MaxNofEvents=fMaxNEventMix;
    const Int_t    MinNofTracks=fMinNofTracksMix;
    const Int_t    nofMBins=fMBinEdges.GetSize()-1;
    Double_t 	 MBinsTemp[nofMBins+1];
    for(int i=0; i<=nofMBins; ++i) MBinsTemp[i] = fMBinEdges.At(i);
    const Int_t    nofZBins=fZBinEdges.GetSize()-1;//5;
    Double_t 	 ZBinsTemp[nofZBins+1];
    for(int i=0; i<=nofZBins; ++i) ZBinsTemp[i] = fZBinEdges.At(i);
    //Create the AliEventPoolManager 
    AliEventPoolManager* poolMgr = new AliEventPoolManager(MaxNofEvents, MinNofTracks, nofMBins, (Double_t*)MBinsTemp, nofZBins, (Double_t*)ZBinsTemp);
    poolMgr->SetTargetValues(MinNofTracks,1.0E-4,1.0);
    correlator->InitEventMixing(poolMgr);
    
    //initialize track worker and add to the output if appropriate
    TString tracksname = Form("tracktrigger_correlation_%.0f_%.0f", fMinTriggerPt, fMaxTriggerPt);
    TString triggertype = "tracks";
    TString triggerinit = Form("minTriggerPt=%.1f maxTriggerPt=%.1f minAssociatedPt=%.1f maxAssociatedPt=%.1f collisiontype=%s triggertype=%s", fMinTriggerPt, fMaxTriggerPt, fMinAssociatedPt, fMaxAssociatedPt,collisiontype.Data(),triggertype.Data());
    AliCorrelation3p* workertracktrigger =new AliCorrelation3p(tracksname, fMBinEdges, fZBinEdges);
    workertracktrigger->SetAcceptanceCut(fAcceptancecut);
    workertracktrigger->SetBinningVersion(fBinVer);
    workertracktrigger->Init(triggerinit);
    correlator->Add(workertracktrigger);
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 0));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 1));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 2));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 3));
    fOutput->Add(workertracktrigger);
    }
  if(!fQA&&fqatask){
    //Create the appropriate ThreeParticleCorrelators and add the used one to be fCorrelator.
    AliThreeParticleCorrelator<AliCorrelation3p_noQA>* correlator=new AliThreeParticleCorrelator<AliCorrelation3p_noQA>;
    fCorrelator=correlator;
    //Initialize QA histograms and add them to fOutput
    InitializeQAhistograms();
    //Intitialize the Multiplicity and ZVertex bins.
    const Int_t    MaxNofEvents=fMaxNEventMix;
    const Int_t    MinNofTracks=fMinNofTracksMix;
    const Int_t    nofMBins=fMBinEdges.GetSize()-1;
    Double_t 	 MBinsTemp[nofMBins+1];
    for(int i=0; i<=nofMBins; ++i) MBinsTemp[i] = fMBinEdges.At(i);
    const Int_t    nofZBins=fZBinEdges.GetSize()-1;//5;
    Double_t 	 ZBinsTemp[nofZBins+1];
    for(int i=0; i<=nofZBins; ++i) ZBinsTemp[i] = fZBinEdges.At(i);
    //Create the AliEventPoolManager 
    AliEventPoolManager* poolMgr = new AliEventPoolManager(MaxNofEvents, MinNofTracks, nofMBins, (Double_t*)MBinsTemp, nofZBins, (Double_t*)ZBinsTemp);
    poolMgr->SetTargetValues(MinNofTracks,1.0E-4,1.0);
    correlator->InitEventMixing(poolMgr);
    
    //initialize track worker and add to the output if appropriate
    TString tracksname = Form("tracktrigger_correlation_%.0f_%.0f", fMinTriggerPt, fMaxTriggerPt);
    TString triggertype = "tracks";
    TString triggerinit = Form("minTriggerPt=%.1f maxTriggerPt=%.1f minAssociatedPt=%.1f maxAssociatedPt=%.1f collisiontype=%s triggertype=%s", fMinTriggerPt, fMaxTriggerPt, fMinAssociatedPt, fMaxAssociatedPt,collisiontype.Data(),triggertype.Data());
    AliCorrelation3p_noQA* workertracktrigger =new AliCorrelation3p_noQA(tracksname, fMBinEdges, fZBinEdges);
    workertracktrigger->SetAcceptanceCut(fAcceptancecut);
    workertracktrigger->SetBinningVersion(fBinVer);
    workertracktrigger->Init(triggerinit);
    correlator->Add(workertracktrigger);
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 0));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 1));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 2));
    fOutput->Add(correlator->GetCorrespondingME(workertracktrigger, 3));
    fOutput->Add(workertracktrigger);
    }
  
  if(fQA){
    InitializeQAhistograms();
  }
}



ThreePartTreeCorrelations::~ThreePartTreeCorrelations()
{
  // destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }
  if(fRunNumberList) delete fRunNumberList;
  fRunNumberList=NULL;
}


void ThreePartTreeCorrelations::Correlate(){
  if(!GetInput())return;
  CreateOutputObjects();
  for(Int_t i=0;i<fTree->GetEntries();i++){
    if(i%10000==0)cout << i << endl;
    fTree->GetEvent(i);
    fRun = fEvent->GetRunNr();
    if(fQA){
      TAxis* runnumberaxis= dynamic_cast<TH1D*>(fOutput->FindObject("EventsperRun"))->GetXaxis();
      if (runnumberaxis){double RunBin = runnumberaxis->FindBin(Form("%i",fRun));fRunFillValue = runnumberaxis->GetBinCenter(RunBin);}   
    }
    GetCentralityAndVertex();
  //   if(!SelectEvent()) return;//events are rejected.
    FillHistogram("centVsZVertex",fCentralityPercentile,fVertex[2]);
    //initialize the period dependent cuts.
    if(fQA){
      //Fill Events/run histogram.
      FillHistogram("EventsperRun", fRunFillValue);
      FillHistogram("NEventsVertex",fRunFillValue,fVertex[2]);
      FillHistogram("NEventsCent",fRunFillValue,fCentralityPercentile);
    }
    //To fill with tracks and pions:
    TObjArray* allrelevantParticles = NULL;
    
    if(!fQA){
      allrelevantParticles = new TObjArray();
      allrelevantParticles->SetOwner();//In order for it to work in the event pool.
    }
    fNTriggers=0.0;//Reset fNTriggers
    fNAssociated=0.0;//Reset fNAssociated
    //Fill all the tracks
    GetTracks(allrelevantParticles, fEvent);
    
    FillHistogram("Ntriggers",fNTriggers);
    FillHistogram("NAssociated",fNAssociated);
    if(fQA){
      FillHistogram("NTriggersperRun",fRunFillValue,fNTriggers);
      FillHistogram("NAssociatedperRun",fRunFillValue,fNAssociated);
    }
    if(fNTriggers>=1)FillHistogram("NAssociatedETriggered",fNAssociated);
    //if fQA the correlations are not build.
    if(!fQA){
      if(!fqatask){dynamic_cast<AliThreeParticleCorrelator<AliCorrelation3p>*>(fCorrelator)->SetEventVzM(fVertex[2],fCentralityPercentile);}
      if(fqatask){dynamic_cast<AliThreeParticleCorrelator<AliCorrelation3p_noQA>*>(fCorrelator)->SetEventVzM(fVertex[2],fCentralityPercentile);}
      //Do the actual correlations.
      if(!((fNTriggers+fNAssociated)==0))fCorrelator->Execute(NULL, allrelevantParticles);//correlate for events that contain at least one trigger or associated.
      else delete allrelevantParticles;
      //Post the output
    }  
//     delete allrelevantParticles;
  }
  
  
  
  if(fOutfile){
    fOutfile->mkdir(fname.Data());
    fOutfile->GetDirectory(fname.Data())->cd();
  }
  if(fOutput)fOutput->Write();
  delete fOutput; fOutfile->Close();
  fInfile->Close();
//   delete fEvent; delete fTree;
}


Int_t ThreePartTreeCorrelations::GetTracks(TObjArray* allrelevantParticles, AliFilteredEvent *pEvent)
{
  Int_t nofTracks = 0;
  Int_t MultBin; Int_t VZbin;
  if(fWeights){
    MultBin = fWeights->GetXaxis()->FindBin(fCentralityPercentile);
    VZbin   = fWeights->GetYaxis()->FindBin(fVertex[2]);
  }
  nofTracks=pEvent->GetNtrks();
  FillHistogram("trackCount",nofTracks);
  for (int i=0; i<nofTracks; i++) {
    Double_t Weight = 1.0;
    AliFilteredTrack* t= dynamic_cast<AliFilteredTrack*>(pEvent->GetTrack(i));
    t->Calculate();
    if (!t) continue;
    if(fWeights){
      if(t->Pt()<2.0){
	Int_t pTbin  = fWeights->GetZaxis()->FindBin(t->Pt());
	Weight = fWeights->GetBinContent(MultBin,VZbin,pTbin);
      }
      else{
	Weight = fWeightshpt->GetBinContent(MultBin,VZbin)*fpTfunction->Eval(t->Pt());
      }
    }
    if(fQA){
      FillHistogram("TracksperRun",fRunFillValue);
    }
    FillHistogram("trackUnselectedPt",t->Pt());
    FillHistogram("trackUnselectedPhi",t->Phi());
    FillHistogram("trackUnselectedTheta",t->Theta());
    if (!IsSelected(t)) continue;
    if(allrelevantParticles){
	AliFilteredTrack * filp = new AliFilteredTrack(*t);
	filp->SetEff(Weight);      
      allrelevantParticles->Add(filp);
    }
    if(fQA){
      FillHistogram("selectedTracksperRun",fRunFillValue);
      FillHistogram("NTracksVertex",fRunFillValue,fVertex[2]);
      FillHistogram("NTracksCent",fRunFillValue,fCentralityPercentile);
    }
    FillHistogram("trackPt",t->Pt());
    if(fQA&&fWeights&&fWeightshpt){
      if(t->Pt()<4.0){
	FillHistogram("Track_Cent_Vertex_lpT",fCentralityPercentile,fVertex[2],t->Pt());
      }
      else{
	FillHistogram("Track_Cent_Vertex_eta",fCentralityPercentile,fVertex[2]);
      }
    }
    FillHistogram("trackPhi",t->Phi());
    FillHistogram("trackTheta",t->Theta());
    if(IsSelectedTrigger(t)){
      fNTriggers+=1;
      FillHistogram("trackTriggerPt",t->Pt());
      FillHistogram("trackTriggerPhi",t->Phi());
      FillHistogram("trackTriggerTheta",t->Theta());
    }
    if(IsSelectedAssociated(t)){
      fNAssociated+=1;
      FillHistogram("trackAssociatedPt",t->Pt());
      FillHistogram("trackAssociatedPhi",t->Phi());
      FillHistogram("trackAssociatedTheta",t->Theta());
    }
  }

  return nofTracks;
}



void ThreePartTreeCorrelations::SetMixingScheme(Int_t MaxNEventMix, Int_t MinNofTracksMix, TArrayD MBinEdges, TArrayD ZBinEdges)
{
  fMaxNEventMix= MaxNEventMix;
  fMinNofTracksMix = MinNofTracksMix;
  for(int i=0; i<MBinEdges.GetSize()-1; ++i)
    if(MBinEdges.At(i) > MBinEdges.At(i+1)) AliFatal("edges are not sorted");
  for(int i=0; i<ZBinEdges.GetSize()-1; ++i)
    if(ZBinEdges.At(i) > ZBinEdges.At(i+1)) AliFatal("edges are not sorted");  
  fMBinEdges = MBinEdges;
  fZBinEdges = ZBinEdges;
  fMaxMult = fMBinEdges.At(fMBinEdges.GetSize()-1);
}

void ThreePartTreeCorrelations::FillHistogram(const char* key, Double_t x)
{
  TH1 * hist = dynamic_cast<TH1*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

void ThreePartTreeCorrelations::FillHistogram(const char* key, Double_t x, Double_t y)
{
  TH2 * hist = dynamic_cast<TH2*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y) ;
  else if(dynamic_cast<TH1*>(fOutput->FindObject(key)))dynamic_cast<TH1*>(fOutput->FindObject(key))->Fill(x,y);
  else AliError(Form("can not find histogram (of instance TH2 or TH1) <%s> ",key)) ;
}

void ThreePartTreeCorrelations::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z)
{
  TH3 * hist = dynamic_cast<TH3*>(fOutput->FindObject(key)) ;
  if(hist)
    hist->Fill(x,y,z) ;
  else if(dynamic_cast<TH2*>(fOutput->FindObject(key)))dynamic_cast<TH2*>(fOutput->FindObject(key))->Fill(x,y,z);
  else AliError(Form("can not find histogram (of instance TH3 or TH2) <%s> ",key)) ;
}

void ThreePartTreeCorrelations::FillHistogram(const char* key, Double_t x, Double_t y, Double_t z,Double_t a, Double_t b)
{
  THnD * hist = dynamic_cast<THnD*>(fOutput->FindObject(key)) ;
  if(hist){
    Double_t s[5] = {x,y,z,a,b};
    hist->Fill(s) ;
  }
  else AliError(Form("can not find histogram (of instance TH3) <%s> ",key)) ;
}

Bool_t ThreePartTreeCorrelations::IsSelected(AliVParticle* p)
{
  //Performs selection cuts for tracks and triggers
  
  if (dynamic_cast<AliFilteredTrack*>(p) && IsSelectedTrack(p)) return IsSelectedTrigger(p)||IsSelectedAssociated(p);
  return kFALSE;
}
Bool_t ThreePartTreeCorrelations::IsSelectedTrack(AliVParticle* p)
{
  if(dynamic_cast<AliFilteredTrack*>(p)->IsGlobalHybrid()&&(fCutMask==0||fCutMask>3))return kTRUE;
  if(dynamic_cast<AliFilteredTrack*>(p)->IsBIT4()&&(fCutMask==1))return kTRUE;
  if(dynamic_cast<AliFilteredTrack*>(p)->IsBIT5()&&(fCutMask==2))return kTRUE;
  if(dynamic_cast<AliFilteredTrack*>(p)->IsBIT6()&&(fCutMask==3))return kTRUE;
  return kFALSE;
}
Bool_t ThreePartTreeCorrelations::IsSelectedTrigger(AliVParticle* p)
{
  if (p->Pt()<=fMinTriggerPt) return kFALSE;
  if (fMaxTriggerPt>fMinTriggerPt && p->Pt()>fMaxTriggerPt) return kFALSE;
  float etatrigger=p->Eta();
  if (etatrigger<=-fAcceptancecut || etatrigger>=fAcceptancecut) return kFALSE;
  return kTRUE;
}

Bool_t ThreePartTreeCorrelations::IsSelectedAssociated(AliVParticle* p)
{
  if (p->Pt()<=fMinAssociatedPt) return kFALSE;
  if (fMaxAssociatedPt>fMinAssociatedPt && p->Pt()>fMaxAssociatedPt) return kFALSE;
  float etaAssociated=p->Eta();
  if (etaAssociated<=-fAcceptancecut || etaAssociated>=fAcceptancecut) return kFALSE;
  return kTRUE;
}
void ThreePartTreeCorrelations::GetCentralityAndVertex()
{
  fCentralityPercentile = fEvent->GetCentrality();
 //Get the primary Vertex
  fVertex[0] = fEvent->GetfVertexX();
  fVertex[1] = fEvent->GetfVertexY();
  fVertex[2] = fEvent->GetfVertexZ();
}
void ThreePartTreeCorrelations::InitializeQAhistograms()
{
  //Function that initializes the QA histograms 
  if (!fOutput) return;
  //QA histograms
  fOutput->Add(new TH1D("trackCount", "trackCount", 1000,  0, 15000));
  fOutput->Add(new TH1D("trackUnselectedPt"   , "trackPt"   , 1000,  0, 20));
  fOutput->Add(new TH1D("trackPt"   			, "trackPt"   				, 1000,  0, 20));
  fOutput->Add(new TH1D("trackAssociatedPt" , "Pt of associated Tracks", 1000, fMinAssociatedPt, fMaxAssociatedPt));
  fOutput->Add(new TH1D("trackTriggerPt" , "Pt of Trigger Tracks", 1000, fMinTriggerPt, fMaxTriggerPt));
  fOutput->Add(new TH1D("trackUnselectedPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackTriggerPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackAssociatedPhi"  , "trackPhi"  ,  180,  0., 2*TMath::Pi()));
  fOutput->Add(new TH1D("trackUnselectedTheta", "trackTheta",  180, 0, TMath::Pi()));
  fOutput->Add(new TH1D("trackTheta", "trackTheta",  180, 0.0, TMath::Pi()));
  fOutput->Add(new TH1D("trackTriggerTheta", "trackTheta",  180, 0.0, TMath::Pi()));
  fOutput->Add(new TH1D("trackAssociatedTheta", "trackTheta",  180, 0.0, TMath::Pi()));
  fOutput->Add(new TH1D("Ntriggers","Number of triggers per event",50,-0.5,49.5));
  fOutput->Add(new TH1D("NAssociated","Number of Associated per event",200,-0.5,199.5));
  fOutput->Add(new TH1D("NAssociatedETriggered","Number of Associated per event that contains a trigger.",200,-0.5,199.5));
  if(fWeights)fOutput->Add(fWeights);
  if(fWeightshpt)fOutput->Add(fWeightshpt);
  if(fpTfunction)fOutput->Add(fpTfunction);

  if(fCollisionType==PbPb)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, 100, 100, -10, 10));
  if(fCollisionType==pp)fOutput->Add(new TH2D("centVsZVertex", "centvszvertex", 100, 0, fMaxNumberOfTracksInPPConsidered, 100, -10, 10));

  if(fQA){
    //QA per run histograms:
    TH1D * eventsperrun 	= new TH1D("EventsperRun", "# Events per Run", fNruns, 0, 1);
    TH1D * TracksperRun 	= new TH1D("TracksperRun", "# tracks per Run", fNruns, 0,1);
    TH1D * selectedTracksperRun = new TH1D("selectedTracksperRun", "# selected tracks per Run", fNruns, 0,1);
    TH2D * NTriggersperRun 	= new TH2D("NTriggersperRun","# triggers per event per Run",fNruns, 0,1,50,-0.5,49.5);
    TH2D * NAssociatedperRun 	= new TH2D("NAssociatedperRun","# associated per event per Run",fNruns, 0,1,100,-0.5,99.5);
    TH2D * NTracksVertex	= new TH2D("NTracksVertex","#selected tracks per run and vertex",fNruns,0,1,100,-10.0,10.0);
    TH2D * NEventsVertex	= new TH2D("NEventsVertex","Events per run and vertex",fNruns,0,1,100,-10.0,10.0);
    TH2D * NTracksCent	   	= new TH2D("NTracksCent","#selected tracks per run and vertex",fNruns,0,1,100,0.0,100.0);
    TH2D * NEventsCent	   	= new TH2D("NEventsCent","Events per run and vertex",fNruns,0,1,100,0.0,100.0);
    
    for(int i=0; i<fNruns; i++){
      TString lable = Form("%i",fRunNumberList[i]);
      eventsperrun->GetXaxis()->SetBinLabel(i+1, lable);
      eventsperrun->GetXaxis()->LabelsOption("v");
      TracksperRun->GetXaxis()->SetBinLabel(i+1, lable);
      TracksperRun->GetXaxis()->LabelsOption("v");
      selectedTracksperRun->GetXaxis()->SetBinLabel(i+1, lable);
      selectedTracksperRun->GetXaxis()->LabelsOption("v");
      NTriggersperRun->GetXaxis()->SetBinLabel(i+1,lable);
      NTriggersperRun->GetXaxis()->LabelsOption("v");
      NAssociatedperRun->GetXaxis()->SetBinLabel(i+1,lable);
      NAssociatedperRun->GetXaxis()->LabelsOption("v");    
      NTracksVertex->GetXaxis()->SetBinLabel(i+1,lable);
      NTracksVertex->GetXaxis()->LabelsOption("v"); 
      NEventsVertex->GetXaxis()->SetBinLabel(i+1,lable);
      NEventsVertex->GetXaxis()->LabelsOption("v"); 
      NTracksCent->GetXaxis()->SetBinLabel(i+1,lable);
      NTracksCent->GetXaxis()->LabelsOption("v"); 
      NEventsCent->GetXaxis()->SetBinLabel(i+1,lable);
      NEventsCent->GetXaxis()->LabelsOption("v"); 
    }
    fOutput->Add(eventsperrun);
    fOutput->Add(TracksperRun);
    fOutput->Add(selectedTracksperRun);
    fOutput->Add(NTriggersperRun);
    fOutput->Add(NAssociatedperRun);
    fOutput->Add(NTracksVertex);
    fOutput->Add(NEventsVertex);
    fOutput->Add(NTracksCent);
    fOutput->Add(NEventsCent);
  }
  if(fQA&&fWeights&&fWeightshpt){
    TH3D * histtrackslpt = (TH3D*)(fWeights->Clone("Track_Cent_Vertex_lpT"));
    histtrackslpt->Reset();
    histtrackslpt->SetTitle("Tracks in Centrality vs Vertex vs pT");
    fOutput->Add(histtrackslpt);
    TH2D * histtrackshpt = (TH2D*)(fWeightshpt->Clone("Track_Cent_Vertex_eta"));
    histtrackshpt->Reset();    
    histtrackshpt->SetTitle("Tracks in Centrality vs Vertex");
    fOutput->Add(histtrackshpt);    
  }
}
void ThreePartTreeCorrelations::MakeRunNumbers()
{
  //Initialize array for run numbers.
  Int_t runnumbersP10b[fNRunsP10b] = {117222, 117220, 117116, 117112, 117109, 117099, 117092, 117063, 117060, 117059, 117053, 117052, 117050, 117048, 116787, 116645, 116643, 116574, 116571, 116562, 116432, 116431, 116429,116403, 116402, 116372, 116360, 116358, 116288, 116102, 116081, 116079, 115521, 115414, 115406, 115401, 115399, 115393, 115369,115345, 115335, 115328,  115327, 115322, 115318, 115312, 115310,  115193, 115186, 115056, 114931, 114930, 114924, 114920, 114918, 114798, 114786};
  Int_t runnumbersP10c[fNRunsP10c] = {121040, 121039, 120829, 120825, 120824, 120823, 120822, 120821, 120820, 120758, 120750, 120741, 120671, 120617, 120616, 120505, 120504, 120503, 120244,120079, 120076, 120073, 120072, 120069, 120067, 119862, 119859, 119856, 119853, 119849, 119846, 119845, 119844, 119842, 119841, 119163, 119161, 119159, 118561, 118560, 118558, 118556, 118518, 118512, 118507, 118506};
  Int_t runnumbersP10d[fNRunsP10d] = {126432, 126425, 126424, 126422, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126350, 126285, 126284, 126283, 126168,126167, 126160, 126158, 126097, 126090, 126088, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843,  125842, 125633, 125632, 125630, 125628, 125296, 125295, 125186, 125156, 125140, 125139, 125134, 125133, 125101, 125100, 125097, 125085, 125083, 125023, 124751, 122375, 122374};
  Int_t runnumbersP10e[fNRunsP10e] = {130850, 130848, 130847, 130844, 130842, 130840, 130834, 130804, 130803, 130802, 130799, 130798, 130795, 130793, 130704, 130696, 130628, 130623, 130621, 130620, 130609, 130608, 130601, 130526, 130524, 130520, 130519, 130517, 130481, 130480, 130479, 130375, 130360, 130358, 130356, 130354, 130343, 130342, 130178, 130172, 130168, 130158, 130157, 130151, 130149, 129983, 129966, 129962, 129961, 129960, 129959, 129744, 129742, 129738, 129736, 129735, 129734, 129729, 129726, 129725, 129723, 129666, 129659, 129653, 129652, 129651, 129650, 129647, 129641, 129639, 129599, 129587, 129586, 129540, 129536, 129528, 129527, 129525, 129524, 129523, 129521, 129520, 129519, 129516, 129515, 129514, 129513, 129512, 129042, 128913, 128855, 128853, 128850, 128843, 128836, 128835, 128834, 128833, 128824, 128823, 128820, 128819, 128778, 128777, 128678, 128677, 128621, 128615, 128611, 128609, 128605, 128596, 128594, 128592, 128590, 128582, 128506, 128505, 128504, 128503, 128498, 128495, 128494, 128486, 128452, 128366}; 
  Int_t runnumbersP10h[fNRunsP10h] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137366, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161, 137135};
  Int_t runnumbersP11a[fNRunsP11a] = {146860, 146859, 146858, 146856, 146824, 146817, 146807, 146806, 146805, 146804, 146803, 146802, 146801, 146748, 146747, 146746, 146402, 146369,146292, 146287, 146282, 146277, 146273, 146272, 146223, 146220, 146208, 146158, 146156, 146153, 146152, 146148, 146147, 146141, 146099, 146079, 146072, 146071, 146027, 146026, 146025, 146024, 146023, 145674, 145455, 145385, 145384, 145383, 145379, 145355, 145354, 145353, 145314, 145300, 145292, 145290, 145289, 145288};
  Int_t runnumbersP11h[fNRunsP11h] = {170593, 170572, 170388, 170387, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 170040, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167988, 167987, 167985, 167920, 167915};

  if (fperiod == ThreePartTreeCorrelations::P10b){
    fNruns = fNRunsP10b;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10b[i];
    //Set the correct collision type
    fCollisionType = ThreePartTreeCorrelations::pp;
  }
  else if (fperiod == ThreePartTreeCorrelations::P10c){
    fNruns = fNRunsP10c;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10c[i];
    //Set the correct collision type
    fCollisionType = ThreePartTreeCorrelations::pp;
  }
  else if (fperiod == ThreePartTreeCorrelations::P10d){
    fNruns = fNRunsP10d;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10d[i];
    //Set the correct collision type
    fCollisionType = ThreePartTreeCorrelations::pp;
  }
  else if (fperiod == ThreePartTreeCorrelations::P10e){
    fNruns = fNRunsP10e;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10e[i];
    //Set the correct collision type
    fCollisionType = ThreePartTreeCorrelations::pp;
  }
  else if (fperiod == ThreePartTreeCorrelations::P11a){
    fNruns = fNRunsP11a;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP11a[i];
    //Set the correct collision type
    fCollisionType = ThreePartTreeCorrelations::pp;
  }
  else if (fperiod ==  ThreePartTreeCorrelations::P10h){
    fNruns = fNRunsP10h;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP10h[i];
    //Set the correct collision type
    fCollisionType = ThreePartTreeCorrelations::PbPb;
  }
  else if (fperiod == ThreePartTreeCorrelations::P11h){
    fNruns = fNRunsP11h;
    fRunNumberList = new Int_t[fNruns];
    for(int i = 0; i<fNruns; i++) fRunNumberList[i] = runnumbersP11h[i];
    //Set the correct collision type
    fCollisionType = ThreePartTreeCorrelations::PbPb;
  }
  if(fCollisionType==pp)for(int i=0;i<fMBinEdges.GetSize();i++){
      fMBinEdges.AddAt(fMaxNumberOfTracksInPPConsidered*fMBinEdges.At(i)/fMBinEdges.At(fMBinEdges.GetSize()-1),i);
    }
  return ;
}
