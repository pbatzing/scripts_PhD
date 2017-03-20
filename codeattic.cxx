//This file contains old versions of code.

//from AliCorrelation3p.h:
  int ProcessEventMixing(const AliCorrelation3p* pCorrME, Bool_t Divide=false);
  int ProcessEventMixingmergefirst(const AliCorrelation3p* pCorrME);
  TH1* AverageHists(Int_t Hist, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin, const char* name = "mergehist") const;  
  TH1D* CollectHists1d(Int_t Hist, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin, const char* name="collection1d") const;
  TH2D* CollectHists2d(Int_t Hist, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin, const char* name="collection2d") const;
  TH3D* CollectHists3d(Int_t Hist, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin, const char* name="collection3d") const;
/ overloaded from TObject: draw histograms
  virtual void Draw(Option_t *option=""){Draw(0,-1,0,-1,"",false,option);} // *MENU*
  void Draw( Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin,const char* codir,Bool_t bDivide=false,Option_t *option="");
  virtual void Drawmergefirst(Option_t *option=""); // *MENU*
  virtual void DrawBinStatistics(Option_t *option="");
  
  /// draw set of histograms
  void DrawHistograms(const TObjArray* histograms, const char* options, Int_t minMbin,Int_t maxMbin, Int_t minZbin,Int_t maxZbin) const;
  /// set properties and draw histogram
  void DrawHistogram(const char* name, TObject* pObj, const char* drawoption, const char* titleX, const char* titleY=NULL, const char* titleZ=NULL) const;
  /// set properties and draw histogram
  void DrawHistogram(TObject* pObj, const char* drawoption, const char* titleX, const char* titleY=NULL, const char* titleZ=NULL) const {
    DrawHistogram(NULL, pObj, drawoption, titleX, titleY, titleZ);
  }
  /// draw delta eta histogram in different modes
  void DrawDeltaEta(const TObjArray* histograms, int index, int canvasno,const char* dir, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin,Bool_t eachbin=kFALSE,Bool_t average=kFALSE) const;
  void DrawDeltaEtamergedfirst(const TObjArray* histograms, int index, int canvasno,const char* dir) const;
  void SetAxisTitles(TObject* pObj, const char* titleX, const char* titleY=NULL, const char* titleZ=NULL) const;
  void CollectHistsFile(const char* name);



//from AliCorrelation3p.cxx
  // ///FROM HERE ON EVERYTHING IS USED FOR DRAWING HISTOGRAMS ONLY:
// //Only some stats.
// void AliCorrelation3p::DrawBinStatistics(Option_t* options)
// {
//   //loops over all bins, and draws statistics for number of events, number of triggers and number of associated.
//   TString delimiter(" ");
//   TString directory("");
//   TStringToken token(options, delimiter);
//   while (token.NextToken()) {
//     const char* key=0;
//     TString argument=token;
//     if (argument.CompareTo("-h")==0 ||
// 	argument.CompareTo("--help")==0) {
//       AliInfo(Form("Draw options:"
// 		   "\n\t  dir=<directory> saves all plots in the given directory"
// 		   ));
//       return;
//     }
//     key="dir=";
//     if (argument.BeginsWith(key)) {
//       argument.ReplaceAll(key, "");
//       directory.Append(argument);
// 
//     }    
//   }
//   
//   TString MultAxis = TString("Multiplicity [");
//   if(fCollisionType==pp)MultAxis.Append("# tracks]");
//   if(fCollisionType==PbPb||fCollisionType==pPb)MultAxis.Append("%]");
// 
//   
//   const double Pii=TMath::Pi();
//   TCanvas* c=NULL;
//   c=new TCanvas("BinStats");
//   c->SetWindowSize(1600,1114);
//   c->Divide(1,2);
//   c->cd(1);DrawHistogram("Distribution of Vz",((TH2D*)fHistograms->At(kcentrvsvz))->ProjectionY(),"E","Vz[cm]");
//   c->cd(2);DrawHistogram("Distribution of Multiplicity",((TH2D*)fHistograms->At(kcentrvsvz))->ProjectionX(),"E",MultAxis.Data());
//   c->Print(Form("%s%s",directory.Data(),"eventstats.eps"));
//   c->Print(Form("%s%s",directory.Data(),"eventstats.png"));
//   c->SaveAs(Form("%s%s",directory.Data(),"eventstats.root"));
//   
//   TH1D* Ntriggerscent = new TH1D("hntriggerscent","Number of triggers per centrality",fMBinEdges.GetSize()-1,fMBinEdges.GetArray());
//   TH1D* Ntriggersvz = new TH1D("hntriggersvz","Number of triggers per vz",fZBinEdges.GetSize()-1,fZBinEdges.GetArray());
//   TH1D* Nassociatedcent = new TH1D("hnassociatedcent","Number of triggers per centrality",fMBinEdges.GetSize()-1,fMBinEdges.GetArray());
//   TH1D* Nassociatedvz = new TH1D("hnassociatedvz","Number of triggers per vz",fZBinEdges.GetSize()-1,fZBinEdges.GetArray());  
//   TH1D* pTTriggerCheck =  new TH1D("hTpTcheck","Trigger pT",100,0.0,(fMaxTriggerPt>fMinTriggerPt?fMaxTriggerPt:fMinTriggerPt));
//   TH1D* pTAssociatedCheck =  new TH1D("hApTcheck","Associated pT",100,0.0,(fMaxTriggerPt>fMinTriggerPt?fMaxTriggerPt:fMinTriggerPt));
//   TH1D* phiTriggerCheck =  new TH1D("hTphicheck","Trigger phi",270,-.5*Pii ,2.5*Pii);
//   TH1D* phiAssociatedCheck =  new TH1D("hAphicheck","Associated phi",270,-.5*Pii ,2.5*Pii);
//   TH1D* etaTriggerCheck =  new TH1D("hTetacheck","Trigger eta",100,-3.0,3.0);
//   TH1D* etaAssociatedCheck =  new TH1D("hAetacheck","Associated eta",100,-3.0,3.0);
//   TH1D* eventsMb = new TH1D("hneventsmbin","Number of events per centrality",fMBinEdges.GetSize()-1,fMBinEdges.GetArray());
//   TH1D* eventsVz = new TH1D("hneventsvz"  ,"Number of events per vz"        ,fZBinEdges.GetSize()-1,fZBinEdges.GetArray());
//   for(int mbin = 0;mbin<(fMBinEdges.GetSize()-1);mbin++){
//     for(int zbin= 0;zbin<(fZBinEdges.GetSize()-1);zbin++){
//       double MBinfillvalue = fMBinEdges.At(mbin)+0.5*(fMBinEdges.At(mbin+1)-fMBinEdges.At(mbin));
//       double ZBinfillvalue = fZBinEdges.At(zbin)+0.5*(fZBinEdges.At(zbin)-fZBinEdges.At(zbin));
//       Ntriggerscent->Fill(MBinfillvalue, ((TH1F*)fHistograms->At(GetNumberHist(kHistTriggerpT,mbin,zbin)))->Integral());
//       Ntriggersvz->Fill(ZBinfillvalue, ((TH1F*)fHistograms->At(GetNumberHist(kHistTriggerpT,mbin,zbin)))->Integral());
//       Nassociatedcent->Fill(MBinfillvalue, ((TH1F*)fHistograms->At(GetNumberHist(kHistAssociatedpT,mbin,zbin)))->Integral());
//       Nassociatedvz->Fill(ZBinfillvalue, ((TH1F*)fHistograms->At(GetNumberHist(kHistAssociatedpT,mbin,zbin)))->Integral());
//       pTTriggerCheck->Add((TH1F*)fHistograms->At(GetNumberHist(kHistTriggerpT,mbin,zbin)));
//       pTAssociatedCheck->Add((TH1F*)fHistograms->At(GetNumberHist(kHistAssociatedpT,mbin,zbin)));
//       phiTriggerCheck->Add((TH1F*)fHistograms->At(GetNumberHist(kHistTriggerPhi,mbin,zbin)));
//       phiAssociatedCheck->Add((TH1F*)fHistograms->At(GetNumberHist(kHistAssociatedPhi,mbin,zbin)));    
//       etaTriggerCheck->Add((TH1F*)fHistograms->At(GetNumberHist(kHistTriggerEta,mbin,zbin)));
//       etaAssociatedCheck->Add((TH1F*)fHistograms->At(GetNumberHist(kHistAssociatedEta,mbin,zbin)));
//       eventsMb->Fill(MBinfillvalue,((TH2F*)fHistograms->At(GetNumberHist(khQAtocheckadressing,mbin,zbin)))->Integral());
//       eventsVz->Fill(ZBinfillvalue,((TH2F*)fHistograms->At(GetNumberHist(khQAtocheckadressing,mbin,zbin)))->Integral());
//     }
//   }
//   c->Clear();
//   c->Divide(2,2);
//   c->cd(1);
//   DrawHistogram("Number of triggers per centrality",Ntriggerscent,"E", MultAxis.Data(),"# triggers");
//   c->cd(2);
//   DrawHistogram("Number of Associated per centrality",Nassociatedcent,"E",MultAxis.Data(),"# Associated");
//   c->cd(3);
//   DrawHistogram("Number of triggers per Vz",Ntriggersvz,"E","Vz [cm]","# triggers");
//   c->cd(4);
//   DrawHistogram("Number of Associated per Vz",Nassociatedvz,"E","Vz [cm]","# Associated");
//   c->Print(Form("%s%s",directory.Data(),"Trackstats.eps"));
//   c->Print(Form("%s%s",directory.Data(),"Trackstats.png"));
//   c->SaveAs(Form("%s%s",directory.Data(),"Trackstats.root"));
//   c->Clear();
//   c->Divide(2,3);
//   c->cd(1);gPad->SetLogy();
//   DrawHistogram("Triggers when summing over all vertex and multiplicity bins.",pTTriggerCheck,"E","pT[GeV/c]","# triggers");
//   c->cd(3);gPad->SetLogy();
//   DrawHistogram("All Triggers found.",(TH1F*)fHistograms->At(kHistpTTriggerallbins),"E","pT[GeV/c]","# triggers");
//   c->cd(2);gPad->SetLogy();
//   DrawHistogram("Associated summed over all vertex and multiplicity",pTAssociatedCheck,"E","pT[GeV/c]","# Associated");
//   c->cd(4);gPad->SetLogy();
//   DrawHistogram("All Associated Found.",(TH1F*)fHistograms->At(kHistpTAssociatedallbins),"E","pT[GeV/c]","# Associated");
//   TH1D* RatioTrig = (TH1D*)pTTriggerCheck->Clone("trigclone");
//   RatioTrig->Divide((TH1F*)fHistograms->At(kHistpTTriggerallbins));
//   c->cd(5);
//   DrawHistogram("Ratio (supposed to be one if there are no events with exactly one trigger and one associated)",RatioTrig,"","pT[GeV/c]");
//   TH1D* RatioAss = (TH1D*)pTAssociatedCheck->Clone("assclone");
//   RatioAss->Divide((TH1F*)fHistograms->At(kHistpTAssociatedallbins));
//   c->cd(6);
//   DrawHistogram("Ratio (supposed to be one if there are no events with exactly one trigger and one associated)",RatioAss,"","pT[GeV/c]");
//   c->Print(Form("%s%s",directory.Data(),"pTcomparison.eps"));
//   c->Print(Form("%s%s",directory.Data(),"pTcomparison.png"));
//   c->SaveAs(Form("%s%s",directory.Data(),"pTcomparison.root"));
//   c->Clear();
//   c->Divide(2,3);
//   c->cd(1);gPad->SetLogy();
//   DrawHistogram("Triggers when summing over all vertex and multiplicity bins.",phiTriggerCheck,"E","#phi [rad]","# triggers");
//   c->cd(3);gPad->SetLogy();
//   DrawHistogram("All Triggers found.",(TH1F*)fHistograms->At(kHistPhiTriggerallbins),"E","#phi [rad]","# triggers");
//   c->cd(2);gPad->SetLogy();
//   DrawHistogram("Associated summed over all vertex and multiplicity",phiAssociatedCheck,"E","#phi [rad]","# Associated");
//   c->cd(4);gPad->SetLogy();
//   DrawHistogram("All Associated Found.",(TH1F*)fHistograms->At(kHistPhiAssociatedallbins),"E","#phi [rad]","# Associated");
//   TH1D* RatioTrigphi = (TH1D*)phiTriggerCheck->Clone("trigclonephi");
//   RatioTrigphi->Divide((TH1F*)fHistograms->At(kHistPhiTriggerallbins));
//   c->cd(5);
//   DrawHistogram("Ratio (supposed to be one if there are no events with exactly one trigger and one associated)",RatioTrigphi,"","#phi [rad]");
//   TH1D* RatioAssphi = (TH1D*)phiAssociatedCheck->Clone("assclonephi");
//   RatioAssphi->Divide((TH1F*)fHistograms->At(kHistPhiAssociatedallbins));
//   c->cd(6);
//   DrawHistogram("Ratio (supposed to be one if there are no events with exactly one trigger and one associated)",RatioAssphi,"","#phi [rad]");
//   c->Print(Form("%s%s",directory.Data(),"phicomparison.eps"));
//   c->Print(Form("%s%s",directory.Data(),"phicomparison.png"));
//   c->SaveAs(Form("%s%s",directory.Data(),"phicomparison.root"));
//   c->Clear();
//   c->Divide(2,3);
//   c->cd(1);gPad->SetLogy();
//   DrawHistogram("Triggers when summing over all vertex and multiplicity bins.",etaTriggerCheck,"E","#eta []","# triggers");
//   c->cd(3);gPad->SetLogy();
//   DrawHistogram("All Triggers found.",(TH1F*)fHistograms->At(kHistEtaTriggerallbins),"E","#eta []","# triggers");
//   c->cd(2);gPad->SetLogy();
//   DrawHistogram("Associated summed over all vertex and multiplicity",etaAssociatedCheck,"E","#eta []","# Associated");
//   c->cd(4);gPad->SetLogy();
//   DrawHistogram("All Associated Found.",(TH1F*)fHistograms->At(kHistEtaAssociatedallbins),"E","#eta []","# Associated");
//   TH1D* RatioTrigeta = (TH1D*)etaTriggerCheck->Clone("trigcloneeta");
//   RatioTrigeta->Divide((TH1F*)fHistograms->At(kHistEtaTriggerallbins));
//   c->cd(5);
//   DrawHistogram("Ratio (supposed to be one if there are no events with exactly one trigger and one associated)",RatioTrigeta,"","#eta []");
//   TH1D* RatioAsseta = (TH1D*)etaAssociatedCheck->Clone("asscloneeta");
//   RatioAsseta->Divide((TH1F*)fHistograms->At(kHistEtaAssociatedallbins));
//   c->cd(6);
//   DrawHistogram("Ratio (supposed to be one if there are no events with exactly one trigger and one associated)",RatioAsseta,"","#eta []");
//   c->Print(Form("%s%s",directory.Data(),"etacomparison.eps"));
//   c->Print(Form("%s%s",directory.Data(),"etacomparison.png"));
//   c->SaveAs(Form("%s%s",directory.Data(),"etacomparison.root"));
//   c->Clear();
//   c->Divide(2,2);
//   c->cd(1);
//   DrawHistogram("Number of Events per multiplicity (collected from bins)",eventsMb,"E",MultAxis.Data(),"#events");
//   c->cd(2);
//   DrawHistogram("Number of Events per Vz (collected from bins)",eventsVz,"E",MultAxis.Data(),"#events");
//   c->cd(3);
//   DrawHistogram("Number of Events per multiplicity.",((TH2D*)fHistograms->At(kcentrvsvzbin))->ProjectionX(),"E","Vz [cm]");
//   c->cd(4);
//   DrawHistogram("Number of Events per Vz",((TH2D*)fHistograms->At(kcentrvsvzbin))->ProjectionY(),"E","Vz[cm]");
//   c->Print(Form("%s%s",directory.Data(),"BinComparison.eps"));
//   c->Print(Form("%s%s",directory.Data(),"BinComparison.png"));
//   c->SaveAs(Form("%s%s",directory.Data(),"BinComparison.root"));
//   c->Clear();
//   delete c;
//   delete Ntriggerscent;
//   delete Ntriggersvz;
//   delete Nassociatedcent;
//   delete Nassociatedvz;
//   delete pTTriggerCheck;
//   delete pTAssociatedCheck;
//   delete phiTriggerCheck;
//   delete phiAssociatedCheck;
//   delete etaTriggerCheck;
//   delete etaAssociatedCheck;
//   delete eventsMb;
//   delete eventsVz;
//  
// }
// ///Merge everything before correcting (as if there were no bins.
// void AliCorrelation3p::Drawmergefirst(Option_t *options)
// {
//   /// overloaded from TObject: draw histograms
//   TString* opbackup = new TString(options);
//   DrawHistograms(fHistograms, Form("%s%s",opbackup->Data()," mergedfirst"),0,fMBinEdges.GetSize()-1,0,fZBinEdges.GetSize()-1);
//   if (fMixedEvent) {
//     AliCorrelation3p correctedmergedfirst;
//     // note: the assignment operator duplicates also the mixed event object
//     // wheras Copy avoids this if the target object does not host a mixed event
//     this->Copy(correctedmergedfirst);
//     TString name(GetName());
//     name+="_correctedmergedfirst";
//     correctedmergedfirst.SetName(name);
//     correctedmergedfirst.ProcessEventMixingmergefirst(fMixedEvent);
//     correctedmergedfirst.Draw(Form("%s%s",opbackup->Data(),"corrected/ mergedfirst"));
//   }
// }
// int AliCorrelation3p::ProcessEventMixingmergefirst(const AliCorrelation3p* pCorrME)
// {
//   /// add histograms from another instance
//   if (!pCorrME || !fHistograms || !pCorrME->fHistograms) return -1;
// 
//   if (pCorrME->fMixedEvent) {
//     ProcessEventMixingmergefirst(pCorrME->fMixedEvent);
//     // AliCorrelation3p correctedME;
//     // pCorrME->Copy(correctedME);
//     // correctedME.ProcessEventMixing(pCorrME->fMixedEvent);
//     // ProcessEventMixing(&correctedME);
//     return 0;
//   }
//   
//  for (int i=khDeltaPhi; i<khPhiPhiDEta; i++) {
//     if (fHistograms->At(GetNumberHist(i,0,0))==NULL || pCorrME->fHistograms->At(GetNumberHist(i,0,0))==NULL) continue;
//     TH1* target=reinterpret_cast<TH1*>(fHistograms->At(GetNumberHist(i,0,0)));
//     TH1* source=reinterpret_cast<TH1*>(pCorrME->fHistograms->At(GetNumberHist(i,0,0)));
//     for(Int_t mbin=0;mbin<fMBinEdges.GetSize()-1;mbin++){
//       for(Int_t zbin=0;zbin<fZBinEdges.GetSize()-1;zbin++){
// 	if (fHistograms->At(GetNumberHist(i,mbin,zbin))==NULL || pCorrME->fHistograms->At(GetNumberHist(i,mbin,zbin))==NULL) continue;
// 	TH1* targettemp=reinterpret_cast<TH1*>(fHistograms->At(GetNumberHist(i,mbin,zbin)));
// 	TH1* sourcetemp=reinterpret_cast<TH1*>(pCorrME->fHistograms->At(GetNumberHist(i,mbin,zbin)));
// 	if (!target || !source) continue;
// 	TString name(fHistograms->At(GetNumberHist(i,mbin,zbin))->GetName());
// 	if (name.CompareTo(targettemp->GetName())!=0) {
// 	  AliWarning(Form("skipping incompatible objects at position %d: %s vs %s", i, source->GetName(), target->GetName()));
// 	  continue;
// 	}
// 	if (source->IsA()!=target->IsA()) {
// 	  AliWarning(Form("skipping incompatible classes at position %d: %s vs %s", i, source->ClassName(), target->ClassName()));
// 	  continue;
// 	}
// // 	if(target->GetEntries()<1){
// // 	  for(int k=0;k<target->GetNbinsX();k++){
// // 	    for(int m=0;m<target->GetNbinsY();m++){
// // 	      targettemp->SetBinContent(k,m,0.0);
// // 	    }
// // 	  }
// // 	  continue;//if there is no statistics, dont do it.
// // 	}
// 	target->Add(targettemp);
// 	source->Add(sourcetemp);
//       }
//     }
//     Double_t scalingfactor=1.0;//scaling factor for eta normalization.
//     if(i>khPhiEta){
//       Double_t inttarget = dynamic_cast<TH2D*>(target)->Integral(1, dynamic_cast<TH2D*>(target)->GetNbinsX(),1,1);
//       Double_t intsource = dynamic_cast<TH2D*>(source)->Integral(1,source->GetNbinsX(),1,1);
//       if(!((inttarget==0)||(intsource==0)))scalingfactor=inttarget/intsource;
//       else scalingfactor =0;
//       scalingfactor=dynamic_cast<TH1D*>(fHistograms->At(kHistTriggerpT))->Integral();
//       source->Scale(scalingfactor);
//     }//Integral along the deltaeta_12/deltaeta axis.
//     target->Divide(source);
//   }
//   return 0;
// }
// void AliCorrelation3p::DrawDeltaEtamergedfirst(const TObjArray* histograms, int index, int canvasno, const char* dir) const
// {
//   /// draw canvas with delta eta histogram in different modes
//   if (histograms) {
//     TString name;
//     int padno=1;
//     name.Form("%s_%d", GetName(), canvasno);
//     TCanvas* c=new TCanvas(name);
//     name.Insert(0,dir);
//     c->SetWindowSize(1600,1114);
//     c->SetTitle(Form("%s: #Delta#Phi vs #Delta#Eta12 of associated particles", GetName()));
//     c->Divide(2,2);
//     padno=1;
//     if(index!= khPhiEta1){
//     c->cd(padno++); DrawHistogram(dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(index,0,0))), "surf2", "#Delta#eta_{12}", "#Delta#Phi", "yield a.u.");
//     c->cd(padno++); DrawHistogram(dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(index,0,0))), "surf3", "#Delta#eta_{12}", "#Delta#Phi", "yield a.u.");
//     }
//     if(index== khPhiEta1){
//     c->cd(padno++); DrawHistogram(dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(index,0,0))), "surf2", "#Delta#eta", "#Delta#Phi", "yield a.u.");
//     c->cd(padno++); DrawHistogram(dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(index,0,0))), "surf3", "#Delta#eta", "#Delta#Phi", "yield a.u.");
//     }
//     TH1* projY=(dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(index,0,0))))->ProjectionY();
//     projY->GetYaxis()->SetRangeUser(0., 1.1*projY->GetBinContent(projY->GetMaximumBin()));
//     c->cd(padno++); DrawHistogram("", projY, "", "#Delta#Phi", "yield a.u.");
//     c->Print(name+".eps");
//     c->Print(name+".png");
//     name+=".root";
//     c->SaveAs(name);
//   }
// }
// //Draw "normal"
// void AliCorrelation3p::Draw(Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin,const char* codir,Bool_t bDivide,Option_t *options)
// {
//   /// overloaded from TObject: draw histograms
//   TString* opbackup = new TString(options);
//   if(minMbin<0)minMbin=0;
//   if(maxMbin<0||(maxMbin>(fMBinEdges.GetSize()-2)))maxMbin=fMBinEdges.GetSize()-2;//negative means max is max
//   if(minMbin>maxMbin)minMbin=maxMbin;//only the bin given by maxMbin
//   if(minZbin<0)minZbin=0;
//   if(maxZbin<0||(maxZbin>(fZBinEdges.GetSize()-2)))maxZbin=fZBinEdges.GetSize()-2;//negative means max is max
//   if(minZbin>maxMbin)minZbin=maxZbin;//only the bin given by maxZbin
// 
//   DrawHistograms(fHistograms, Form("%s",opbackup->Data()),minMbin,maxMbin,minZbin,maxZbin);
// 
//   if (fMixedEvent) {
//     AliCorrelation3p corrected;
//     // note: the assignment operator duplicates also the mixed event object
//     // wheras Copy avoids this if the target object does not host a mixed event
//     this->Copy(corrected);
//     TString name(GetName());
//     name+="_corrected";
//     corrected.SetName(name);
//     corrected.ProcessEventMixing(fMixedEvent, bDivide);
//     TFile f1("corrected.root","RECREATE");
//     f1.cd();
//     corrected.Write();
//     f1.Close();
//     if(TString(codir).CompareTo("")==0)DrawHistograms(corrected.fHistograms,Form("%s%s",opbackup->Data(),"../corrected/ corrected"),minMbin,maxMbin,minZbin,maxZbin);
//     else DrawHistograms(corrected.fHistograms,Form("%s%s%s",opbackup->Data(),codir,"/ corrected"),minMbin,maxMbin,minZbin,maxZbin);
//   }
// }
// int AliCorrelation3p::ProcessEventMixing(const AliCorrelation3p* pCorrME, Bool_t bDivide)
// {
//   /// add histograms from another instance
//   
//   if (!pCorrME || !fHistograms || !pCorrME->fHistograms) return -1;
// 
//   if (pCorrME->fMixedEvent) {
//     ProcessEventMixing(pCorrME->fMixedEvent,bDivide);
//     return 0;
//   }
//   
//   
//   Double_t scalingfactorsig[fMBinEdges.GetSize()-1][fZBinEdges.GetSize()-1];//First collect the number of triggers filled with in each bin.
//   Double_t scalingfactormixed[fMBinEdges.GetSize()-1][fZBinEdges.GetSize()-1];
// //   Double_t scalingfactorintegral[fMBinEdges.GetSize()-1][fZBinEdges.GetSize()-1];
// //   Double_t temp1=0;
// //   Double_t temp2=0;
// //   Double_t pii = TMath::Pi();
//   for(int mb=0;mb<fMBinEdges.GetSize()-1;mb++){
//     for(int zb=0;zb<fZBinEdges.GetSize()-1;zb++){
//       scalingfactorsig[mb][zb]= dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistNTriggers,mb,zb)))->Integral();
//       scalingfactormixed[mb][zb] = dynamic_cast<TH1D*>(pCorrME->fHistograms->At(GetNumberHist(kHistNTriggers,mb,zb)))->Integral();
// //       temp1=0;
// //       temp2=0;
// //       TH3D* threedhist = dynamic_cast<TH3D*>(fHistograms->At(GetNumberHist(khPhiPhiDEta,mb,zb)));
// //       TH3D* threedhistmix = dynamic_cast<TH3D*>(pCorrME->fHistograms->At(GetNumberHist(khPhiPhiDEta,mb,zb)));
// //       temp1 += threedhist->Integral(1,threedhist->GetXaxis()->GetNbins()+1,threedhist->GetYaxis()->FindBin(pii),threedhist->GetYaxis()->GetNbins()+1,threedhist->GetZaxis()->FindBin(0.5*pii-0.2),threedhist->GetZaxis()->FindBin(0.5*pii+0.2));
// //       temp1 += threedhist->Integral(1,threedhist->GetXaxis()->GetNbins()+1,threedhist->GetYaxis()->FindBin(0.5*pii-0.2),threedhist->GetYaxis()->FindBin(0.5*pii+0.2),threedhist->GetZaxis()->FindBin(pii),threedhist->GetZaxis()->GetNbins()+1);
// //       temp2 += threedhistmix->Integral(1,threedhistmix->GetXaxis()->GetNbins()+1,threedhistmix->GetYaxis()->FindBin(pii),threedhistmix->GetYaxis()->GetNbins()+1,threedhistmix->GetZaxis()->FindBin(0.5*pii-0.2),threedhistmix->GetZaxis()->FindBin(0.5*pii+0.2));
// //       temp2 += threedhistmix->Integral(1,threedhistmix->GetXaxis()->GetNbins()+1,threedhistmix->GetYaxis()->FindBin(0.5*pii-0.2),threedhistmix->GetYaxis()->FindBin(0.5*pii+0.2),threedhistmix->GetZaxis()->FindBin(1.0*pii),threedhistmix->GetZaxis()->GetNbins()+1);      
// //       if(temp2!=0)scalingfactorintegral[mb][zb] = temp1/temp2;
// //       else scalingfactorintegral[mb][zb] = 0;
//     }
//   }
// 
//   for(Int_t mbin=0;mbin<fMBinEdges.GetSize()-1;mbin++){
//     for(Int_t zbin=0;zbin<fZBinEdges.GetSize()-1;zbin++){
//       if(fHistograms->At(GetNumberHist(khPhiPhiDEta,mbin,zbin))==NULL||pCorrME->fHistograms->At(GetNumberHist(khPhiPhiDEta,mbin,zbin))==NULL) continue;
//       TH3D* target = dynamic_cast<TH3D*>(fHistograms->At(GetNumberHist(khPhiPhiDEta,mbin,zbin)));
//       TH3D* source = dynamic_cast<TH3D*>(pCorrME->fHistograms->At(GetNumberHist(khPhiPhiDEta,mbin,zbin))->Clone(Form("%s%i%i","mixingclone3d",mbin,zbin)));
//       target->Sumw2();
//       source->Sumw2();
//       if (!target||!source) continue;
//       TString name(pCorrME->fHistograms->At(GetNumberHist(khPhiPhiDEta,mbin,zbin))->GetName());
//       if(name.CompareTo(target->GetName())!=0){
// 	AliWarning(Form("skipping incompatible objects at position %d: %s vs %s", khPhiPhiDEta, source->GetName(), target->GetName()));
// 	continue;
//       }
//       if(source->IsA()!=target->IsA()){
// 	AliWarning(Form("skipping incompatible classes at position %d: %s vs %s", khPhiPhiDEta, source->ClassName(), target->ClassName()));
// 	continue;
//       }
//       if(target->GetEntries()<1||scalingfactormixed[mbin][zbin]==0||scalingfactorsig[mbin][zbin]==0){//either nothing in the signal, or zero scaling factor.
// 	for(int k=0; k<=target->GetNbinsX()+1;k++){
// 	  for(int m=0;m<=target->GetNbinsY()+1;m++){
// 	    for(int n=0; n<=target->GetNbinsZ()+1;n++){
// 	      target->SetBinContent(k,m,n,0.0);
// 	    }
// 	  }
// 	}
// 	continue;//if there is no statistics, do nothing.
//       }
//       if(!bDivide) {
// 	target->Add(source,-scalingfactorsig[mbin][zbin]/scalingfactormixed[mbin][zbin]);//Corrected Yield.
// 	target->Scale(1.0/scalingfactorsig[mbin][zbin]);
//       }
//       else{
// 	source->Scale(scalingfactorsig[mbin][zbin]/scalingfactormixed[mbin][zbin]);
//  	target->Divide(source);
//       }
//     }
//   }
//   
// //   TH3D* DPhiDPhiDEta12s = this->CollectHists3d(khPhiPhiDEta,0,fMBinEdges.GetSize()-2,0,fZBinEdges.GetSize()-2);
// //   TCanvas* c=NULL;
// //   c=new TCanvas("correctedhistogram");
// //   c->SetWindowSize(1600,1114);
// //   c->Divide(2,2);
// //   int padno = 1;
// //   TH2D* Temphist;
// //   TH2D* Temphist2;
// //   TH2D* Temphist3;
// //   TH2D* Temphist4;
// //   c->cd(padno++);
// //   Temphist = this->slice(DPhiDPhiDEta12s,"yz",1,DPhiDPhiDEta12s->GetNbinsX());
// //   DrawHistogram("#Delta#Phi for Associated 1 vs #Delta#Phi for Associated 2",Temphist, "surf3", "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
// //   c->cd(padno++);
// //   Temphist2 = this->slice(DPhiDPhiDEta12s,"yz",DPhiDPhiDEta12s->GetXaxis()->FindBin(-0.4),DPhiDPhiDEta12s->GetXaxis()->FindBin(0.4));
// //   DrawHistogram("#Delta#Phi for A1 vs #Delta#Phi for A2 (#Delta#eta_{12} <= 0.4)",Temphist2, "surf3", "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
// //   c->cd(padno++);
// //   Temphist3 = this->slice(DPhiDPhiDEta12s,"yz",DPhiDPhiDEta12s->GetXaxis()->FindBin(-1),DPhiDPhiDEta12s->GetXaxis()->FindBin(-0.4));
// //   Temphist3->Add(this->slice(DPhiDPhiDEta12s,"yz",DPhiDPhiDEta12s->GetXaxis()->FindBin(0.4),DPhiDPhiDEta12s->GetXaxis()->FindBin(1.0)));
// //   DrawHistogram("#Delta#Phi for A1 vs #Delta#Phi for A2 (0.4 < #Delta#eta_{12} <= 1.0)",Temphist3, "surf3", "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
// //   c->cd(padno++);
// //   Temphist4 = this->slice(DPhiDPhiDEta12s,"yz",1,DPhiDPhiDEta12s->GetXaxis()->FindBin(-1));
// //   Temphist4->Add(this->slice(DPhiDPhiDEta12s,"yz",DPhiDPhiDEta12s->GetXaxis()->FindBin(1),DPhiDPhiDEta12s->GetNbinsX()));
// //   DrawHistogram("#Delta#Phi for A1 vs #Delta#Phi (#Delta#eta_{12} > 1.0)",Temphist4, "surf3", "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
// //   c->Print("correctedhistogram.eps");
// //   c->Print("correctedhistogram.png");
// //   c->SaveAs("correctedhistogram.root");
// //     
//   for (int i=khDeltaPhi; i<khPhiPhiDEta; i++) {
//     for(Int_t mbin=0;mbin<fMBinEdges.GetSize()-1;mbin++){
//       for(Int_t zbin=0;zbin<fZBinEdges.GetSize()-1;zbin++){
// 	if (fHistograms->At(GetNumberHist(i,mbin,zbin))==NULL || pCorrME->fHistograms->At(GetNumberHist(i,mbin,zbin))==NULL) continue;
// 	TH2D* target=dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(i,mbin,zbin)));
// 	TH2D* source=dynamic_cast<TH2D*>(pCorrME->fHistograms->At(GetNumberHist(i,mbin,zbin))->Clone(Form("%s%i%i%i","mixingclone2d",mbin,zbin,i)));
// 	target->Sumw2();
// 	source->Sumw2();
// 	if (!target || !source) continue;
// 	TString name(pCorrME->fHistograms->At(GetNumberHist(i,mbin,zbin))->GetName());
// 	if (name.CompareTo(target->GetName())!=0) {
// 	  AliWarning(Form("skipping incompatible objects at position %d: %s vs %s", i, source->GetName(), target->GetName()));
// 	  continue;
// 	}
// 	if (source->IsA()!=target->IsA()) {
// 	  AliWarning(Form("skipping incompatible classes at position %d: %s vs %s", i, source->ClassName(), target->ClassName()));
// 	  continue;
// 	}
// 
// 	if(target->GetEntries()<1||scalingfactormixed[mbin][zbin]==0||scalingfactorsig[mbin][zbin]==0){
// 	  for(int k=0;k<=target->GetNbinsX()+1;k++){
// 	    for(int m=0;m<=target->GetNbinsY()+1;m++){
// 	      target->SetBinContent(k,m,0.0);
// 	    }
// 	  }
// 	  continue;//if there is no statistics, dont do it.
// 	}
// 	if(!bDivide){
// 	  target->Add(source,-scalingfactorsig[mbin][zbin]/scalingfactormixed[mbin][zbin]);
// 	  target->Scale(1.0/scalingfactorsig[mbin][zbin]);
// 	}
// 	else {
// 	  source->Scale(scalingfactorsig[mbin][zbin]/scalingfactormixed[mbin][zbin]);
//  	  target->Divide(source);
// 	}
//       }
//     }
//   }
//   return 0;
// }
// void AliCorrelation3p::DrawDeltaEta(const TObjArray* histograms, int index, int canvasno, const char* dir, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin,Bool_t eachbin,Bool_t average) const
// {
//   /// draw canvas with delta eta histogram in different modes
//   if (histograms) {
//     TString name;
//     int padno=1;
//     name.Form("%s_%d", GetName(), canvasno);
//     TCanvas* c=new TCanvas(name);
//     name.Insert(0,dir);
//     c->SetWindowSize(1600,1114);
//     c->SetTitle(Form("%s: #Delta#Phi vs #Delta#Eta12 of associated particles", GetName()));
//     c->Divide(2,2);
//     padno=1;
//     TString zaxisunit;
//     TString xaxisunit=TString("#Delta#eta_{12}");
//     if(index==khPhiEta1){zaxisunit=TString("#associated per trigger");xaxisunit=TString("#Delta#eta");}
// //     if(index==khPhiEtaScaled)zaxisunit=TString("#associated per trigger");
//     if(index!=khPhiEta1)zaxisunit=TString("#pairs per trigger");
//     
//     TH2D* Correlation;
//     if(eachbin==kFALSE){
//       if(!average) Correlation =      CollectHists2d(index,minMbin,maxMbin,minZbin,maxZbin,"deltaetadeltaphihistogram2d");
//       if( average) Correlation = (TH2D*)AverageHists(index,minMbin,maxMbin,minZbin,maxZbin,"deltaetadeltaphihistogram2d");
//       c->cd(padno++); DrawHistogram(Correlation, "surf2", xaxisunit.Data(), "#Delta#Phi",zaxisunit.Data());
//       c->cd(padno++); DrawHistogram(Correlation, "surf3", xaxisunit.Data(), "#Delta#Phi",zaxisunit.Data());
//       TH1D* projY  = Correlation->ProjectionY();
//       projY->GetYaxis()->SetRangeUser(0., 1.1*projY->GetBinContent(projY->GetMaximumBin()));
//       c->cd(padno++); DrawHistogram(Form("%s%s","Integral over ",xaxisunit.Data()), projY, "E", "#Delta#Phi", zaxisunit.Data());
//       
//       Int_t binpih = Correlation->GetYaxis()->FindBin(TMath::Pi()+0.2);
//       Int_t binpil = Correlation->GetYaxis()->FindBin(TMath::Pi()-0.2);
//       TH1* projX   = Correlation->ProjectionX(Form("%s%s",Correlation->GetName(),"_px"),binpil,binpih);
//       projX->GetYaxis()->SetRangeUser(0., 1.1*projX->GetBinContent(projX->GetMaximumBin()));
//       c->cd(padno++); DrawHistogram("Integral over a slice of the away side peak around #pi#pm 0.2", projX, "E", xaxisunit.Data(), zaxisunit.Data());
//       c->Print(name+".eps");
//       c->Print(name+".png");
//       c->SaveAs(Form("%s%s",name.Data(),".root"));
//        
//       name.Form("%s_%d", GetName(), canvasno);
//       name.Append("fromTH3D");
//       c=new TCanvas(name);
//       name.Insert(0,dir);
//       c->SetWindowSize(1600,1114);
//       c->SetTitle(Form("%s: #Delta#Phi vs #Delta#Eta12 of associated particles", GetName()));
//       c->Divide(2,2);
//       padno=1;
//       TH3D* DPhiDPhiDEta12;
//       if(!average) DPhiDPhiDEta12 =      CollectHists3d(khPhiPhiDEta,minMbin,maxMbin,minZbin,maxZbin,"deltaetadeltaphihistogram3d");
//       if( average) DPhiDPhiDEta12 = (TH3D*)AverageHists(khPhiPhiDEta,minMbin,maxMbin,minZbin,maxZbin,"deltaetadeltaphihistogram3d");
//       TH2D* Localhist;
//       if(index==khPhiEta){
// 	Localhist = slice(DPhiDPhiDEta12,"yx",1,DPhiDPhiDEta12->GetNbinsZ(),"NoPhiCut",average);	
//       }
//       else if(index==khPhiEta_phicut1){
// 	    Localhist = DeltaEtaCut(DPhiDPhiDEta12,"lesspi2","PhiCut1",average);
//       }
//       else if(index==khPhiEta_phicut2){
// 	    Localhist = DeltaEtaCut(DPhiDPhiDEta12,"lesspi4","PhiCut2",average);
//       }
//       else if(index==khPhiEta_sameside){
// 	    Localhist = DeltaEtaCut(DPhiDPhiDEta12,"sameside","OnTheSameSide",average);
//       }
//       if(index != khPhiEta1){
// 	c->cd(padno++); DrawHistogram(Localhist, "surf2", xaxisunit.Data(), "#Delta#Phi",zaxisunit.Data());
// 	c->cd(padno++); DrawHistogram(Localhist, "surf3", xaxisunit.Data(), "#Delta#Phi",zaxisunit.Data());
// 	TH1* projYf3d=Localhist->ProjectionY();
// 	projYf3d->GetYaxis()->SetRangeUser(0., 1.1*projYf3d->GetBinContent(projYf3d->GetMaximumBin()));
// 	c->cd(padno++); DrawHistogram(Form("%s%s","Integral over ",xaxisunit.Data()), projYf3d, "E", "#Delta#Phi", zaxisunit.Data());
// 	Int_t binpihf3d = Localhist->GetYaxis()->FindBin(TMath::Pi()+0.1);
// 	Int_t binpilf3d = Localhist->GetYaxis()->FindBin(TMath::Pi()-0.1);
// 
// 	TH1* projXf3d=Localhist->ProjectionX("_px",binpilf3d,binpihf3d);
// 	projXf3d->GetYaxis()->SetRangeUser(0., 1.1*projXf3d->GetBinContent(projXf3d->GetMaximumBin()));
// 	c->cd(padno++); DrawHistogram("Integral over a slice of the away side peak around #pi#pm 0.1", projXf3d, "E", xaxisunit.Data(), zaxisunit.Data());
//         c->Print(name+".eps");
//         c->Print(name+".png");
//         c->SaveAs(Form("%s%s",name.Data(),".root"));
//     }
//     return;
//   }
//     if(eachbin==kTRUE){
//       TString tempname;
//       for(int mbin = minMbin;mbin<maxMbin;mbin++){
// 	for(int zbin= minZbin;zbin<maxZbin;zbin++){
// 	  tempname.Form("%s%i%i",name.Data(),mbin,zbin);
// 	  padno=1;
// 	  if(index != khPhiEta1){
// 	    c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)), "surf2", "#Delta#eta_{12}", "#Delta#Phi", zaxisunit.Data());
// 	    c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)), "surf3", "#Delta#eta_{12}", "#Delta#Phi", zaxisunit.Data());
// 	  }
// 	  if(index == khPhiEta1){
// 	    c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)), "surf2", "#Delta#eta", "#Delta#Phi", zaxisunit.Data());
// 	    c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)), "surf3", "#Delta#eta", "#Delta#Phi", zaxisunit.Data());
// 	  }
// 	  TH1* projY=((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)))->ProjectionY();
// 	  projY->GetYaxis()->SetRangeUser(0., 1.1*projY->GetBinContent(projY->GetMaximumBin()));
// 	  c->cd(padno++); DrawHistogram("", projY, "E", "#Delta#Phi", zaxisunit.Data());
// 	  Int_t binpih = ((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)))->GetYaxis()->FindBin(TMath::Pi()+0.2);
// 	  Int_t binpil = ((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)))->GetYaxis()->FindBin(TMath::Pi()-0.2);
// 	  TH1* projX   = ((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)))->ProjectionX(Form("%s%s",((TH2D*)fHistograms->At(GetNumberHist(index,mbin,zbin)))->GetName(),"_px"),binpil,binpih);
// 	  projX->GetYaxis()->SetRangeUser(0., 1.1*projX->GetBinContent(projX->GetMaximumBin()));
// 	  c->cd(padno++); DrawHistogram("Integral over a slice of the away side peak around #pi#pm 0.01", projX, "E", xaxisunit.Data(), zaxisunit.Data());
// 	  c->Print(tempname+".png");
// 	  c->Print(tempname+".eps");
// 	  tempname+=".root";
// 	  c->SaveAs(tempname);
// 	  
// // 	  name.Form("%s_%d", GetName(), canvasno);
// // 	  name.Append("fromTH3D");
// 	  tempname.Form("%s%i%i%s",name.Data(),mbin,zbin,"fromTH3D");
// 	  c=new TCanvas(tempname);
// // 	  tempname.Insert(0,dir);
// 	  c->SetWindowSize(1600,1114);
// 	  c->SetTitle(Form("%s: #Delta#Phi vs #Delta#Eta12 of associated particles", GetName()));
// 	  c->Divide(2,2);
// 	  padno=1;
// 	  TH3D* DPhiDPhiDEta12;
// 	  DPhiDPhiDEta12 = (TH3D*)fHistograms->At(GetNumberHist(khPhiPhiDEta,mbin,zbin));
// 	  TH2D* Localhist;
// 	  if(index==khPhiEta){
// 	    Localhist = slice(DPhiDPhiDEta12,"yx",1,DPhiDPhiDEta12->GetNbinsZ(),"NoPhiCut",average);	
// 	  }
// 	  else if(index==khPhiEta_phicut1){
// 	    Localhist = DeltaEtaCut(DPhiDPhiDEta12,"lesspi2","PhiCut1",average);
// 	  }
// 	  else if(index==khPhiEta_phicut2){
// 	    Localhist = DeltaEtaCut(DPhiDPhiDEta12,"lesspi4","PhiCut2",average);
// 	  }
// 	  else if(index==khPhiEta_sameside){
// 	    Localhist = DeltaEtaCut(DPhiDPhiDEta12,"sameside","OnTheSameSide",average);
// 	  }
// 	  if(index != khPhiEta1){
// 	    c->cd(padno++); DrawHistogram(Localhist, "surf2", xaxisunit.Data(), "#Delta#Phi",zaxisunit.Data());
// 	    c->cd(padno++); DrawHistogram(Localhist, "surf3", xaxisunit.Data(), "#Delta#Phi",zaxisunit.Data());
// 	    TH1* projYf3d=Localhist->ProjectionY();
// 	    projYf3d->GetYaxis()->SetRangeUser(0., 1.1*projYf3d->GetBinContent(projYf3d->GetMaximumBin()));
// 	    c->cd(padno++); DrawHistogram(Form("%s%s","Integral over ",xaxisunit.Data()), projYf3d, "E", "#Delta#Phi", zaxisunit.Data());
// 	    Int_t binpihf3d = Localhist->GetYaxis()->FindBin(TMath::Pi()+0.1);
// 	    Int_t binpilf3d = Localhist->GetYaxis()->FindBin(TMath::Pi()-0.1);
// 
// 	    TH1* projXf3d=Localhist->ProjectionX("_px",binpilf3d,binpihf3d);
// 	    projXf3d->GetYaxis()->SetRangeUser(0., 1.1*projXf3d->GetBinContent(projXf3d->GetMaximumBin()));
// 	    c->cd(padno++); DrawHistogram("Integral over a slice of the away side peak around #pi#pm 0.1", projXf3d, "E", xaxisunit.Data(), zaxisunit.Data());
// 	    c->Print(tempname+".eps");
// 	    c->Print(tempname+".png");
// 	    c->SaveAs(Form("%s%s",tempname.Data(),".root"));
// 	  }
// 	}
//       }
//     }
//   }
// }
// 
// void AliCorrelation3p::DrawHistograms(const TObjArray* histograms, const char* options , Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin) const
// {
//   /// draw a set of histograms
//   bool bDrawAll=true;
//   bool bDrawControl=false;
//   bool bDrawPhiPhi=false;
//   bool bDrawEtaPhi=false;
//   bool allreadymerged =false;
//   bool eachbin = false;
//   bool average = false;
//   TString delimiter(" ");
//   TString directory("");
//   TStringToken token(options, delimiter);
//   while (token.NextToken()) {
//     const char* key=0;
//     TString argument=token;
// 
//     if (argument.CompareTo("-h")==0 ||
// 	argument.CompareTo("--help")==0) {
//       AliInfo(Form("Draw options:"
// 		   "\n\t  control   - draw control histograms"
// 		   "\n\t  phiphi    - draw dPhi-dPhi distributions"
// 		   "\n\t  etaphi    - draw dEta-dPhi distributions"
// 		   "\n\t  etaphi    - draw dEta-dPhi distributions"
// 		   "\n\t by default, all histograms are plotted"
// 		   "\n\t  name=<regexp> plot objects matching regular expression"
// 		   "\n\t  skip=<regexp> skip objects matching regular expression"
// 		   "\n\t  dir=<directory> saves all plots in the given directory"
// 		   "\n\t  mergedfirst draws only the first hists."
// 		   "\n\t  corrected  - average over mult-vz bins, standard is sum."
// 		   ));
//       return;
//     }
//     key="control";
//     if (argument.CompareTo(key)==0) {
//       bDrawAll=false;
//       bDrawControl=true;
//       continue;
//     }
//     key="phiphi";
//     if (argument.CompareTo(key)==0) {
//       bDrawAll=false;
//       bDrawPhiPhi=true;
//       continue;
//     }
//     key="etaphi";
//     if (argument.CompareTo(key)==0) {
//       bDrawAll=false;
//       bDrawEtaPhi=true;
//       continue;
//     }
//     key="name=";
//     if (argument.BeginsWith(key)) {
//       argument.ReplaceAll(key, "");
//       TRegexp re(argument, true);
//       TString oname(GetName());
//       if (!oname.Contains(re))
// 	return;
//     }
//     key="skip=";
//     if (argument.BeginsWith(key)) {
//       argument.ReplaceAll(key, "");
//       TRegexp re(argument, true);
//       TString oname(GetName());
//       if (oname.Contains(re))
// 	return;
//     }
//     key="dir=";
//     if (argument.BeginsWith(key)) {
//       argument.ReplaceAll(key, "");
//       directory.Append(argument);
//     }    
//     key="mergedfirst";
//     if (argument.BeginsWith(key)) {
//       argument.ReplaceAll(key, "");
//       allreadymerged = true;
//     }    
//     key="draweachbin";
//     if(argument.BeginsWith(key)){
//       argument.ReplaceAll(key,"");
//       eachbin = true;
//       bDrawAll=false;
//       bDrawEtaPhi=true;
//     }
//     key = "corrected";
//     if(argument.CompareTo(key)==0){
//       average = true;
//     }
//   }
//   if (histograms) {
//     TCanvas* c=NULL;
//     TString name;
//     int canvasno=1;
//     int padno=1;
//     const char* drawoption="surf3";
//     const char* d2ddrawoption="E";
//     TH3D* DPhiDPhiDEta12;
//     
//     if(!average) DPhiDPhiDEta12 =      CollectHists3d(khPhiPhiDEta,minMbin,maxMbin,minZbin,maxZbin);
//     if(average)  DPhiDPhiDEta12 = (TH3D*)AverageHists(khPhiPhiDEta,minMbin,maxMbin,minZbin,maxZbin,"PhiPhiEta12");//If it is corrected, it is an average.
// 
//     c=new TCanvas("summed3dhist");
//     DrawHistogram("#Delta #Phi_1 vs #Delta #Phi_2 vs #Delta #eta_{12}",DPhiDPhiDEta12,"","#Delta #eta_{12}","#Delta#Phi_1","#Delta#Phi_{2}");
//     c->SaveAs(Form("%s%s",directory.Data(),"/3dhist.root"));
//     
//     
//     name.Form("%s_%d", GetName(), canvasno++);
//     if (bDrawAll || bDrawControl) {
//     c=new TCanvas(name);
//     name.Insert(0,directory);
//     c->SetWindowSize(1600,1114);
//     //c->SetTitle(Form("%s: particle properties", GetName()));//meaningless, gets not shown
//     c->Divide(3,3);
//     padno=1;
//     if(!allreadymerged){
//       c->cd(padno++);gPad->SetLogy(1);DrawHistogram(CollectHists1d(kHistpT,minMbin,maxMbin,minZbin,maxZbin)           , d2ddrawoption, "p_{T} [GeV/c]");
//       c->cd(padno++);gPad->SetLogy(0);DrawHistogram(CollectHists1d(kHistPhi,minMbin,maxMbin,minZbin,maxZbin)          , d2ddrawoption, "#Phi [rad]");
//       c->cd(padno++);gPad->SetLogy(0);DrawHistogram(CollectHists1d(kHistEta,minMbin,maxMbin,minZbin,maxZbin)          , d2ddrawoption, "#eta []");
//       c->cd(padno++);gPad->SetLogy(1);DrawHistogram(CollectHists1d(kHistTriggerpT,minMbin,maxMbin,minZbin,maxZbin)    , d2ddrawoption, "p_{T} [GeV/c]");
//       c->cd(padno++);gPad->SetLogy(0);DrawHistogram(CollectHists1d(kHistTriggerPhi,minMbin,maxMbin,minZbin,maxZbin)   , d2ddrawoption, "#Phi [rad]");
//       c->cd(padno++);gPad->SetLogy(0);DrawHistogram(CollectHists1d(kHistTriggerEta,minMbin,maxMbin,minZbin,maxZbin)   , d2ddrawoption, "#eta []");
//       c->cd(padno++);gPad->SetLogy(1);DrawHistogram(CollectHists1d(kHistAssociatedpT,minMbin,maxMbin,minZbin,maxZbin) , d2ddrawoption, "p_{T} [GeV/c]");
//       c->cd(padno++);gPad->SetLogy(0);DrawHistogram(CollectHists1d(kHistAssociatedPhi,minMbin,maxMbin,minZbin,maxZbin), d2ddrawoption, "#Phi [rad]");
//       c->cd(padno++);gPad->SetLogy(0);DrawHistogram(CollectHists1d(kHistAssociatedEta,minMbin,maxMbin,minZbin,maxZbin), d2ddrawoption, "#eta []");
//     }
//     if(allreadymerged){
//       c->cd(padno++);gPad->SetLogy(1); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistpT,0,0))           , d2ddrawoption, "p_{T}  [GeV/c]");
//       c->cd(padno++);gPad->SetLogy(0); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistPhi,0,0))          , d2ddrawoption, "#Phi (rad)");
//       c->cd(padno++);gPad->SetLogy(0); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistEta,0,0))          , d2ddrawoption, "#eta");
//       c->cd(padno++);gPad->SetLogy(1); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistTriggerpT,0,0))    , d2ddrawoption, "p_{T}  [GeV/c]");
//       c->cd(padno++);gPad->SetLogy(0); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistTriggerPhi,0,0))   , d2ddrawoption, "#Phi (rad)");
//       c->cd(padno++);gPad->SetLogy(0); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistTriggerEta,0,0))   , d2ddrawoption, "#eta");
//       c->cd(padno++);gPad->SetLogy(1); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistAssociatedpT,0,0)) , d2ddrawoption, "p_{T}  [GeV/c]");
//       c->cd(padno++);gPad->SetLogy(0); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistAssociatedPhi,0,0)), d2ddrawoption, "#Phi (rad)");
//       c->cd(padno++);gPad->SetLogy(0); DrawHistogram((TH1D*)fHistograms->At(GetNumberHist(kHistAssociatedEta,0,0)), d2ddrawoption, "#eta");
//     }
//     c->Print(name+".eps");
//     c->Print(name+".png");
//     name+=".root";
//     c->SaveAs(name);
//     }
// 
//     name.Form("%s_%d", GetName(), canvasno++);
//     if (bDrawAll || bDrawPhiPhi) {
//     c=new TCanvas(name);
//     name.Insert(0,directory);
//     c->SetWindowSize(1600,1114);
//     c->SetTitle(Form("%s: #Delta#Phi for associated particles", GetName()));
//     c->Divide(2,2);
//     padno=1;
//     if(!allreadymerged){
//       c->cd(padno++); 
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi     ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi1"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(  AverageHists(khDeltaPhi     ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi1"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); 
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi_near,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi2"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(  AverageHists(khDeltaPhi_near,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi2"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); 
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi_mid ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi3"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(  AverageHists(khDeltaPhi_mid ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi3"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); 
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi_far ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi4"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(  AverageHists(khDeltaPhi_far ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhi4"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       
//     }
//     if(allreadymerged){
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi_near,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi_mid,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi_far,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//     }
//     c->Print(name+".eps");
//     c->Print(name+".png");
//     c->SaveAs(Form("%s%s",name.Data(),".root"));
//     //Same from the TH3D
//     name.Append("fromTH3D");
//     c=new TCanvas(name);
//     c->SetWindowSize(1600,1114);
//     c->Divide(2,2);
//     padno=1;
//     TH2D* Temphist;
//     TH2D* Temphist2;
//     TH2D* Temphist3;
//     TH2D* Temphist4;
//     c->cd(padno++);
//     Temphist = slice(DPhiDPhiDEta12,"yz",1,DPhiDPhiDEta12->GetNbinsX(),"phiphislice1",average);
//     DrawHistogram("#Delta#Phi for Associated 1 vs #Delta#Phi for Associated 2",Temphist, drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//     c->cd(padno++);
//     Temphist2 = slice(DPhiDPhiDEta12,"yz",DPhiDPhiDEta12->GetXaxis()->FindBin(-0.4),DPhiDPhiDEta12->GetXaxis()->FindBin(0.4),"phiphislice2",average);
//     DrawHistogram("#Delta#Phi for A1 vs #Delta#Phi for A2 (#Delta#eta_{12} <= 0.4)",Temphist2, drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//     c->cd(padno++);
//     Temphist3 = slice(DPhiDPhiDEta12,"yz",DPhiDPhiDEta12->GetXaxis()->FindBin(-1),DPhiDPhiDEta12->GetXaxis()->FindBin(-0.4),"phiphislice3_1",average);
//     Temphist3->Add(slice(DPhiDPhiDEta12,"yz",DPhiDPhiDEta12->GetXaxis()->FindBin(0.4),DPhiDPhiDEta12->GetXaxis()->FindBin(1.0),"phiphislice3_2",average));
//     DrawHistogram("#Delta#Phi for A1 vs #Delta#Phi for A2 (0.4 < #Delta#eta_{12} <= 1.0)",Temphist3, drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//     c->cd(padno++);
//     Temphist4 = slice(DPhiDPhiDEta12,"yz",1,DPhiDPhiDEta12->GetXaxis()->FindBin(-1),"phiphislice4_1",average);
//     Temphist4->Add(slice(DPhiDPhiDEta12,"yz",DPhiDPhiDEta12->GetXaxis()->FindBin(1),DPhiDPhiDEta12->GetNbinsX(),"phiphislice4_2",average));
//     DrawHistogram("#Delta#Phi for A1 vs #Delta#Phi (#Delta#eta_{12} > 1.0)",Temphist4, drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//     c->Print(name+".eps");
//     c->Print(name+".png");
//     name+=".root";
//     c->SaveAs(name);
//     }
// 
//     name.Form("%s_%d", GetName(), canvasno++);
//     if (bDrawAll || bDrawPhiPhi) {
//     c=new TCanvas(name);
//     name.Insert(0,directory);
//     c->SetWindowSize(1600,1114);
//     c->SetTitle(Form("%s: #Delta#Phi for associated particles with charge selection", GetName()));
//     c->Divide(2,2);
//     padno=1;    
//     if(!allreadymerged){
//       c->cd(padno++); 
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi_ass_like             ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich1"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(  AverageHists(khDeltaPhi_ass_like             ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich1"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++);
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi_ass_unlike           ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich2"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(  AverageHists(khDeltaPhi_ass_unlike           ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich2"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); 
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi_ass_like__trig_unlike,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich3"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(  AverageHists(khDeltaPhi_ass_like__trig_unlike,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich3"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); 
//       if(!average)DrawHistogram(CollectHists2d(khDeltaPhi_all_like             ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich4"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       if( average)DrawHistogram(CollectHists2d(khDeltaPhi_all_like             ,minMbin,maxMbin,minZbin,maxZbin,"PhiPhich4"), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       
//     }
//     if(allreadymerged){
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi_ass_like,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi_ass_unlike,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi_ass_like__trig_unlike,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//       c->cd(padno++); DrawHistogram((TH2D*)fHistograms->At(GetNumberHist(khDeltaPhi_all_like,0,0)), drawoption, "#Delta#Phi_{1}", "#Delta#Phi_{2}", "yield a.u.");
//     }
//     c->Print(name+".eps");
//     c->Print(name+".png");
//     name+=".root";
//     c->SaveAs(name);
//     }
// 
//     if ((bDrawAll || bDrawEtaPhi)&&!allreadymerged){ DrawDeltaEta(histograms, khPhiEta        ,  canvasno, directory.Data(),minMbin,maxMbin,minZbin,maxZbin,eachbin,average); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&& allreadymerged){ DrawDeltaEtamergedfirst(histograms, khPhiEta        , canvasno, directory.Data()); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&&!allreadymerged){ DrawDeltaEta(histograms, khPhiEta1       ,  canvasno, directory.Data(),minMbin,maxMbin,minZbin,maxZbin,eachbin,average); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&& allreadymerged){ DrawDeltaEtamergedfirst(histograms, khPhiEta1       , canvasno, directory.Data()); canvasno++;}
// //     if ((bDrawAll || bDrawEtaPhi)&&!allreadymerged){ DrawDeltaEta(histograms, khPhiEtaScaled        , canvasno, directory.Data(),minMbin,maxMbin,minZbin,maxZbin,eachbin); canvasno++;}
// //     if ((bDrawAll || bDrawEtaPhi)&& allreadymerged){ DrawDeltaEtamergedfirst(histograms, khPhiEtaScaled        , canvasno, directory.Data()); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&&!allreadymerged){ DrawDeltaEta(histograms, khPhiEta_phicut1,  canvasno, directory.Data(),minMbin,maxMbin,minZbin,maxZbin,eachbin,average); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&& allreadymerged){ DrawDeltaEtamergedfirst(histograms, khPhiEta_phicut1, canvasno, directory.Data()); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&&!allreadymerged){ DrawDeltaEta(histograms, khPhiEta_phicut2,  canvasno, directory.Data(),minMbin,maxMbin,minZbin,maxZbin,eachbin,average); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&& allreadymerged){ DrawDeltaEtamergedfirst(histograms, khPhiEta_phicut2, canvasno, directory.Data()); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&&!allreadymerged){ DrawDeltaEta(histograms, khPhiEta_sameside, canvasno, directory.Data(),minMbin,maxMbin,minZbin,maxZbin,eachbin,average); canvasno++;}
//     if ((bDrawAll || bDrawEtaPhi)&& allreadymerged){ DrawDeltaEtamergedfirst(histograms, khPhiEta_sameside, canvasno, directory.Data()); canvasno++;}
//   }
// }
/*
void AliCorrelation3p::DrawHistogram(const char* title, TObject* pObj, const char* drawoption, const char* titleX, const char* titleY, const char* titleZ) const
{
  /// set properties and draw histogram
  if (!pObj) return;
  TNamed* pNamed=dynamic_cast<TNamed*>(pObj);
  if (!pNamed) return;
  if (title) pNamed->SetTitle(title);
  SetAxisTitles(pObj, titleX, titleY, titleZ);
  dynamic_cast<TH1*>(pObj)->SetTitleSize(0.045,"xyz");
  dynamic_cast<TH1*>(pObj)->SetStats(0);

  pObj->Draw(drawoption);
}*/
//Functions used to collect histograms and set axis.
// void AliCorrelation3p::SetAxisTitles(TObject* pObj, const char* titleX, const char* titleY, const char* titleZ) const
// {
//   // set axis titles of histograms
//   TH1* pHist=dynamic_cast<TH1*>(pObj);
//   if (!pHist) return;
//   if (pHist->GetXaxis() && titleX) pHist->GetXaxis()->SetTitle(titleX);
//   if (pHist->GetYaxis() && titleY) pHist->GetYaxis()->SetTitle(titleY);
//   if (pHist->GetZaxis() && titleZ) pHist->GetZaxis()->SetTitle(titleZ);
// }
// TH1D* AliCorrelation3p::CollectHists1d(Int_t Hist, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin, const char* name) const
// {
//   if(minMbin>(fMBinEdges.GetSize()-1)||maxMbin>(fMBinEdges.GetSize()-1)||minMbin<0||maxMbin<0) {cout << "Error parsing the Multiplicity bins between "<<minMbin<<" and "<<maxMbin<<endl;return NULL;}
//   if(minZbin>(fZBinEdges.GetSize()-1)||maxZbin>(fZBinEdges.GetSize()-1)||minZbin<0||maxZbin<0) {cout << "Error parsing the Vertex bins between "<<minZbin<<" and "<<maxZbin<<endl;return NULL;}
//   TH1D* temphist = (TH1D*)fHistograms->At(GetNumberHist(Hist,minMbin,minZbin));
//   TH1D* Collection = (TH1D*)temphist->Clone(name);
//   for(Int_t i=minMbin;i<maxMbin;i++){
//     for(Int_t j=minZbin;j<maxZbin;j++){
//       if(!i==0&&!j==0) Collection->Add((TH1D*)fHistograms->At(GetNumberHist(Hist,i,j)));
//     }
//   }
//   return Collection;
// }
// TH2D* AliCorrelation3p::CollectHists2d(Int_t Hist, Int_t minMbin, Int_t maxMbin, Int_t minZbin, Int_t maxZbin, const char* name) const
// {
//   if(minMbin>(fMBinEdges.GetSize()-1)||maxMbin>(fMBinEdges.GetSize()-1)||minMbin<0||maxMbin<0) {cout << "Error parsing the Multiplicity bins between "<<minMbin<<" and "<<maxMbin<<endl;return NULL;}
//   if(minZbin>(fZBinEdges.GetSize()-1)||maxZbin>(fZBinEdges.GetSize()-1)||minZbin<0||maxZbin<0) {cout << "Error parsing the Vertex bins between "<<minZbin<<" and "<<maxZbin<<endl;return NULL;}
//   TH2D* temphist = (TH2D*)fHistograms->At(GetNumberHist(Hist,minMbin,minZbin));
//   TH2D* Collection = (TH2D*)temphist->Clone(name);
//   for(Int_t i=minMbin;i<maxMbin;i++){
//     for(Int_t j=minZbin;j<maxZbin;j++){
//       if(!i==0&&!j==0) Collection->Add((TH2D*)fHistograms->At(GetNumberHist(Hist,i,j)));
//     }
//   }
//   return Collection;
// }
// 
// TH3D* AliCorrelation3p::CollectHists3d(Int_t Hist, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin, const char* name) const
// {
//   if(minMbin>(fMBinEdges.GetSize()-1)||maxMbin>(fMBinEdges.GetSize()-1)||minMbin<0||maxMbin<0) {cout << "Error parsing the Multiplicity bins between "<<minMbin<<" and "<<maxMbin<<endl;return NULL;}
//   if(minZbin>(fZBinEdges.GetSize()-1)||maxZbin>(fZBinEdges.GetSize()-1)||minZbin<0||maxZbin<0) {cout << "Error parsing the Vertex bins between "<<minZbin<<" and "<<maxZbin<<endl;return NULL;}
//   TH3D* temphist;
//   temphist = (TH3D*)fHistograms->At(GetNumberHist(Hist,minMbin,minZbin));
//   TH3D* Collection = (TH3D*)temphist->Clone(name);
//   for(Int_t i=minMbin;i<=maxMbin;i++){
//     for(Int_t j=minZbin;j<=maxZbin;j++){
//       if(!i==0&&!j==0){ 
// 	Collection->Add((TH3D*)fHistograms->At(GetNumberHist(Hist,i,j)));
//       }
//     }
//   }
//   return Collection;
// }
// TH1* AliCorrelation3p::AverageHists(Int_t Hist, Int_t minMbin, Int_t maxMbin,Int_t minZbin,Int_t maxZbin,const char* name) const
// {//Averages over multiplicity and vertex bins.
//   TH1* empty = NULL;
//   TH1D* temphist1d;
//   TH2D* temphist2d;
//   TH3D* temphist3d;
//   temphist1d = dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(Hist,minMbin,minZbin)));
//   temphist2d = dynamic_cast<TH2D*>(fHistograms->At(GetNumberHist(Hist,minMbin,minZbin)));
//   temphist3d = dynamic_cast<TH3D*>(fHistograms->At(GetNumberHist(Hist,minMbin,minZbin)));
//   Double_t BinContent = 0.0;
//   Double_t Normalization = 0.0;
//   Double_t locerr=0.0;  
//   if(temphist1d){//histogram is of 1d type
//     TH1D* Collection = (TH1D*)temphist1d->Clone(name);
//     for(int xbin=0; xbin<=Collection->GetNbinsX()+1;xbin++){
//       for(int ybin=0; ybin<=Collection->GetNbinsY()+1;ybin++){
// 	BinContent = 0;
// 	Normalization=0;
// 	locerr = 0;
// 	for(Int_t i=minMbin;i<=maxMbin;i++){
// 	  for(Int_t j=minZbin;j<=maxZbin;j++){
// 	    temphist1d = (TH1D*)fHistograms->At(GetNumberHist(Hist,i,j));
// 	    locerr = temphist1d->GetBinError(xbin);
// 	    if(locerr!=0)BinContent += temphist1d->GetBinContent(xbin)/(locerr*locerr);
// 	    if(locerr!=0)Normalization += 1/(locerr*locerr);
// 	  }
// 	}
// 	if(Normalization!=0)Collection->SetBinContent(xbin, BinContent/Normalization);
// 	else Collection->SetBinContent(xbin,0.0);
// 	if(Normalization!=0)Collection->SetBinError(xbin,1/TMath::Sqrt(Normalization));
// 	else Collection->SetBinError(xbin,0);
//       }
//     }
//     return Collection;
//   }
//   if(temphist2d){//histogram is of 2d type
//     TH2D* Collection = (TH2D*)temphist2d->Clone(name);
//     for(int xbin=0; xbin<=Collection->GetNbinsX()+1;xbin++){
//       for(int ybin=0; ybin<=Collection->GetNbinsY()+1;ybin++){
// 	BinContent = 0;
// 	Normalization=0;
// 	locerr = 0;
// 	for(Int_t i=minMbin;i<=maxMbin;i++){
// 	  for(Int_t j=minZbin;j<=maxZbin;j++){
// 	    temphist2d = (TH2D*)fHistograms->At(GetNumberHist(Hist,i,j));
// 	    locerr = temphist2d->GetBinError(xbin,ybin);
// 	    if(locerr!=0)BinContent += temphist2d->GetBinContent(xbin,ybin)/(locerr*locerr);
// 	    if(locerr!=0)Normalization += 1/(locerr*locerr);
// 	  }
// 	}
// 	if(Normalization!=0)Collection->SetBinContent(xbin, ybin, BinContent/Normalization);
// 	else Collection->SetBinContent(xbin,ybin,0);
// 	if(Normalization!=0)Collection->SetBinError(xbin,ybin,1.0/TMath::Sqrt(Normalization));
// 	else Collection->SetBinError(xbin,ybin,0.0);
//       }
//     }
//     return Collection;
//   }
//   if(temphist3d){
//     TH3D* Collection = (TH3D*)temphist3d->Clone(name);
//     for(int xbin=0;xbin<=Collection->GetNbinsX()+1;xbin++){
//       for(int ybin=0;ybin<=Collection->GetNbinsY()+1;ybin++){
// 	for(int zbin=0;zbin<=Collection->GetNbinsZ()+1;zbin++){
// 	  BinContent = 0;
// 	  Normalization=0;
// 	  for(Int_t i=minMbin;i<=maxMbin;i++){
// 	    for(Int_t j=minZbin;j<=maxZbin;j++){
// 	      temphist3d = (TH3D*)fHistograms->At(GetNumberHist(Hist,i,j));
// 	      locerr = temphist3d->GetBinError(xbin,ybin,zbin);
// 	      if(locerr!=0)
// 	      {
// 		BinContent += temphist3d->GetBinContent(xbin,ybin,zbin)/(locerr*locerr);
// 		Normalization += 1/(locerr*locerr);
// 	      }
// 	    }
// 	  }
// 	  if(Normalization!=0)Collection->SetBinContent(xbin, ybin,zbin, BinContent/Normalization);
// 	  else Collection->SetBinContent(xbin,ybin,zbin,0);
// 	  if(Normalization!=0)Collection->SetBinError(xbin,ybin,zbin,1.0/TMath::Sqrt(Normalization));
// 	  else Collection->SetBinError(xbin,ybin,zbin,0);
// 	}
//       }
//     }
//     return Collection;
//   }
//   return empty;
// }

// void AliCorrelation3p::CollectHistsFile(const char* name)
// {
// //   CreateWeights();
//   Double_t scalingfactorsig;//First collect the number of triggers filled with in each bin.
//   for(int mb=0;mb<fMBinEdges.GetSize()-1;mb++){
//     for(int zb=0;zb<fZBinEdges.GetSize()-1;zb++){
//       scalingfactorsig += dynamic_cast<TH1D*>(fHistograms->At(GetNumberHist(kHistNTriggers,mb,zb)))->Integral();
//     }
//   }
//   TFile outfile(name,"RECREATE");
//   outfile.cd();
//   for (int i=khDeltaPhi;i<khPhiPhiDEta;i++){
//     TString name1 = TString(fHistograms->At(GetNumberHist(i,0,0))->GetName());
//     TH2D* target = CollectHists2d(i,0,fMBinEdges.GetSize()-2,0,fZBinEdges.GetSize()-2);
//     target->SetName(name1.Data());
//     target->Write();
//     target->SetName(Form("%s%s",name1.Data(),"divbyntriggers"));
//     target->Scale(1.0/scalingfactorsig);
//     target->Write();
//   }
//   TString name2 = TString(fHistograms->At(GetNumberHist(khPhiPhiDEta,0,0))->GetName());
//   TH3D* target2 = CollectHists3d(khPhiPhiDEta,0,fMBinEdges.GetSize()-2,0,fZBinEdges.GetSize()-2);
//   target2->SetName(name2.Data());
//   target2->Write();
//   target2->SetName(Form("%s%s",name2.Data(),"divbyntriggers"));
//   target2->Scale(1.0/scalingfactorsig);
//   target2->Write();
//   
//   
//   outfile.Close();
// }

  
  
int AliCorrelation3p::Fill(const AliVParticle* ptrigger, const AliVParticle* p1, const AliVParticle* p2, const int weight)
{
  /// fill histograms from particles, fills each histogram exactly once.
  if (!ptrigger || !p1 || !p2) return -EINVAL;
  const double Pii=TMath::Pi();
  HistFill(GetNumberHist(kHistNassoc,fMBin,fVzBin),weight);
  
  // phi difference associated 1 to trigger particle
  Double_t DeltaPhi1 = ptrigger->Phi() - p1->Phi();
  while(DeltaPhi1<-0.5*Pii||DeltaPhi1>1.5*Pii){
    if (DeltaPhi1<-0.5*Pii) DeltaPhi1 += 2*Pii;
    if (DeltaPhi1>1.5*Pii)  DeltaPhi1 -= 2*Pii;
  }
  // phi difference associated 2 to trigger particle
  Double_t DeltaPhi2 = ptrigger->Phi() - p2->Phi();
  while(DeltaPhi2<-0.5||DeltaPhi2>1.5*Pii){
    if (DeltaPhi2<-0.5*Pii) DeltaPhi2 += 2*Pii;
    if (DeltaPhi2>1.5*Pii)  DeltaPhi2 -= 2*Pii;
  }
  // phi difference associated particles
  Double_t DeltaPhi12 = p1->Phi() - p2->Phi();
  if (DeltaPhi12<0.) DeltaPhi12*=-1.;

  // eta difference
  Double_t DeltaEta1  = ptrigger->Eta() - p1->Eta();
  Double_t DeltaEta2  = ptrigger->Eta() - p2->Eta();
  Double_t DeltaEta12 = p1      ->Eta() - p2->Eta();
  
//   if (TMath::Abs(DeltaEta1)<0.02  ||
//       TMath::Abs(DeltaEta2)<0.02  ||
//       TMath::Abs(DeltaEta12)<0.02) {
//     // avoid pair efficiency effects
//     // return 0;
//   }

//   if (fInvMassCut) {
//     TObjArray particles(3);
// 
//     // inv mass between trigger and either associated particle
//     // unfortunately, TObject does not allow const objects, but want to use the
//     // AliAnalysisCuts interface and parameter needs to be TObject type
//     particles.Add(const_cast<AliVParticle*>(ptrigger));
//     particles.Add(const_cast<AliVParticle*>(p1));
//     particles.Add(const_cast<AliVParticle*>(p2));
//     if (!fInvMassCut->IsSelected(&particles)){
//       particles.Clear();
//       return 0;
//     }
// 
//     // inv mass between associated particles
//     particles.Clear();
//     particles.Add(const_cast<AliVParticle*>(p1));
//     particles.Add(const_cast<AliVParticle*>(p2));
//     if (!fInvMassCut->IsSelected(&particles)){
//       particles.Clear();
//       return 0;
//     }
//   }
// 

  HistFill(GetNumberHist(khDeltaPhi,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
  if (p1->Charge()*p2->Charge()>0) {
    HistFill(GetNumberHist(khDeltaPhi_ass_like,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
    if (ptrigger->Charge()*p1->Charge()>0) {
      HistFill(GetNumberHist(khDeltaPhi_all_like,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
    } else {
      HistFill(GetNumberHist(khDeltaPhi_ass_like__trig_unlike,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
    }
  } else {
    HistFill(GetNumberHist(khDeltaPhi_ass_unlike,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
  }

  HistFill(GetNumberHist(khPhiEta,fMBin,fVzBin),DeltaEta12,DeltaPhi1);
  HistFill(GetNumberHist(khPhiPhiDEta,fMBin,fVzBin),DeltaEta12,DeltaPhi1,DeltaPhi2);
  if(weight>1)  HistFill(GetNumberHist(khPhiPhiDEtaScaled,fMBin,fVzBin),DeltaEta12,DeltaPhi1,DeltaPhi2,1.0/(weight-1));

  if (DeltaPhi12<fhPhiEtaDeltaPhi12Cut2) {
    HistFill(GetNumberHist(khPhiEta_phicut2,fMBin,fVzBin),DeltaEta12,DeltaPhi1);
  } 
  if (DeltaPhi12<fhPhiEtaDeltaPhi12Cut1) {
    HistFill(GetNumberHist(khPhiEta_phicut1,fMBin,fVzBin),DeltaEta12,DeltaPhi1);
  }
  if ((DeltaPhi1<0.5*Pii && DeltaPhi2<0.5*Pii )|| (DeltaPhi1>0.5*Pii &&DeltaPhi2>0.5*Pii)){
    HistFill(GetNumberHist(khPhiEta_sameside,fMBin,fVzBin),DeltaEta12,DeltaPhi1);
  }
  if (DeltaEta1 <0.0) DeltaEta1  *= -1.0;
  if (DeltaEta2 <0.0) DeltaEta2  *= -1.0;
  if (DeltaEta12<0.0) DeltaEta12 *= -1.0;
  if (DeltaEta12>1.0) {
    HistFill(GetNumberHist(khDeltaPhi_far,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
  } else if ((DeltaEta12>0.4)&&(DeltaEta12<=1.0)) {
    HistFill(GetNumberHist(khDeltaPhi_mid,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
  } else if (DeltaEta12<=0.4) {
    HistFill(GetNumberHist(khDeltaPhi_near,fMBin,fVzBin),DeltaPhi1,DeltaPhi2);
  }
  return 0;
}
