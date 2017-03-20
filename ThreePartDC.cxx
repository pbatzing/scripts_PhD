
#include "ThreeParticleorrectionFunctions.h"
#include "AliCorrelation3p.h"
#include "AliCorrelation3p_noQA.h"
#include "TSystem.h"
#include "THashList.h"
#include "TFileMerger.h"

//Macro to draw all relevant histograms in a folder structure and then collect and correct the resulting histograms.


void Draw(const char* options = ""){
  TString option = TString(options);
  //Possible options: 48 or 816 , Same, META or METrigger
  if(option.CompareTo("")==0){cout << "Please provide options to draw in the second argument. Possible options:\n" 
				   << "48, 816 or 416  			- decide the range of the trigger between 4-8 GeV/c or 8-16 GeV/c or both.\n"
				   << "Same,META or METrigger  	- decide to draw either Same event, META or METrigger historgrams.\n"
				   << "fakecor                  - fakes META, META2 and METrigger no matter what is in the file.\n"
				   << endl;return;}
  TFile * ffile = TFile::Open("AnalysisResults.root", "READ");
  if(!ffile) return;  
  TList * folderlist = ffile->GetListOfKeys();
  bool fakecor = false;
  if(option.Contains("fakecor"))fakecor = true;
  int remade = 0;
  //Loop over the folders:
  for(int i = 0;i<folderlist->GetEntries();i++){
    if(!TString(folderlist->At(i)->GetName()).BeginsWith("ThreePart"))continue;
    if(TString(folderlist->At(i)->GetName()).BeginsWith("ThreePart")&&((TString(folderlist->At(i)->GetName()).Contains("4_8")&&option.Contains("48")))||(TString(folderlist->At(i)->GetName()).Contains("8_16")&&option.Contains("816")||option.Contains("416"))){
      cout << folderlist->At(i)->GetName()<<endl;
      const char* foldername = folderlist->At(i)->GetName();
      TDirectory * folderdir = ffile->GetDirectory(foldername);
      TList * Containerlist = folderdir->GetListOfKeys();
      //Loop over the Correlation objects:
      for(int j=0;j<Containerlist->GetEntries();j++){
	if(TString(Containerlist->At(j)->GetName()).BeginsWith("ThreePart")){
	  const char* containername = Containerlist->At(j)->GetName();
	  cout << containername<<endl;
	  TList * containterdir = dynamic_cast<TList*>(folderdir->Get(containername));
	  //find out if there are META types:
	  bool isMETA = false;
	  bool isMETA2 = false;
	  bool isMETrigger = false;
	  for(int l = 0;l<containterdir->GetEntries();l++){	
	   if( TString(containterdir->At(l)->GetName()).Contains("META")&&!TString(containterdir->At(l)->GetName()).Contains("META2"))isMETA=true;
	   if( TString(containterdir->At(l)->GetName()).Contains("META2"))isMETA2=true;
	   if( TString(containterdir->At(l)->GetName()).Contains("METrigger"))isMETrigger=true;
	  }
	  if(!(isMETA&&isMETA2&&isMETrigger)&&!option.Contains("Same"))fakecor = true;//if there is missing anything, create them except in Same
	  
	  //Loop over the objects in the task:
	  for(int k = 0;k<containterdir->GetEntries();k++){	
	    if(TString(containterdir->At(k)->GetName()).Contains("trigger_correlation")){
	      //Same:
	      if(!TString(containterdir->At(k)->GetName()).Contains("ME")&&(option.Contains("Same")||(fakecor&&option.Contains("MET")))){
		AliCorrelation3p* signaltrack;
		AliCorrelation3p_noQA* signaltracknoQA;
		signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));
		signaltracknoQA = dynamic_cast<AliCorrelation3p_noQA*>(containterdir->At(k));
		if(option.Contains("Same")&&!option.Contains("MET")){
		  if(signaltrack){
		    if(j==0&&i==0&&remade ==0){signaltrack->MakeResultsFile(Form("%s",containername),true);remade =1;}
		    else signaltrack->MakeResultsFile(Form("%s",containername));
		  }
		  if(signaltracknoQA){
		    if(j==0&&i==0&&remade ==0){signaltracknoQA->MakeResultsFile(Form("%s",containername),true);remade =1;}
		    else signaltracknoQA->MakeResultsFile(Form("%s",containername));
		  }
		}
		if(option.Contains("META")&&!option.Contains("META2")){
		  if(signaltrack){
		      signaltrack->MakeResultsFile(Form("%s/META",containername),false,true);
		    }
		}
		if(option.Contains("META2")){
		  if(signaltrack){
		      signaltrack->MakeResultsFile(Form("%s/META2",containername),false,true);
		    }
		}		
		if(option.Contains("METrigger")){
		  if(signaltrack){
		      signaltrack->MakeResultsFile(Form("%s/META2",containername),false,true);
		    }
		}
	      }
	      //Meta:
	      if(TString(containterdir->At(k)->GetName()).Contains("META")&&option.Contains("META")&&!fakecor){
		AliCorrelation3p* signaltrack;
		AliCorrelation3p_noQA* signaltracknoQA;
		signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));//->Clone(Form("%sclone",containterdir->At(k)->GetName())));
		signaltracknoQA = dynamic_cast<AliCorrelation3p_noQA*>(containterdir->At(k));
		if(signaltrack){
		  if(TString(containterdir->At(k)->GetName()).Contains("META2"))	signaltrack->MakeResultsFile(Form("%s/META2",containername));
		  else 									signaltrack->MakeResultsFile(Form("%s/META",containername));
		}	
		if(signaltracknoQA){
		  if(TString(containterdir->At(k)->GetName()).Contains("META2"))	signaltracknoQA->MakeResultsFile(Form("%s/META2",containername));
		  else 									signaltracknoQA->MakeResultsFile(Form("%s/META",containername));
		}
	      }

	      //METrigger:
	      if(TString(containterdir->At(k)->GetName()).Contains("METrigger")&&option.Contains("METrigger")&&!fakecor){
		AliCorrelation3p* signaltrack;
		AliCorrelation3p_noQA* signaltracknoQA;
		signaltrack = dynamic_cast<AliCorrelation3p*>(containterdir->At(k));//->Clone(Form("%sclone",containterdir->At(k)->GetName())));
		signaltracknoQA = dynamic_cast<AliCorrelation3p_noQA*>(containterdir->At(k));
		if(signaltrack){
		  signaltrack->MakeResultsFile(Form("%s/METrigger",containername));
		}
		if(signaltracknoQA){
		  signaltracknoQA->MakeResultsFile(Form("%s/METrigger",containername));
		}
	      }
	    }
	  }//End Loop Objects in the task
	}
      }//End Loop Correlationobjects
    }
  }//End folder loop
  

  ffile->Close();
  delete ffile;
}

void Collect(const char* options = ""){
  TString option = TString(options);
//   {cout << "Please provide options to collect in the second argument. Possible options:\n" 
//        << "DivFirst  	- Collect the Same/Mixed first, then divide by the total number of triggers.\n"
//        << "pp or PbPb   - Changes the folder structures to accomodate for the collision system.\n"
//        << "folder=dir   - Adds dir to the folder structure.\n"
//        << "vertexcut=cut- only sums bins with a vertex>=cut.\n"
//        << "debug        - prints more output.\n"
//        << "010		- sums one centrality bin from 0 to 10 percent.\n"
//        << endl;}
  bool isPbPb = true;
  if(option.Contains("pp")) isPbPb = false;
  bool divfirst = true;
  if(!option.Contains("DivFirst")) divfirst = false; 
  bool debug = false;
  if(option.Contains("debug")) debug = true;
  bool b010 = false;
  if(option.Contains("010")) b010=true;
  TString folders = TString("");
  Double_t vertexcut = 10.0;
  if((option.Contains(" folder=")&&!option.Contains(" folder= "))||(option.Contains(" vertexcut=")&&!option.Contains(" vertexcut= "))){
    TObjArray * optionarray =option.Tokenize(" ");
    for(int i =0;i<optionarray->GetEntries();i++){
      if(dynamic_cast<TObjString*>(optionarray->At(i))->GetString().Contains("folder=")){
	TObjArray * folderarray = dynamic_cast<TObjString*>(optionarray->At(i))->GetString().Tokenize("=");
	folders.Append(dynamic_cast<TObjString*>(folderarray->At(1))->GetString());
	delete folderarray;
      }
      if(dynamic_cast<TObjString*>(optionarray->At(i))->GetString().Contains("vertexcut=")){
	TObjArray * vertexarray = dynamic_cast<TObjString*>(optionarray->At(i))->GetString().Tokenize("=");
	vertexcut = dynamic_cast<TObjString*>(vertexarray->At(1))->GetString().Atof();
	delete vertexarray;
      }      
    }
    delete optionarray;
  }
  if(isPbPb){
    if(divfirst){
      if(folders.CompareTo("")!=0){
	cout << "Collecting same/mixed centrality bins for all vertexbins with a vertex smaller then " << vertexcut << " and dividing by the total number of triggers after that into subfolder "<<folders.Data() << "."<<endl;
      }
      if(folders.CompareTo("")==0){
	cout << "Collecting same/mixed centrality bins for all vertexbins with a vertex smaller then " << vertexcut << " and dividing by the total number of triggers after that."<<endl;
      }
    }
    if(!divfirst){
      if(folders.CompareTo("")!=0){
	cout << "Collecting same and mixed in centrality bins for all vertexbins with a vertex smaller then " << vertexcut << " and calculating same/(ntriggers*mixed) after that into subfolder "<<folders.Data() << "."<<endl;
      }
      if(folders.CompareTo("")==0){
	cout << "Collecting same/(ntriggers*mixed) centrality bins for all vertexbins with a vertex smaller then " << vertexcut << " and calculating same/(ntriggers*mixed) after that into subfolder "<<folders.Data() << "."<<endl;
      }
    }
    
  }
  else{
//     cout << "pp needs more implementation."<<endl; 
  }

  TFile * Infile =  TFile::Open("results.root","READ");
  TFile * Outfile = TFile::Open("CollectedResults.root","UPDATE");
  TList * folderlist = Infile->GetListOfKeys();  
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    if(!TString(folder).Contains("ThreePartTracks"))continue;
    TObjArray * multdirlist = new TObjArray(6);

    TObjArray * namelist = TString(folder).Tokenize("_");
    TString foldername = TString("Correlation_Trigger_") + dynamic_cast<TObjString*>(namelist->At(1))->GetString() + TString("_") + dynamic_cast<TObjString*>(namelist->At(2))->GetString();
    if(folders.CompareTo(""))foldername.Append(folders);
    if(!Outfile->GetDirectory(foldername.Data())) Outfile->mkdir(foldername.Data());
    TString Subfoldername = TString("Associated_") + dynamic_cast<TObjString*>(namelist->At(3))->GetString() + TString("_") + dynamic_cast<TObjString*>(namelist->At(4))->GetString().Remove(1,8);
    if(Subfoldername.CompareTo("Associated_8_11")==0)Subfoldername = TString("Associated_8_16");
    TDirectory* dir= Outfile->GetDirectory(foldername.Data())->GetDirectory(Subfoldername.Data());
    if(!dir)dir = Outfile->GetDirectory(foldername.Data())->mkdir(Subfoldername.Data());

    if(!isPbPb){
      BinDirs * divsame = new BinDirs(dir, "" ,true);
      multdirlist->Add(divsame);
    }
    //List of directories for multiplicity bins:
    TList * directories = GetMZDirectories(Infile->GetDirectory(folder),vertexcut,debug);
    if(isPbPb){
      //Go through the list and find the multiplicity/centrality binning in order to sum over them and add them to the TObjArray:
      TString  s = TString("");
//       cout<< directories->GetEntries()<<endl;
      for(int i=0; i<directories->GetEntries();i++){
	TString totbin = TString(dynamic_cast<BinDirs*>(directories->At(i))->Same()->GetName());
	if(totbin.Contains("Z")){
	  if(totbin.Contains("BinM(0.00)->(5.00)")&&b010)totbin.ReplaceAll("BinM(0.00)->(5.00)","BinM(0.00)->(10.00)");
	  if(totbin.Contains("BinM(5.00)->(10.00)")&&b010) continue;
	  TString * bin = new TString(totbin.Tokenize("Z")->At(0)->GetName());
	  if(bin->CompareTo(s.Data())!=0&&bin->Contains("BinM")){
	    dir->mkdir(bin->Data());
	    BinDirs * dirs = new BinDirs(dir,bin->Data(),true);
	    multdirlist->Add(dirs);
// 	    cout << bin->Data()<<endl;
	    s.Clear();
	    s.Append(bin->Data());
	  }
	}
      }
    }
    //Get Tokens for all histograms:
    TStringToken histtokensbinstats = GetHistTokens(Infile->GetDirectory(folder)->GetDirectory("bin_stats"));
    while(histtokensbinstats.NextToken()){CollectHistbinstats(histtokensbinstats.Data(),directories,multdirlist,isPbPb);}
    TStringToken histtokens = GetHistTokens(dynamic_cast<BinDirs*>(directories->At(1))->SameDir("divided"));
    while(histtokens.NextToken()){CollectHist(histtokens.Data(),directories,multdirlist,false,divfirst,isPbPb);}
  }
  
  Infile->Close();
  Outfile->Close();
}

void Correct(const char* options = ""){
  TString option = TString(options);
  bool isPbPb = true;
  bool independentMETAscale = true;
  if(option.Contains("pp")) isPbPb = false;
  if(option.Contains("META12")) independentMETAscale = false;
  
  //open the file in UPDATE mode:
  TFile * rfile = TFile::Open("CollectedResults.root","UPDATE");
  TList * folderlist = rfile->GetListOfKeys();
  for(int j = 0;j<folderlist->GetEntries();j++){
    if(!TString(folderlist->At(j)->GetName()).Contains("Correlation"))continue;
    TDirectory * folderdir = dynamic_cast<TDirectory*>(dynamic_cast<TKey*>(folderlist->At(j))->ReadObj());
    //Subfolderlist (which associated:
    TList * subfolderlist = folderdir->GetListOfKeys();
    for(int k = 0; k<subfolderlist->GetEntries();k++){
      TDirectory * subfolderdir = dynamic_cast<TDirectory*>(dynamic_cast<TKey*>(subfolderlist->At(k))->ReadObj());
      //Find a list over all directories beginning with BinM and create a BinDirs object for each:
      TList * dirs = subfolderdir->GetListOfKeys();
      TObjArray * BinDirAr = new TObjArray();
      for(int i = 0;i<dirs->GetEntries();i++)
      {
	if(TString(dirs->At(i)->GetName()).BeginsWith("BinM")&&isPbPb){
	  BinDirs * dir = new BinDirs(subfolderdir,dirs->At(i)->GetName(),false);
	  BinDirAr->Add(dir);
	}
	if(TString(dirs->At(i)->GetName()).CompareTo("divided")==0&&!isPbPb){
	  BinDirs * dir = new BinDirs(subfolderdir->GetDirectory(""),subfolderdir->GetDirectory("META"),subfolderdir->GetDirectory("META2"),subfolderdir->GetDirectory("METrigger"),false);
	  BinDirAr->Add(dir);
	}
      }
      for(int i=0;i<BinDirAr->GetEntriesFast();i++){
        if(!isPbPb)	Correctpp(	dynamic_cast<BinDirs*>(BinDirAr->At(i)),independentMETAscale);
        if(isPbPb)	CorrectPbPb(	dynamic_cast<BinDirs*>(BinDirAr->At(i)),independentMETAscale);
      }
    }
    
    
  }
  rfile->Close();
}

void Runs11h(const char* options,TString workingdir){
  //macro to compare the different runs in LHC11h
  TString option = TString(options);
  TFileMerger* fileMerger = new TFileMerger();
  if(option.Contains("merge")){ 
    fileMerger->OutputFile("AnalysisResults.root");
  }
  TFile * Outfile = TFile::Open("LHC11hRuns.root","RECREATE");
  TDirectory * All = Outfile->mkdir("All");
  TDirectory * Cent05 = Outfile->mkdir("Cent05");
  TDirectory * Cent510 = Outfile->mkdir("Cent510");
  TDirectory * Cent1020 = Outfile->mkdir("Cent1020");
  TDirectory * Cent2040 = Outfile->mkdir("Cent2040");
  TDirectory * Cent4060 = Outfile->mkdir("Cent4060");
  TDirectory * Cent6080 = Outfile->mkdir("Cent6080");
  TDirectory * Cent8090 = Outfile->mkdir("Cent8090");
  TDirectory * Vertex_10_5 = Outfile->mkdir("Vertex_10_5");
  TDirectory * Vertex_5_2 = Outfile->mkdir("Vertex_5_2");
  TDirectory * Vertex_20 = Outfile->mkdir("Vertex_20");
  TDirectory * Vertex02 = Outfile->mkdir("Vertex02");
  TDirectory * Vertex25 = Outfile->mkdir("Vertex25");
  TDirectory * Vertex510 = Outfile->mkdir("Vertex510");

  Int_t GoodrunnumbersP11h[108] = {170593,170572,170388,170387,170315,170313,170312,170311,170309,170308,170306,170270,170269,170268,170230,170228,170207,170204,170203,170193,170163,170159,170155,170091,170089,170088,170085,170084,170083,170081,170040,170027,169965,169923,169859,169858,169855,169846,169838,169837,169835,169591,169590,169588,169587,169586,169557,169555,169554,169553,169550,169515,169512,169506,169504,169498,169475,169420,169419,169418,169417,169415,169411,169238,169167,169160,169156,169148,169145,169144,169138,169099,169094,169091,169045,169044,169040,169035,168992,168988,168826,168777,168514,168512,168511,168467,168464,168460,168458,168362,168361,168342,168341,168325,168322,168311,168310,168115,168108,168107,168105,168076,168069,167988,167987,167985,167920,167915};
  Int_t AllrunnumbersP11h[119]  = {167902,167903,167915,167920,167985,167987,167988,168066,168068,168069,168076,168104,168105,168107,168108,168115,168212,168310,168311,168322,168325,168341,168342,168361,168362,168458,168460,168461,168464,168467,168511,168512,168514,168777,168826,168984,168988,168992,169035,169040,169044,169045,169091,169094,169099,169138,169143,169144,169145,169148,169156,169160,169167,169238,169411,169415,169417,169418,169419,169420,169475,169498,169504,169506,169512,169515,169550,169553,169554,169555,169557,169586,169587,169588,169590,169591,169835,169837,169838,169846,169855,169858,169859,169923,169956,169965,170027,170036,170040,170081,170083,170084,170085,170088,170089,170091,170155,170159,170163,170193,170203,170204,170207,170228,170230,170268,170269,170270,170306,170308,170309,170311,170312,170313,170315,170387,170388,170572,170593};
  
  TH1D * eventsperrun 	= new TH1D("EventsperRun", "# Events per Run", 119, 0, 1);
  TH1D * nTriggers48 	= new TH1D("nTriggers48", "# Triggers (4-8) per Run", 119, 0, 1);
  TH1D * nTriggers816 	= new TH1D("nTriggers816", "# Triggers (4-8) per Run", 119, 0, 1);
  TH1D * nAss23 	= new TH1D("nAss23", "# Associated (2-3) per Run", 119, 0, 1);
  TH1D * nAss34 	= new TH1D("nAss34", "# Associated (3-4) per Run", 119, 0, 1);
  TH1D * nAss46 	= new TH1D("nAss46", "# Associated (4-6) per Run", 119, 0, 1);
  TH1D * nAss68 	= new TH1D("nAss68", "# Associated (6-8) per Run", 119, 0, 1);
  
  Double_t VZbins[6] = {-10.0,-5.0,-2.0,2.0,5.0,10};

  TH2D * eventsperrunCent05 	= new TH2D("EventsperRunCent05", "# Events per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nTriggers48Cent05 	= new TH2D("nTriggers48Cent05", "# Triggers (4-8) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nTriggers816Cent05 	= new TH2D("nTriggers816Cent05", "# Triggers (4-8) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss23Cent05 		= new TH2D("nAss23Cent05", "# Associated (2-3) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss34Cent05 		= new TH2D("nAss34Cent05", "# Associated (3-4) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss46Cent05 		= new TH2D("nAss46Cent05", "# Associated (4-6) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss68Cent05	 	= new TH2D("nAss68Cent05", "# Associated (6-8) per Run against vertex", 119, 0, 1,5,VZbins);
 
  TH2D * eventsperrunCent510 	= new TH2D("EventsperRunCent510", "# Events per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nTriggers48Cent510 	= new TH2D("nTriggers48Cent510", "# Triggers (4-8) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nTriggers816Cent510 	= new TH2D("nTriggers816Cent510", "# Triggers (4-8) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss23Cent510 		= new TH2D("nAss23Cent510", "# Associated (2-3) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss34Cent510 		= new TH2D("nAss34Cent510", "# Associated (3-4) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss46Cent510 		= new TH2D("nAss46Cent510", "# Associated (4-6) per Run against vertex", 119, 0, 1,5,VZbins);
  TH2D * nAss68Cent510	 	= new TH2D("nAss68Cent510", "# Associated (6-8) per Run against vertex", 119, 0, 1,5,VZbins);
  
  for(int i=0; i<119; i++){
//     if(i >=50) continue;
    Int_t color = 2;
    Int_t goodrunsintex = -1;
    for (int j=0;j<108;j++){if(AllrunnumbersP11h[i]==GoodrunnumbersP11h[j])color=4;goodrunsintex = j+1;}
    TString lable = Form("#color[%i]{%i}",color,AllrunnumbersP11h[i]);
    eventsperrun->GetXaxis()->SetBinLabel(i+1, lable);
    eventsperrun->GetXaxis()->LabelsOption("v");
    nTriggers48->GetXaxis()->SetBinLabel(i+1, lable);
    nTriggers48->GetXaxis()->LabelsOption("v");
    nTriggers816->GetXaxis()->SetBinLabel(i+1, lable);
    nTriggers816->GetXaxis()->LabelsOption("v");
    nAss23->GetXaxis()->SetBinLabel(i+1, lable);
    nAss23->GetXaxis()->LabelsOption("v");
    nAss34->GetXaxis()->SetBinLabel(i+1, lable);
    nAss34->GetXaxis()->LabelsOption("v");
    nAss46->GetXaxis()->SetBinLabel(i+1, lable);
    nAss46->GetXaxis()->LabelsOption("v");
    nAss68->GetXaxis()->SetBinLabel(i+1, lable);
    nAss68->GetXaxis()->LabelsOption("v");    
    eventsperrunCent05->GetXaxis()->SetBinLabel(i+1, lable);
    eventsperrunCent05->GetXaxis()->LabelsOption("v");
    nTriggers48Cent05->GetXaxis()->SetBinLabel(i+1, lable);
    nTriggers48Cent05->GetXaxis()->LabelsOption("v");
    nTriggers816Cent05->GetXaxis()->SetBinLabel(i+1, lable);
    nTriggers816Cent05->GetXaxis()->LabelsOption("v");
    nAss23Cent05->GetXaxis()->SetBinLabel(i+1, lable);
    nAss23Cent05->GetXaxis()->LabelsOption("v");
    nAss34Cent05->GetXaxis()->SetBinLabel(i+1, lable);
    nAss34Cent05->GetXaxis()->LabelsOption("v");
    nAss46Cent05->GetXaxis()->SetBinLabel(i+1, lable);
    nAss46Cent05->GetXaxis()->LabelsOption("v");
    nAss68Cent05->GetXaxis()->SetBinLabel(i+1, lable);
    nAss68Cent05->GetXaxis()->LabelsOption("v"); 
    eventsperrunCent510->GetXaxis()->SetBinLabel(i+1, lable);
    eventsperrunCent510->GetXaxis()->LabelsOption("v");
    nTriggers48Cent510->GetXaxis()->SetBinLabel(i+1, lable);
    nTriggers48Cent510->GetXaxis()->LabelsOption("v");
    nTriggers816Cent510->GetXaxis()->SetBinLabel(i+1, lable);
    nTriggers816Cent510->GetXaxis()->LabelsOption("v");
    nAss23Cent510->GetXaxis()->SetBinLabel(i+1, lable);
    nAss23Cent510->GetXaxis()->LabelsOption("v");
    nAss34Cent510->GetXaxis()->SetBinLabel(i+1, lable);
    nAss34Cent510->GetXaxis()->LabelsOption("v");
    nAss46Cent510->GetXaxis()->SetBinLabel(i+1, lable);
    nAss46Cent510->GetXaxis()->LabelsOption("v");
    nAss68Cent510->GetXaxis()->SetBinLabel(i+1, lable);
    nAss68Cent510->GetXaxis()->LabelsOption("v"); 
    
    TString nowdir = workingdir + TString(Form("%i/AnalysisResults.root",AllrunnumbersP11h[i]));
    if( gSystem->AccessPathName(nowdir.Data())||AllrunnumbersP11h[i]==000000 ) {
      cout << "File: " << nowdir.Data() << " not there!!!" <<  endl;
      continue;
      }
    TFile * Infilenow = TFile::Open(nowdir.Data(),"READ");
    if(Infilenow->IsOpen()){
      if(option.Contains("mergeall")){    
	printf("adding %s: ", nowdir.Data());
	fileMerger->AddFile(nowdir.Data());
      }
      else if(option.Contains("mergegood")&&goodrunsintex!=-1){
	printf("adding %s: ", nowdir.Data());
	fileMerger->AddFile(nowdir.Data());
      }
      cout << "File: " << nowdir.Data() << " exists!!!" <<  endl;
      TDirectory* dir1 = Infilenow->GetDirectory("ThreePartTracksPbPb_4_8");
      THashList * container11 = dynamic_cast<THashList*>(dir1->Get("ThreePartTracksPbPb_4_8_2_3Coutput1"));	
      THashList * container21 = dynamic_cast<THashList*>(dir1->Get("ThreePartTracksPbPb_4_8_3_4Coutput1"));
      THashList * container31 = dynamic_cast<THashList*>(dir1->Get("ThreePartTracksPbPb_4_8_4_6Coutput1"));
      THashList * container41 = dynamic_cast<THashList*>(dir1->Get("ThreePartTracksPbPb_4_8_6_8Coutput1"));
      TDirectory* dir2 = Infilenow->GetDirectory("ThreePartTracksPbPb_8_16");
      THashList * container12 = dynamic_cast<THashList*>(dir2->Get("ThreePartTracksPbPb_8_16_2_3Coutput1"));	
      int nevents = dynamic_cast<TH2D*>( container11->FindObject("centVsZVertex"))->Integral();
      eventsperrun->SetBinContent(i+1,nevents);
      int ntriggers48 = dynamic_cast<TH1D*>( container11->FindObject("trackTriggerPt"))->Integral();
      nTriggers48->SetBinContent(i+1,ntriggers48);
      int ntriggers816 = dynamic_cast<TH1D*>( container12->FindObject("trackTriggerPt"))->Integral();
      nTriggers816->SetBinContent(i+1,ntriggers816);
      int nass23 =dynamic_cast<TH1D*>( container11->FindObject("trackAssociatedPt"))->Integral();
      nAss23->SetBinContent(i+1,nass23);
      int nass34 =dynamic_cast<TH1D*>( container21->FindObject("trackAssociatedPt"))->Integral();
      nAss34->SetBinContent(i+1,nass34);
      int nass46 =dynamic_cast<TH1D*>( container31->FindObject("trackAssociatedPt"))->Integral();
      nAss46->SetBinContent(i+1,nass46);
      int nass68 =dynamic_cast<TH1D*>( container41->FindObject("trackAssociatedPt"))->Integral();
      nAss68->SetBinContent(i+1,nass68);
      
      //Cent05, Mbin = 0
      TH1D* tmphist;
      AliCorrelation3p * tmpCor1 = dynamic_cast<AliCorrelation3p*>(container11->FindObject("tracktrigger_correlation_4_8"));
      AliCorrelation3p * tmpCor2 = dynamic_cast<AliCorrelation3p*>(container21->FindObject("tracktrigger_correlation_4_8"));
      AliCorrelation3p * tmpCor3 = dynamic_cast<AliCorrelation3p*>(container31->FindObject("tracktrigger_correlation_4_8"));
      AliCorrelation3p * tmpCor4 = dynamic_cast<AliCorrelation3p*>(container41->FindObject("tracktrigger_correlation_4_8"));
      AliCorrelation3p * tmpCor5 = dynamic_cast<AliCorrelation3p*>(container12->FindObject("tracktrigger_correlation_8_16"));

      for(int Vz = 0;Vz <5; Vz++){
	tmphist = dynamic_cast<TH1D*>(tmpCor1->GetHistogram(AliCorrelation3p::khQAtocheckadressing,0,Vz,"tmp"));
	if(tmphist){
	  eventsperrunCent05->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor1->GetHistogram(AliCorrelation3p::kHistTriggerpT,0,Vz,"tmp"));
	if(tmphist){
	  nTriggers48Cent05->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor5->GetHistogram(AliCorrelation3p::kHistTriggerpT,0,Vz,"tmp"));
	if(tmphist){
	  nTriggers816Cent05->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor1->GetHistogram(AliCorrelation3p::kHistAssociatedpT,0,Vz,"tmp"));
	if(tmphist){
	  nAss23Cent05->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor2->GetHistogram(AliCorrelation3p::kHistAssociatedpT,0,Vz,"tmp"));
	if(tmphist){
	  nAss34Cent05->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor3->GetHistogram(AliCorrelation3p::kHistAssociatedpT,0,Vz,"tmp"));
	if(tmphist){
	  nAss46Cent05->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor4->GetHistogram(AliCorrelation3p::kHistAssociatedpT,0,Vz,"tmp"));
	if(tmphist){
	  nAss68Cent05->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}      
	
	tmphist = dynamic_cast<TH1D*>(tmpCor1->GetHistogram(AliCorrelation3p::khQAtocheckadressing,1,Vz,"tmp"));
	if(tmphist){
	  eventsperrunCent510->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor1->GetHistogram(AliCorrelation3p::kHistTriggerpT,1,Vz,"tmp"));
	if(tmphist){
	  nTriggers48Cent510->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor5->GetHistogram(AliCorrelation3p::kHistTriggerpT,1,Vz,"tmp"));
	if(tmphist){
	  nTriggers816Cent510->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor1->GetHistogram(AliCorrelation3p::kHistAssociatedpT,1,Vz,"tmp"));
	if(tmphist){
	  nAss23Cent510->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor2->GetHistogram(AliCorrelation3p::kHistAssociatedpT,1,Vz,"tmp"));
	if(tmphist){
	  nAss34Cent510->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor3->GetHistogram(AliCorrelation3p::kHistAssociatedpT,1,Vz,"tmp"));
	if(tmphist){
	  nAss46Cent510->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}
	tmphist = dynamic_cast<TH1D*>(tmpCor4->GetHistogram(AliCorrelation3p::kHistAssociatedpT,1,Vz,"tmp"));
	if(tmphist){
	  nAss68Cent510->SetBinContent(i+1,Vz+1,tmphist->Integral());
	  delete tmphist;
	}      
      }
      tmpCor1 = NULL;tmpCor2 = NULL;tmpCor3 = NULL;tmpCor4 = NULL;tmpCor5 = NULL;
    
      Infilenow->Close();
      delete container11;delete container21;delete container31;delete container41;delete container12;
    }

  }
  if(option.Contains("merge")){fileMerger->Merge();  delete fileMerger;}
  
  All->cd();
  eventsperrun->Sumw2();
  eventsperrun->Write();
  nTriggers48->Sumw2();
  nTriggers48->Write();
  nTriggers48->Divide(eventsperrun);
  nTriggers48->SetTitle("# Triggers (4-8) per event");
  nTriggers48->Write("nTriggers48pEvent");
  canvasruns("nTriggers48pEventC","#Triggers 4GeV/c-8GeV/c per Event.",nTriggers48,All);
  nTriggers816->Sumw2();
  nTriggers816->Write();
  nTriggers816->Divide(eventsperrun);
  nTriggers816->SetTitle("# Triggers (8-16) per event");
  nTriggers816->Write("nTriggers816pEvent");
  canvasruns("nTriggers816pEventC","#Triggers 8GeV/c-16GeV/c per Event.",nTriggers816,All);
  nAss23->Sumw2();
  nAss23->Write();
  nAss23->Divide(eventsperrun);
  nAss23->SetTitle("# Associated (2-3) per event");
  nAss23->Write("nAss23pEvent");
  canvasruns("nAss23pEventC","#Associated 2GeV/c-3GeV/c per Event.",nAss23,All);
  nAss34->Sumw2();
  nAss34->Write();
  nAss34->Divide(eventsperrun);
  nAss34->SetTitle("# Associated (3-4) per event");
  nAss34->Write("nAss34pEvent");
  canvasruns("nAss34pEventC","#Associated 3GeV/c-4GeV/c per Event.",nAss34,All);
  nAss46->Sumw2();
  nAss46->Write();
  nAss46->Divide(eventsperrun);
  nAss46->SetTitle("# Associated (4-6) per event");
  nAss46->Write("nAss46pEvent");
  canvasruns("nAss46pEventC","#Associated 4GeV/c-6GeV/c per Event.",nAss46,All);
  nAss46->Sumw2();
  nAss68->Write();
  nAss68->Divide(eventsperrun);
  nAss68->SetTitle("# Associated (6-8) per event");
  nAss68->Write("nAss68pEvent");  
  canvasruns("nAss68pEventC","#Associated 6GeV/c-8GeV/c per Event.",nAss68,All);

  Cent05->cd();
  eventsperrunCent05->Sumw2();
  eventsperrunCent05->Write();
  nTriggers48Cent05->Sumw2();
  nTriggers48Cent05->Write();
  nTriggers48Cent05->Divide(eventsperrunCent05);
  nTriggers48Cent05->Write("nTriggers48Cent05pEvent");
  canvasruns("nTriggers48pEventC","#Triggers 4GeV/c-8GeV/c per Event.",nTriggers816Cent05,Cent05);
  nTriggers816Cent05->Sumw2();
  nTriggers816Cent05->Write();
  nTriggers816Cent05->Divide(eventsperrunCent05);
  nTriggers816Cent05->Write("nTriggers816Cent05pEvent");
  canvasruns("nTriggers816Cent05pEventC","#Triggers 8GeV/c-16GeV/c per Event.",nTriggers816Cent05,Cent05);
  nAss23Cent05->Sumw2();
  nAss23Cent05->Write();
  nAss23Cent05->Divide(eventsperrunCent05);
  nAss23Cent05->Write("nAss23Cent05pEvent");  
  canvasruns("nAss23Cent05pEventC","#Associated 2GeV/c-3GeV/c per Event.",nAss23Cent05,Cent05);
  nAss34Cent05->Sumw2();
  nAss34Cent05->Write();
  nAss34Cent05->Divide(eventsperrunCent05);
  nAss34Cent05->Write("nAss34Cent05pEvent");   
  canvasruns("nAss34Cent05pEventC","#Associated 3GeV/c-4GeV/c per Event.",nAss34Cent05,Cent05);
  nAss46Cent05->Sumw2();
  nAss46Cent05->Write();
  nAss46Cent05->Divide(eventsperrunCent05);
  nAss46Cent05->Write("nAss46Cent05pEvent"); 
  canvasruns("nAss46Cent05pEventC","#Associated 4GeV/c-6GeV/c per Event.",nAss46Cent05,Cent05);
  nAss68Cent05->Sumw2();
  nAss68Cent05->Write();
  nAss68Cent05->Divide(eventsperrunCent05);
  nAss68Cent05->Write("nAss68Cent05pEvent"); 
  canvasruns("nAss68Cent05pEventC","#Associated 6GeV/c-8GeV/c per Event.",nAss68Cent05,Cent05);
 
  Cent510->cd();
  eventsperrunCent510->Sumw2();
  eventsperrunCent510->Write();
  nTriggers48Cent510->Sumw2();
  nTriggers48Cent510->Write();
  nTriggers48Cent510->Divide(eventsperrunCent510);
  nTriggers48Cent510->Write("nTriggers48Cent510pEvent"); 
  canvasruns("nTriggers48Cent510pEventC","#Triggers 4GeV/c-8GeV/c per Event.",nTriggers48Cent510,Cent510);
  nTriggers816Cent510->Sumw2();
  nTriggers816Cent510->Write();
  nTriggers816Cent510->Divide(eventsperrunCent510);
  nTriggers816Cent510->Write("nTriggers816Cent510pEvent");  
  canvasruns("nTriggers816Cent510pEventC","#Triggers 8GeV/c-16GeV/c per Event.",nTriggers816Cent510,Cent510);
  nAss23Cent510->Sumw2();
  nAss23Cent510->Write();
  nAss23Cent510->Divide(eventsperrunCent510);
  nAss23Cent510->Write("nAss23Cent510pEvent");  
  canvasruns("nAss23Cent510pEventC","#Associated 2GeV/c-3GeV/c per Event.",nAss23Cent510,Cent510);
  nAss34Cent510->Sumw2();
  nAss34Cent510->Write();
  nAss34Cent510->Divide(eventsperrunCent510);
  nAss34Cent510->Write("nAss34Cent510pEvent"); 
  canvasruns("nAss34Cent510pEventC","#Associated 3GeV/c-4GeV/c per Event.",nAss34Cent510,Cent510);
  nAss46Cent510->Sumw2();
  nAss46Cent510->Write();
  nAss46Cent510->Divide(eventsperrunCent510);
  nAss46Cent510->Write("nAss46Cent510pEvent"); 
  canvasruns("nAss46Cent510pEventC","#Associated 4GeV/c-6GeV/c per Event.",nAss46Cent510,Cent510);
  nAss68Cent510->Sumw2();
  nAss68Cent510->Write();
  nAss68Cent510->Divide(eventsperrunCent510);
  nAss68Cent510->Write("nAss68Cent510pEvent"); 
  canvasruns("nAss68Cent510pEventC","#Associated 6GeV/c-8GeV/c per Event.",nAss68Cent510,Cent510);

  Outfile->Close();
}

void pTBins(const char* options){
  //Function to draw comparisons between pT bins:
  TString option = TString(options);
  //open the file in UPDATE mode:
  TFile * rfile = TFile::Open("CollectedResults.root","READ");
  TFile * ofile = TFile::Open("Compare.root","RECREATE");
  TList * folderlist = rfile->GetListOfKeys();
  for(int j = 0;j<folderlist->GetEntries();j++){
    TDirectory * folderdir = dynamic_cast<TDirectory*>(dynamic_cast<TKey*>(folderlist->At(j))->ReadObj());
    if(!TString(folderdir->GetName()).BeginsWith("Correlation_Trigger_")) continue;
    TDirectory* pTdir =  resultsdirectory(ofile->GetDirectory(""),folderdir->GetName());
    TList * subfolderlist = folderdir->GetListOfKeys();
    pTdir->cd();
    TCanvas* canvas05 =  canvaspTbins("Canvas_0_5_perc_cent","Canvas, 0 to 5 % centrality",subfolderlist);
    TCanvas* canvas510 =  canvaspTbins("Canvas_5_10_perc_cent","Canvas, 5 to 10 % centrality",subfolderlist);
    TCanvas* canvas1020 =  canvaspTbins("Canvas_10_20_perc_cent","Canvas, 10 to 20 % centrality",subfolderlist);
    TCanvas* canvas2040 =  canvaspTbins("Canvas_20_40_perc_cent","Canvas, 20 to 40 % centrality",subfolderlist);
    
    
    for(int k = 0; k<subfolderlist->GetEntries();k++){
      TDirectory * subfolderdir = dynamic_cast<TDirectory*>(dynamic_cast<TKey*>(subfolderlist->At(k))->ReadObj());
      TDirectory * subpTdir = resultsdirectory(pTdir,subfolderdir->GetName());
      
      TCanvas * can = new TCanvas("tmp");
      can->Divide(2,2);
      TH2D* hist2dsame = PrepareHist(subfolderdir,"BinM(0.00)->(5.00)/divided","DPhi_1_DEta_12_SameSide","C=0% - 5%, divided","ZYAM"); 
      can->cd(1);
      hist2dsame->Draw("colz");
      can->cd(2);
      TH2D* hist2dmeta =PrepareHist(subfolderdir,"META/BinM(0.00)->(5.00)/divided","DPhi_1_DEta_12_SameSide","C=0% - 5%, META","ZYAM"); 
      hist2dmeta->Scale(dynamic_cast<TParameter<double>*>(subfolderdir->GetDirectory("BinM(0.00)->(5.00)")->GetDirectory("iteration1")->Get("METASCALE"))->GetVal());
      hist2dmeta->Draw("colz");
      can->cd(3);
      TH2D* hist2dmetrigger =PrepareHist(subfolderdir,"METrigger/BinM(0.00)->(5.00)/divided","DPhi_1_DEta_12_SameSide","C=0% - 5%, METRIGGER","ZYAM"); 
      hist2dmetrigger->Scale(dynamic_cast<TParameter<double>*>(subfolderdir->GetDirectory("BinM(0.00)->(5.00)")->GetDirectory("iteration1")->Get("METRIGGERSCALE"))->GetVal());
      hist2dmetrigger->Draw("colz");      
      can->cd(4);      
      TH2D* hist2dsamemmeta = dynamic_cast<TH2D*>(hist2dsame->Clone("hist2dsamemmeta"));
      hist2dsamemmeta->Add(hist2dmeta,-1.0);      
      hist2dsamemmeta->Add(hist2dmetrigger,-1.0);      
      hist2dsamemmeta->SetTitle("Same-META-METRIGGER");
      hist2dsamemmeta->Draw("colz");
      subpTdir->cd();
      can->Write("C0_5");
      delete hist2dsamemmeta;
      delete can;
      can = new TCanvas("tmp");
      can->Divide(2,2);
      can->cd(1);
      hist2dsame = PrepareHist(subfolderdir,"BinM(5.00)->(10.00)/divided","DPhi_1_DEta_12_SameSide","C=5% - 10%, divided","ZYAM"); 
      hist2dsame->Draw("colz");
      can->cd(2);
      hist2dmeta =PrepareHist(subfolderdir,"META/BinM(5.00)->(10.00)/divided","DPhi_1_DEta_12_SameSide","C=5% - 10%, META","ZYAM"); 
      hist2dmeta->Scale(dynamic_cast<TParameter<double>*>(subfolderdir->GetDirectory("BinM(5.00)->(10.00)")->GetDirectory("iteration1")->Get("METASCALE"))->GetVal());
      hist2dmeta->Draw("colz");
      can->cd(3);
      hist2dmetrigger =PrepareHist(subfolderdir,"METrigger/BinM(5.00)->(10.00)/divided","DPhi_1_DEta_12_SameSide","C=5% - 10%, METRIGGER","ZYAM"); 
      hist2dmetrigger->Scale(dynamic_cast<TParameter<double>*>(subfolderdir->GetDirectory("BinM(5.00)->(10.00)")->GetDirectory("iteration1")->Get("METRIGGERSCALE"))->GetVal());
      hist2dmetrigger->Draw("colz");      
      can->cd(4);      
      hist2dsamemmeta = dynamic_cast<TH2D*>(hist2dsame->Clone("hist2dsamemmeta"));
      hist2dsamemmeta->Add(hist2dmeta,-1.0);      
      hist2dsamemmeta->Add(hist2dmetrigger,-1.0);      
      hist2dsamemmeta->SetTitle("Same-META-METRIGGER");
      hist2dsamemmeta->Draw("colz");
      subpTdir->cd();
      can->Write("C5_10");
      delete hist2dsamemmeta;
      delete can;
      
      
      TH2D* SameSide = dynamic_cast<TH2D*>(subfolderdir->GetDirectory("BinM(0.00)->(5.00)")->GetDirectory("iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      SameSide->GetXaxis()->SetTitleOffset(1.0);
      SameSide->GetYaxis()->SetTitleOffset(1.0);
      SameSide->GetXaxis()->SetTitleSize(0.04);
      SameSide->GetYaxis()->SetTitleSize(0.04);
      RemovePlateau(SameSide->Integral(SameSide->GetXaxis()->FindBin(-1.0),SameSide->GetXaxis()->FindBin(1.0),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5))/(SameSide->GetXaxis()->FindBin(1.0)-SameSide->GetXaxis()->FindBin(-1.0)),SameSide);
      canvas05->cd(k+1);
      SameSide->Draw("colz");
      
      SameSide = dynamic_cast<TH2D*>(subfolderdir->GetDirectory("BinM(5.00)->(10.00)")->GetDirectory("iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      SameSide->GetXaxis()->SetTitleOffset(1.0);
      SameSide->GetYaxis()->SetTitleOffset(1.0);
      SameSide->GetXaxis()->SetTitleSize(0.04);
      SameSide->GetYaxis()->SetTitleSize(0.04);
      RemovePlateau(SameSide->Integral(SameSide->GetXaxis()->FindBin(-1.0),SameSide->GetXaxis()->FindBin(1.0),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5))/(SameSide->GetXaxis()->FindBin(1.0)-SameSide->GetXaxis()->FindBin(-1.0)),SameSide);
      canvas510->cd(k+1);
      SameSide->Draw("colz");
      
      SameSide = dynamic_cast<TH2D*>(subfolderdir->GetDirectory("BinM(10.00)->(20.00)")->GetDirectory("iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      SameSide->GetXaxis()->SetTitleOffset(1.0);
      SameSide->GetYaxis()->SetTitleOffset(1.0);
      SameSide->GetXaxis()->SetTitleSize(0.04);
      SameSide->GetYaxis()->SetTitleSize(0.04);
      RemovePlateau(SameSide->Integral(SameSide->GetXaxis()->FindBin(-1.0),SameSide->GetXaxis()->FindBin(1.0),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5))/(SameSide->GetXaxis()->FindBin(1.0)-SameSide->GetXaxis()->FindBin(-1.0)),SameSide);
      canvas1020->cd(k+1);
      
      SameSide->Draw("colz");      
      SameSide = dynamic_cast<TH2D*>(subfolderdir->GetDirectory("BinM(20.00)->(40.00)")->GetDirectory("iteration1")->Get("DPhi_1_DEta_12_SameSide"));
      SameSide->GetXaxis()->SetTitleOffset(1.0);
      SameSide->GetYaxis()->SetTitleOffset(1.0);
      SameSide->GetXaxis()->SetTitleSize(0.04);
      SameSide->GetYaxis()->SetTitleSize(0.04);
      RemovePlateau(SameSide->Integral(SameSide->GetXaxis()->FindBin(-1.0),SameSide->GetXaxis()->FindBin(1.0),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5),SameSide->GetYaxis()->FindBin(TMath::Pi()*0.5))/(SameSide->GetXaxis()->FindBin(1.0)-SameSide->GetXaxis()->FindBin(-1.0)),SameSide);
      canvas2040->cd(k+1);
      SameSide->Draw("colz"); 
      
    }
    pTdir->cd();
    canvas05->Write();
    canvas510->Write();
    canvas1020->Write();
    canvas2040->Write();
    delete canvas05;delete canvas510; delete canvas1020; delete canvas2040;
  }
  rfile->Close();
  ofile->Close();
}

void yield(const char* options){
  bool ispp = false;
  bool isPbPb = false;
  bool etasubstracted = false;
  TString delimiter(" ");
  TStringToken token(options, delimiter);
  while (token.NextToken()) {
    TString argument=token;
    if (argument.CompareTo("-h")==0 ||argument.CompareTo("--help")==0) {
      cout<<Form(" options:"
		    "\n\t  pp   - uses pp fitting ranges to substract the correlated background."
		    "\n\t  PbPb     - uses PbPb fitting ranges to substract the correlated background."
		    )<<endl;
      return;
      }
    if(argument.CompareTo("pp")==0){
      ispp = true;
      isPbPb = false;
      
    }
    if(argument.CompareTo("PbPb")==0){
      isPbPb = true;
      ispp = false;
   }
   if(argument.CompareTo("Eta")==0) etasubstracted = true;
  }
  if(!(ispp||isPbPb))return;
  
  TFile * rfile = TFile::Open("CollectedResults.root","UPDATE");
  TObjArray * dirarray = new TObjArray();
  TObjArray * yieldarray = new TObjArray();
  TList * folderlist = rfile->GetListOfKeys();
  for(int l = 0;l<folderlist->GetEntries();l++){
    const char* folder = folderlist->At(l)->GetName();  
    if(!TString(folder).Contains("Correlation_Trigger_"))continue;
    TDirectory * folderdir = rfile->GetDirectory(folder);
    TList * dirs = folderdir->GetListOfKeys();
    for(int i = 0;i<dirs->GetEntries();i++){
      if(!TString(dirs->At(i)->GetName()).BeginsWith("Associated"))continue;
      if(ispp){
	if(!etasubstracted){
	  unsigned int k=1;
	  while(k<2){
	    TDirectory * tmp1 = dynamic_cast<TDirectory*>(folderdir->Get(Form("%s/iteration%u",dirs->At(i)->GetName(),k)));	  
	    if(!tmp1){k = 11;continue;}
	    TDirectory * tmp2 = dynamic_cast<TDirectory*>(folderdir->Get(Form("%s/iteration%u",dirs->At(i)->GetName(),k+1)));
	    if(tmp1&&!tmp2){
	      dirarray->Add(tmp1);
	      yieldarray->Add(resultsdirectory(folderdir->GetDirectory(Form("%s",dirs->At(i)->GetName())),"yield"));
	      k = 11;
	    }
	    if(tmp1&&tmp2){k +=1;}
	  }
	}
	else{
	 TDirectory * tmp1 = dynamic_cast<TDirectory*>(folderdir->Get(Form("%s/Etasubstracted",dirs->At(i)->GetName())));
	 if(tmp1){
	   dirarray->Add(tmp1);
	   yieldarray->Add(resultsdirectory(folderdir->GetDirectory(Form("%s/Etasubstracted",dirs->At(i)->GetName())),"yield"));
	 }
	}
      }

      if(isPbPb){
	TDirectory* subfolderdir= folderdir->GetDirectory(dirs->At(i)->GetName());
	TList * subfolderlist = subfolderdir->GetListOfKeys();
	for(int j = 0;j<subfolderlist->GetEntries();j++){
	  if(TString(subfolderlist->At(j)->GetName()).BeginsWith("Bin")){
	    if(!etasubstracted){
	      unsigned int k=1;
// 	      cout << subfolderlist->At(j)->GetName()<<endl;
	      while(k<2){
		TDirectory * tmp1 = dynamic_cast<TDirectory*>(subfolderdir->Get(Form("%s/iteration%u",subfolderlist->At(j)->GetName(),k)));
		if(!tmp1){k = 11;continue;}
		TDirectory * tmp2 = dynamic_cast<TDirectory*>(subfolderdir->Get(Form("%s/iteration%u",subfolderlist->At(j)->GetName(),k+1)));
		if(tmp1&&!tmp2){
		  dirarray->Add(tmp1);
		  yieldarray->Add(resultsdirectory(subfolderdir->GetDirectory(Form("%s",subfolderlist->At(j)->GetName())),"yield"));
		  k = 11;
		}
		if(tmp1&&tmp2){k +=1;}
	      }
	    }
	    else{
	      TDirectory * tmp1 = dynamic_cast<TDirectory*>(subfolderdir->Get(Form("%s/Etasubstracted",subfolderlist->At(j)->GetName())));
	      if(tmp1){
		dirarray->Add(tmp1);
		yieldarray->Add(resultsdirectory(subfolderdir->GetDirectory(Form("%s/Etasubstracted",subfolderlist->At(j)->GetName())),"yield"));
	      }
	    }
	  }
	}
      }
    }
  }
  
  Double_t etalimit = 0.0;
  if(ispp) etalimit = 1.4;
  if(isPbPb) etalimit = 1.0;
  TDirectory * dir;
  TDirectory * yielddir;
  for(int i=0;i<dirarray->GetEntriesFast();i++){
    dir = dynamic_cast<TDirectory*>(dirarray->At(i));
//     dir->pwd();
    yielddir = dynamic_cast<TDirectory*>(yieldarray->At(i));
    extractbinyield(dir,yielddir,etalimit,options);
  }
  rfile->Close();
  delete yieldarray;delete dirarray;
}

void Present(const char* options){
  TString option = TString(options); 
  bool is11a = false;
  bool is10d = false;
  bool is10h = false;
  bool is11h = false;
  if(option.Contains("11a")) is11a = true;
  if(option.Contains("10d")) is10d = true;
  if(option.Contains("10h")) is10h = true;
  if(option.Contains("11h")) is11h = true;
  TString mode = TString("");
  if(option.Contains(" mode=")&&!option.Contains(" mode= ")){
    TObjArray * optionarray =option.Tokenize(" ");
    for(int i =0;i<optionarray->GetEntries();i++){
      if(dynamic_cast<TObjString*>(optionarray->At(i))->GetString().Contains("mode=")){
	TObjArray * folderarray = dynamic_cast<TObjString*>(optionarray->At(i))->GetString().Tokenize("=");
	mode.Append(dynamic_cast<TObjString*>(folderarray->At(1))->GetString());
	delete folderarray;
      }
    }
    delete optionarray;
  }
  else{
    cout << "Please provide a mode in the options on the form 'mode=$mode' where $mode can be divided, iteration#,META,META2,METrigger or yield."<<endl;
  }  
  
  if((is11a && is10h)||(is11a&&is11h)||is10h&&is11h||(is10d&&is11a)||(is10d&&is11h)||(is10d&&is10h))return;
  if(!(is11a||is10d||is10h||is11h)) return;
  
  TString * path = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/results/compare/results.root");
  TString * inpath = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  if(is10d)inpath->Append("LHC10d/Train/CollectedResults.root");
  if(is11a)inpath->Append("LHC11a/Train/CollectedResults.root");
  if(is10h)inpath->Append("LHC10h/Train/CollectedResults.root");
  if(is11h)inpath->Append("LHC11h/Train/CollectedResults.root");
  
  TFile * infile =  TFile::Open(inpath->Data(),"READ");
  TFile * Outfile =  TFile::Open(path->Data(),"UPDATE");
  TDirectory * outdir = NULL;

  if(is10d){
    outdir = resultsdirectory(Outfile->GetDirectory(""),"LHC10d");
  }
  if(is10h){
    outdir = resultsdirectory(Outfile->GetDirectory(""),"LHC10h");
  }  
  if(is11a){
    outdir = resultsdirectory(Outfile->GetDirectory(""),"LHC11a");
  }
  if(is11h){
    outdir = resultsdirectory(Outfile->GetDirectory(""),"LHC11h");
  }  
  if(is11h||is10h)PresentPbPb(infile,outdir,mode);  
  if(is11a||is10d)Presentpp(infile,outdir,mode);  
  infile->Close();Outfile->Close();
}

void pTBinsPer(const char* options){
  TString * path = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/results/compare/results.root");
  TString * inpath10h = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/LHC10h/TrainFine/CollectedResults.root");
  TString * inpath11h = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/LHC11h/TrainFine/CollectedResults.root");
  TString * inpath11a = new TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/LHC11a/TrainFine/CollectedResults.root");
  TFile * infile10h =  TFile::Open(inpath10h->Data(),"READ");
  TFile * infile11h =  TFile::Open(inpath11h->Data(),"READ");
  TFile * infile11a =  TFile::Open(inpath11a->Data(),"READ");
  TFile * Outfile =  TFile::Open(path->Data(),"UPDATE");
  
  TDirectory * outdir = resultsdirectory(Outfile->GetDirectory(""),"pTbins");
  TObjArray  * arr10h = new TObjArray();
  TObjArray  * arr11h = new TObjArray();
  TObjArray  * arr11a = new TObjArray();
  TObjArray  * arr10hsvert = new TObjArray();
  TObjArray  * arr11hsvert = new TObjArray();
  TObjArray  * arr11asvert = new TObjArray();
    
  arr10h->Add(infile10h->GetDirectory("Correlation_Trigger_8_16/Associated_2_3/BinM(60.00)->(90.00)"));
  arr10h->Add(infile10h->GetDirectory("Correlation_Trigger_8_16/Associated_3_4/BinM(60.00)->(90.00)"));  
  arr10h->Add(infile10h->GetDirectory("Correlation_Trigger_8_16/Associated_4_6/BinM(60.00)->(90.00)"));    
  arr10h->Add(infile10h->GetDirectory("Correlation_Trigger_8_16/Associated_6_8/BinM(60.00)->(90.00)"));  
  arr10hsvert->Add(infile10h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_2_3/BinM(60.00)->(90.00)"));
  arr10hsvert->Add(infile10h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_3_4/BinM(60.00)->(90.00)"));  
  arr10hsvert->Add(infile10h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_4_6/BinM(60.00)->(90.00)"));    
  arr10hsvert->Add(infile10h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_6_8/BinM(60.00)->(90.00)"));   
  
  arr11h->Add(infile11h->GetDirectory("Correlation_Trigger_8_16/Associated_2_3/BinM(0.00)->(5.00)"));
  arr11h->Add(infile11h->GetDirectory("Correlation_Trigger_8_16/Associated_3_4/BinM(0.00)->(5.00)"));  
  arr11h->Add(infile11h->GetDirectory("Correlation_Trigger_8_16/Associated_4_6/BinM(0.00)->(5.00)"));    
  arr11h->Add(infile11h->GetDirectory("Correlation_Trigger_8_16/Associated_6_8/BinM(0.00)->(5.00)"));  
  arr11hsvert->Add(infile11h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_2_3/BinM(0.00)->(5.00)"));
  arr11hsvert->Add(infile11h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_3_4/BinM(0.00)->(5.00)"));  
  arr11hsvert->Add(infile11h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_4_6/BinM(0.00)->(5.00)"));    
  arr11hsvert->Add(infile11h->GetDirectory("Correlation_Trigger_8_16sVert/Associated_6_8/BinM(0.00)->(5.00)"));   

  arr11a->Add(infile11a->GetDirectory("Correlation_Trigger_8_16/Associated_2_3"));
  arr11a->Add(infile11a->GetDirectory("Correlation_Trigger_8_16/Associated_3_4"));  
  arr11a->Add(infile11a->GetDirectory("Correlation_Trigger_8_16/Associated_4_6"));    
  arr11a->Add(infile11a->GetDirectory("Correlation_Trigger_8_16/Associated_6_8"));  
  arr11asvert->Add(infile11a->GetDirectory("Correlation_Trigger_8_16sVert/Associated_2_3"));
  arr11asvert->Add(infile11a->GetDirectory("Correlation_Trigger_8_16sVert/Associated_3_4"));  
  arr11asvert->Add(infile11a->GetDirectory("Correlation_Trigger_8_16sVert/Associated_4_6"));    
  arr11asvert->Add(infile11a->GetDirectory("Correlation_Trigger_8_16sVert/Associated_6_8"));   
  
  for(int i=0; i<arr10h->GetEntries();i++){
    TStringToken tokens("divided mixed_event iteration1"," ");
    TDirectory* dir10h   = dynamic_cast<TDirectory*>(arr10h->At(i));
    TDirectory* dir10hvz = dynamic_cast<TDirectory*>(arr10hsvert->At(i));
    TDirectory* dir11h   = dynamic_cast<TDirectory*>(arr11h->At(i));
    TDirectory* dir11hvz = dynamic_cast<TDirectory*>(arr11hsvert->At(i));
    TDirectory* dir11a   = dynamic_cast<TDirectory*>(arr11a->At(i));
    TDirectory* dir11avz = dynamic_cast<TDirectory*>(arr11asvert->At(i));
    TDirectory* ptdir = outdir;
    if(i==0)ptdir = resultsdirectory(outdir,"8_16_2_3");
    if(i==1)ptdir = resultsdirectory(outdir,"8_16_3_4");
    if(i==2)ptdir = resultsdirectory(outdir,"8_16_4_6");
    if(i==3)ptdir = resultsdirectory(outdir,"8_16_6_8");
    
    while(tokens.NextToken()){
      TCanvas* canall  = new TCanvas("allbins");
      canall->Divide(2,2);
      TCanvas* canvert = new TCanvas("vertbins");
      canvert->Divide(2,2);
      TCanvas* candiv  = new TCanvas("divbins");
      candiv->Divide(2,2);
      TDirectory* outsub = resultsdirectory(ptdir,tokens);
      TH2D * PhiEta12ss_10h = dynamic_cast<TH2D*>(dir10h->GetDirectory(tokens)->Get("DPhi_1_DEta_12_SameSide"));
      TH2D * PhiEta12ss_10hvz = dynamic_cast<TH2D*>(dir10hvz->GetDirectory(tokens)->Get("DPhi_1_DEta_12_SameSide"));
      TH2D * PhiEta12ss_11h = dynamic_cast<TH2D*>(dir11h->GetDirectory(tokens)->Get("DPhi_1_DEta_12_SameSide"));
      TH2D * PhiEta12ss_11hvz = dynamic_cast<TH2D*>(dir11hvz->GetDirectory(tokens)->Get("DPhi_1_DEta_12_SameSide"));      
      TH2D * PhiEta12ss_11a = dynamic_cast<TH2D*>(dir11a->GetDirectory(tokens)->Get("DPhi_1_DEta_12_SameSide"));
      TH2D * PhiEta12ss_11avz = dynamic_cast<TH2D*>(dir11avz->GetDirectory(tokens)->Get("DPhi_1_DEta_12_SameSide"));        
      
      TH2D* PhiEta12ss_10h_c = dynamic_cast<TH2D*>(PhiEta12ss_10h->Clone("DPhi_1_DEta_12_SameSide_C"));
      TH2D* PhiEta12ss_10hvz_c = dynamic_cast<TH2D*>(PhiEta12ss_10hvz->Clone("DPhi_1_DEta_12_SameSide_C2"));
      
      canall->cd(1);
      PhiEta12ss_11h->SetTitle("0%-5% (LHC11h)");
      PhiEta12ss_11h->Draw("colz");
      canall->cd(2);
      PhiEta12ss_10h->SetTitle("60%-90% (LHC10h)");      
      PhiEta12ss_10h->Draw("colz");
      canall->cd(3);
      PhiEta12ss_11a->SetTitle("minimum bias pp (LHC11a)");      
      PhiEta12ss_11a->Draw("colz"); 
      canall->cd(4);
      PhiEta12ss_10h_c->SetTitle("60%-90% PbPb/pp");      
      PhiEta12ss_10h_c->Divide(PhiEta12ss_11a);
      PhiEta12ss_10h_c->Draw("colz");
      outsub->cd();
      canall->Write();
      delete canall;
      canvert->cd(1);
      PhiEta12ss_11hvz->SetTitle("0%-5% (LHC11h)");
      PhiEta12ss_11hvz->Draw("colz");
      canvert->cd(2);
      PhiEta12ss_10hvz->SetTitle("60%-90% (LHC10h)");      
      PhiEta12ss_10hvz->Draw("colz");
      canvert->cd(3);
      PhiEta12ss_11avz->SetTitle("minimum bias pp (LHC11a)");      
      PhiEta12ss_11avz->Draw("colz"); 
      canvert->cd(4);
      PhiEta12ss_10hvz_c->SetTitle("60%-90% PbPb/pp");      
      PhiEta12ss_10hvz_c->Divide(PhiEta12ss_11avz);
      PhiEta12ss_10hvz_c->Draw("colz");
      outsub->cd();
      canvert->Write();
      delete canvert;
      candiv->cd(1);
      PhiEta12ss_11hvz->SetTitle("0%-5% (LHC11h)");
      PhiEta12ss_11hvz->Divide(PhiEta12ss_11h);
      PhiEta12ss_11hvz->Draw("colz");
      candiv->cd(2);
      PhiEta12ss_10hvz->SetTitle("60%-90% (LHC10h)");      
      PhiEta12ss_10hvz->Divide(PhiEta12ss_10h);
      PhiEta12ss_10hvz->Draw("colz");
      candiv->cd(3);
      PhiEta12ss_11avz->SetTitle("minimum bias pp (LHC11a)");      
      PhiEta12ss_11avz->Divide(PhiEta12ss_11a);
      PhiEta12ss_11avz->Draw("colz"); 
      outsub->cd();
      candiv->Write();      
      delete candiv;
      delete PhiEta12ss_10h;delete PhiEta12ss_10hvz;
      delete PhiEta12ss_11a;delete PhiEta12ss_11avz;
      delete PhiEta12ss_11h;delete PhiEta12ss_11hvz;
      delete PhiEta12ss_10h_c;
    }
  }
  infile10h->Close();
  infile11a->Close();
  infile11h->Close();
  Outfile->Close();
}


void ComparisonFile(const char* options){
  //This takes input from subfolder /pp and /PbPb to compare centrality bins, PbPb with pp and pT bins with each other.
  TString workingdir 	= TString(gSystem->WorkingDirectory());
  TFile * Infile10d  = NULL;
  TFile * Infile11a = NULL;
  TFile * Infile10h = NULL;
  TFile * Infile11h = NULL;
  if(!gSystem->AccessPathName(Form("%s/LHC10d/CollectedResults.root",workingdir.Data()))) Infile10d 	= TFile::Open(Form("%s/LHC10d/CollectedResults.root",workingdir.Data()),"READ");
  if(!gSystem->AccessPathName(Form("%s/LHC11a/CollectedResults.root",workingdir.Data()))) Infile11a 	= TFile::Open(Form("%s/LHC11a/CollectedResults.root",workingdir.Data()),"READ");
  if(!gSystem->AccessPathName(Form("%s/LHC10h/CollectedResults.root",workingdir.Data()))) Infile10h 	= TFile::Open(Form("%s/LHC10h/CollectedResults.root",workingdir.Data()),"READ");
  if(!gSystem->AccessPathName(Form("%s/LHC11h/CollectedResults.root",workingdir.Data()))) Infile11h 	= TFile::Open(Form("%s/LHC11h/CollectedResults.root",workingdir.Data()),"READ");
  TFile * Outfile 	= TFile::Open("Compare.root","RECREATE");
  if(!Infile10d||!Infile11a||!Infile11h||!Infile10h) return; // We want to compare all of them.


//collect all bins for pp
  TObjArray * LHC10dfolderarray = new TObjArray();
  TList * folderlist = Infile10d->GetListOfKeys();  
  
  int nptbins10d = 0;
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    if(!TString(folder).Contains("Correlation_Trigger"))continue;
    TDirectory * subdir = Infile10d->GetDirectory(folder);
    TList * subfolderlist = subdir->GetListOfKeys();
    for(int k = 0;k<subfolderlist->GetEntries();k++){
      const char* subfolder = subfolderlist->At(k)->GetName();
      if(!TString(subfolder).Contains("Associated"))continue;
      nptbins10d += 1;
      TDirectory * subsubdir = subdir->GetDirectory(subfolder);
      subsubdir->SetTitle(Form("ppLHC10d%s_%s",folder,subfolder));
      LHC10dfolderarray->Add(subsubdir);
    }
  }
  TObjArray * LHC11afolderarray = new TObjArray();
  folderlist = Infile11a->GetListOfKeys();  
  int nptbins11a = 0;
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    if(!TString(folder).Contains("Correlation_Trigger"))continue;
    TDirectory * subdir = Infile11a->GetDirectory(folder);
    TList * subfolderlist = subdir->GetListOfKeys();
    for(int k = 0;k<subfolderlist->GetEntries();k++){
      const char* subfolder = subfolderlist->At(k)->GetName();
      if(!TString(subfolder).Contains("Associated"))continue;
      nptbins11a += 1;
      TDirectory * subsubdir = subdir->GetDirectory(subfolder);
      subsubdir->SetTitle(Form("ppLHC11a%s_%s",folder,subfolder));
      LHC11afolderarray->Add(subsubdir);
    }
  }
  //collect all bins for PbPb
  TObjArray * LHC10hfolderarray = new TObjArray();
  int nptbins10h = 0;
  folderlist = Infile10h->GetListOfKeys();  
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    if(!TString(folder).Contains("Correlation_Trigger"))continue;
    TDirectory * subdir = Infile10h->GetDirectory(folder);
    TList * subfolderlist = subdir->GetListOfKeys();
    for(int k = 0;k<subfolderlist->GetEntries();k++){
      const char* subfolder = subfolderlist->At(k)->GetName();
      if(!TString(subfolder).Contains("Associated"))continue;
      nptbins10h += 1;
      TDirectory * subsubdir = subdir->GetDirectory(subfolder);
      subsubdir->SetTitle(Form("PbPbLHC10h%s_%s",folder,subfolder));
      LHC10hfolderarray->Add(subsubdir);
    }
  }  
  TObjArray * LHC11hfolderarray = new TObjArray();
  int nptbins11h = 0;
  folderlist = Infile11h->GetListOfKeys();  
  for(int j = 0;j<folderlist->GetEntries();j++){
    const char* folder = folderlist->At(j)->GetName();
    if(!TString(folder).Contains("Correlation_Trigger"))continue;
    TDirectory * subdir = Infile11h->GetDirectory(folder);
    TList * subfolderlist = subdir->GetListOfKeys();
    for(int k = 0;k<subfolderlist->GetEntries();k++){
      const char* subfolder = subfolderlist->At(k)->GetName();
      if(!TString(subfolder).Contains("Associated"))continue;
      nptbins11h += 1;
      TDirectory * subsubdir = subdir->GetDirectory(subfolder);
      subsubdir->SetTitle(Form("PbPbLHC11h%s_%s",folder,subfolder));
      LHC11hfolderarray->Add(subsubdir);
    }
  }  
    
  //Draw all pt bins for pp on one canvas:
  TDirectory * outdirpppt = Outfile->mkdir("PP_LHC10d_pT_bins");
  TString options2 = TString(Form("pp_pTbins_%i",nptbins10d));
  Compare(LHC10dfolderarray,outdirpppt,"DPHIDPHI#DETA_12#DETA_12ss#yield#nearpeak#awaypeak",options2,"");
  options2.Clear();
  TDirectory * outdirpppt2 = Outfile->mkdir("PP_LHC11a_pT_bins");
  options2 = TString(Form("pp_pTbins_%i",nptbins11a));
  Compare(LHC11afolderarray,outdirpppt2,"DPHIDPHI#DETA_12#DETA_12ss#yield#nearpeak#awaypeak",options2,"");  
  

  //Draw all pt bins for PbPb on one canvas:  
  options2.Clear();
  options2.Append(Form("PbPb_pTbins_%i",nptbins10h));
  TDirectory * outdirPbPbpt = Outfile->mkdir("PbPb_LHC10h_pT_bins");
  Compare(LHC10hfolderarray,outdirPbPbpt,"DPHIDPHI#DETA_12#DETA_12ss#yield#nearpeak#awaypeak",options2,"");

  //Draw all pt bins for PbPb on one canvas:  
  options2.Clear();
  options2.Append(Form("PbPb_pTbins_%i",nptbins11h));
  TDirectory * outdirPbPbpt2 = Outfile->mkdir("PbPb_LHC11h_pT_bins");
  Compare(LHC11hfolderarray,outdirPbPbpt2,"DPHIDPHI#DETA_12#DETA_12ss#yield#nearpeak#awaypeak",options2,"");

  //to compare yields we need both arrays:
  TObjArray * comparearray = new TObjArray();
  for(int i = 0; i<LHC10dfolderarray->GetEntries();i++){
    TDirectory* dir = dynamic_cast<TDirectory*>(LHC10dfolderarray->At(i));
    comparearray->Add(dir);
  }

  for(int i = 0; i<LHC11hfolderarray->GetEntries();i++){
    TDirectory* dir = dynamic_cast<TDirectory*>(LHC11hfolderarray->At(i));
    comparearray->Add(dir);
  }
  for(int i = 0; i<LHC10hfolderarray->GetEntries();i++){
    TDirectory* dir = dynamic_cast<TDirectory*>(LHC10hfolderarray->At(i));
    comparearray->Add(dir);
  }  
  for(int i = 0; i<LHC11afolderarray->GetEntries();i++){
    TDirectory* dir = dynamic_cast<TDirectory*>(LHC11afolderarray->At(i));
    comparearray->Add(dir);
  }
  //and now compare them:
  options2.Clear();
  TDirectory * outdiryields = Outfile->mkdir("Peak_shape_and_yield");
  CompareYields(comparearray,outdiryields,options2);

  TDirectory * outdirmethods = Outfile->mkdir("Comparisons");
//   PlotMethods(comparearray,outdirmethods,"all");
  Outfile->Close();
  TFile * Outfiletest 	= TFile::Open("Compare.root","UPDATE");
  TDirectory * outdirplots = Outfiletest->mkdir("pp_PbPb_together");
  PlotCollectedValues(Outfiletest->GetDirectory("Peak_shape_and_yield"),outdirplots,options);
  

  Outfiletest->Close();

  Infile10d->Close();
  Infile11a->Close();
  Infile10h->Close();
  Infile11h->Close();

}

void Sign(const char* options){
  //signs a root folder in the CollectedResults with whatever is in options.
  TFile * rfile = TFile::Open("CollectedResults.root","UPDATE");
  TPaveText * Production;
  Production = new TPaveText(.05,.1,.95,.8);
  Production->AddText(options);
  Production->Write();
  delete Production;
  rfile->Close();
}

void ExamplePlots(const char* options){
  //Prints a couple of example plots in the appropriate folders.
  TString DropBoxFolder = TString("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/illustrations/");
  TString BaseFolder    = TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  TString file11h = TString(Form("%sresults/Physics/leading/LHC11h/CollectedResults.root",BaseFolder.Data()));
  TString file10d = TString(Form("%sresults/Physics/leading/LHC10d/CollectedResults.root",BaseFolder.Data()));
  TString file11a = TString(Form("%sresults/Physics/leading/LHC11a/CollectedResults.root",BaseFolder.Data()));
  TString file10h = TString(Form("%sresults/Physics/leading/LHC10h/CollectedResults.root",BaseFolder.Data()));  
  
  
  Double_t zmax;
  Double_t zmin;
  Double_t sMETA;
  Double_t sMETA2;
  Double_t sMETrigger;
  TString PhiPhi = TString("DPhi_1_DPHI_2")		;
  TString PhiEta = TString("DPhi_1_DEta_12_SameSide")	;
  TString TAdir;
  
  
  TAdir.Append("Correlation_Trigger_8_16/Associated_4_6");
  Plot2dHist(zmax,zmin,file11h.Data(),PhiPhi.Data(),Form("%s/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sDiv_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","Correlation function in PbPb LHC11h with C=0%-10%.");
  Plot2dHist(sMETA,sMETA2,sMETrigger,file11h.Data(),PhiPhi.Data(),Form("%s/BinM(0.00)->(10.00)/iteration1",TAdir.Data()),Form("%sSubs_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","Fully correlated correlation function.");  
  Plot2dHist(file11h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sMETA_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","META1 distribution.",sMETA);
  Plot2dHist(file11h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META2/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sMETA2_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","META2 distribution.",sMETA2);
  Plot2dHist(file11h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/METrigger/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sMETrigger_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","METrigger distribution.",sMETrigger);
  Plot2dHist(zmax,zmin,file11h.Data(),PhiEta.Data(),Form("%s/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sDivss_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","Correlation function in PbPb LHC11h with C=0%-10%.");
  Plot2dHist(sMETA,sMETA2,sMETrigger,file11h.Data(),PhiEta.Data(),Form("%s/BinM(0.00)->(10.00)/iteration1",TAdir.Data()),Form("%sSubsss_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","Fully correlated correlation function.");  
  Plot2dHist(file11h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sMETAss_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","META1 distribution.",sMETA);
  Plot2dHist(file11h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META2/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sMETA2ss_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","META2 distribution.",sMETA2);
  Plot2dHist(file11h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/METrigger/BinM(0.00)->(10.00)/divided",TAdir.Data()),Form("%sMETriggerss_11h_0-10.eps",DropBoxFolder.Data()),"surf3","","","METrigger distribution.",sMETrigger);
    
  
  Plot2dHist(zmax,zmin,file10d.Data(),PhiPhi.Data(),Form("%s/divided",TAdir.Data()),Form("%sDiv_10d.eps",DropBoxFolder.Data()),"surf3","","","Correlation function in pp LHC10d.");
  Plot2dHist(sMETA,sMETA2,sMETrigger,file10d.Data(),PhiPhi.Data(),Form("%s/iteration1",TAdir.Data()),Form("%sSubs_10d.eps",DropBoxFolder.Data()),"surf3","","","Fully correlated correlation function.");  
  Plot2dHist(file10d.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META/divided",TAdir.Data()),Form("%sMETA_10d.eps",DropBoxFolder.Data()),"surf3","","","META1 distribution.",sMETA);
  Plot2dHist(file10d.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META2/divided",TAdir.Data()),Form("%sMETA2_10d.eps",DropBoxFolder.Data()),"surf3","","","META2 distribution.",sMETA2);
  Plot2dHist(file10d.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/METrigger/divided",TAdir.Data()),Form("%sMETrigger_10d.eps",DropBoxFolder.Data()),"surf3","","","METrigger distribution.",sMETrigger);
  Plot2dHist(zmax,zmin,file10d.Data(),PhiEta.Data(),Form("%s/divided",TAdir.Data()),Form("%sDivss_10d.eps",DropBoxFolder.Data()),"surf3","","","Correlation function in pp LHC10d.");
  Plot2dHist(sMETA,sMETA2,sMETrigger,file10d.Data(),PhiEta.Data(),Form("%s/iteration1",TAdir.Data()),Form("%sSubsss_10d.eps",DropBoxFolder.Data()),"surf3","","","Fully correlated correlation function.");  
  Plot2dHist(file10d.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META/divided",TAdir.Data()),Form("%sMETAss_10d.eps",DropBoxFolder.Data()),"surf3","","","META1 distribution.",sMETA);
  Plot2dHist(file10d.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META2/divided",TAdir.Data()),Form("%sMETA2ss_10d.eps",DropBoxFolder.Data()),"surf3","","","META2 distribution.",sMETA2);
  Plot2dHist(file10d.Data(),zmax,zmin,PhiEta.Data(),Form("%s/METrigger/divided",TAdir.Data()),Form("%sMETriggerss_10d.eps",DropBoxFolder.Data()),"surf3","","","METrigger distribution.",sMETrigger);
  
  TAdir.Clear();
  
  
  TString trig; 
  TString ass ; 
  TString centeps;
  TString centfolder;
  TString Folder11h = TString("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/PbPbLHC11hCorrelation/");
  TString Folder10d = TString("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/ppLHC10dCorrelation/");
  TString Folder11a = TString("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/ppLHC11aCorrelation/");
  TString Folder10h = TString("/home/paulbatzing/Dropbox/uio/PhD/analysisnote/figures/PbPbLHC10hCorrelation/");  
  for( int t = 1 ; t <= 3 ; t++){
    if(t==1)trig.Append("4_8");
    if(t==2)trig.Append("8_16");
    if(t==3)trig.Append("16_50");
    for(int a =1; a <= 6 ; a++){
      if(a==1)ass.Append("2_3");
      if(a==2)ass.Append("3_4");
      if(a==3)ass.Append("4_6");
      if(a==4)ass.Append("6_8");      
      if(a==5&&t>1)ass.Append("8_16");
      if(a==6&&t>2)ass.Append("16_51");
      if(ass.CompareTo("")){
	TAdir.Append(Form("Correlation_Trigger_%s/Associated_%s",trig.Data(),ass.Data()));
	//LHC10d
	Plot2dHist(zmax,zmin,file10d.Data(),PhiPhi.Data(),Form("%s/divided",TAdir.Data()),Form("%sDiv/PhiPhi_T%s_A%s_0-10.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(sMETA,sMETA2,sMETrigger,file10d.Data(),PhiPhi.Data(),Form("%s/iteration1",TAdir.Data()),Form("%sITERATION1/PhiPhiT%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file10d.Data(),PhiPhi.Data(),Form("%s/META/divided",TAdir.Data()),Form("%sMETA/PhiPhi_T_%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file10d.Data(),PhiPhi.Data(),Form("%s/META2/divided",TAdir.Data()),Form("%sMETA2/PhiPhi_T_%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file10d.Data(),PhiPhi.Data(),Form("%s/METrigger/divided",TAdir.Data()),Form("%sMETRIGGER/PhiPhi_T_%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");

	Plot2dHist(zmax,zmin,file10d.Data(),PhiEta.Data(),Form("%s/divided",TAdir.Data()),Form("%sDiv/PhiEtass_T%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(sMETA,sMETA2,sMETrigger,file10d.Data(),PhiEta.Data(),Form("%s/iteration1",TAdir.Data()),Form("%sITERATION1/PhiEtassT%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file10d.Data(),PhiEta.Data(),Form("%s/META/divided",TAdir.Data()),Form("%sMETA/PhiEtass_T_%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file10d.Data(),PhiEta.Data(),Form("%s/META2/divided",TAdir.Data()),Form("%sMETA2/PhiEtass_T_%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file10d.Data(),PhiEta.Data(),Form("%s/METrigger/divided",TAdir.Data()),Form("%sMETRIGGER/PhiEtass_T_%s_A%s.eps",Folder10d.Data(),trig.Data(),ass.Data()),"colz","","","");
	//LHC11a
	Plot2dHist(zmax,zmin,file11a.Data(),PhiPhi.Data(),Form("%s/divided",TAdir.Data()),Form("%sDiv/PhiPhi_T%s_A%s_0-10.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(sMETA,sMETA2,sMETrigger,file11a.Data(),PhiPhi.Data(),Form("%s/iteration1",TAdir.Data()),Form("%sITERATION1/PhiPhiT%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file11a.Data(),PhiPhi.Data(),Form("%s/META/divided",TAdir.Data()),Form("%sMETA/PhiPhi_T_%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file11a.Data(),PhiPhi.Data(),Form("%s/META2/divided",TAdir.Data()),Form("%sMETA2/PhiPhi_T_%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file11a.Data(),PhiPhi.Data(),Form("%s/METrigger/divided",TAdir.Data()),Form("%sMETRIGGER/PhiPhi_T_%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");

	Plot2dHist(zmax,zmin,file11a.Data(),PhiEta.Data(),Form("%s/divided",TAdir.Data()),Form("%sDiv/PhiEtass_T%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(sMETA,sMETA2,sMETrigger,file11a.Data(),PhiEta.Data(),Form("%s/iteration1",TAdir.Data()),Form("%sITERATION1/PhiEtassT%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file11a.Data(),PhiEta.Data(),Form("%s/META/divided",TAdir.Data()),Form("%sMETA/PhiEtass_T_%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file11a.Data(),PhiEta.Data(),Form("%s/META2/divided",TAdir.Data()),Form("%sMETA2/PhiEtass_T_%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");
	Plot2dHist(zmax,zmin,file11a.Data(),PhiEta.Data(),Form("%s/METrigger/divided",TAdir.Data()),Form("%sMETRIGGER/PhiEtass_T_%s_A%s.eps",Folder11a.Data(),trig.Data(),ass.Data()),"colz","","","");

	for(int cent = 1;cent<=5;cent++){
	  if(cent==1){centeps.Append("0-10");centfolder.Append("0.00)->(10.00");}
	  if(cent==2){centeps.Append("10-20");centfolder.Append("10.00)->(20.00");}	    
	  if(cent==3){centeps.Append("20-40");centfolder.Append("20.00)->(40.00");} 
	  if(cent==4){centeps.Append("40-60");centfolder.Append("40.00)->(60.00");}
	  if(cent==5){centeps.Append("60-90");centfolder.Append("60.00)->(90.00");}	  
	  
	  cout << "";
	  //LHC11h
	  Plot2dHist(zmax,zmin,file11h.Data(),PhiPhi.Data(),Form("%s/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sDiv/PhiPhi_T%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(sMETA,sMETA2,sMETrigger,file11h.Data(),PhiPhi.Data(),Form("%s/BinM(%s)/iteration1",TAdir.Data(),centfolder.Data()),Form("%sITERATION1/PhiPhiT%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(file11h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA/PhiPhi_T_%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA);
	  Plot2dHist(file11h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META2/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA2/PhiPhi_T_%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA2);
	  Plot2dHist(file11h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/METrigger/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETRIGGER/PhiPhi_T_%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETrigger);
	  Plot2dHist(zmax,zmin,file11h.Data(),PhiEta.Data(),Form("%s/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sDiv/PhiEtass_T%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(sMETA,sMETA2,sMETrigger,file11h.Data(),PhiEta.Data(),Form("%s/BinM(%s)/iteration1",TAdir.Data(),centfolder.Data()),Form("%sITERATION1/PhiEtassT%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(file11h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA/PhiEtass_T_%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA);
	  Plot2dHist(file11h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META2/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA2/PhiEtass_T_%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA2);
	  Plot2dHist(file11h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/METrigger/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETRIGGER/PhiEtass_T_%s_A%s_%s.eps",Folder11h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETrigger);
	  //LHC10h
	  Plot2dHist(zmax,zmin,file10h.Data(),PhiPhi.Data(),Form("%s/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sDiv/PhiPhi_T%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(sMETA,sMETA2,sMETrigger,file10h.Data(),PhiPhi.Data(),Form("%s/BinM(%s)/iteration1",TAdir.Data(),centfolder.Data()),Form("%sITERATION1/PhiPhiT%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(file10h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA/PhiPhi_T_%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA);
	  Plot2dHist(file10h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/META2/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA2/PhiPhi_T_%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA2);
	  Plot2dHist(file10h.Data(),zmax,zmin,PhiPhi.Data(),Form("%s/METrigger/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETRIGGER/PhiPhi_T_%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETrigger);
	  Plot2dHist(zmax,zmin,file10h.Data(),PhiEta.Data(),Form("%s/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sDiv/PhiEtass_T%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(sMETA,sMETA2,sMETrigger,file10h.Data(),PhiEta.Data(),Form("%s/BinM(%s)/iteration1",TAdir.Data(),centfolder.Data()),Form("%sITERATION1/PhiEtassT%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()));
	  Plot2dHist(file10h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA/PhiEtass_T_%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA);
	  Plot2dHist(file10h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/META2/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETA2/PhiEtass_T_%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETA2);
	  Plot2dHist(file10h.Data(),zmax,zmin,PhiEta.Data(),Form("%s/METrigger/BinM(%s)/divided",TAdir.Data(),centfolder.Data()),Form("%sMETRIGGER/PhiEtass_T_%s_A%s_%s.eps",Folder10h.Data(),trig.Data(),ass.Data(),centeps.Data()),"colz","","",Form("%s",TAdir.Data()),sMETrigger);
	


	  centeps.Clear();
	  centfolder.Clear();
	}
	TAdir.Clear();	
      }
      ass.Clear(); 
    }
    trig.Clear();
  }
  
  
  //Event and Track statistics:
  
  
  //Some Trickery to get the correct statistics:
  TString Eventstats10d1 = TString(Form("%s/LHC10d/TreeRunning/leading/t48/a34/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats10d2 = TString(Form("%s/LHC10d/TreeRunning/leading/t816/a23/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats10d3 = TString(Form("%s/LHC10d/TreeRunning/leading/t16/a16/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats10h1 = TString(Form("%s/LHC10h/TreeRunning/leading/t48/a34/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats10h2 = TString(Form("%s/LHC10h/TreeRunning/leading/t816/a23/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats10h3 = TString(Form("%s/LHC10h/TreeRunning/leading/t16/a16/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats11a1 = TString(Form("%s/LHC11a/TreeRunning/leading/t48/a34/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats11a2 = TString(Form("%s/LHC11a/TreeRunning/leading/t816/a23/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats11a3 = TString(Form("%s/LHC11a/TreeRunning/leading/t16/a16/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats11h1 = TString(Form("%s/LHC11h/TreeRunning/leading/t48/a34/AnalysisResults.root",BaseFolder.Data()));
  TString Eventstats11h2 = TString(Form("%s/LHC11h/TreeRunning/leading/t816/a23/AnalysisResults.root",BaseFolder.Data()));  
  TString Eventstats11h3 = TString(Form("%s/LHC11h/TreeRunning/leading/t16/a16/AnalysisResults.root",BaseFolder.Data()));
  
  ETStats(Eventstats10d1,Eventstats10d2,Eventstats10d3,false,Form("%sLHC10d",DropBoxFolder.Data()));
  ETStats(Eventstats10h1,Eventstats10h2,Eventstats10h3,true,Form("%sLHC10h",DropBoxFolder.Data()));
  ETStats(Eventstats11a1,Eventstats11a2,Eventstats11a3,false,Form("%sLHC11a",DropBoxFolder.Data()));
  ETStats(Eventstats11h1,Eventstats11h2,Eventstats11h3,true,Form("%sLHC11h",DropBoxFolder.Data()));
  
  
}