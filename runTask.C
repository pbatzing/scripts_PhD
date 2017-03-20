

#if defined(__CINT__) && !defined(__MAKECINT__)
///////////////////////////////////////////////////////////////////////////////////////////////////
//
// environment
//
int macroVerbosity=0;
const char* defaultAnalysisName="myanalysis";
const char* libraryDependencies=
  "libANALYSIS.so "
  "libANALYSISalice.so "
  "libSTEERBase.so "
  "libESD.so "
  "libPWGLFforward2.so "
  "libAOD.so "
  "libPWGCFCorrelationsThreePart.so "
  ;

TString BuildJDLName(const char* analysisName);
AliAnalysisManager* LoadAnalysisManager(const char* filename);
AliInputEventHandler* LoadCustomInputHandler(const char* name);
TChain* CreateChain(const char* filename);
void ErrorConfigurationFile(const char* fileName);

void runTask(	     const char* mode,
		     const char* input,
		     const char* tasknames=NULL,
		     const char* analysisName=defaultAnalysisName,
		     const char* arguments="",
		     int nevents=-1,
		     int firstevent =-1,
		     const char* odirin=""
		     )
{
  const char* gridDataDir=NULL;
  const char* dataPattern=NULL;
  const char* friendDataPattern=NULL;
  bool bRunLocal = strcmp(mode,"local")==0;
  int mergeMode=0;
  TString odir = TString(odirin);
  if ((strcmp(mode, "merge")==0 && (mergeMode=1)>0) ||
      (strcmp(mode, "collect")==0 && (mergeMode=2)>0)) {
    mode="terminate";
    odir=input;
    input="0";
  }
  //
  // argument settings
  //
  bool bRunAnalysis=true;
  bool bDefaultOutput=true;
  bool bMergeOnGrid=mergeMode==2?false:true; // default true in all cases but 'collect'
  bool mcData=false;
  int nTestFiles=2;
  int nMaxInputFiles=100;
  TString strArguments(arguments);
  TString mergeDirs;
  TString strCustomInputHandler;
  TObjArray* tokens=strArguments.Tokenize(" ");
  
   if (tokens) {
    for (int iToken=0; iToken<tokens->GetEntriesFast(); iToken++) {
      TObject* token=tokens->At(iToken);
      if (!token) continue;
      TString arg=token->GetName();
      const char* key=NULL;

      if (arg.CompareTo("--help")==0 ||
          arg.CompareTo("-h")==0 ||
          arg.CompareTo("help")==0 ||
          arg.CompareTo("options")==0) {
	// printing help when called without arguments
	run_single_task();
	return;
      }
      key="--mcData";
      if (arg.CompareTo(key)==0) {
	// this is an argument to the macro, don't propagate it further to tasks
	// switch indicates that the input data is mc data, the run numbers have
	// a different format in real data
	// NOTE: not to be confused with option 'mc' which is propagated to tasks
	// and switches processing and output modes inside tasks
	mcData=true;
	continue;
      }
      key="--merge=";
      if (arg.BeginsWith(key)) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	arg.ReplaceAll(key, "");
	if (arg.CompareTo("local")==0) {
	  // download all files and merge locally
	  bMergeOnGrid=false;
	} else if (arg.CompareTo("collect")==0) {
	  // download the output of merging on Grid
	  // macro must have been called in mode "terminate" with option
	  // --merge=grid and the correct grid working directory
	  bMergeOnGrid=false;
	} else if (arg.CompareTo("grid")==0) {
	  // merge output on grid,  the correct grid working directory
	  // must be provided
	  bMergeOnGrid=true;
	}
	continue;
      }
      key="--nTestFiles=";
      if (arg.BeginsWith(key)) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	arg.ReplaceAll(key, "");
	nTestFiles=arg.Atoi();
	continue;
      }
      key="--noDefaultOutput";
      if (arg.CompareTo(key)==0) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	bDefaultOutput=false;
	continue;
      }
      key="--stopBeforeRunning";
      if (arg.CompareTo(key)==0) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	bRunAnalysis=false;
	continue;
      }
      key="--maxInputFiles=";
      if (arg.BeginsWith(key)) {
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
	arg.ReplaceAll(key, "");
	nMaxInputFiles=arg.Atoi();
	continue;
      }
      key="--InputHandler=";
      if (arg.BeginsWith(key)) {
	arg.ReplaceAll(key, "");
	strCustomInputHandler=arg;
      }
      if (!arg.BeginsWith("-") && mergeMode>0) {
	// treat as subdirectories in the remote grid dir
	mergeDirs+=" "; mergeDirs+=arg;
	// this is an argument to the macro, don't propagate it further to tasks
	strArguments.ReplaceAll(arg, "");
      }
    }
    delete tokens;
  }
  
  //
  // load task classes and find and load all dependencies
  //
  //Test if the include path is already added, and add if not:
  if( !TString(gSystem->GetIncludePath()).Contains("-I$ROOTSYS/include"))gSystem->AddIncludePath("-I$ROOTSYS/include");
  if( !TString(gSystem->GetIncludePath()).Contains("-I$ALICE_PHYSICS/include"))gSystem->AddIncludePath("-I$ALICE_PHYSICS/include");
  if( !TString(gSystem->GetIncludePath()).Contains("-I$ALICE_ROOT/include"))gSystem->AddIncludePath("-I$ALICE_ROOT/include");
  TString libraries=libraryDependencies;
  TObjArray* pTokens=libraries.Tokenize(" ");
  TString buf;
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast(); i++) {
      TString library=pTokens->At(i)->GetName();
      if (!library.EndsWith(".so")) {
	cerr << "libraries need to have ending '.so' in order to be correctly processed by alien plugin, please correct library name '" << library << "'" << endl;
      }
      buf = pTokens->At(i)->GetName();
      buf.ReplaceAll(".so", "");
      if (gSystem->Load(buf)==0) {
	cout << "loading " << pTokens->At(i)->GetName() << endl;
      }
    }
    delete pTokens;
  }
  TString taskNames=tasknames;
  TString addTaskMacros="";
  TString delimiter(" ");
  TStringToken taskNameTokens(taskNames, delimiter);
  while (taskNameTokens.NextToken()) {
    TString taskSource(taskNameTokens);
    TString taskHeader(taskNameTokens);
    bool bIsAddTask=false;
    if (taskSource.EndsWith(".C")) {
      // suppose that's an 'AddTask' macro
      taskHeader="";
      bIsAddTask=true;
      addTaskMacros+=" "; addTaskMacros+=taskSource;
    }
  }
  //
  // grid defaults
  //
  const char* gridConfigFile="grid-config.C";
  TString strGridConfigFile=gridConfigFile;
  if (gSystem->AccessPathName(strGridConfigFile)!=0) {
    strGridConfigFile.Prepend("/");
    strGridConfigFile.Prepend(gSystem->Getenv("HOME"));
    if (gSystem->AccessPathName(strGridConfigFile)!=0) {
      if (!bRunLocal && mergeMode==0) {
	ErrorConfigurationFile(gridConfigFile);
	return;
      }
      strGridConfigFile="";
    }
  }

  // load the grid configuration file if not merging and not running locally
  if (strGridConfigFile.IsNull()==0 && !bRunLocal) {
    cout << "loading grid configuration from file '" << strGridConfigFile << "':" << endl;
    gROOT->LoadMacro(strGridConfigFile);
    cout << " alienAliPhysicsVersion   =" << alienAliPhysicsVersion <<endl;
    cout << " defaultGridDataDir       =" << defaultGridDataDir  << endl;
    cout << " defaultDataPattern       =" << defaultDataPattern  << endl;
    cout << " defaultFriendDataPattern =" << defaultFriendDataPattern  << endl;

    if (gridDataDir==NULL) gridDataDir=defaultGridDataDir;
    if (dataPattern==NULL) dataPattern=defaultDataPattern;
    if (friendDataPattern==NULL) friendDataPattern=defaultFriendDataPattern;
  } else if (bRunLocal) {
    if (dataPattern==NULL) {
      // thats a very crude logic, I guess it can fail in some special cases
      TString strin=input;
      if (strin.Contains("AOD"))
	dataPattern="AOD";
      else if (strin.Contains("dstTree"))
	dataPattern="dstTree";
    }
  }

  if(!bRunLocal) {
    // Connect to AliEn
    TGrid::Connect("alien://");
  }
  
  //
  // make the analysis manager
  //
  AliAnalysisManager *pManager=NULL;
  pManager=new AliAnalysisManager("AnalysisManager");
  if (!pManager) {
    cerr << "failed to create AnalysisManager" << endl;
    return;
  }
  AliInputEventHandler *pInputHandler = NULL;
  TString strDataPattern(dataPattern);
  if (!strCustomInputHandler.IsNull()) {
    pInputHandler=LoadCustomInputHandler(strCustomInputHandler);
  }
  else if (strDataPattern.Contains("AOD")) pInputHandler=new AliAODInputHandler;
  else if (strDataPattern.Contains("dstTree")) pInputHandler=new AliFilteredEventInputHandler;
  else if (mergeMode>0) pInputHandler=new AliInputEventHandler; // a default handler, not used in the end
  else {
    cerr << "can not determine input type from data pattern '" << strDataPattern << "'" << endl;
    return;
  }
  if (!pInputHandler) {
    cerr << "failed to created input handler" << endl;
    return;
  }
  pManager->SetInputEventHandler(pInputHandler);  
  pManager->SetNSysInfo(1000);

  TString ofile=Form("%s.root", analysisName);

  //
  // init for local or GRID analysis
  //
  AliAnalysisAlien *alienHandler = NULL; // for grid analysis
  TChain *chain=NULL; // for local analysis
  TString strInput=input;
  if (bRunLocal) {
    ///////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // local analysis
    //
    if (strInput.BeginsWith("alien://")) {
      // file on Grid -> connect to AliEn
      TGrid::Connect("alien://");
    }
    if (strInput.EndsWith(".root") && gSystem->AccessPathName(strInput)==0) {
      // open a local file
      chain=CreateChain(strInput.Data());
    }
    if (chain) {
      // nothing to do here, just forward to the original
      // functionality if the chain was not created already
    } else 
      if(strInput.EndsWith("AliAOD.root")){
      // single local AOD file
      chain = new TChain("aodTree"); 
      chain->Add(strInput);
    } else if(strInput.EndsWith("AOD")){
      // fetch aod tree from the setup macro
      const char* aodTreeName="aodTree";
      if (gDirectory!=NULL) {
	TObject* chainObject=gDirectory->FindObject(aodTreeName);
	if (chainObject) {
	  chain=dynamic_cast<TChain*>(chainObject);
	}
      }
      if (!chain) {
	::Error("run_single_task", Form("failed to fetch aod tree object from setup; the chain with name '%s' has to be created before calling this macro", aodTreeName));
	return;
      }
    } else if(strInput.EndsWith("dstTree.root")){
      // single local AOD file
      chain = new TChain("DstTree"); 
      chain->Add(strInput);
    } else {
      ::Error("run_single_task", Form("invalid input"));
      return;
    }
  } else {
    //
    // grid analysis
    //
    bool bSetRun=true;
    if (!strInput.IsDigit()) {
      // support for external macros specifying the the runs to be
      // analyzed
      // the input is expected to be an external plugin with name 'input'
      // and all run numbers being set
      TObject* pObj=gDirectory->FindObject(input);
      if (pObj) alienHandler=dynamic_cast<AliAnalysisAlien*>(pObj);
      if (!alienHandler) {
	::Error("run_single_task", Form("can not find plugin of name '%s', please setup alien handler with name and run numbers before calling this macro", input));
	return;
      }
      bSetRun=false;
    } else {
      alienHandler=new AliAnalysisAlien();
    }
    if (!alienHandler) {
      ::Error("run_single_task", Form("failed to create alien handler"));
      return;
    }

    // do not check for copying to grid (CLOSE_SE)
    alienHandler->SetCheckCopy(kFALSE);

    // Set the run mode (can be "full", "test", "offline", "submit" or "terminate")
    alienHandler->SetRunMode(mode);

    // number of files in test mode configurable via argument '--nTestFiles='
    if(mode=="test") alienHandler->SetNtestFiles(nTestFiles);
  
    // check the versions available on alien with the command 'packages'
    if (mergeMode==0) {
    alienHandler->SetAliPhysicsVersion(alienAliPhysicsVersion);
    }

    // using only default output
    // the alien plugin automatically recognizes all output files associated to output
    // containers, all files are treated in the standard output and added to the
    // root-archieve.root, which also seems to be needed for merging on Grid
    // see further down for using non-default output
    alienHandler->SetDefaultOutputs(bDefaultOutput);

    // data alien directory
    alienHandler->SetGridDataDir(gridDataDir);
  
    // Set data search pattern
    alienHandler->SetDataPattern(dataPattern);
    alienHandler->SetFriendChainName(friendDataPattern);

 
    if (bSetRun) {
      // only set if input is a run number
      if (!mcData && !strInput.BeginsWith("000"))
	alienHandler->SetRunPrefix("000");   // real data

      alienHandler->AddRunNumber(input);
    }

    if (mergeMode>0) {
      // the merge and collect modes have been added to simplify merging on grid
      // the treatment of arguments are a bit different in order to reduce list
      // of required arguments
      TString delimiter(" ");
      if (mergeDirs.IsNull()) mergeDirs="000";
      TStringToken dir(mergeDirs, delimiter);
      while (dir.NextToken()) {
	alienHandler->AddDataFile(dir.Data());
      }
      // use the specified directory names rather than a counter
      alienHandler->SetOutputToRunNo();
    }

    // define working and output directories
    TDatime dt;
    if(odir.IsNull())
      odir=(Form("gridwork/%04d-%02d-%02d_%02d-%02d", dt.GetYear(), dt.GetMonth(), dt.GetDay(), dt.GetHour(), dt.GetMinute()));
//     cout << odir << endl;
    alienHandler->SetGridWorkingDir(odir); // relative to $HOME
    alienHandler->SetGridOutputDir("output");   // relative to working dir

    // Optionally set a name for the generated analysis macro (default MyAnalysis.C)
    TString macroName; macroName.Form("run_%s.C",analysisName); macroName.ReplaceAll("-","_");
    alienHandler->SetAnalysisMacro(macroName);
    //alienHandler->SetExecutable("comparison.sh");
    alienHandler->SetExecutable(Form("run_%s.sh",analysisName));
    alienHandler->SetSplitMaxInputFileNumber(nMaxInputFiles);
  
    // Optionally resubmit threshold.
    alienHandler->SetMasterResubmitThreshold(90); // in %
    alienHandler->SetTTL(86400);// in sec
  
    // Optionally set input format (default xml-single)
    alienHandler->SetInputFormat("xml-single");
 
    // Optionally modify the name of the generated JDL (default analysis.jdl)
    alienHandler->SetJDLName(BuildJDLName(analysisName).Data());
 
    // Optionally modify job price (default 1)
    alienHandler->SetPrice(1);
  
    // Optionally modify split mode (default 'se')
    alienHandler->SetSplitMode("se");
  
    // configure merging on grid,
    // argument '--merge=collect' sets 'false' for fetching the merged output
    alienHandler->SetMergeViaJDL(bMergeOnGrid); 

    alienHandler->SetOneStageMerging(kFALSE);
    alienHandler->SetMaxMergeStages(4);
    alienHandler->SetMaxMergeFiles(10);
  }

  // Connect plugin to the analysis manager
  if (alienHandler) {
    if(TString(defaultDataPattern).Contains("dstTree.root"))alienHandler->SetTreeName("DstTree");
    alienHandler->SetDropToShell(false);
    pManager->SetGridHandler(alienHandler);
  }
  
  //
  // create task from the name, create output container, connect slots
  //
  TObjArray* taskMacroTokens=addTaskMacros.Tokenize(" ");
  if (taskMacroTokens) {
    for (int iTaskMacroToken=0; iTaskMacroToken<taskMacroTokens->GetEntriesFast(); iTaskMacroToken++) {
      TString taskSource= taskMacroTokens->At(iTaskMacroToken)->GetName();

      taskSource+="+g";
      TString configuration;
      if(!strArguments.Contains("file=")) configuration+=Form(" file=%s",ofile.Data()); 
      if(!strArguments.Contains("name=")) configuration+=Form(" name=%s",analysisName); 
      configuration+=" "; configuration+=strArguments.Data();
      if (gDirectory) gDirectory->Add(new TNamed("run_single_task_configuration", configuration.Data()));
      cout << "here "<< taskMacroTokens->At(iTaskMacroToken)->GetName()<<endl;
      gROOT->Macro(taskMacroTokens->At(iTaskMacroToken)->GetName());
      cout << "there"<<endl;
    }
    delete taskMacroTokens;
  }

  if (!bDefaultOutput) {
    // fetch all output files from the output containers
    TString ofiles;
    TIter nextcontainer(pManager->GetContainers());
    TObject* objContainer=NULL;
    while ((objContainer=nextcontainer())!=NULL) {
      AliAnalysisDataContainer* container=dynamic_cast<AliAnalysisDataContainer*>(objContainer);
      if (!container) continue;
      ofiles+=container->GetFileName();
      ofiles+=" ";
    }

    alienHandler->SetOutputFiles(ofiles);
    // Optionally define the files to be archived.
    alienHandler->SetOutputArchive("log_archive.zip:stdout,stderr");
  }
  //
  // run
  //

  if (!pManager->InitAnalysis()) {
    cerr << "failed to initialize analysis" << endl;
    return;
  }
  if (nevents<0) nevents=1000000000;
  if (firstevent<0) firstevent =1;
  pManager->PrintStatus();
  if (!bRunAnalysis) return;
  if (bRunLocal) {
    pManager->StartAnalysis("local", chain, nevents,firstevent);
  } else {
    pManager->StartAnalysis("grid", nevents,firstevent);
  }
  
}

TString BuildJDLName(const char* analysisName)
{
  TString jdlname = TString(Form("run_%s.jdl",(analysisName!=NULL?analysisName:"analysis")));
  return jdlname;
}

TChain* CreateChain(const char* filename)
{
  // create input TChain object with tree name derived from input file
  TChain* chain=NULL;
  TFile* file=TFile::Open(filename);
  if (!file || file->IsZombie()) {
    return NULL;
  }

  TList* keys=file->GetListOfKeys();
  if (!keys || keys->GetEntries()==0) {
    cerr << "can not find any keys in file " << filename << endl;
    return NULL;
  }

  TObject* pObj=NULL;
  TObject* pKey=NULL;
  TIter nextkey(keys);
  while (pKey=nextkey()) {
    file->GetObject(pKey->GetName(), pObj);
    if (!pObj || pObj->IsA()!=TTree::Class()) continue;
    TChain* chain = new TChain(pObj->GetName()); 
    chain->Add(filename);
    break;
  }
  file->Close();
  delete file;
  if (chain) {
    cout << "created chain " << chain->GetName() << endl;
  }
  return chain;
}

void ErrorConfigurationFile(const char* fileName) {
  cout << endl;
  cout << "/// -------------------------------------------------------------------" << endl;
  cout << "/// Warning: can not find configuration file '" << fileName << "'" << endl;
  cout << "/// please create a configuration file in either local or HOME directory, or in" << endl;
  cout << "/// specified location. Below is an example, fill in your preferred defaults." << endl;
  cout << "/// -------------------------------------------------------------------" << endl;
  cout << endl;
  cout << "const char* alienAliPhysicsVersion=\"v5-01-Rev-29\";" << endl;
  cout << "const char* defaultGridDataDir=\"/alice/data/2011/LHC11f\";" << endl;
  cout << "const char* defaultDataPattern=\"*ESDs.root\";" << endl;
  cout << "const char* defaultFriendDataPattern=\"\";" << endl;
  cout << "{} // note this empty body";
  cout << endl;
}
#endif
