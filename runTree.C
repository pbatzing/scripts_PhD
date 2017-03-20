#include "TString.h"
#include "TSystem.h"
#include "THashList.h"
#include <iostream>
#include <fstream>

void runTree(const char* mode="", const char* dataset="", const char* directory = "Tree", double nevensperjob = 250000, int jobn = -1){
  //macro that runs the run-single-task macro with minimal input.
  if(TString("").CompareTo(mode)==0||(TString("").CompareTo(dataset)==0&&TString("compile").CompareTo(mode)!=0)){cout << "please provide the mode and run number."<<endl;return;}
  if(TString("local").CompareTo(mode)!=0&&TString("goodruns").CompareTo(mode)!=0&&TString("test").CompareTo(mode)!=0&&TString("full").CompareTo(mode)!=0&&TString("merge").CompareTo(mode)!=0&&TString("collect").CompareTo(mode)!=0&&TString("compile").CompareTo(mode)){cout << "The supported modes are local, test, full, merge and collect."<<endl;return;}
  TString basedir=TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  TString Directory = TString(directory);
  gROOT->LoadMacro(Form("%srunTask.C",basedir.Data()));
  
  ofstream outputFile;
  ofstream outputFile2;
  if(TString("goodruns").CompareTo(mode)==0){
    outputFile.open("goodrunsposition.txt");
    outputFile2.open("goodruns.txt");
    gSystem->GetFromPipe("ln -s /home/paulbatzing/alice/paul/macros/runMergeSet.C");
    gSystem->GetFromPipe("ln -s /home/paulbatzing/alice/paul/macros/MergeSet.C");    
  }
  TString runningbasedir =TString(gSystem->pwd());
  gSystem->MakeDirectory(Form("%s/runs",runningbasedir.Data()));
  
  Int_t maxnevents = 0;

  
  if(TString(dataset).BeginsWith("LHC10h")){
    Int_t runnumbersP10h[93] = {139510, 139507, 139505, 139503, 139465, 139438, 139437, 139360, 139329, 139328, 139314, 139310, 139309, 139173, 139107, 139105, 139038, 139037, 139036, 139029, 139028, 138872, 138871, 138870, 138837, 138732, 138730, 138666, 138662, 138653, 138652, 138638, 138624, 138621, 138583, 138582, 138579, 138578, 138534, 138469, 138442, 138439, 138438, 138396, 138364, 138275, 138225, 138201, 138197, 138192, 138190, 137848, 137844, 137752, 137751, 137724, 137722, 137718, 137704, 137693, 137692, 137691, 137686, 137685, 137639, 137638, 137608, 137595, 137549, 137546, 137544, 137541, 137539, 137531, 137530, 137443, 137441, 137440, 137439, 137434, 137432, 137431, 137430, 137243, 137236, 137235, 137232, 137231, 137230, 137162, 137161};    
    for(int i=0; i<91;i++){
      int nevents = 0;
      filerun = TFile::Open(Form("%sLHC10h/%i/AnalysisResults.root",basedir.Data(),runnumbersP10h[i]),"READ");
      THashList* list = dynamic_cast<THashList*>( filerun->GetDirectory("ThreePartBuildTreePbPb_1_100")->Get("ThreePartBuildTreePbPb_1_100_0_0Coutput1"));
      TH1* events = dynamic_cast<TH1*>( list->FindObject("Eventafterselection"));
      nevents = events->GetEntries();
      filerun->Close();
      if(nevents>maxnevents)maxnevents=nevents;
      cout << "Run "<<runnumbersP10h[i]<< " has "<<nevents << " events, needing "<< nevents/nevensperjob<<" jobs."<<endl;
    }
    int njobs =  maxnevents/nevensperjob;
    if((maxnevents/nevensperjob - njobs)>0.2)njobs+=1;
    cout << "The maximum number of events in a run is: "<< maxnevents<<" needing " << njobs << " jobs"<< endl;    
            
    for(int i=0;i<njobs;i++){
      int  firsteventnr = i*nevensperjob;
      if(i ==(njobs-1)){nevensperjob= maxnevents-firsteventnr;}
      if((jobn>0)&&(i!=(jobn-1)))continue;
      TString rundir = Form("%s/runs/%02i",runningbasedir.Data(), i+1);
      gSystem->MakeDirectory(rundir.Data());
      gSystem->ChangeDirectory(rundir.Data());
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0){
	gSystem->GetFromPipe("find ./ -type f -delete").Data();
	gSystem->GetFromPipe("find ./ -type l -delete").Data();
	gSystem->GetFromPipe("ln -s ../../TaskConfig.C").Data();
	gSystem->GetFromPipe("ln -s ../../grid-config.C").Data();
      }
      if(TString("goodruns").CompareTo(mode)==0){
	outputFile  << Form("%02i",i+1) <<" "<<rundir.Data()<<"/AnalysisResults.root"<<endl;
	outputFile2 << Form("%02i",i+1) <<endl;	

      }
      
      const char* addtaskmacro = Form("%sAddTaskThreePartTracksPbPbTree.C",basedir.Data());
      const char* taskname = Form("%sp%s_%i",dataset,directory,i+1);
      const char* griddir = Form("%s/p%s/%i",dataset,Directory.Data(),i+1);
//       if local, do nothing, execute:
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)runTask(mode,"000000001",addtaskmacro, taskname, "000 --maxInputFiles=1",nevensperjob,firsteventnr,griddir);
      if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)runTask(mode,griddir,addtaskmacro,taskname,"000");
    }
    gSystem->ChangeDirectory(runningbasedir.Data());
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0){
      outputFile.close();
      outputFile2.close();
    }
    if(TString("goodruns").CompareTo(mode)==0){
      outputFile.close();
      outputFile2.close();
    }
    
  }  
  if(TString(dataset).BeginsWith("LHC11h")){
    Int_t runnumbersP11h[110] = {170593, 170572, 170388, 1703872 , 1703871, 170315, 170313, 170312, 170311, 170309, 170308, 170306, 170270, 170269, 170268, 170230, 170228, 170207, 170204, 170203, 170193, 170163, 170159, 170155, 170091, 170089, 170088, 170085, 170084, 170083, 170081, 1700401, 1700402, 170027, 169965, 169923, 169859, 169858, 169855, 169846, 169838, 169837, 169835, 169591, 169590, 169588, 169587, 169586, 169557, 169555, 169554, 169553, 169550, 169515, 169512, 169506, 169504, 169498, 169475, 169420, 169419, 169418, 169417, 169415, 169411, 169238, 169167, 169160, 169156, 169148, 169145, 169144, 169138, 169099, 169094, 169091, 169045, 169044, 169040, 169035, 168992, 168988, 168826, 168777, 168514, 168512, 168511, 168467, 168464, 168460, 168458, 168362, 168361, 168342, 168341, 168325, 168322, 168311, 168310, 168115, 168108, 168107, 168105, 168076, 168069, 167988, 167987, 167985, 167920, 167915};
    if(TString("local").CompareTo(mode)==0){
      const char* addtaskmacro = Form("%sAddTaskThreePartTracksPbPbTree.C",basedir.Data());
      const char* taskname = Form("%sp%s",dataset,directory);
      cout << taskname<<endl;
      runTask(mode,"dstTree.root",addtaskmacro,taskname,"",100000,1);
      return;
    }

    for(int i=0; i<110;i++){
      int nevents = 0;
      TFile *filerun;
      if(runnumbersP11h[i]!=1703872&&runnumbersP11h[i]!=1703871&&runnumbersP11h[i]!=1700401&&runnumbersP11h[i]!=1700402){
	filerun = TFile::Open(Form("%sLHC11h/%i/AnalysisResults.root",basedir.Data(),runnumbersP11h[i]),"READ");
	THashList* list = dynamic_cast<THashList*>( filerun->GetDirectory("ThreePartBuildTreePbPb_1_100")->Get("ThreePartBuildTreePbPb_1_100_0_0Coutput1"));
	TH1* events = dynamic_cast<TH1*>( list->FindObject("Eventafterselection"));
	nevents = events->GetEntries();
	filerun->Close();
      }
      else{
	if(runnumbersP11h[i]==1703871)filerun = TFile::Open(Form("%sLHC11h/170387/merge1/AnalysisResults.root",basedir.Data()),"READ");
	if(runnumbersP11h[i]==1703872)filerun = TFile::Open(Form("%sLHC11h/170387/merge2/AnalysisResults.root",basedir.Data()),"READ");
	if(runnumbersP11h[i]!=1700401)filerun = TFile::Open(Form("%sLHC11h/170040/merge1/AnalysisResults.root",basedir.Data()),"READ");
	if(runnumbersP11h[i]!=1700402)filerun = TFile::Open(Form("%sLHC11h/170040/merge2/AnalysisResults.root",basedir.Data()),"READ");
	THashList* list = dynamic_cast<THashList*>( filerun->GetDirectory("ThreePartBuildTreePbPb_1_100")->Get("ThreePartBuildTreePbPb_1_100_0_0Coutput1"));
	TH1* events = dynamic_cast<TH1*>( list->FindObject("Eventafterselection"));
	nevents = events->GetEntries();
	filerun->Close();
      }
      if(nevents>maxnevents)maxnevents=nevents;
      cout << "Run "<<runnumbersP11h[i]<< " has "<<nevents << " events, needing "<< nevents/nevensperjob<<" jobs."<<endl;
    }
    int njobs =  maxnevents/nevensperjob;
    if((maxnevents/nevensperjob - njobs)>0.2)njobs+=1;
    cout << "The maximum number of events in a run is: "<< maxnevents<<" needing " << njobs << " jobs"<< endl;

    
    for(int i=0;i<njobs;i++){
      int  firsteventnr = i*nevensperjob;
      if(i ==(njobs-1)){nevensperjob= maxnevents-firsteventnr;}
      if((jobn>0)&&(i!=(jobn-1)))continue;
      TString rundir = Form("%s/runs/%02i",runningbasedir.Data(), i+1);
      gSystem->MakeDirectory(rundir.Data());
      gSystem->ChangeDirectory(rundir.Data());
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0){
	gSystem->GetFromPipe("find ./ -type f -delete").Data();
	gSystem->GetFromPipe("find ./ -type l -delete").Data();
	gSystem->GetFromPipe("ln -s ../../TaskConfig.C").Data();
	gSystem->GetFromPipe("ln -s ../../grid-config.C").Data();
      }
      if(TString("goodruns").CompareTo(mode)==0){
	outputFile  << Form("%02i",i+1) <<" "<<rundir.Data()<<"/AnalysisResults.root"<<endl;
	outputFile2 << Form("%02i",i+1) <<endl;	

      }
      
      const char* addtaskmacro = Form("%sAddTaskThreePartTracksPbPbTree.C",basedir.Data());
      const char* taskname = Form("%sp%s_%i",dataset,directory,i+1);
      const char* griddir = Form("%s/p%s/%i",dataset,Directory.Data(),i+1);
//       if local, do nothing, execute:
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)runTask(mode,"000000001",addtaskmacro, taskname, "000 --maxInputFiles=1",nevensperjob,firsteventnr,griddir);
      if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)runTask(mode,griddir,addtaskmacro,taskname,"000");
    }
    gSystem->ChangeDirectory(runningbasedir.Data());
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0){
      outputFile.close();
      outputFile2.close();
    }
    if(TString("goodruns").CompareTo(mode)==0){
      outputFile.close();
      outputFile2.close();
    }
  }
  if(TString(dataset).BeginsWith("LHC11a")){
    Int_t runnumbersP11a[16] = {146746, 146747, 146748, 146801 , 146802, 146803, 146804, 146805, 146806, 146807, 146817, 146824, 146856, 146858, 146859, 146860};
      cout << nevensperjob<<endl;
    nevensperjob *= 10;
      cout << nevensperjob<<endl;
    for(int i=0; i<16;i++){
      int nevents = 0;
      TFile *filerun;
      filerun = TFile::Open(Form("%sLHC11a/%i/AnalysisResults.root",basedir.Data(),runnumbersP11a[i]),"READ");
      THashList* list = dynamic_cast<THashList*>( filerun->GetDirectory("ThreePartBuildTree_1_100")->Get("ThreePartBuildTree_1_100_0_0Coutput1"));
      TH1* events = dynamic_cast<TH1*>( list->FindObject("Eventafterselection"));
      nevents = events->GetEntries();
      filerun->Close();
      if(nevents>maxnevents)maxnevents=nevents;
      cout << "Run "<<runnumbersP11a[i]<< " has "<<nevents << " events, needing "<< nevents/nevensperjob<<" jobs."<<endl;
      
    }
    int njobs =  maxnevents/nevensperjob;
    if((maxnevents/nevensperjob - njobs)>0.2)njobs+=1;
    cout << "The maximum number of events in a run is: "<< maxnevents<<" needing " << njobs << " jobs"<< endl;
    
    for(int i=0;i<njobs;i++){
      int  firsteventnr = i*nevensperjob;
      if(i ==(njobs-1)){nevensperjob= maxnevents-firsteventnr;}
      if((jobn>0)&&(i!=(jobn-1)))continue;
      TString rundir = Form("%s/runs/%02i",runningbasedir.Data(), i+1);
      gSystem->MakeDirectory(rundir.Data());
      gSystem->ChangeDirectory(rundir.Data());
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0){
	gSystem->GetFromPipe("find ./ -type f -delete").Data();
	gSystem->GetFromPipe("find ./ -type l -delete").Data();
	gSystem->GetFromPipe("ln -s ../../TaskConfig.C").Data();
	gSystem->GetFromPipe("ln -s ../../grid-config.C").Data();
      }
      if(TString("goodruns").CompareTo(mode)==0){
	outputFile  << Form("%02i",i+1) <<" "<<rundir.Data()<<"/AnalysisResults.root"<<endl;
	outputFile2 << Form("%02i",i+1) <<endl;	
      }
      
      const char* addtaskmacro = Form("%sAddTaskThreePartTracksTree.C",basedir.Data());
      const char* taskname = Form("%sp%s_%i",dataset,directory,i+1);
      const char* griddir = Form("%s/p%s/%i",dataset,Directory.Data(),i+1);
  //       if local, do nothing, execute:
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)runTask(mode,"000000001",addtaskmacro, taskname, "000 --maxInputFiles=1",nevensperjob,firsteventnr,griddir);
      if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)runTask(mode,griddir,addtaskmacro,taskname,"000");
    }
    gSystem->ChangeDirectory(runningbasedir.Data());
    if(TString("goodruns").CompareTo(mode)==0){
      outputFile.close();
      outputFile2.close();
    }
  }
  if(TString(dataset).BeginsWith("LHC10d")){
    Int_t runnumbersP11a[62] = {126432, 126425, 126424, 126422, 126409, 126408, 126407, 126406, 126405, 126404, 126403, 126359, 126352, 126351, 126350, 126285, 126284, 126283, 126168, 126167, 126160, 126158, 126097, 
      126090, 126088, 126082, 126081, 126078, 126073, 126008, 126007, 126004, 125855, 125851, 125850, 125849, 125848, 125847, 125844, 125843, 125842, 125633, 125632, 125630, 125628, 125296, 125295, 125186, 125156, 125140, 
      125139, 125134, 125133, 125101, 125100, 125097, 125085, 125083, 125023, 124751, 122375, 122374};
    nevensperjob = 3000000;
    
    if(TString("local").CompareTo(mode)==0){
      const char* addtaskmacro = Form("%sAddTaskThreePartTracksTree.C",basedir.Data());
      const char* taskname = Form("%sp%s",dataset,directory);
      cout << taskname<<endl;
      runTask(mode,"dstTree.root",addtaskmacro,taskname,"",100000,1);
      return;
    }
    
    for(int i=0; i<62;i++){
      int nevents = 0;
      TFile *filerun;
      filerun = TFile::Open(Form("%sLHC10d/%i/AnalysisResults.root",basedir.Data(),runnumbersP11a[i]),"READ");
      THashList* list = dynamic_cast<THashList*>( filerun->GetDirectory("ThreePartBuildTree_1_100")->Get("ThreePartBuildTree_1_100_0_0Coutput1"));
      TH1* events = dynamic_cast<TH1*>( list->FindObject("Eventafterselection"));
      nevents = events->GetEntries();
      filerun->Close();
      if(nevents>maxnevents)maxnevents=nevents;
      cout << "Run "<<runnumbersP11a[i]<< " has "<<nevents << " events, needing "<< nevents/nevensperjob<<" jobs."<<endl;
      
    }
    int njobs =  maxnevents/nevensperjob;
    if((maxnevents/nevensperjob - njobs)>0.2)njobs+=1;
    cout << "The maximum number of events in a run is: "<< maxnevents<<" needing " << njobs << " jobs"<< endl;
    for(int i=0;i<njobs;i++){
      int  firsteventnr = i*nevensperjob;
      if(i ==(njobs-1)){nevensperjob= maxnevents-firsteventnr;}
      if((jobn>0)&&(i!=(jobn-1)))continue;
      TString rundir = Form("%s/runs/%02i",runningbasedir.Data(), i+1);
      gSystem->MakeDirectory(rundir.Data());
      gSystem->ChangeDirectory(rundir.Data());
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0){
	gSystem->GetFromPipe("find ./ -type f -delete").Data();
	gSystem->GetFromPipe("find ./ -type l -delete").Data();
	gSystem->GetFromPipe("ln -s ../../TaskConfig.C").Data();
	gSystem->GetFromPipe("ln -s ../../grid-config.C").Data();
      }
      if(TString("goodruns").CompareTo(mode)==0){
	outputFile  << Form("%02i",i+1) <<" "<<rundir.Data()<<"/AnalysisResults.root"<<endl;
	outputFile2 << Form("%02i",i+1) <<endl;	

      }
      
      const char* addtaskmacro = Form("%sAddTaskThreePartTracksTree.C",basedir.Data());
      const char* taskname = Form("%sp%s_%i",dataset,directory,i+1);
      const char* griddir = Form("%s/p%s/%i",dataset,Directory.Data(),i+1);
  //       if local, do nothing, execute:
      if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)runTask(mode,"000000001",addtaskmacro, taskname, "000 --maxInputFiles=1",nevensperjob,firsteventnr,griddir);
      if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)runTask(mode,griddir,addtaskmacro,taskname,"000");
    }
    gSystem->ChangeDirectory(runningbasedir.Data());
    if(TString("goodruns").CompareTo(mode)==0){
      outputFile.close();
      outputFile2.close();
    }
  }
}
