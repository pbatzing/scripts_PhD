#include "TString.h"
#include "iostream"

void runan(const char* mode="", const char* run="", Double_t TriggerpTmin = 4.0,Double_t TriggerpTmax = 8.0, Double_t AssociatedpTmin = 2.0, Double_t AssociatedpTmax = 4.0){
  //macro that runs the run-single-task macro with minimal input.
  if(TString("").CompareTo(mode)==0||(TString("").CompareTo(run)==0&&TString("compile").CompareTo(mode)!=0)){cout << "please provide the mode and run number."<<endl;return;}
  if(TString("local").CompareTo(mode)!=0&&TString("test").CompareTo(mode)!=0&&TString("full").CompareTo(mode)!=0&&TString("merge").CompareTo(mode)!=0&&TString("collect").CompareTo(mode)!=0&&TString("compile").CompareTo(mode)){cout << "The supported modes are local, test, full, merge and collect."<<endl;return;}
  TString basedir=TString("/home/paulbatzing/alice/paul/ThreeParticle/correlation3p/");
  gROOT->LoadMacro(Form("%sdefinerunnumbers.C",basedir.Data()));
  cout << run<<endl;
  TString directory = definerunnumbers(run);
  TString runmc =TString("");
  if(TString(run).Contains("mc"))runmc.Append(TString(run).Tokenize("m")->At(0)->GetName());
  cout << directory.Data()<<endl;
  if(directory.CompareTo("undefined")==0&&TString("compile").CompareTo(mode)!=0){cout<<"The run number was not found. Please edit definerunnumbers.C to include it."<<endl;return;}
  gROOT->LoadMacro(Form("%ssetupPbPb.C",basedir.Data()));
  setupPbPb();
  gROOT->LoadMacro(Form("%srun-single-task.C",basedir.Data()));
  gROOT->LoadMacro(Form("%srunTask.C",basedir.Data()));
  
  if(TString("compile").CompareTo(mode)==0){run_single_task(mode,"AddTaskThreePartTracks.C");return;}
  if(directory.BeginsWith("LHC10b")||directory.BeginsWith("LHC10c")||directory.BeginsWith("LHC10e")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode, run, Form("%sAddTaskThreePartTracks.C",basedir.Data()), Form("p%s",run),"000",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sAddTaskThreePartTracks.C",basedir.Data()), Form("p%s",run),"000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode, "AliAOD.root", Form("%sAddTaskThreePartTracks.C",basedir.Data()), Form("test%s",run));
  }
  if(directory.BeginsWith("LHC10h")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode, run, Form("%sAddTaskThreePartBuildTreePbPb.C",basedir.Data()), Form("p%s",run), "000", directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sAddTaskThreePartBuildTreePbPb.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode, "AliAOD.root", Form("%sAddTaskThreePartBuildTreePbPb.C",basedir.Data()), Form("test%s",run));
  }  
  if(directory.BeginsWith("LHC11h")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode, run,  Form("%sAddTaskThreePartBuildTreePbPb.C",basedir.Data()), Form("p%s",run), "000", directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sAddTaskThreePartBuildTreePbPb.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode, "AliAOD.root", Form("%sAddTaskThreePartTracksPbPb.C",basedir.Data()), Form("test%s",run));
  }  
  if(directory.BeginsWith("LHC11a/")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode, run,  Form("%sAddTaskThreePartBuildTree.C",basedir.Data()), Form("p%s",run), "000", directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sAddTaskThreePartBuildTree.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode, "AliAOD.root", Form("%sAddTaskThreePartTracks.C",basedir.Data()), Form("test%s",run));
  }    
  if(directory.BeginsWith("LHC10d")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode, run,  Form("%sAddTaskThreePartBuildTree.C",basedir.Data()), Form("p%s",run), "000", directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sAddTaskThreePartBuildTree.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode, "AliAOD.root", Form("%sAddTaskThreePartTracks.C",basedir.Data()), Form("test%s",run));
  }    
  if(directory.BeginsWith("LHC14j4d")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode, runmc.Data(),  Form("%sMCproductions/LHC14j4d/AddTaskThreePartBuildTree_14j4d.C",basedir.Data()), Form("p%s",run), "000 --mcData", directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC14j4d/AddTaskThreePartBuildTree_14j4d.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)runTask(mode, "dstTree.root", Form("%sMCproductions/LHC14j4d/AddTaskThreePartTrackEfficiencies_14j4d.C",basedir.Data()), Form("test%s",run));
  }    
  if(directory.BeginsWith("LHC11a10abis")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC11a10a/AddTaskThreePartBuildTreePbPb_11a10a.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC11a10a/AddTaskThreePartBuildTreePbPb_11a10a.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0){run_single_task(mode, "dstTree.root", Form("%sMCproductions/LHC11a10a/AddTaskThreePartTrackEfficienciesPbPb_11a10a.C",basedir.Data()));
    }
  }
  if(directory.BeginsWith("LHC10f6a")||directory.BeginsWith("LHC12f1a")||directory.BeginsWith("LHC12f1b")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12af1a/AddTaskThreePartBuildTree_12f1a.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12af1a/AddTaskThreePartBuildTree_12f1a.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)runTask(mode, "dstTree.root", Form("%sMCproductions/LHC12af1a/AddTaskThreePartTrackEfficienciespp_12f1ab.C",basedir.Data()), Form("test%s",run));
    
  }
  if(directory.BeginsWith("LHC12a17a")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17a/AddTaskThreePartTrackEfficienciesPbPb_12a_17a.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17a/AddTaskThreePartTrackEfficienciesPbPb_12a_17a.C",basedir.Data()), Form("p%s",run), "000");
  }
  if(directory.BeginsWith("LHC12a17b")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17b/AddTaskThreePartTrackEfficienciesPbPb_12a_17b.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17b/AddTaskThreePartTrackEfficienciesPbPb_12a_17b.C",basedir.Data()), Form("p%s",run), "000");
  }
  if(directory.BeginsWith("LHC12a17c")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17c/AddTaskThreePartTrackEfficienciesPbPb_12a_17c.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17c/AddTaskThreePartTrackEfficienciesPbPb_12a_17c.C",basedir.Data()), Form("p%s",run), "000");
  }
   if(directory.BeginsWith("LHC12a17d")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17d/AddTaskThreePartTrackEfficienciesPbPb_12a_17d.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17d/AddTaskThreePartTrackEfficienciesPbPb_12a_17d.C",basedir.Data()), Form("p%s",run), "000");
  }
   if(directory.BeginsWith("LHC12a17e")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17e/AddTaskThreePartTrackEfficienciesPbPb_12a_17e.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17e/AddTaskThreePartTrackEfficienciesPbPb_12a_17e.C",basedir.Data()), Form("p%s",run), "000");
  }
   if(directory.BeginsWith("LHC12a17f")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17f/AddTaskThreePartTrackEfficienciesPbPb_12a_17f.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17f/AddTaskThreePartTrackEfficienciesPbPb_12a_17f.C",basedir.Data()), Form("p%s",run), "000");
  }  																	
   if(directory.BeginsWith("LHC12a17g")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17g/AddTaskThreePartBuildTreePbPb_12a_17g.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17g/AddTaskThreePartBuildTreePbPb_12a_17g.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode,"dstTree.root",Form("%sMCproductions/LHC12a17g/AddTaskThreePartTrackEfficienciesPbPb_12a_17g.C",basedir.Data()));
     
  }  
  if(directory.BeginsWith("LHC12a17h")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17h/AddTaskThreePartBuildTreePbPb_12a_17h.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17h/AddTaskThreePartBuildTreePbPb_12a_17h.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode,"dstTree.root",Form("%sMCproductions/LHC12a17h/AddTaskThreePartTrackEfficienciesPbPb_12a_17h.C",basedir.Data()));
  }     
  if(directory.BeginsWith("LHC12a17i")){
    if(TString("test").CompareTo(mode)==0||TString("full").CompareTo(mode)==0)run_single_task(mode,runmc.Data(),Form("%sMCproductions/LHC12a17i/AddTaskThreePartBuildTreePbPb_12a_17i.C",basedir.Data()),Form("p%s",run),"000 --mcData",directory.Data());
    if(TString("merge").CompareTo(mode)==0||TString("collect").CompareTo(mode)==0)run_single_task(mode, directory.Data(), Form("%sMCproductions/LHC12a17i/AddTaskThreePartBuildTreePbPb_12a_17i.C",basedir.Data()), Form("p%s",run), "000");
    if(TString("local").CompareTo(mode)==0)run_single_task(mode,"dstTree.root",Form("%sMCproductions/LHC12a17i/AddTaskThreePartTrackEfficienciesPbPb_12a_17i.C",basedir.Data()));
  }  
}
