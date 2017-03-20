const char* libraryDependencies=        
  "libANALYSIS.so "
  "libANALYSISalice.so "
  "libSTEERBase.so "
  "libESD.so "
  "libPWGLFforward2.so "
  "libAOD.so "
  "libPWGCFCorrelationsThreePart.so "
  ;
void setupPbPb()
{
  TString libraries(libraryDependencies);
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include $ALICE_PHYSICS/lib");
  TObjArray* pTokens=libraries.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntriesFast(); i++) {
      if (gSystem->Load(pTokens->At(i)->GetName())==0) {
        cout << "loading " << pTokens->At(i)->GetName() << endl;
      }
    }
    delete pTokens;
  }
  if (gDirectory) gDirectory->Add(new TNamed("analysis_libraries", libraryDependencies));

}
