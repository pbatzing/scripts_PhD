###Shell script that loops over all periods and corrects and gets yields:
ALICE_WORK_DIR=$HOME/alice/build
echo $ALICE_WORK_DIR
eval "`alienv shell-helper`"
alienv load AliPhysics/latest-code-release
cd ~/alice/paul/ThreeParticle/correlation3p
rm ThreePartDC_cxx.so ThreePartDC_cxx.d
rm ThreeParticleCorrectionFunctions_cxx.so ThreeParticleCorrectionFunctions_cxx.d
aliroot -b -q -l runThreePartDC.C'("Compile","")'