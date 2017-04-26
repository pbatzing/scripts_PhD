###Shell script that loops over all periods and corrects and gets yields:
ALICE_WORK_DIR=$HOME/alice/build
echo $ALICE_WORK_DIR
eval "`alienv shell-helper`"
alienv load AliPhysics/latest-code-release
clear

for type in alltriggers leading bit4leading bit5leading bit6leading bit56leading
do
  echo "Drawing and collecting all periods using type "$type"."
  for periods in LHC10d LHC11a
  do
    echo "Starting with period "$periods"."
    cd ~/alice/paul/ThreeParticle/correlation3p/$periods/TreeRunning/$type
    #draw everything:
    for trigger in t48 t816 t16
    do
      cd $trigger
      for associated in a23 a34 a46 a68 a816 a16
	do
	  if cd $associated
	    then
	      rm results.root
	      rm CollectedResults.root
	      aliroot -b -q -l runThreePartDC.C'("Draw","416 Same ")'
	      aliroot -b -q -l runThreePartDC.C'("Draw","416 META ")'
	      aliroot -b -q -l runThreePartDC.C'("Draw","416 METrigger ")'
	      #and collect without vertexcut:
	      
	      aliroot -b -q -l runThreePartDC.C'("Collect","pp DivFirst vertexcut=10.0")'
	      cd ..
	  fi
	done
      cd ..
    done
    cd ~/alice/paul/ThreeParticle/correlation3p/results/Physics
    cp 'goodrunsposition_'$periods'_'$type'.txt' goodrunsposition.txt
    aliroot -b -q -l runMergeSet.C
    mv AnalysisResults.root $type/$periods/Backup/CollectedResults.root
    echo moved to $type/$periods'/Backup/CollectedResults.root'

    cd ~/alice/paul/ThreeParticle/correlation3p/$periods/TreeRunning/$type
    for trigger in t48 t816 t16
    do
      cd $trigger
      for associated in a23 a34 a46 a68 a816 a16
	do
	  if cd $associated
	    then
	      rm CollectedResults.root
	      #and collect without vertexcut:
	      aliroot -b -q -l runThreePartDC.C'("Collect","pp DivFirst vertexcut=7.5")'
	      cd ..
	  fi
	done
      cd ..
    done
    cd ~/alice/paul/ThreeParticle/correlation3p/results/Physics
    cp 'goodrunsposition_'$periods'_'$type'.txt' goodrunsposition.txt
    aliroot -b -q -l runMergeSet.C
    mv AnalysisResults.root $type'_vertexcut'/$periods'/Backup/CollectedResults.root'
    echo moved to $type'_vertexcut'/$periods'/Backup/CollectedResults.root'

  done
done


