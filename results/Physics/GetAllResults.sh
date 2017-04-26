###Shell script that loops over all periods and corrects and gets yields:
ALICE_WORK_DIR=$HOME/alice/build
echo $ALICE_WORK_DIR
eval "`alienv shell-helper`"
alienv load AliPhysics/latest-code-release
clear
for type in alltriggers alltriggers_vertexcut leading leading_vertexcut bit4leading bit4leading_vertexcut bit5leading bit5leading_vertexcut bit6leading bit6leading_vertexcut bit56leading bit56leading_vertexcut
do
  echo "starting on production of type "$type"."
  if cd $type  
    then
      for ppproduction in LHC10d LHC11a 
      do 
	echo "Period "$ppproduction" in production "$type"."
	if cd $ppproduction
	  then
	    #test if the aliroot script is here, link it if not
	    [[ -f runThreePartDC.C ]] || ln -s ../../../../runThreePartDC.C
	    cp Backup/CollectedResults.root .
	    aliroot -b -q -l runThreePartDC.C'("Correct","pp")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","pp")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","pp Eta")'
	    cd ..
	fi
      done
      for PbPbproduction in LHC10h LHC11h
      do
	echo "Period "$PbPbproduction" in production "$type"."
	if cd $PbPbproduction
	  then
	    #test if the aliroot script is here, link it if not
	    [[ -f runThreePartDC.C ]] || ln -s ../../../../runThreePartDC.C
	    cp Backup/CollectedResults.root .
	    aliroot -b -q -l runThreePartDC.C'("Correct","PbPb")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","PbPb")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","PbPb Eta")'
	    cd ..
	fi
      done
      aliroot -b -q -l runThreePartDC.C'("Compare","'$type'")'
      cd ..;
  fi
done
for type in leading_vertexcut
do
  echo "starting on production of type "$type"."
  if cd $type  
    then
      for ppproduction in LHC10d LHC11a 
      do 
	echo "Period "$ppproduction" in production "$type"."
	if cd $ppproduction
	  then
	    #test if the aliroot script is here, link it if not
	    [[ -f runThreePartDC.C ]] || ln -s ../../../../runThreePartDC.C
	    cp Backup/CollectedResults.root .
	    aliroot -b -q -l runThreePartDC.C'("Correct","pp")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","pp")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","pp Eta")'
	    cd ..
	fi
      done
      for PbPbproduction in LHC10h LHC11h
      do
	echo "Period "$PbPbproduction" in production "$type"."
	if cd $PbPbproduction
	  then
	    #test if the aliroot script is here, link it if not
	    [[ -f runThreePartDC.C ]] || ln -s ../../../../runThreePartDC.C
	    cp Backup/CollectedResults.root .
	    aliroot -b -q -l runThreePartDC.C'("Correct","PbPb")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","PbPb")'
	    aliroot -b -q -l runThreePartDC.C'("Yield","PbPb Eta")'
	    cd ..
	fi
      done
      aliroot -b -q -l runThreePartDC.C'("Compare","'$type'")'
      aliroot -b -q -l runThreePartDC.C'("ExamplePlots")'
      cd ..;
  fi
done