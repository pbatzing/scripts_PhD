###Shell script that loops over all periods and corrects and gets yields:
ALICE_WORK_DIR=$HOME/alice/build
echo $ALICE_WORK_DIR
eval "`alienv shell-helper`"
alienv load AliPhysics/latest-code-release
clear
for type in alltriggers alltriggers_vertexcut leading #bit4leading bit4leading_vertexcut bit5leading bit5leading_vertexcut bit6leading bit6leading_vertexcut bit56leading bit56leading_vertexcut
do
  echo "Plotting production of type "$type"."
  if cd $type  
    then
      aliroot -b -q -l runThreePartDC.C'("Compare","'$type'")'
      cd ..;
  fi
done

for type in leading_vertexcut
do
  echo "Plotting production of type "$type"."
  if cd $type  
    then
      aliroot -b -q -l runThreePartDC.C'("Compare","'$type'")'
      aliroot -b -q -l runThreePartDC.C'("ExamplePlots")'
      cd ..;
  fi
done