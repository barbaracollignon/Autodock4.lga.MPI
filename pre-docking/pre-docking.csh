#!/bin/csh -f
##################################################################
# Task-Parallel MPI-based Ligand-Protein Docking                 #
# EXPERIMENTAL PROJECT using Autodock4.lga.MPI                   #
# Oct 23 2010                                                    #
# Author: Barbara Collignon                                      #
##################################################################

#########################
# PRE-DOCKING PROCEDURE # 
#########################

#AutoDockTool (ADT) has to be installed  
#echo "INSTALL AUTODOCKTOOL1.5.2"
#cd mgltools_i86Linux2_1.5.2
#echo "In" `pwd`
#source install.csh >& ../ADTinstall.out
#rehash
#source initMGLtools.csh
#setenv PYTHONPATH `pwd`/MGLToolsPckgs
#cd ../
#echo " "

set dir="MY_PROTEIN"
echo "PROJECT: $dir"
echo " "
echo "PRE-DOCKING PROCEDURE"

#Prepare the Database of Ligands
#The database should have been previously "washed" 

cd $dir
echo "In" `pwd`
echo "I) PDBQT structure files for Flex Lig and Rigid Lig are ALREADY PREPARED"

#Prepare Receptor and Grids
foreach rec ("myprotein")
   set i=1
   echo " "

   #Prepare Receptor
   cd $rec 
   echo "II) Prepare the PDBQT structure file for the receptor: $rec"
   echo "    In `pwd`"
 pythonsh $HOME/PROJECT/mgltools_i86Linux2_1.5.2/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_receptor4.py -r $rec\_nolig.pdb -A hydrogens -o $rec.pdbqt 
   # Add the charge on Zn, Fe... as Autodock4 does not
  cat $rec.pdbqt | sed 's/0.000 Zn/2.000 Zn/' > tmp.dat
  mv tmp.dat $rec.pdbqt
#  cat $rec.pdbqt | sed 's/0.000 K/2.000 K/' > tmp.dat
#  mv tmp.dat $rec.pdbqt
 
   echo " "

   #prepare GPF Input for the grid calculation
   echo "III) Prepare the GPF parameter file for the grid calculation"
   echo "     In `pwd`"
 
   echo "     Also:"
   echo "          The ADT GUI can be used to vizualize/adjust the grid and its center"
 pythonsh $HOME/PROJECT/mgltools_i86Linux2_1.5.2/MGLToolsPckgs/AutoDockTools/Utilities24/prepare_gpf4.py -r $rec.pdbqt -d ../database_rig/ -p gridcenter="20.606,-7.987,-4.594" -p npts="70,70,70" -o $rec.gpf >& $HOME/PROJECT/prepare_gridinput.err
   
   
   echo " "

   #prepare Grids
   echo "     1) Calculate ASCII Grids using Autogrid4"
   rm *map*
  ##cp  /ccs/home/bci/PROJECT/autodock4.2.1/autodocksuite-4.2.1/src/autogrid-4.2.1/autogrid4 .
   ../../autogrid4/autogrid4 -p $rec.gpf -l $rec.glg 


   echo " " 
 
   echo "     2) Prepare Binary Maps"
   echo "        Processing ./make_binary_static.exe -m map.dpf"
   echo "     "
   echo "        map.dpf contains the name of the ASCII files to be converted into binaries"
   echo "        Be aware of the map.dpf file format"
   cp ../map_hdac.dpf .
   ./make_binary_static.exe -m map_hdac.dpf > map.out
#    ./make_binary_malloc_newt.exe -m map.dpf > map.out
   echo " "

   echo "     3) Prepare Receptor Conformers from EN.NMA or MD and corresponding Grids"
   set bool="no"

   #Prepare receptor conformers from EN.NMA or MD
   if("$bool" == "yes") then
      echo "Prepare receptor conformers from EN.NMA or MD"
      tar xvzf MRC*gz

      foreach rec_dup (`ls model*pdb`)
         echo "Prepare Rec duplicate:$rec_dup"
         set dup=`basename $rec_dup .pdb`
         mkdir $dup
         mv $rec_dup $dup
         cd $dup

         prepare_receptor4.py -r $dup.pdb -o $dup.pdbqt
         cat $dup.pdbqt | sed 's/0.000 Zn/2.000 Zn/' > tmp.dat
         mv tmp.dat $dup.pdbqt
         cat $dup.pdbqt | sed 's/0.000 Fe/2.000 Fe/' > tmp.dat
         mv tmp.dat $dup.pdbqt

         #prepare GPF Inputs
         echo "Prepare GPF"
         #cp /home/barbara/AUTODOCK/DUD/testlig.pdbqt .
         #prepare_gpf4.py -r $dup.pdbqt -d  /home/barbara/mol_database/database_rig/ -p gridcenter="85.243 56.199 46.071" -p npts="70,70,70" -o $dup.gpf
         #prepare_gpf4.py -r $dup.pdbqt -d  /home/barbara/mol_database/database_rig -p gridcenter="1.465 1.788 -7.455" -p npts="70,70,70" -o $dup.gpf

         #prepare Grids
         echo "Prepare Grids using Autogrid4"
         ../../../autogrid4/autogrid4 -p $dup.gpf -l $dup.glg

         echo "Prepare Binary Maps"
         # map.dpf contains the name of the ASCII files to be converted into binaries
         # Be aware of the map.dpf file format
        ./make_binary_static.exe -m map.dpf

         cd ../
      end
 
   else if("$bool" == "no") then
   echo "        There are no protein conformers: No more grid calculations"
   endif 

cd ../
end

echo " "

#Prepare Input Parameter Files For Autodock4.lga.MPI
echo "IV) Prepare the 2 DPF input parameter files for Autodock4.lga.MPI"
echo "    In `pwd`"
foreach db (rig flex)

   echo " "
   echo "    In database_$db"
          
   cd database_$db/

   echo "    1) Prepare the common parameter file: dock_par.dpf"
   echo "       The user can use the template file from Collignon et al., 2010"
   echo "       `ls -al ../../dock_par.dpf`"
   echo "       Also:"  
   echo "            A ligand reference can be specified through the keyword: RMSREF"
   echo "            For flexible protein docking it is recommended to generate Multiple Receptor Conformations (MRC)"

   # Make ligand list
   echo "    3) the list of ligands with their specific parameters: list_par.dpf ALREADY PREPARED"


cd $HOME/PROJECT/
echo " "
echo "VI) Copy the input data for the docking procedure"
echo "    Set the path to the actual working directory"
setenv path_to_working_directory /tmp/work/bci/AUTODOCK/MPI/
mkdir /tmp/work/bci/AUTODOCK/
mkdir /tmp/work/bci/AUTODOCK/MPI/
mkdir /tmp/work/bci/AUTODOCK/MPI/$dir/
echo "    Copy Input Data in $path_to_working_directory"
./script_copy.csh -v $path_to_working_directory
echo " "
echo "Now you can execute the parallel docking procedure by using the following command:" 
echo "~/PROJECT> qsub launch_docking.csh"
