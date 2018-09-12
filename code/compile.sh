#!/bin/bash

# Simple shortcut for compiling my programs.
# Very very basic script.
# Accept as argument the name of the program
# (that is the name of the source without .f90

rootlib=/hpc/hub_oudenaarden/aalemany/bcn
root=/hpc/hub_oudenaarden/aalemany/bcn

# now a list of shorthand names for the libraries
nrtype=$rootlib/lib/nrtype.f90
nrutil=$rootlib/lib/nrutil.f90
ran_state=$rootlib/lib/ran_state.f90
nr=$rootlib/lib/nr.f90
amFileUtil=$rootlib/lib/amFileUtil.f90
amFluctTh=$rootlib/lib/amFluctTh.f90
amPhysConst=$rootlib/lib/amPhysConst.f90
amProbDistrib=$rootlib/lib/amProbDistrib.f90
HairpinParameters=$rootlib/lib/HairpinParameters.f90
aaElasticModel=$rootlib/lib/aaElasticModel.f90
aaBootstrapMethods=$rootlib/lib/aaBootstrapMethods.f90
aaFDCanalysis=$rootlib/lib/aaFDCanalysis.f90

if [ $# -ne 1 ] 
then
	echo "Need one and only one argument"
	exit
else
  file=${1%.f*}
# Separa trajectòries
    if [[ $file = "DataSplitPulling" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "DataSplitPaused" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "DataSplitFastSlow" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "DataSplitMichelle" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "DataSplitTemperaturaPulling" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $file.f90 -o $root/bin/$file
# Troba velocitat en nm/s
    elif [[ $file = "findRate" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $nr $file.f90 -o $root/bin/$file
    elif [[ $file = "findRatev2" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $nr $file.f90 -o $root/bin/$file
    elif [[ $file = "findFeedback" ]]
      then gfortran $nrtype  $nrutil $amFileUtil $nr $file.f90 -o $root/bin/$file
# Alienació de trajectòries.
    elif [[ $file = "fitAlign" ]]
	then gfortran $nrtype  $nrutil $nr $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "normalizeEventsFit" ]]
        then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "dumbAlign_302" ]]
	then gfortran $nrtype $nrutil $nr $amFileUtil $file.f90 -o $root/bin/$file 
    elif [[ $file = "dumbAlignCycles" ]]
	then gfortran $nrtype $nrutil $nr $amFileUtil $file.f90 -o $root/bin/$file 
    elif [[ $file = "DumbAlignFJ" ]]
	then gfortran $nrtype $nrutil $nr $amFileUtil $file.f90 -o $root/bin/$file 
    elif [[ $file = "dumbAlignEventsNorm" ]]
	then gfortran $nrtype $nrutil $nr $amFileUtil $file.f90 -o $root/bin/$file 
    elif [[ $file = "normalizeEventsDumb" ]]
        then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "normalizeEventsDumbFJ" ]]
        then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
#     elif [[ $file = "distingirUF" ]]
#         then gfortran $nrtype  $nrutil $nr $amFileUtil $file.f90 -o $root/bin/$file
# Detecció de forces de ruptura
    elif [[ $file = "findRuptureForce" ]]
	then gfortran $nrtype  $nrutil $nr $file.f90 -o $root/bin/$file
# Anàlisi de forces de ruptura
    elif [[ $file = "ForceAnalysis" ]]
        then gfortran $nrtype $nrutil $nr $amFileUtil $amProbDistrib $file.f90 -o $root/bin/$file
    elif [[ $file = "ForceAnalysisWLCFJC" ]]
        then gfortran $nrtype $nrutil $nr $amFileUtil $amProbDistrib  $file.f90 -o $root/bin/$file
    elif [[ $file = "findTransitionk" ]]
      then gfortran $nrtype $nrutil $amFileUtil $ran_state $nr $aaElasticModel $file.f90 -o $root/bin/$file
    elif [[ $file = "ssResponse" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
# Càlcul de treball
    elif [[ $file = "computeWork" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file 
    elif [[ $file = "computeFilteredWork" ]]
      then gfortran $nrtype $nrutil $amFileUtil $ran_state $aaBootstrapMethods $file.f90 -o $root/bin/$file 
    elif [[ $file = "DeltaWork" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file 
# Avaluació de paràmetres elàstics de les handles
    elif [[ $file = "keff" ]]
      then gfortran $nrtype $nrutil $nr $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "keff2" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "keff3" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "keffAveraging" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "kbkxEvaluation" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "kssDNAEvaluation" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "chi2_WLC" ]]
      then gfortran $nrtype $nrutil $amFileUtil $aaElasticModel $file.f90 -o $root/bin/$file
# Anàlisi del treball
    elif [[ $file = "MyWorkAnalysisWLCFJC" ]]
      then gfortran $nrtype $nrutil $nr $amFileUtil $amProbDistrib $amFluctTh $file.f90 -o $root/bin/$file
    elif [[ $file = "BarnaseWorkAnalysisWLCFJC" ]]
      then gfortran $nrtype $nrutil $nr $amFileUtil $amProbDistrib $amFluctTh $file.f90 -o $root/bin/$file
    elif [[ $file = "HandlesWork" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file     
    elif [[ $file = "Stretching" ]]
	then gfortran $nrtype $nrutil $nr $aaElasticModel $file.f90 -o $root/bin/$file
    elif [[ $file = "wdis_m_2estats" ]]
	then gfortran $file.f -o $root/bin/$file
# Analisi feedback
    elif [[ $file = "feedbackAnalysis" ]]
	then gfortran $nrtype $nrutil $amFileUtil $aaFDCanalysis $file.f90 -o $root/bin/$file
    elif [[ $file = "ForcesWorksFeedback" ]]
	then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    elif [[ $file = "ForcesWorksFeedback2" ]]
	then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
# Altres
    elif [[ $file = "ModelEmpiricTIB" ]]
	then gfortran $nrtype $nrutil $file.f90 -o $root/bin/$file
    elif [[ $file = "ModelManning" ]]
	then gfortran $nrtype $nrutil $file.f90 -o $root/bin/$file
    elif [[ $file = "varf_bin" ]]
      then gfortran $nrtype $nrutil $amFileUtil $ran_state $nr $file.f90 -o $root/bin/$file 
    elif [[ $file = "avTrajPow2" ]]
      then gfortran $nrtype $nrutil $amFileUtil $file.f90 -o $root/bin/$file
    else 
	gfortran $nrtype  $nrutil $file.f90 -o $root/bin/$file
    fi
fi

