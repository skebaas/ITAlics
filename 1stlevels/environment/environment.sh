#! /bin/bash

#default path
export PATH=/usr/local/software/anaconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/skebaas/go/bin

#Set up fsl
#export FSLDIR=/usr/local/software/fsl-5
export FSLDIR=/usr/local/software/fsl-5/
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=$FSLDIR/bin:$PATH

export FSL_DIR=$FSLDIR

export PATH=$PATH:/usr/local/software/gold/afni/bin
export topdir=/usr/local/software/gold

export PATH=$PATH:$topdir/bin

export MATLABPATH=/usr/local/software/gold/matlab/defults:/usr/local/software/gold/matlab/embarc-1.0:/usr/local/software/gold/matlab/embarc-2.0:/usr/local/software/gold/matlab/ExploreDTI:/usr/local/software/gold/matlab/nifti:/usr/local/software/gold/matlab/PPPI:/usr/local/software/gold/matlab/scripts:/usr/local/software/gold/matlab/spm8:/usr/local/software/gold/matlab/tools:/usr/local/MATLAB/R2015a/toolbox:/usr/local/MATLAB/R2015a/toolbox/images/images
