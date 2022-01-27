#! /bin/bash

#default path
export PATH=/usr/local/software/anaconda3/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin:/home/skebaas/go/bin

#Set up fsl
#export FSLDIR=/usr/local/software/fsl-5
export FSLDIR=/usr/local/software/fsl-5/
source $FSLDIR/etc/fslconf/fsl.sh
export PATH=$FSLDIR/bin:$PATH

export FSL_DIR=$FSLDIR

export PATH=$PATH:/usr/local/software/ITAlics/afni/bin
export topdir=/usr/local/software/ITAlics

export PATH=$PATH:$topdir/bin
export MATLABPATH=$topdir/matlab/defults:$topdir/matlab/embarc-2.0:$topdir/matlab/PPPI:$topdir/matlab/spm8

