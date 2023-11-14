Updated to incorporate fmriprep_20.2.6 to fix bug with mask header https://github.com/nipreps/fmriprep/issues/2507. 

Fmriprep container was build utilizing singularity's --sandbox flag (1). This allows us to make edits to include tedana in the pipeline.

Edited the multi_echo interface to allow for tedana to be used instead of t2smap (2)

Removed skullstripping second pass step done before tedana (3)(a) (3)(b)



references:
(1) singularity build --sandbox fmriprep_20.2.6 docker://nipreps/fmriprep:20.2.6 \n
(2)edited file at /usr/local/miniconda/lib/python3.7/site-packages/fmriprep/interfaces/multiecho.py to change '_cmd' from 't2smap' to 'tedana' and edit outputs in def _list_outputs(self) \n
(3)\n
  \t(a)To correct for this issue, we needed to find a way to remove the skullstrip_second_pass for only this one section.  We want to interfere with the code as little as possible.  I achieved this by creating a custom workflow in /usr/local/miniconda/lib/python3.7/site-packages/fmriprep/workflows/bold/custom_utils.py` script can be found in attachments as 'custom_utils.py' \n
  \t(b) After creating the custom_utils.py it had to be added into the main bold workflow /usr/local/miniconda/lib/python3.7/site-packages/fmriprep/workflows/bold/base.py at line 414

