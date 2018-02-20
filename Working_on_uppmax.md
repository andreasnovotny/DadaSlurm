# Working on Uppmax - Instructions in summary.
Andreas Novotny
2018-02-19

https://github.com/andreasnovotny



### 1. Login

Login using the secure shell command:

`ssh *USERNAME*@rackham.uppmax.uu.se`

### 2. Upload and download files in uppmax

*Always run this commands from your lapdop terminal*


secure copy
`scp FROM: TO:`

Download file:
`scp *USERNAME*@rackham.uppmax.uu.se:/path/to/file.ext /local/path`

Uppload file:
`scp /local/path/to/file.ext *USERNAME*@rackham.uppmax.uu.se:/path/`

Directories:
Add the -r option:
`scp -r FROM: TO:`

### 3. The file system at Uppmax
`cd /path/to/dir/`

##### Your home directory

`cd ~`
or
`cd /home/*USERNAME*`

Only visible for the user. Dectory is backed up! Do NOT run any analysis here!

##### Project directories

`cd /proj/*PROJECT-NAME*`

-Directory is backed up! Do NOT run any analysis here - or you will run out of memory! Directory is public! Do NOT store nonpublic data here.

`cd /proj/*PROJECT-NAME*/private`

Dectory is backed up!Do NOT run any analysis here - or you will run out of memory! Directory is private! **USE THIS DIRECTORY FOR BACK UP STORAGE:**
	of raw data files; scripts; final final final results ect.

`cd /proj/*PROJECT-NAME*/nobackup`

Directory is public! Do NOT store nonpublic data here.

`cd /proj/*PROJECT-NAME*/nobackup/private`

Directory is only available for project members
Directory is NOT backed up **USE THIS DIRECORY FOR ANALYSIS:** temporary storage of analysis results ect.


### 4. The module system att uppmax:
Make modules available with:
`module load *MODULE*`

Important examples:

To Access all installed R packages:
`module load R_packages/3.4.3`

To Access all installed bioinformatics software:
`module load bioinfo-tools`
AND!
`module load *MODULE*`



### 5. The SLURM system for launcing big jobs:

Run analyses from Bash file scripts (myscript.sh)
Use the following template for SURM commands:

---
`#!/bin/bash -l`

`#SBATCH -A snic2017-x-xxx` <- the paying project

`#SBATCH -p core` <- core or node partition?

`#SBATCH -n 16` <- Number of cores

`#SBATCH -t 12:00:00` <- time (hh:mm:ss)

`#SBATCH -J dada2_pipeline` <- name the job

`module load python3` <- load needed modules

`FUNCTIONS.ect...`

---

Launch the script by using:
`sbatch /path/to/myscript.sh`

Stop process by:
`scancel *JOBID*`

Monitor jobs by:
`jobinfo -u *USERNAME*`

---
More information can be found at: http://uppmax.uu.se/support/user-guides/
