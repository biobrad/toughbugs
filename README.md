# toughbugs
A Tool for running popular Antimicrobial Resistance Software

toughbugs is a shell script that works in ubuntu shell, not 100% sure of posix compatibility so results may vary.
I have included a conda .yml file to create the environment.

*** warning: The script still needs a fair bit of work, it doesn't support arguments at this stage and has operator required options. ***

Until I add arguments to the script, it will need to be manually updated to provide the location of the databases. 

To run, it requires the fasta and fastq.gz files for your sequences to be in the folder where you run the .sh file from. 
It is a requirement that the conda environment is also active.

With multiple databases and tools available to examine short read sequences for antimicrobial determinants, it became necessary to build a pipeline that utilises all of these platforms into a single, easy to execute pipeline.

Toughbugs requires a few databases to be downloaded and updated, namely ariba, Amrfinder and CARD. after this has been done, all that is needed is a fasta file, fastq.gz files and the pipeline will run from the environment in which the files are located.

