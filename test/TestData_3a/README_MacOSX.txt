If you are running this test suite on the Macintosh OSX operating system, make sure you set the 
BOWTIE_INDEXES environment variable to contain the path to the bowtie indexes. For example, if 
the indexes built by bowtie-build are in the directory /Users/sgivan/Documents/bowtie/indexes, 
set the environment variable like this:

export BOWTIE_INDEXES='/Users/sgivan/Documents/bowtie/indexes'

Also, make sure all the associated software are in your path (bowtie, tophat, cufflinks, 
fastq_quality_trimmer, fastq_quality_filter). When defining the PATH environment variable, 
be sure to fully qualify the path from the root directory:

export PATH="/Users/sgivan/bin:$PATH"

Using the shortcut "~sgivan" in the path definition seems to cause problems.

