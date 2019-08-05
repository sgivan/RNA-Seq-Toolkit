#!/usr/bin/env python

# This is a script to take a yaml file and run the RNA-Seq Toolkit.
# It will probably stop at the creation of a DESeq2 Rscript that can
# be submitted to the cluster.
#
import os, sys, argparse, shutil
import re, yaml, subprocess, time

# start with command line options
argparser = argparse.ArgumentParser(description="Parse tab-delimited file")

# infile is the yaml file
argparser.add_argument("--infile", type=str,  help="yaml file to parse", default="file")
argparser.add_argument("--verbose", action="store_true", help="verbose messages to terminal")

args = argparser.parse_args()

if args.verbose: print 'input file: %(filename)s' % { "filename": args.infile }

f_yaml = open(args.infile, 'r') # open infile

# create yaml object
#config = yaml.load(f_yaml)
config = yaml.load(f_yaml, Loader=yaml.FullLoader)
curdir=os.getcwd()
if args.verbose: print "current working directory: '%(workdir)s'" % { "workdir": curdir }

os.environ['PATH']=config['rst_path'] + "/bin" + ":" + os.environ['PATH']

clength=len(config['input']['control'])
elength=len(config['input']['experimental'])

if args.verbose:
    print "dumping config"
    print yaml.dump(config)
#
# define some functions
#
# instead of composing this twice
# create a function to build the file structure
#
def create_file_struct(sample_number, fileset, config, curdir):

    number_of_seq_files=len(fileset)
    if number_of_seq_files == 2:
        if args.verbose: print "\tworking with paired-end data"

        if config['paired']:
            if args.verbose: print "\tthis confirms configuration file"
        else:
            if args.verbose: print "\tthis conflicts with configuration file\n\tplease revise\n\texiting now"
            sys.exit(1)

        try:
            os.chdir(config['working_datadir'])
        except:
            print "can't chdir to '%(dirname)s'" % { 'dirname': config['working_datadir'] }
            print "exiting now"
            sys.exit(2)

        wdir=os.getcwd()

        try:
            os.mkdir("Sample_" + str(sample_number))
        except:
            print "can't create directory 'Sample_%(dirdigit)s' in '%(workdirname)s'" % { "dirdigit": sample_number, "workdirname": config['working_datadir'] }
            print "exiting now"
            sys.exit(3)

        if args.verbose: print "\tdirectory 'Sample_%(dirdigit)s' created in '%(workdirname)s'" % { "dirdigit": sample_number, "workdirname": config['working_datadir'] }

        try:
            os.chdir("Sample_" + str(sample_number))
        except:
            print OSError

        if args.verbose: print "\tnow in directory " + os.getcwd()


        lcnt=0
        for sfile in fileset:
            target=os.path.join(curdir,config['original_datadir'],sfile)
            if args.verbose: print "\tcreating symlink called '%(linkname)s' pointing to '%(linktarget)s'" % { "linkname": sfile, "linktarget": target }
            try:
                os.symlink(target,sfile)
            except OSError as e:
                print e.errno
                print e.filename
                print e.strerror

            lcnt += 1
            if args.verbose: print "\tcreating symlink called 'set%(linkname)s.fq' pointing to '%(linktarget)s'" % { "linkname": lcnt, "linktarget": sfile }
            try:
                os.symlink(sfile, "set" + str(lcnt) + ".fq")
            except OSError as e:
                print e.errno
                print e.filename
                print e.strerr

        os.chdir(curdir)
        wdir=os.getcwd()
        print "\tnow in directory" + wdir

    else:
        if args.verbose: print "\tworking with non-paired-end data"

def monitor_cluster_jobs(jobs):
#   not perfect, but it works
    if args.verbose:
        print "cluster jobs to monitor:"
        print jobs

    wait=1
    while wait:
        wait=0

        print "looping through job id's"
        for jobid in jobs:
            if args.verbose:
                print "checking job #%(id)s" % { 'id': jobid }

            rtn=0
            try:
                rtn=subprocess.check_output("tracejob -a -l -m " + jobid, shell=True)
                if args.verbose:
                    print "rtn value: '%(val)s'" % { 'val': rtn }

            except subprocess.CalledProcessError as e:
                print "can't call tracejob with jobid %(jobid)i: %(ecode)i" % { "jobid": jobid, "ecode": e.returncode }
                sys.exit(15)
  
            if re.search(r'Exit_status',rtn):
                if args.verbose: print "job " + jobid + " finished"
                continue
            else:
                wait=1
#                if args.verbose: print "waiting for job " + jobid
                print "waiting for job " + jobid
                time.sleep(5)
                break

def check_compressed_files():

    if args.verbose:
        print "checking to see if %(dirname)s/%(morepath)s is compressed" % {
                'dirname': config['working_alignment_dir'],
                'morepath': "/Sample_1/read_1",
                }
    rtn = ''
    try:
        rtn = subprocess.check_output("file -L " + curdir + "/" + config['working_alignment_dir'] + "/Sample_1/read_1", shell=True)
    except subprocess.CalledProcessError as e:
        print "can't call file to test if file is compressed"
        sys.exit(16)

    if args.verbose:
        print "rtn: " + rtn

    compressed = 0
    if re.search(r'gzip', rtn):
        if args.verbose: "input files are compressed"
        compressed = 1
    else:
        if args.verbose: "input files are not compressed"
        compressed = -1

    return compressed
    
            
#
# end of functions
#

#
# Begin running RNA-Seq Toolkit pipeline
# based on configuration in yaml file
#
jobs=[]
filemap={}
filemap['control']=[]
filemap['experimental']=[]

if 'setup_files' in config.keys() and config['setup_files'] != False:
    # create directory to contain
    # copies of input files
    # so we keep original files and only work with copies
    print "\nsetting up file structure for input files\n"
    if not os.access(config['working_datadir'], os.F_OK):
        if args.verbose: print "creating new directory to place renamed files: %(newdir)s." % { "newdir": config['working_datadir'] }
        if os.access(os.path.split(config['working_datadir'])[0], os.W_OK):
            os.mkdir(config['working_datadir'])

        else:
            print "can't create %(newdir)s" % { "newdir": config['working_datadir'] }
            sys.exit(4)
    else:
        print "%(newdir)s already exists. Will not overwrite -- please rename or move the diretory." % { "newdir": config['working_datadir'] }
        print "Exiting now."
        sys.exit(5)

#
# end of working directory section
#
    sample_number=0

#    clength=len(config['input']['control'])
    if args.verbose: print "control replicates: %(clen)i" % { "clen": clength }

    for i in config['input']['control']:
        sample_number += 1
        if args.verbose: print "\n\tcontrol - %(repname)s will be given symbolic name 'Sample_%(sint)s'" % { "repname": i, "sint": sample_number }
        create_file_struct(sample_number, config['input']['control'][i], config, curdir)
        filemap['control'].append([i, 'Sample_' + str(sample_number)])

    if 'experimental' in config['input']:
#        elength=len(config['input']['experimental'])

        if args.verbose: print "experimental data sets: %(length)i" % { "length": elength }
#        efiles={}
        for i in config['input']['experimental']:
#            efiles[i]=[]
#
#   This cut corresponds to the sample replicate. There can be any number of sample replicates. which will have either a single file (non-PE) or a pair of files (Paired End)
#
#            if args.verbose: print "sample replicates in set %(eset)s: %(filenames)s" % { "eset": i, "filenames": config['input']['experimental'][i] }
#            number_of_reps=len(config['input']['experimental'][i])
#            if args.verbose: print "number of replicates: %(numseqs)i." % { "numseqs": number_of_reps }

#            for j in config['input']['experimental'][i]:

             sample_number += 1
#                if args.verbose: print "\n\t%(setname)s - %(repname)s will be given symbolic name 'Sample_%(sint)s'" % { "setname": i, "repname": j, "sint": sample_number }
             if args.verbose: print "\n\t%(setname)s - %(repname)s will be given symbolic name 'Sample_%(sint)s'" % { "setname": "experimental", "repname": i, "sint": sample_number }

#                create_file_struct(sample_number, config['input']['experimental'][i][j], config, curdir)
             create_file_struct(sample_number, config['input']['experimental'][i], config, curdir)
#                efiles[i].append([j, 'Sample_' + str(sample_number)])

#             filemap['experimental']=efiles
             filemap['experimental'].append([i, "Sample_" + str(sample_number)])
#             if args.verbose: print "experimental data sets: %(length)i" % { "length": elength }

#
#   print out yaml file containing map of original files to standardized files
#
    mapfile = file('filemap.yaml', 'w')
    yaml.dump(filemap, mapfile)
    mapfile.close()

# End of creating yaml file map

    if args.verbose: print str(sample_number) + ' samples'
    print "Input file setup finished."

# 
# do preprocessing
# typically, this is removing reads that match rRNA genes and/or mt genome
#
if 'preprocess' in config.keys() and config['preprocess'] != False:
    print '\n\n pre-process the input data.'

#
# create working directory
# warn if it's already present
# #
    try:
        os.mkdir(config['working_alignment_dir'])
    except OSError as e:
        print "can't create directory '%(dirname)s'." % { "dirname": config['working_alignment_dir'] }
        print e.errno
        print e

#
# cd into working directory
# exit if error
#
    try:
        os.chdir(config['working_alignment_dir'])
    except OSError as e:
        print "can't chdir into '%(dirname)s'." % { "dirname": config['working_alignment_dir'] }
        print e.errno
        print e
        sys.exit(50)

# create symlinks to index directories
#
    print "creating symlnks to preprocess and alignment index files in " + config['working_alignment_dir']

# Create symlink to directory containing
# bowtie1 index files of filter sequences.
# Filter sequences are the ones we don't
# want in the input files
# First, check if it already exists. If it does, exit.
    if os.access('index.preprocess', os.F_OK):
        print "Will not overwrite current 'index.preprocess' symlink.\nPlease remove it."
        sys.exit(6)

# Create symlink
# fail if error
    try:
        os.symlink(config['filter_datadir'], 'index.preprocess')
    except OSError as e:
        print "can't create index.preprocess symlink pointing to '%(dirname)s.'" % { 'dirname': config['filter_datadir'] }
        print e.errno
        print e.filename
        sys.exit(60)

#
# Check if symlink already exists to alignment indices.
# These are the STAR indices of the reference transcriptome.
# Exit if it already exists.
    if os.access('index.align', os.F_OK):
        print "Will now overwrite current 'index.align' symlink.\nPlease remove it."
        sys.exit(7)

#
# Create symlink to STAR index files.
# Exit if error.
    try:
        os.symlink(config['index_datadir'], 'index.align')
    except OSError as e:
        print "can't create index.align symlink pointing to '%(dirname)s.'" % { 'dirname': config['index_datadir'] }
        print e.errno
        print e.filename
        sys.exit(70)
#        print e.strerr

# Check if a current symlink exists named "index".
# Exit if it already exists.
    if os.access('index', os.F_OK):
        print "Will now overwrite current 'index' symlink.\nPlease remove it."
        sys.exit(8)

#
# Create the "index" symlink.
# At this point, point the index symlink to the filter indices.
# Exit if error.
    try:
        os.symlink('index.preprocess', 'index')
    except OSError as e:
        print "can't create index symlink pointing to index.preproces"
        print e.errno
        print e.filename
        sys.exit(80)
#        print e.strerr

# Now that directory structure is created,
# move on to running pre-processing routines.
    if args.verbose: print "preprocess and align symlinks created in " + curdir

    setup_script=os.path.join(config['rst_path'], 'bin', 'setup.sh')
    if args.verbose: print "calling setup script '%(scriptname)s'." % { "scriptname": setup_script }
    try:
#        out=subprocess.check_call(setup_script, shell=True)
        dpth = os.path.join("..", config['working_datadir'])
        print "will symlink Sample_* directories in %(dirpath)s." % { "dirpath": dpth }
        out=subprocess.check_call(setup_script + " " + dpth, shell=True)
    except subprocess.CalledProcessError as e:
        print "call to symlink failed"
        print "error code: %(ecode)i" % { "ecode": e.returncode }
        sys.exit(9)

    if args.verbose: print "running RST preprocessing routines"

    rst_script=os.path.join(config['rst_path'], 'bin', 'RNAseq_process_data.sh')

    out=""
    scriptargs = " --preprocess_only --notrim --submit --threads " + str(config['procs']) + " --queue " + str(config['jobQ'])
    if config['seq_compressed'] == True:
        scriptargs += " --gzip"

    scriptargs += " Sample_*"

    if args.verbose:
        print "scriptargs: '%(argstring)s'" % { 'argstring': scriptargs }

    try:
        out=subprocess.check_output(rst_script + " " + scriptargs, shell=True)
    except subprocess.CalledProcessError as e:
        print "call to %(rst)s failed with error code %(ecode)i" % { "rst": rst_script, "ecode": e.returncode }
        sys.exit(10)
    if args.verbose:
        print "script output: '%(output)s'" % { 'output': out }

    os.chdir(curdir)
#   the preprocessing scripts will decompress the fastq files,
#   so set the seq_compressed config setting to False
    config['seq_compressed'] = False

# Monitor cluster jobs and wait for them to finish
# before exiting ending this section of the pipeline
    slurmjobs = []
    for line in str.splitlines(out):
        match = re.search("master\.cm\.cluster", line)
        if match:
            words = line.split(".")
            id = words[0]
            slurmjobs.append(id)

    monitor_cluster_jobs(slurmjobs)


if 'align' in config.keys() and config['align'] != False:
    print "\n\n setting up alignment directory"

#
# We may or may not be working with a pre-made directory structure.
# So, we must check if the expected directories, symlinks and files already exist.
# In general, don't die if they already exist (as we do in the setup section, above)
#
    try:
        os.mkdir(config['working_alignment_dir'])
    except OSError as e:
        print "can't create directory '%(dirname)s'." % { "dirname": config['working_alignment_dir'] }
        print e.errno
        print e
#        sys.exit()

    try:
        os.chdir(config['working_alignment_dir'])
    except OSError as e:
        print "can't chdir into '%(dirname)s'." % { "dirname": config['working_alignment_dir'] }
        print e.errno
        print e
        sys.exit(90)


    print "creating symlnks to alignment index files in " + config['working_alignment_dir']

    if os.access('index.align', os.F_OK):
        print "Will not overwrite current 'index.align' symlink."
    else:

        try:
            os.symlink(config['index_datadir'], 'index.align')
        except OSError as e:
            print "can't create index.align symlink pointing to '%(dirname)s.'" % { 'dirname': config['index_datadir'] }
            print e.errno
            print e.filename

    if os.access('index', os.F_OK):
#        print "Will not overwrite current 'index' symlink in the %(dirname)s directory.\nPlease remove it." % { 'dirname': config['working_alignment_dir'] }
        print "An 'index' link currently exists in %(dirname)s.\nIt will be removed and replaced with a new one." % { 'dirname': config['working_alignment_dir'] }

        try:
            os.remove('index')
        except OSError as e:
            print "can't remove index symlink in directory %(dirname)s.\nPlease remove it manually." % { 'dirname': config['working_alignment_dir'] }
            print e.errno
            print e
            sys.exit(91)

    try:
        os.symlink('index.align', 'index')
    except OSError as e:
        print "can't create index symlink pointing to index.align"
        print e.errno
        print e.filename
        sys.exit(92)

    if args.verbose: print "align symlink created in " + curdir

    setup_script=os.path.join(config['rst_path'], 'bin', 'setup.sh')
    if args.verbose: print "calling setup script '%(scriptname)s'." % { "scriptname": setup_script }
    try:
        dpth = os.path.join("..", config['working_datadir'])
        print "will symlink Sample_* directories in %(dirpath)s." % { "dirpath": dpth }
        out=subprocess.check_call(setup_script + " " + dpth, shell=True)
    except subprocess.CalledProcessError as e:
        print "call to symlink failed"
        print "error code: %(ecode)i" % { "ecode": e.returncode }
#        sys.exit(9)

    os.chdir(curdir)

    if args.verbose: print "\n\naligning data to reference genome sequence"

    os.chdir(config['working_alignment_dir'])

    rst_script=os.path.join(config['rst_path'], 'bin', 'RNAseq_process_data.sh')

    if args.verbose: print "checking to see if input files are compressed with gzip"
    chk = check_compressed_files()
    if chk > 0:
        rst_script = rst_script + " --gzip"
        
    try:
        subprocess.check_call(
                rst_script + " --partial --submit --threads " + str(config['procs']) + \
                " --queue " + str(config['jobQ']) + \
                " --min_intron_length " + str(config['min_intron_length']) + " --max_intron_length " + str(config['max_intron_length']) + \
                " Sample_*",
                shell=True
                )
    except subprocess.CalledProcessError as e:
        print "can't call %(scriptname)s: %(ecode)i" % { "scriptname": rst_script, "ecode": e.returncode }
        sys.exit(13)

#
#   I need to monitor jobs here and not continue until last job is finished
#
    print('''

    The sequence alignment jobs have been submitted to the cluster. You can monitor their progress using
    the qstat command:
    qstat -u $USER

    ''')

    time.sleep(5)

if 'diff_expression' in config.keys() and config['diff_expression'] != False:

    os.chdir(curdir)

    filemapfile = file('filemap.yaml', 'r')
    filemap=yaml.load(filemapfile, Loader=yaml.FullLoader)

    files = []

    keys = ['control','experimental']
    for stype in keys:
#        print "stype %s: %s" % (stype, filemap.get(stype))
        for sname in filemap[stype]:
#            print "sname: '%s'" % sname[1]
            files.append(sname[1])

    print "files: '%s'" % files

    try:
        os.mkdir(config['working_DEA_dir'])
    except OSError as e:
        print "can't create directory '%(dirname)s'." % { "dirname": config['working_DEA_dir'] }
        print e.errno
        print e
        sys.exit()

    try:
        os.chdir(config['working_DEA_dir'])
    except OSError as e:
        print "can't chdir into '%(dirname)s'." % { "dirname": config['working_DEA_dir'] }
        print e.errno
        print e

#    try:
#        os.mkdir('DEA')
#    except OSError as e:
#        print "can't create DEA directory: %(ecode)s" % { "ecode": e.strerror }
#        sys.exit(13)
#
#    os.chdir('DEA')

    for filename in files:
        fmatch=re.match("Sample_", filename)
        if (fmatch):
            print "Creating symlink to '%(filename)s'" % { 'filename': filename }
#            os.symlink("../align/" + str(filename), str(filename))
            os.symlink("../" + config['working_alignment_dir'] + "/" + str(filename), str(filename))
            os.chdir(filename)
#            rst_script=os.path.join(config['rst_path'], 'bin', 'STAR_merge_gene_counts.py')
#
#            if config['paired']:
#                try:
#                    subprocess.check_call([rst_script])
#                except OSError as e:
#                    print "can't run %(scriptname)s" % { 'scriptname': rst_script }
#                    sys.exit(14)
#            else:
#                try:
#                    subprocess.check_call([rst_script, '--seonly'])
#                except OSError as e:
#                    print "can't run %(scriptname)s --seonly" % { 'scriptname': rst_script }

#            os.chdir(os.path.join(curdir, 'DEA'))
            os.chdir(os.path.join(curdir, config['working_DEA_dir']))
#
#   copy & run make_gene_cnts_per_sample.sh script from rst directory to curdir
#
#    shutil.copyfile(os.path.join(config['rst_path'], 'bin', 'make_gene_cnts_per_sample.sh'), 'make_gene_cnts_per_sample.sh')
#
#    try:
#        #subprocess.check_call(['sh', 'make_gene_cnts_per_sample.sh'])
#        subprocess.check_call("sh make_gene_cnts_per_sample.sh", shell=True)
#    except OSError as e:
#        print "can't run make_gene_cnts_per_sample.sh: %(estring)s" % { 'estring': e.strerror }
#        sys.exit(15)
#
#   run join_gene_cnts_opt.sh script from rst directory 
#   this will create the gene exp matrix file
#
#    rst_script = os.path.join(config['rst_path'], 'bin', 'join_gene_cnts_opt.sh')
#
#    print "clength: %(clength)i, elength: %(elength)i." % { "clength": clength, "elength": elength }
#
#    try:
#        subprocess.check_call(rst_script + " " + "-c " + str(clength) + " -e " + str(clength + 1) + " -E " + str(clength + elength), shell=True)
#    except OSError as e:
#        print "can't run join statement: %s" % e.strerror
#
#
#    for j in config['input']['experimental'][i]:
#        jlength=length(j)

    datafilename='C_v_E.txt'
    deseq2_script=os.path.join(config["rst_path"], "bin", "create_DESeq2_cmd_batch_file.py")
    strand = config['strand'] + 2
    try:
        subprocess.check_call(deseq2_script + " " + "--numberOfControls " + str(clength) + \
                " --numberOfExperimentals " + str(elength) + " --datafile " + datafilename + \
                " --org " + config['org'] + " --gProfilerkey " + config['gProfilerkey'] + \
                " --dbkey " + config['dbkey'] + " --strand " + str(strand) + " --aligndir " + config['working_alignment_dir'] + \
                " > DESeq2.Rscript", \
                shell=True)
    except OSError as e:
        print "can't run %(scriptname)s : %(errorstr)s" % { "scriptname": deseq2_script, "errorstr": e.strerror }

    print("""

    A file named DESeq2.Rscript has been created in the DEA directory (path = %(DEAdir)s/DESeq2.Rscript).

    Before you can proceed, please make sure all alignment jobs have finished. Use this command:

    qstat -u $USER

    This file can be directly submitted to a PBS cluster using this command:
    qsub DESeq2.Rscript
    To run that specific command, you must cd into the DEA directory:
    cd %(DEAdir)s
    For the analysis to run to completion, you must have R in your path and several R
    libraries need to be installed:
    tidyverse
    edgeR
    DESeq2
    ggplot2
    dplyr
    vegan
    pheatmap
    gprofiler2

    Finally, the organism-specific R annotation database for your organism of interest must be installed.
    You can find a list of databases here: https://www.bioconductor.org/packages/release/BiocViews.html#___Organism.
    Search that table for databases that start with "org." There are also some settings in the YAML configuration
    file pertinent to this functionality.

    In the VARI BBC, you can accomplish all of this, except for possibly the organism annotation database, by loading
    the bbc R module:
    module load bbc/R/R-3.6.0

    """) % { 'DEAdir': config['working_DEA_dir'] }

