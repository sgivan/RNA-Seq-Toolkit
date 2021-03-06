#!/usr/bin/env python

import os, sys, argparse, shutil
import re, yaml, subprocess, time

# start with command line options
argparser = argparse.ArgumentParser(description="Parse tab-delimited file")

argparser.add_argument("--infile", type=str,  help="file to parse", default="file")
argparser.add_argument("--verbose", action="store_true", help="verbose messages to terminal")

args = argparser.parse_args()

if args.verbose: print 'input file: %(filename)s' % { "filename": args.infile }

f_yaml = open(args.infile, 'r') # open infile

config = yaml.load(f_yaml)
curdir=os.getcwd()
if args.verbose: print "current working directory: '%(workdir)s'" % { "workdir": curdir }

os.environ['PATH']=config['rst_path'] + "/bin" + ":" + os.environ['PATH']

clength=len(config['input']['control'])

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
            print "\tthis conflicts with configuration file\n\tplease revise\n\texiting now"
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

def monitor_slurm_jobs(slurmjobs):
#   not perfect, but it works
    if args.verbose:
        print "slurm jobs to monitor:"
        print slurmjobs

    wait=1
    while wait:
        wait=0

        for jobid in slurmjobs:

            rtn=0
            try:
                rtn=subprocess.check_output("squeue -o %T --noheader -j " + jobid, shell=True)
            except subprocess.CalledProcessError as e:
                print "can't call squeue with jobid %(jobid)i: %(ecode)i" % { "jobid": jobid, "ecode": e.returncode }
                sys.exit(15)
  
            if rtn:
                wait=1
                if args.verbose: print "waiting for job " + jobid
                time.sleep(5)
            
#
# end of functions
#

#
# create new working directory
# fail if the directory already exists
#
slurmjobs=[]
filemap={}
filemap['control']=[]
filemap['experimental']=[]
if config['setup_files']:
    print "\nsetting up file structure for input files\n"
    if not os.access(config['working_datadir'], os.F_OK):
        if args.verbose: print "creating new directory to place renamed files: %(newdir)s." % { "newdir": config['working_datadir'] }
        if os.access(os.path.split(config['working_datadir'])[0], os.W_OK):
            os.mkdir(config['working_datadir'])

        else:
            print "can't create the directory"
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
        elength=len(config['input']['experimental'])

        if args.verbose: print "experimental data sets: %(length)i" % { "length": elength }
        efiles={}
        for i in config['input']['experimental']:
            efiles[i]=[]
#
#   This cut corresponds to the sample replicate. There can be any number of sample replicates. which will have either a single file (non-PE) or a pair of files (Paired End)
#
            if args.verbose: print "sample replicates in set %(eset)s: %(filenames)s" % { "eset": i, "filenames": config['input']['experimental'][i] }
            number_of_reps=len(config['input']['experimental'][i])
            if args.verbose: print "number of replicates: %(numseqs)i." % { "numseqs": number_of_reps }

            for j in config['input']['experimental'][i]:

                sample_number += 1
                if args.verbose: print "\n\t%(setname)s - %(repname)s will be given symbolic name 'Sample_%(sint)s'" % { "setname": i, "repname": j, "sint": sample_number }

                create_file_struct(sample_number, config['input']['experimental'][i][j], config, curdir)
                efiles[i].append([j, 'Sample_' + str(sample_number)])

        filemap['experimental']=efiles

#
#   print out yaml file containing map of original files to standardized files
#
    mapfile = file('filemap.yaml', 'w')
    yaml.dump(filemap, mapfile)
    mapfile.close()

# End of yaml file map

    if args.verbose: print str(sample_number) + ' samples'
    print "Input file setup finished."

if config['preprocess']:
    print '\n\n pre-process the input data.'

    try:
        os.mkdir(config['working_alignment_dir'])
    except OSError as e:
        print "can't create directory '%(dirname)s'." % { "dirname": config['working_alignment_dir'] }
        print e.errno
        print e

    try:
        os.chdir(config['working_alignment_dir'])
    except OSError as e:
        print "can't chdir into '%(dirname)s'." % { "dirname": config['working_alignment_dir'] }
        print e.errno
        print e


    print "creating symlnks to preprocess and alignment index files in " + config['working_alignment_dir']

    if os.access('index.preprocess', os.F_OK):
        print "Will not overwrite current 'index.preprocess' symlink.\nPlease remove it."
        sys.exit(6)

    try:
        os.symlink(config['filter_datadir'], 'index.preprocess')
    except OSError as e:
        print "can't create index.preprocess symlink pointing to '%(dirname)s.'" % { 'dirname': config['filter_datadir'] }
        print e.errno
        print e.filename
#        print e.strerr

    if os.access('index.align', os.F_OK):
        print "Will now overwrite current 'index.align' symlink.\nPlease remove it."
        sys.exit(7)

    try:
        os.symlink(config['index_datadir'], 'index.align')
    except OSError as e:
        print "can't create index.align symlink pointing to '%(dirname)s.'" % { 'dirname': config['index_datadir'] }
        print e.errno
        print e.filename
#        print e.strerr

    if os.access('index', os.F_OK):
        print "Will now overwrite current 'index' symlink.\nPlease remove it."
        sys.exit(8)

    try:
        os.symlink('index.preprocess', 'index')
    except OSError as e:
        print "can't create index symlink pointing to index.preproces"
        print e.errno
        print e.filename
#        print e.strerr

    if args.verbose: print "preprocess and align symlinks created in " + curdir

    setup_script=os.path.join(config['rst_path'], 'bin', 'setup.sh')
    if args.verbose: print "calling setup script '%(scriptname)s'." % { "scriptname": setup_script }
    try:
        out=subprocess.check_call(setup_script, shell=True)
    except subprocess.CalledProcessError as e:
        print "call to symlink failed"
        print "error code: %(ecode)i" % { "ecode": e.returncode }
        sys.exit(9)

    if args.verbose: print "running RST preprocessing routines"

    rst_script=os.path.join(config['rst_path'], 'bin', 'RNAseq_process_data.sh')

    out=""
    try:
        out=subprocess.check_output(rst_script + " --preprocess_only --submit --threads " + str(config['threads']) + " Sample_*", shell=True)
    except subprocess.CalledProcessError as e:
        print "call to %(rst)s failed with error code %(ecode)i" % { "rst": rst_script, "ecode": e.returncode }
        sys.exit(10)
        
    os.chdir(curdir)

    for line in str.splitlines(out):
        match = re.match("OUTPUT", line)
        if match:
            words = str.split(line)
            id = words[-1]
            slurmjobs.append(id)

    if config['align']: monitor_slurm_jobs(slurmjobs)

if config['align']:
    if args.verbose: print "\n\naligning data to reference genome sequence"

    os.chdir(config['working_alignment_dir'])

    try:
        os.remove('index')
    except OSError as e:
        print "can't remove index symlink: %(ecode)i" % { "ecode": e.errno }
        sys.exit(11)

    try:
        os.symlink('index.align', 'index')
    except OSError as e:
        print "can't create index symlink to index.align: %(ecode)i" % { "ecode": e.errno }
        sys.exit(12)

    rst_script=os.path.join(config['rst_path'], 'bin', 'RNAseq_process_data.sh')
    try:
        subprocess.check_call(rst_script + " --partial --submit --threads " + str(config['threads']) + " Sample_*", shell=True)
    except subprocess.CalledProcessError as e:
        print "can't call %(scriptname)s: %(ecode)i" % { "scriptname": rst_script, "ecode": e.returncode }
        sys.exit(13)

if config['diff_expression']:

    filemapfile = file('filemap.yaml', 'r')
    filemap=yaml.load(filemap_file)

    try:
        os.mkdir('DEA')
    except OSError as e:
        print "can't create DEA directory: %(ecode)s" % { "ecode": e.strerror }
        sys.exit(13)

    os.chdir('DEA')

    files=os.listdir('../align')

    for filename in files:
        fmatch=re.match("Sample_", filename)
        if (fmatch):
            print "Creating symlink to '%(filename)s'" % { 'filename': filename }
            os.symlink("../align/" + str(filename), str(filename))
            os.chdir(filename)
            rst_script=os.path.join(config['rst_path'], 'bin', 'STAR_merge_gene_counts.py')

            if config['paired']:
                try:
                    subprocess.check_call([rst_script])
                except OSError as e:
                    print "can't run %(scriptname)s" % { 'scriptname': rst_script }
                    sys.exit(14)
            else:
                try:
                    subprocess.check_call([rst_script, '--seonly'])
                except OSError as e:
                    print "can't run %(scriptname)s --seonly" % { 'scriptname': rst_script }

            os.chdir(os.path.join(curdir, 'DEA'))
#
#   copy & run make_gene_cnts_per_sample.sh script from rst directory to curdir
#
    shutil.copyfile(os.path.join(config['rst_path'], 'bin', 'make_gene_cnts_per_sample.sh'), 'make_gene_cnts_per_sample.sh')

    try:
        subprocess.check_call(['sh', 'make_gene_cnts_per_sample.sh'])
    except OSError as e:
        print "can't run make_gene_cnts_per_sample.sh: %(estring)s" % { 'estring': e.strerror }
        sys.exit(15)
#
#   copy & run join_gene_cnts.sh script from rst directory to curdir
#
    shutil.copyfile(os.path.join(config['rst_path'], 'bin', 'join_gene_cnts.sh'), 'join_gene_cnts.sh')

    try:
        subprocess.check_call(['sh', 'join_gene_cnts.sh'])
    except OSError as e:
        print "can't run join_gene_cnts.sh: %(estring)s" % { 'estring': e.strerror }
        sys.exit(16)

    for j in config['input']['experimental'][i]:
        jlength=length(j)

        datafilename=
        deseq2_script=os.path.join(config['rst_path'], 'bin', 'create_DESeq2_cmd_sbatch_file.py')
        subprocess.check_call([deseq2_script, '--numberOfControls=' + elength, '--numberOfExperimentals=' + jlength, '--datafile', '--prefix')


