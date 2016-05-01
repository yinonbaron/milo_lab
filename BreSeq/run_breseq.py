#This script should be executed from within a dir contating seperate sub-dirs for each sample (fastq files)
#Sample sub-dir can contain more than one fastq file that belong to this sample

import sys, os, shutil, time


#True for -p parameters (population mode), fasle for isolated clone mode (ommits -p)
POP = True

#path to breseq + command
s = '/home/ronm/breseq-0.27.1-Linux-x86_64/bin/breseq -j 4 '

if POP:
	s = s + '-p '

#Path to refrences files here
BWRef = '/home/ronm/NGS/Refs/CP009273_1.gb'
plasmidRef = '/home/ronm/NGS/Refs/F1_plasmid.gb'

#add refrence files to command
s = s + '-r %s -r %s -o ' % (BWRef,plasmidRef)

#get current directory
mypath =  os.getcwd()

#output directory will be /Outputs from current dir
dest = mypath +'/Outputs'



print "Analyzing directories in %s... " % mypath
print 'The output directory will be %s ' % dest
print '\n'

# go over all .gz files in dir, unzip if required
allfiles = [ f for f in os.listdir(mypath) if os.path.isfile(os.path.join(mypath,f)) ]
for f in allfiles:
    if f.endswith(".gz"):
        currentdir = mypath+'/'
        filename = currentdir + f
        com = 'gunzip -d %s' %filename
        print 'unzipping %s...\n' % filename
        os.system(com)

# go over all sub-dirs in current dir
alldirs = [ f for f in os.listdir(mypath) if os.path.isdir(os.path.join(mypath,f)) ]
for d in alldirs:
	real_dir = 0
	print 'Processing directory %s : ' %d
	gotopath = mypath+'/'+d
	os.chdir(gotopath)
	#For each sample make its own output sub-dir within /Outputs
	destpath = dest +'/' + d + ' '
	#add to command the output destenation after the -o parameter.
	com = s + destpath
	# go over all files in sub-dir looking for fastq
	allfiles = [ f for f in os.listdir(gotopath) if os.path.isfile(os.path.join(gotopath,f)) ]
	onlyfiles = []
	for f in allfiles:
	#if fastq is found, add it to the list of fastq files to be analyzed
		if f.endswith(".fastq"):
			real_dir =1
                	onlyfiles.append(f)
	com = com + ' '.join(onlyfiles)
	print com
	print '\n'
	#if at least one fastq file was found in this dir --> run com and copy/rename summary to /SumOfRes dir
	if real_dir :
		os.system(com)
		sourcedir = dest +'/' + d  + '/output'
		destdir = mypath  + '/SumOfRes/' + d + '/'

		if not os.path.exists(destdir):
	    		os.makedirs(destdir)
		src_files = os.listdir(sourcedir)
		for file_name in src_files:
			oldFileName = file_name
			if file_name.endswith(".gd"):
				newFileName = d + '.gd'
			elif file_name.startswith("index"):
				newFileName = d + '.html'
			else:
				newFileName = oldFileName

			full_new_file_name = os.path.join(sourcedir, newFileName)
			full_old_file_name = os.path.join(sourcedir, oldFileName)
			os.rename(full_old_file_name, full_new_file_name)
	    		if (os.path.isfile(full_new_file_name)):
				shutil.copy(full_new_file_name, destdir)











