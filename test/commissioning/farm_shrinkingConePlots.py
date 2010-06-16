import os
import string
import subprocess

from shrinkingConePlots import distributions
from shrinkingConePlots import efficiencies
from shrinkingConePlots import efficiency_versus

# Based on submitToBatch by Christian

# get name of directory in which config files will be created;
workingDirectory = os.getcwd()
submissionDirectory = os.path.join(workingDirectory,  "lxbatch")

if not os.path.isdir(submissionDirectory):
    os.mkdir(submissionDirectory)

# compose name of modified config file including the replacements
script = string.Template("""#!/bin/csh
limit vmem unlim
cd $workingDir
setenv SCRAM_ARCH slc5_ia32_gcc434
eval `scramv1 runtime -csh`
source root526.csh
python shrinkingConePlots.py $myArgs
""")


# Plot all distributions
for distribution in distributions.keys():
    print "Building script for plot", distribution
    script_file = os.path.join(submissionDirectory, 'plot_dist_%s.csh' % distribution)
    scf = open(script_file,"w")
    scf.write(script.substitute( 
        {'workingDir' : workingDirectory, 
         'myArgs' : '-p %s' % distribution
        }
    ))
    scf.close()

    # make shell script executable
    os.chmod(script_file,0744)
    
    logFile = script_file.replace(".csh", ".out")
    jobName = distribution
    bsubCommand = 'bsub -q 1nh ' + ' -J ' + jobName + ' -L /bin/csh -eo ' + logFile + ' -oo ' + logFile
    bsubCommand += ' < ' + script_file
    subprocess.call(bsubCommand, shell = True)

for eff in efficiencies.keys():
    for eff_vs in efficiency_versus.keys():
        distribution = "_".join([eff, eff_vs])
        print "Building script for plot", distribution
        script_file = os.path.join(submissionDirectory, 'plot_dist_%s.csh' % distribution)
        scf = open(script_file,"w")
        scf.write(script.substitute( 
            {'workingDir' : workingDirectory, 
             'myArgs' : '-p %s -v %s' % (eff, eff_vs),
            }
        ))
        scf.close()

        # make shell script executable
        os.chmod(script_file,0744)
        
        logFile = script_file.replace(".csh", ".out")
        jobName = distribution
        bsubCommand = 'bsub -q 1nh ' + ' -J ' + jobName + ' -L /bin/csh -eo ' + logFile + ' -oo ' + logFile
        bsubCommand += ' < ' + script_file
        subprocess.call(bsubCommand, shell = True)
