import TauAnalysis.TauIdEfficiency.tools.castor_mirror as mirror

# Import the samples into this namespace
from samples import minbias_mc_pythia6, minbias_mc_pythia8, qcd_mc_pythia6, qcd_mc_pythia8, ztautau_mc, data

mirror.LOCAL_DIRECTORY = "/tmp/tau_commissioning_friis"

#_sample_list = [ data, ztautau_mc, minbias_mc_pythia6, minbias_mc_pythia8, qcd_mc_pythia6, qcd_mc_pythia8 ]
_sample_list = [ data, ztautau_mc, qcd_mc_pythia6, qcd_mc_pythia8 ]

# Update the samples to use any existing local files
for sample in _sample_list:
    mirror.update_sample_to_use_local_files(sample)

# Script usage - copy all files 
if __name__ == "__main__":
    print "Copying CASTOR files to local area:", mirror.LOCAL_DIRECTORY
    # 20 concurrent rfcp jobs
    mirror.mirror_samples(_sample_list, max_jobs = 20)




