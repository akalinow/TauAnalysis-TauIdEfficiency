import TauAnalysis.TauIdEfficiency.tools.castor_mirror as mirror

# Import the samples into this namespace
from samples import minbias_mc_pythia6, minbias_mc_pythia8, qcd_mc_pythia6, qcd_mc_pythia6_recoTrackDowngrade, qcd_mc_pythia8, qcd_mc_pythia8_recoTrackDowngrade, ztautau_mc, zllPU156bx_mc, data, data_runs132440to133802, data_runs135821to141887, data_runs141950to144114

#mirror.LOCAL_DIRECTORY = "/tmp/tau_commissioning_friis"
mirror.LOCAL_DIRECTORY = "/tmp/tau_commissioning_veelken"

#_sample_list = [ data, ztautau_mc, minbias_mc_pythia6, minbias_mc_pythia8, qcd_mc_pythia6, qcd_mc_pythia8 ]
#_sample_list = [ data, ztautau_mc, qcd_mc_pythia6, qcd_mc_pythia8 ]
#_sample_list = [ data, qcd_mc_pythia8 ]
#_sample_list = [ ztautau_mc ]
#_sample_list = [ data, qcd_mc_pythia8, ztautau_mc ]
#_sample_list = [ data_runs132440to133802, data_runs135821to141887, data_runs141950to144114, data ]
_sample_list = [ data ]
#_sample_list = [ zllPU156bx_mc ]

# Update the samples to use any existing local files
for sample in _sample_list:
    mirror.update_sample_to_use_local_files(sample)

# Script usage - copy all files 
if __name__ == "__main__":
    print "Copying CASTOR files to local area:", mirror.LOCAL_DIRECTORY
    # 20 concurrent rfcp jobs
    mirror.mirror_samples(_sample_list, max_jobs = 20)




