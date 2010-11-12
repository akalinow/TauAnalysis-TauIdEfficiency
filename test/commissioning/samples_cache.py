import TauAnalysis.TauIdEfficiency.tools.castor_mirror as mirror

# Import the samples into this namespace
from samples import minbias_mc_pythia6, minbias_mc_pythia8, \
        qcd_mc_pythia6, qcd_mc_pythia8, \
        ztautau_mc, zttPU156bx_mc, zllPU156bx_mc, \
        data, data_runs132440to133802, data_runs135821to141887, \
        data_runs141950to144114

#mirror.LOCAL_DIRECTORY = "/tmp/tau_commissioning_friis"
mirror.LOCAL_DIRECTORY = "/tmp/tau_commissioning_veelken"

class SampleWrapper(object):
    " Wrapper class to automatically update a sample when it is accessed "
    def __init__(self, wrapped):
        self.wrapped = wrapped
        self.mirrored = False
    def __getattr__(self, attribute):
        if not self.mirrored:
            print "samples_cache is updating", self.wrapped.name
            mirror.update_sample_to_use_local_files(self.wrapped)
            self.mirrored = True
        return getattr(self.wrapped, attribute)
    #def __setattr__(self, attribute):
    #    setattr(self.wrapped, attribute)

# Wrap all of the samples in the current namespace
current_objects = globals().keys()
for name in current_objects:
    the_object = globals()[name]
    # Check if it is one of the samples
    if hasattr(the_object, 'build_ntuple_manager'):
        print "samples_cache is wrapping:", name
        globals()[name] = SampleWrapper(the_object)


# Script usage - copy all files
#_sample_list = [ data, ztautau_mc, minbias_mc_pythia6, minbias_mc_pythia8, qcd_mc_pythia6, qcd_mc_pythia8 ]
#_sample_list = [ data, ztautau_mc, qcd_mc_pythia6, qcd_mc_pythia8 ]
#_sample_list = [ data, qcd_mc_pythia8 ]
#_sample_list = [ ztautau_mc ]
#_sample_list = [ data, qcd_mc_pythia8, ztautau_mc ]
#_sample_list = [ data_runs132440to133802, data_runs135821to141887, data_runs141950to144114, data ]
#_sample_list = [ data ]
#_sample_list = [ zllPU156bx_mc ]
_sample_list = [ zttPU156bx_mc ]
#_sample_list = [ zttPU156bx_mc, data_runs132440to133802, data_runs135821to141887, data_runs141950to144114, data ]
if __name__ == "__main__":
    print "Copying CASTOR files to local area:", mirror.LOCAL_DIRECTORY
    # 20 concurrent rfcp jobs
    mirror.mirror_samples(_sample_list, max_jobs = 20)




