import TauAnalysis.TauIdEfficiency.tools.castor_mirror as mirror

#--------------------------------------------------------------------------------
# Script usage: python samples_cache.py
#              to copy all files referenced by sample_list to local disk
#--------------------------------------------------------------------------------

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


#--------------------------------------------------------------------------------
# define QCD di-jet samples
#_sample_list = [ qcddijet_mc, data_dijet ]
#_sample_list = [ qcddijet_mc ]
#_sample_list = [ data_dijet ]
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define QCD muon-enriched samples
#_sample_list = [ ppmux_mc, data_ppmux ]
#_sample_list = [ ppmux_mc ]
#_sample_list = [ data_ppmux ]
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define W + jets samples
#_sample_list = [ wmunu_mc, data_wjets ]
#_sample_list = [ wmunu_mc ]
#_sample_list = [ data_wjets ]
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define Z --> tau+ tau- samples for tau id. efficiency plots
#_sample_list = [ ztautau_mc, zttPU156bx_mc ]
#_sample_list = [ ztautau_mc ]
_sample_list = [ zttPU156bx_mc ]
#--------------------------------------------------------------------------------

if __name__ == "__main__":
    print "Copying CASTOR files to local area:", mirror.LOCAL_DIRECTORY
    # run up to 20 rfcp jobs concurrently
    mirror.mirror_samples(_sample_list, max_jobs = 20)




