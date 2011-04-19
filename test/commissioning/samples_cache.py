import TauAnalysis.TauIdEfficiency.tools.castor_mirror2 as mirror

#--------------------------------------------------------------------------------
# Script usage: python samples_cache.py
#              to copy all files referenced by sample_list to local disk
#--------------------------------------------------------------------------------

# Import the samples into this namespace
from samples import ztautau_mc, qcddijet_mc, ppmux15_mc, wjets_mc, \
        data_dijet, data_ppmux, data_wjets

mirror.LOCAL_DIRECTORY = "/tmp/abdollah/tau_fakerate_ntuples" # on lxplus430
##mirror.LOCAL_DIRECTORY = "/data2/friis/tau_fakerate_ntuples"

class SampleWrapper(object):
    " Wrapper class to automatically update a sample when it is accessed "
    def __init__(self, wrapped):
        self.wrapped = wrapped
        self.mirrored = False
    def __getattr__(self, attribute):
        if not self.mirrored:
            print "samples_cache is updating", self.wrapped.name
            mirror.update_sample_to_use_local_files(self.wrapped,
                                                    ##only_local=True)
                                                    only_local=False)
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

_sample_list = []

#--------------------------------------------------------------------------------
# define QCD muon-enriched samples
_sample_list.extend([ ppmux15_mc, data_ppmux ])
#_sample_list.append(ppmux_mc)
#_sample_list.append(data_ppmux)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define W + jets samples
_sample_list.extend([ wjets_mc, data_wjets ])
#_sample_list.append(wjets_mc)
#_sample_list.append(data_wjets)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define Z --> tau+ tau- samples for tau id. efficiency plots
_sample_list.extend([ ztautau_mc ])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define QCD di-jet samples
_sample_list.extend([ qcddijet_mc, data_dijet ])
#_sample_list.append(qcddijet_mc)
#_sample_list.append(data_dijet)
_sample_list = [
    ztautau_mc,
    qcddijet_mc, data_dijet,
    wjets_mc, data_wjets,
    ppmux15_mc, data_ppmux
]
#--------------------------------------------------------------------------------

if __name__ == "__main__":
    print "Copying CASTOR files to local area:", mirror.LOCAL_DIRECTORY
    # run up to 50 rfcp jobs concurrently
    mirror.mirror_samples(_sample_list, max_jobs = 50)
