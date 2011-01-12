import TauAnalysis.TauIdEfficiency.tools.castor_mirror as mirror

#--------------------------------------------------------------------------------
# Script usage: python samples_cache.py
#              to copy all files referenced by sample_list to local disk
#--------------------------------------------------------------------------------

# Import the samples into this namespace
from samples import ztautau_mc, zttPU156bx_mc, ztautau_test,\
        qcddijet_mc, ppmux10_mc, ppmux10PU156bx_mc, ppmux15_mc, ppmux15PU156bx_mc, wmunu_mc, wmunuPU156bx_mc, \
        data_dijet, data_ppmux, data_wjets
        #data_dijet_runs132440to135802, data_dijet_runs135821to141887, data_dijet_runs141950to144114, data_dijet, \
        #data_ppmux_runs132440to145761, data_ppmux_runs145762_147116, data_ppmux_runs147117_149442, data_ppmux, \
        #data_wjets_runs132440to145761, data_wjets_runs145762_147116, data_wjets_runs147117_149442, data_wjets

#mirror.LOCAL_DIRECTORY = "/tmp/tau_commissioning_friis"
mirror.LOCAL_DIRECTORY = "/data1/friis/tau_fakerate_ntuples"

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
_sample_list.extend([ ppmux10_mc, ppmux10PU156bx_mc, ppmux15_mc, ppmux15PU156bx_mc, data_ppmux ])
#_sample_list.append(ppmux_mc)
#_sample_list.append(ppmuxPU156bx_mc)
#_sample_list.append(data_ppmux)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define W + jets samples
_sample_list.extend([ wmunu_mc, wmunuPU156bx_mc, data_wjets ])
#_sample_list.append(wmunu_mc)
#_sample_list.append(wmunuPU156bx_mc)
#_sample_list.append(data_wjets)
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define Z --> tau+ tau- samples for tau id. efficiency plots
_sample_list.extend([ ztautau_mc, zttPU156bx_mc, ztautau_test ])
#_sample_list.extend([ zttPU156bx_mc ])
#--------------------------------------------------------------------------------

#--------------------------------------------------------------------------------
# define QCD di-jet samples
_sample_list.extend([ qcddijet_mc, data_dijet ])
#_sample_list.append(qcddijet_mc)
#_sample_list.append(data_dijet)
#--------------------------------------------------------------------------------

if __name__ == "__main__":
    print "Copying CASTOR files to local area:", mirror.LOCAL_DIRECTORY
    # run up to 20 rfcp jobs concurrently
    mirror.mirror_samples(_sample_list, max_jobs = 30)
