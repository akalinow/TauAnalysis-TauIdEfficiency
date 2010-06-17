from ROOT import gROOT
gROOT.SetBatch(True)
import ROOT
from TauAnalysis.TauIdEfficiency.ntauples.TauNtupleManager \
        import TauNtupleManager

'''
Classes to support managing different samples with 
different integrated luminosities.

Author: Evan K. Friis (UC Davis)

'''

class NtupleSample(object):
    ''' 
    NtupleSample

    Holds a set of root files corresponding to a known luminosity.

    '''
    def __init__(self, name, int_lumi, scaleFactor=1.0, prescale=1.0, files=[], directory=""):
        self.name = name
        self.int_lumi = int_lumi
        self.scaleFactor = scaleFactor
        self.prescale = prescale
        self.files = ["".join([directory, file]) for file in files]
        self.events = None
        # Check if root supports ttree caching
        self.ttree_cache = False
        if ROOT.gROOT.GetVersionCode() >= 334336:
            self.ttree_cache = True

    def effective_luminosity(self):
        ''' Effective integrated luminosity, given prescale

        returns (int_lumi)/(prescale)
        '''
        return float(self.int_lumi)/float(self.prescale)

    def build_events(self):
        # Build TChain of events
        self.events = ROOT.TChain("Events")
        for file in self.files:
            self.events.AddFile(file)
        # Turn on read-ahead caching, described in
        # http://root.cern.ch/drupal/content/spin-little-disk-spin
        # (adds ~100 MB memory consumption)
        if self.ttree_cache:
            print "Enabling TTree Cache"
            self.events.SetCacheSize(10000000)
            self.events.AddBranchToCache("*")

    def get_events(self):
        if not self.events:
            self.build_events()
        # If we can, restart the cache learning phase
        if self.ttree_cache:
            pass
            self.events.SetCacheSize(0) # reset
            self.events.SetCacheSize(10000000)
            self.events.AddBranchToCache("*")
        return self.events

    def set_selection(self, selection=None, selection_file=None):
        ''' Retrive a TEntryList from the given file if it exists '''

        # Save old path
        old_path = ROOT.gDirectory.GetPath()

        if selection is None:
            self.get_events().SetEntryList(None)
            return

        # Check if we already got it
        print "Setting %s as selection for %s" % (selection, self.name)
        name_in_file = "%s_%s" % (self.name, selection.hash_string())
        file = ROOT.TFile.Open(selection_file, "UPDATE")
        file.cd()
        entry_list = file.Get(name_in_file)
        # Build if necessary
        if not entry_list:
            # Build entry list
            print "Selection list not found in file, building..."
            self.get_events().Draw(">>%s" % name_in_file, str(selection), 
                                   "entrylist")
            entry_list = ROOT.gDirectory.Get(name_in_file)
            entry_list.Write()

        print "Got entry list with %i entries" % entry_list.GetN()
        self.get_events().SetEntryList(entry_list)
        #file.Close()
        ROOT.gDirectory.cd(old_path)

    def norm_factor_for_lumi(self, target_int_lumi):
        ''' Return weight need to scale sample to target luminosity '''
        return target_int_lumi*self.scaleFactor/self.effective_luminosity()

    def build_ntuple_manager(self, name):
        ''' Build an ntuple manager corresponding to [name] in the associated files'''
        return TauNtupleManager(self.get_events(), name)

class NtupleSampleCollection(object):
    '''
    NtupleSampleCollection

    Holds a collection a set of NtupleSample collections.
    There are two luminosity modes, 'add', and 'merge'.
    In add mode, each subsample is concatenated together,
    and the effective luminosity is the sum of all the sub samples.
    In the merge case, the effective luminosity is taken from the sub sample
    with the lowest integrated luminosity.

    '''
    def __init__(self, name, subsamples=[], mode='add'):
        self.name = name
        self.subsamples = subsamples
        self.mode = mode

    def set_selection(self, selection=None, selection_file=None):
        ''' Set a static selection on this sample using a TEntrylists '''
        for subsample in self.subsamples:
            subsample.set_selection(selection, selection_file)

    def build_ntuple_manager(self, name):
        ''' Build an ntuple manager corresponding to [name] in the associated files'''
        if not self.subsamples:
            raise ValueError, "cannot retrieve ntuple, no subsamples defined"
        return self.subsamples[0].build_ntuple_manager(name)  

    def effective_luminosity(self):
        ''' Get effective int. luminosity for this sample colleciton '''
        if self.mode=='add':
            return sum(subsample.effective_luminosity() 
                       for subsample in self.subsamples)
        elif self.mode=='merge':
            # Take the lowest integrated luminsotiy one.  
            # (always scale down, never up)
            return min(subsample.effective_luminosity() 
                       for subsample in self.subsamples)

    def norm_factor_for_lumi(self, target_int_lumi):
        return target_int_lumi*1.0/self.effective_luminosity()

    def events_and_weights(self, target_lumi=None):
        ''' Yields a list of event TTree sources with the appropriate weights 

        The entire collection is weighted to correspond to an integrated
        luminosity of [target_lumi].
        '''
        # If no argument is supplied, don't reweight.
        if target_lumi is None:
            target_lumi = self.effective_luminosity()

        overall_norm_factor = self.norm_factor_for_lumi(target_lumi)

        # If we are adding up the subsamples, just scale each by the overall
        # factor of the collection
        if self.mode=='add':
            for subsample in self.subsamples:
                yield (subsample.get_events(), overall_norm_factor)
        # Otherwise, scale each subsample to the corresponding integrated
        # luminosity of the whole sample, then apply the correction factor.
        elif self.mode=='merge':
            for subsample in self.subsamples:
                yield (subsample.get_events(), 
                       subsample.norm_factor_for_lumi(self.effective_luminosity())*
                       overall_norm_factor)

if __name__ == "__main__":
    import unittest
    print "Running tests." 
    class TestSampleManagers(unittest.TestCase):
        def setUp(self):
            self.sample1 = NtupleSample(
                'sample1', 
                10.0, prescale=16.0, files = [])
            self.sample2 = NtupleSample(
                'sample2',
                20.0, prescale=1.0, files=[])
            self.sample3 = NtupleSample(
                'sample3',
                30.0, prescale=2.0, files = [])

            self.merger = NtupleSampleCollection(
                'merged',
                subsamples=[self.sample1, 
                            self.sample2,
                            self.sample3],
                mode='merge')
            self.adder = NtupleSampleCollection(
                'add',
                subsamples=[self.sample1,
                            self.sample2,
                            self.sample3],
                mode='add')

        def testNtupleSample(self):
            # test prescale
            self.assertEqual(self.sample1.effective_luminosity(), 10.0/16.0)
            self.assertEqual(self.sample1.norm_factor_for_lumi(10.0), 16.0)

        def testCollection(self):
            # merger should take the lowest as eff. int. lumi
            self.assertEqual(
                self.merger.effective_luminosity(),
                10.0/16.0)
            self.assertEqual(
                self.adder.effective_luminosity(),
                10.0/16.0 + 20.0 + 30.0/2.0)

            # When added, each subsample's integrated luminosity
            # should add up to the total
            self.assertEqual(
                sum(weight for event, weight in 
                    self.adder.events_and_weights(5.0)), 
                3*5.0/(10/16.0 + 20 + 30/2.0)
            )

            self.assertEqual(
                sum(weight for event, weight in 
                    self.merger.events_and_weights(5.0)), 
                5.0/(10/16.0) + 5.0/20 + 5.0/(30/2.0)
            )
    unittest.main()


