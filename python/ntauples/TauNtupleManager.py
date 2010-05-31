import ROOT
import TauAnalysis.TauIdEfficiency.ntauples.helpers as helpers
from TauAnalysis.TauIdEfficiency.ntauples.TauNtuple import TauNtuple

class TauNtupleManager(object):
    def __init__(self, events, ntuple_name):
        self.ntuple_name = ntuple_name
        # make sure that any aliases are set correctly for a TChain
        if isinstance(events, ROOT.TChain):
            helpers.copy_aliases(events)

        # Get a list of each alias in the TTree 
        aliases = [ alias.GetName().split('#') 
                   for alias in events.GetListOfAliases() ]

        # Find collections used for this ntuple
        collections = []
        for alias in helpers.filter_aliases(aliases, ntuple_name):
            collections.append(alias[1])

        self.ntuples = {}
        # Make a TauNtuple for each collection
        for collection in collections:
            variables = []
            for collection_aliases in helpers.filter_aliases(
                aliases, ntuple_name, collection):
                variables.append(collection_aliases[2])
            # Make a new ntuple
            self.ntuples[collection] = \
                TauNtuple(ntuple_name, collection, variables)

    def get_ntuple(self, ntuple):
        " Retrive the specfied ntuple "
        if ntuple in self.ntuples:
            return self.ntuples[ntuple]
        else: 
            raise KeyError, "Attempt to get non-existent ntuple", ntuple

    def __repr__(self):
        output = ""
        output += "== TauNtupleManager\n"
        output +=  "===  Available ntuples for label %s\n" % self.ntuple_name
        for ntuple in self.ntuples.keys():
            output += "==== %s\n" % ntuple
        return output

if __name__ == "__main__":
    # Example of use
    import os
    cmssw_base = os.environ['CMSSW_BASE']
    example_file = os.path.join(
        cmssw_base, 'src', 'TauAnalysis/TauIdEfficiency/test/example_ntuple.root')
    file = ROOT.TFile.Open(example_file, "READ")
    events = file.Get("Events")

    manager = TauNtupleManager(events, "exampleNtuple")
    # Get list of available collections
    print manager

    # Get one of the collections
    allLayer1Taus = manager.get_ntuple('allLayer1Taus')
    print allLayer1Taus
     
    # Make an expression or a selection
    print allLayer1Taus.expr('$pt') < 5



