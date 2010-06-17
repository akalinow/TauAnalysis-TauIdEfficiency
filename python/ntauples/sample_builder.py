from __future__ import with_statement
from TauAnalysis.TauIdEfficiency.ntauples.SampleManager \
        import NtupleSample, NtupleSampleCollection
import sys

'''

sample_builder

Build NtupleSampleCollections given a list of sub-datasets
and a JSON file that maps the desired datasets to collections
of files and corresponding integrated luminsoity.

Author: Evan K. Friis (UC Davis)

$Id: sample_builder.py,v 1.6 2010/06/17 20:52:44 friis Exp $

'''

def build_sample(lumifile, sample_name, mode, take_every=1, datasets=[]):
    ''' Combine [datasets] into a NtupleSampleCollection

    The combine method ([mode]) can be either "add" (i.e.  concatenate weeks of
    data) or "merge" (i.e. merge pt hat bins of MC QCD).  The [lumifile] should
    be a JSON file produced by the scripts/lumiCalc.py utility.

    The take_every parameter indicates what subset of files to take (for 
    prototyping).  take_every=10 would take every tenth file.  The scaleFactor 
    is increased to account for the missing files.

    '''

    output=None
    with open(lumifile, 'r') as lumi_map_file:
        # Retrive information about the different sub-samples
        lumi_map = None
        if sys.version_info >= (2, 6):
            import json
            lumi_map = json.load(lumi_map_file)
        else:
            lumi_map = eval(lumi_map_file.read())
        # Retrieve all the desired datasets
        datasets_to_add = []
        for dataset in datasets:
            if dataset not in lumi_map:
                raise KeyError, "dataset: %s not found in luminosity directory"\
                        % dataset
            dataset_info = lumi_map[dataset]
            # Add this sample

            allEvents = -1 # allEvents not set will default to scaleFactor == 1
            if dataset_info.get('allEvents') is not None: allEvents = dataset_info['allEvents']

            # Add every nth file
            files_to_add = [ 
                file for index, file in enumerate(dataset_info['files']) 
                if index % take_every == 0]

            datasets_to_add.append(
                NtupleSample(dataset, int_lumi=dataset_info['int_lumi'], allEvents=allEvents,
                             files=files_to_add, directory=dataset_info['directory'], prescale=1.0)
            )
        # Create the merged dataset
        output = NtupleSampleCollection(sample_name, subsamples=datasets_to_add, 
                                        mode = mode)
    return output
