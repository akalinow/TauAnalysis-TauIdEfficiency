import json
from TauAnalysis.TauIdEfficiency.ntauples.SampleManager \
        import NtupleSample, NtupleSampleCollection

'''

sample_builder

Build NtupleSampleCollections given a list of sub-datasets
and a JSON file that maps the desired datasets to collections
of files and corresponding integrated luminsoity.

Author: Evan K. Friis (UC Davis)

$Id: $

'''

def build_sample(lumifile, sample_name, mode, *datasets):
    ''' Combine [datasets] into a NtupleSampleCollection

    The combine method ([mode]) can be either "add" (i.e.  concatenate weeks of
    data) or "merge" (i.e. merge pt hat bins of MC QCD).  The [lumifile] should
    be a JSON file produced by the scripts/lumiCalc.py utility.

    '''

    output=None
    with open(lumifile, 'r') as lumi_map_file:
        # Retrive information about the different sub-samples
        lumi_map = json.load(lumi_map_file)
        # Retrieve all the desired datasets
        datasets_to_add = []
        for dataset in datasets:
            if dataset not in lumi_map:
                raise KeyError, "dataset: %s not found in luminosity directory"\
                        % dataset
            dataset_info = lumi_map[dataset]
            # Add this sample
            datasets_to_add.append( 
                NtupleSample(dataset, int_lumi=dataset_info['int_lumi'],
                             files=dataset_info['files'], prescale=1.0)
            )
        # Create the merged dataset
        output = NtupleSampleCollection(sample_name, subsamples=datasets_to_add, 
                                        mode = mode)
    return output
