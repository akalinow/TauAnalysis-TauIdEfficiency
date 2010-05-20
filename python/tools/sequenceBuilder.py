import FWCore.ParameterSet.Config as cms

'''

tools.sequenceBuilder.py

Create pat::Taus from a specified collection, match
them to the tag and probe jets, and create ntuples
of the resulting collections.

Author: Evan K. Friis (UC Davis)

'''

##import TauAnalysis.TauIdEfficiency.patConfiguration.tauProductionPrototypes as tauProto
##import TauAnalysis.TauIdEfficiency.patConfiguration.matchingPrototypes as matchProto

def buildTauSequence(
    process, 
    collectionName = "",
    producerProtoType = None,
    cleanerProtoType = None,
    cleanerOverlapCheckerSource = "pfJetsTagAndProbes",
    matchedCollections = []):
    '''
    blah
    '''

    # Build basic pat::Taus from input collection
    myPatTaus = cleanerProtoType.clone()
    
    patTauProductionName = collectionName + "BuildPatTaus"

    # Insert pat production into process
    setattr(process, patTauProductionName, myPatTaus) 

    # First step in the sequence is to build the pat taus
    outputSequence = cms.Sequence(myPatTaus)

    # Now we split the pat::Taus into three categories
    # * matched to tag jet
    # * matched to highest pt probe jet
    # * matched to second highest pt probe jet

    myCleaner = cleanerProtoType.clone()
    myCleaner.src = cms.InputTag(patTauProductionName)

    # Set the overlap checker to use our input for each of the overlap matchers
    overlaps = myCleaner.checkOverlaps.parameterNames_()
    for overlap in overlaps:
        overlapPSet = getattr(myCleaner.checkOverlaps, overlap)
        overlapPSet.src.setModuleLabel(cleanerOverlapCheckerSource)

    # Add the overlap checker to the process & sequence
    cleanerName = collectionName + "WithMatches"
    setattr(process, cleanerName, myCleaner)
    outputSequence += myCleaner

    # Now for each of the desired matching collections, produce a collection of
    # pat::Taus associated with that matching (matching was done by cleaner)
    finalOutputNames = []
    for name, cut in matchedCollections:
        selectorName = collectionName + name
        finalOutputNames.append(selectorName)
        mySelector = cms.EDProducer("PATTauSelector",
            src = cms.InputTag(cleanerName),
            cut = cms.string(cut)
        )
        # Add to process
        setattr(process, selectorName, mySelector)
        outputSequence += mySelector
    # We now have a collection of taus corresponding to each of the different
    # matching collections.  Return the full sequence, and a list of the output
    # collections
    return outputSequence

def buildDijetTauSequence(process, **kwargs):
    ''' Build a sequence for the dijet fake rate measurement method '''
    return buildTauSequence(
        process,
        matchedCollections = [
            ('TagJet', 'hasOverlaps("TagJet")'),
            ('HighestPtProbe', 'hasOverlaps("HighestPtProbe")'),
            ('SecondHighestPtProbe', 'hasOverlaps("SecondHighestPtProbe")'),
        ], **kwargs)









