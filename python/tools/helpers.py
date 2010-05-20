def addToSequence(process, sequence, toAdd):
    ''' Add toAdd to sequence in process 
    
    If sequence does not exist in process, it will be created.

    '''
    if hasattr(process, sequence):
        my_seq = getattr(process, sequence)
        my_seq += (toAdd)
    else:
        setattr(process, sequence, cms.Sequence(toAdd))
