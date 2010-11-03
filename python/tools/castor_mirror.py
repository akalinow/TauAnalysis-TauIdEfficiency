import os
import subprocess
import re
import time
import itertools
import TauAnalysis.TauIdEfficiency.tools.castor as castor

# Where to store the files
LOCAL_DIRECTORY = "/tmp/tau_commissioning"

# How to identify castor files
_CASTOR_MATCHER = re.compile("/castor/cern.ch")

def is_on_castor(file):
    return file.find("rfio") != -1

def clean_name(file):
    return file.replace("rfio:", "")

def unixtime_from_timestamp(time_str):
    #Tue Jun  1 09:47:28 2010
    format = "%a %b %d %H:%M:%S %Y"
    parsed = time.strptime(time_str, format)
    return parsed

def local_version(castor_file):
    ''' Map a castor file to a local file '''
    return _CASTOR_MATCHER.sub(LOCAL_DIRECTORY, castor_file)

def local_version_current(castor_file):
    ''' Check if the local copy of [castor_file] exists and is up to date '''
    local_file = local_version(castor_file)
    if not os.path.exists(local_file):
        return False
    local_stat = os.stat(local_file)
    # Get last mod time of local file
    #local_mtime = time.ctime(local_stat.st_mtime)
    local_mtime = local_stat.st_mtime
    local_size = local_stat.st_size
    castor_stat = castor.rfstat(castor_file)
    castor_size = int(castor_stat["Size (bytes)"])
    castor_mtime = time.mktime(
        unixtime_from_timestamp(castor_stat["Last modify"]))
    # Check sizes are same
    if local_size != castor_size:
        print "Local copy of", castor_file, " is the wrong size: %i != %i"% (local_size, castor_size)
        return False
    # Check local file is newer
    if local_mtime < castor_mtime:
        print "Local copy of", castor_file, " is outdated %i > %i"% (local_mtime, castor_mtime)
        return False
    return True


def mirror_file(castor_file):
    ''' Initiate copy a castor file to the local directory.

    Returns a dictionary containing the background subprocess 
    and information about the file.
    '''
    # Get file path for new file
    local_path = local_version(castor_file)
    local_dirname = os.path.dirname(local_path)
    # Create the directories if necessary
    if not os.path.isdir(local_dirname):
        print "Creating local directory %s" % local_dirname
        os.makedirs(local_dirname)
    command = ['rfcp', castor_file, local_path]
    print "Requesting %s" % (castor_file, )
    subproc = subprocess.Popen(command, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE)
    return { "proc" : subproc,
            "pid" : subproc.pid,
            "dest" : local_path,
            "source" : castor_file,
            "ret" : None }

def wait_for_completion(current_jobs, finished_jobs, max_running_jobs):
    " Block until there are at most [max_running_jobs] running "
    while len(current_jobs.keys()) > max_running_jobs:
        time.sleep(1.0)
        # Check if any processes are done
        pids = current_jobs.keys()
        for pid in pids:
            proc_info = current_jobs[pid]
            proc_info['proc'].poll()
            # Update return status
            proc_info['ret'] = proc_info['proc'].returncode
            # Check if processes has terminated
            if proc_info['ret'] is not None:
                print "Copy %s has terminated with return code: %i"\
                        % (proc_info['source'], proc_info['ret'])
                if proc_info['ret'] == 28:
                    raise IOError, " target directory is out of space!  No more files can be copied"
                # Move process to finished jobs
                finished_jobs[pid] = proc_info
                # Remove from running jobs list
                del current_jobs[pid]


def stage_files(castor_files):
    ''' Request that castor stage the relevant files '''
    for castor_file in castor_files:
        # Don't care about keeping track of output
        try: 
            stager = subprocess.Popen(['stager_get', '-M', castor_file], 
                                      stdout=subprocess.PIPE, 
                                      stderr=subprocess.PIPE)
            stager.wait()
        except OSError: # File doesn't exist
            if not stage_files._complained_once:
                print "stager_get command doesn't exist - won't prestage files."
                stage_files._complained_once = True
            pass

stage_files._complained_once = False

def mirror_files(castor_files, max_jobs=20):
    ''' Copy [castor_files] to local disk.

    max_jobs specifies maximum number of concurrent copies
    '''
    # Request all copies be staged
    castor_files = list(castor_files)
    if not len(castor_files):
        return
    print "About to copy %i files" % len(castor_files)
    print "Requesting that CASTOR stage the files about to be copied."
    stage_files(castor_files)

    current_jobs = {}
    finished_jobs = {}
    try:
        for file in castor_files:
            # Check if we can submit this job
            if len(current_jobs.keys()) > max_jobs:
                wait_for_completion(current_jobs, finished_jobs, max_jobs)
            # submit it
            new_proc = mirror_file(file)
            current_jobs[new_proc['pid']] = new_proc

        # Wait for all jobs to finish up
        print "All jobs submitted.  Wait for final ones to finsh."
        wait_for_completion(current_jobs, finished_jobs, 0)
    finally:
        # terminate any running jobs
        for pid, proc in current_jobs.iteritems():
            print "Terminating %s => %s" % (proc['source'], proc['dest'])
            if proc['proc'].poll() is None:
                proc['proc'].terminate()
            proc['ret'] = proc['proc'].returncode
            finished_jobs[pid] = proc
        print " ===== Summary ===== "
        good_jobs = [job for job in finished_jobs.values() if job['ret'] == 0]
        bad_jobs = [job for job in finished_jobs.values() if job['ret'] != 0]
        for proc_info in good_jobs:
            print "%s => %s OK" % (proc_info['source'], proc_info['dest'])
        for proc_info in bad_jobs:
            print "%s => %s FAILED!" % (proc_info['source'], proc_info['dest'])

def needs_local_copy(castor_files, verbose=False):
    ''' Yield files from castor files that aren't current on local disk'''
    for file in castor_files:
        if not local_version_current(file):
            yield file

def expand_file_list(fileEntries):
     for fileEntry in fileEntries:
         if fileEntry.find("*") != -1:
             for file in castor.nsls(clean_name(fileEntry)):
                 yield "rfio:" + file
         else:
             yield fileEntry

def castor_files_in(sample):
    ''' Copy the files associated with this sample to the local drive '''
    for subsample in sample.subsamples:
        for file in expand_file_list(subsample.files):
            if is_on_castor(file):
                yield clean_name(file)

def mirror_sample(sample, max_jobs=20):
    ''' Mirror files associated witha  sample'''
    mirror_files(needs_local_copy(castor_files_in(sample)), max_jobs)

def mirror_samples(samples, max_jobs=20):
    ''' Mirror all the files associated to a list of samples '''
    sample_files = [needs_local_copy(castor_files_in(sample))
                    for sample in samples]
    mirror_files(itertools.chain(*sample_files), max_jobs)

def update_sample_to_use_local_files(sample, tiny_mode=False):
    ''' Update a sample to use any available local files '''
    if tiny_mode:
        print "Warning, tiny mode is on! Only taking 1 file"
    for subsample in sample.subsamples:
        count = 0
        new_file_list = []
        local_count = 0
        for input_file in expand_file_list(subsample.files):
            count += 1            
            if tiny_mode:
                if count > 1:
                    break                
            if is_on_castor(input_file):
                # strip rfio:
                clean_file = clean_name(input_file)
                # Check if it is cached locally
                if local_version_current(clean_file):
                    local_count += 1
                    input_file = "file:%s" % local_version(clean_file)
                else:
                    # Request CASTOR stage this file
                    stage_files([clean_file])                    
            new_file_list.append(input_file)
        # Update the files for this sample
        subsample.files = new_file_list
        print "Subsample %s has (%i/%i) files cached locally" % (
            subsample.name, local_count, len(subsample.files))

