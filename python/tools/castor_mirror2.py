import os
import shlex
import sys
import subprocess
import threading
import Queue
import re
import time
import itertools
import TauAnalysis.Configuration.tools.castor as castor

# Where to store the files
LOCAL_DIRECTORY = [
    #"/tmp/abdollah/tau_fakerate_ntuples",
    #"/tmp/tau_commissioning",
    "/data2/friis/tau_fakerate_ntuples"
    #"/data2/veelken/CMSSW_4_1_x/ntuples/TauIdEffMeas/2011Jun10"
]

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

def local_version(castor_file, local_directory = LOCAL_DIRECTORY):
    ''' Map a castor file to a local file '''
    # Return the first directory where it is already local
    # Otherwise return last restul
    result = None
    for dir in local_directory:
        local_file = _CASTOR_MATCHER.sub(dir, castor_file)
        if os.path.exists(local_file):
            return local_file
        result = local_file
    return result

def local_version_current(castor_file, local_directory = LOCAL_DIRECTORY):
    ''' Check if the local copy of [castor_file] exists and is up to date '''
    local_file = local_version(castor_file, local_directory)
    if not os.path.exists(local_file):
        return False
    local_stat = os.stat(local_file)
    # Get last mod time of local file
    #local_mtime = time.ctime(local_stat.st_mtime)
    local_mtime = time.localtime(local_stat.st_mtime)
    local_size = local_stat.st_size
    # This call is memoized
    castor_stat = list(castor.nslsl(castor_file))[0]
    castor_size = castor_stat["size"]
    #castor_mtime = time.mktime(
    #    unixtime_from_timestamp(castor_stat["Last modify"]))
    castor_mtime = castor_stat['time']
    #print local_mtime, castor_mtime
    # Check sizes are same
    if local_size != castor_size:
        print "Local copy of", castor_file, " is the wrong size: %i != %i"% (local_size, castor_size)
        return False
    # Check local file is newer
    if local_mtime < castor_mtime:
        print "Local copy of", castor_file, " is outdated!"
        print "local:", time.asctime(local_mtime), \
                "castor:", time.asctime(castor_mtime)
        return False
    return True

class CopyWorker(threading.Thread):
    def __init__(self, work_queue, results_queue, local_directory, skipWaiting=False):
        super(CopyWorker, self).__init__()
        self.work_queue = work_queue
        self.results_queue = results_queue
        self.local_directory = local_directory
        self.skipWaiting = skipWaiting
    def run(self):
        while True:
            try:
                castor_file = self.work_queue.get()
                if castor_file is None:
                    break
                self.mirror(castor_file, self.local_directory, self.skipWaiting)
            except KeyboardInterrupt:
                sys.exit(2)
            finally:
                self.work_queue.task_done()

    def mirror(self, castor_file, local_directory, skipWaiting=False):
        ''' Initiate copy a castor file to the local directory.

        Returns a dictionary containing the background subprocess
        and information about the file.
        '''
        # Get file path for new file
        local_path = local_version(castor_file, local_directory)
        local_dirname = os.path.dirname(local_path)
        # Create the directories if necessary
        if not os.path.isdir(local_dirname):
            print "Creating local directory %s" % local_dirname
            os.makedirs(local_dirname, 0777)
        command = ['nice', '--adjustment=15', 'rfcp', castor_file, local_path]
        # CV: check if users run certain interactive processes on local machine;
        #     postpone starting 'rfcp' command in case they do
        jobsToWaitFor = {
            'squires' : [
                'gdb',
                'root',
                ##'screen',
                'xemacs',
                'vi',
                ##'vim'
            ],
            'vasquez' : [
                'root',
                'xemacs'
            ],
            'veelken' : [
                ##'root',
                ##'xemacs'
            ]
        }
        waitFor = 300 # time to wait in seconds
        waitForJobToFinish = True
        while waitForJobToFinish:
            waitForJobToFinish = False
            args = shlex.split('ps aux')
            running_jobs = subprocess.Popen(args, stdout = subprocess.PIPE)
            running_jobs.wait()
            running_jobs = running_jobs.stdout.readlines()
            for running_job in running_jobs:
                columns = running_job.split()
                if len(columns) < 11:
                    continue
                running_job_user = columns[0]
                running_job_executable = columns[10]
                for user, executables in jobsToWaitFor.items():
                    for executable in executables:
                        if user == running_job_user and executable == running_job_executable:
                            print "User %s is running %s --> waiting for %i seconds." % (user, executable, waitFor)
                            waitForJobToFinish = True
            current_time = time.localtime()
	    if current_time.tm_hour >= 22 or current_time.tm_hour <= 8:
	        print "Don't sleep during nighttime!"
                waitForJobToFinish = False
            if skipWaiting:
                print "Don't wait, it's urgent!"
                waitForJobToFinish = False
            if waitForJobToFinish:
                time.sleep(waitFor)        
        print "Requesting %s -> %s" % (castor_file, local_path)
        result = subprocess.call(command)
        #result = 0
        self.results_queue.put((castor_file, result))

def print_results(results_queue):
    while True:
        try:
            results = results_queue.get()
            if results is None:
                break
            castor_file, result = results
            print " finished (%s) file: %s" % (result, castor_file)
        finally:
            results_queue.task_done()

def mirror_files(castor_files, local_directory = LOCAL_DIRECTORY, max_jobs=10, skipWaiting=False):
    ''' Copy [castor_files] to local disk.

    max_jobs specifies maximum number of concurrent copies
    '''
    work_queue = Queue.Queue()
    results_queue = Queue.Queue()
    print "Building %i worker threads" % max_jobs
    workers = [CopyWorker(work_queue, results_queue, local_directory, skipWaiting) for i in range(max_jobs)]
    # Start all our workers
    map(CopyWorker.start, workers)
    results_worker = threading.Thread(
        target = lambda: print_results(results_queue))
    results_worker.start()
    # add our copy jobs to the queue
    for castor_file in castor_files:
        work_queue.put(castor_file)
    # Add poison pills
    for worker in workers:
        work_queue.put(None)
    work_queue.join()
    results_queue.put(None)
    results_queue.join()

def needs_local_copy(castor_files, local_directory = LOCAL_DIRECTORY, verbose=False):
    ''' Yield files from castor files that aren't current on local disk'''
    for file in castor_files:
        if not local_version_current(file, local_directory):
            yield file

def expand_file_list(fileEntries):
     for fileEntry in fileEntries:
         if fileEntry.find("*") != -1:
             for file in castor.nslsl(clean_name(fileEntry)):
                 yield "rfio:" + file['path']
         else:
             yield fileEntry

def castor_files_in(sample):
    ''' Copy the files associated with this sample to the local drive '''
    for subsample in sample.subsamples:
        for file in expand_file_list(subsample.files):
            if is_on_castor(file):
                yield clean_name(file)

def mirror_sample(sample, local_directory = LOCAL_DIRECTORY, max_jobs=10):
    ''' Mirror files associated witha  sample'''
    mirror_files(needs_local_copy(castor_files_in(sample), local_directory), local_directory, max_jobs)

def mirror_samples(samples, local_directory = LOCAL_DIRECTORY, max_jobs=10):
    ''' Mirror all the files associated to a list of samples '''
    sample_files = [needs_local_copy(castor_files_in(sample), local_directory)
                    for sample in samples]
    mirror_files(itertools.chain(*sample_files), local_directory, max_jobs)

def update_sample_to_use_local_files(sample, local_directory = LOCAL_DIRECTORY, tiny_mode=False, only_local=False):
    ''' Update a sample to use any available local files '''
    if tiny_mode:
        print "Warning, tiny mode is on! Only taking 1 file"
    for subsample in sample.subsamples:
        count = 0
        new_file_list = []
        local_count = 0
        skipped_count = 0
        for input_file in expand_file_list(subsample.files):
            #print input_file
            count += 1
            if tiny_mode:
                if count > 1:
                    break
            if is_on_castor(input_file):
                #print "is on castor"
                # strip rfio:
                clean_file = clean_name(input_file)
                # Check if it is cached locally
                if local_version_current(clean_file, local_directory):
                    #print "is current"
                    local_count += 1
                    input_file = "file:%s" % local_version(clean_file, local_directory)
                else:
                    if only_local:
                        skipped_count += 1
                        continue
                    # Request CASTOR stage this file
                    #stage_files([clean_file])
            new_file_list.append(input_file)
        # Update the files for this sample
        subsample.files = new_file_list
        print "Subsample %s has (%i/%i) files cached locally" % (
            subsample.name, local_count, len(subsample.files)+skipped_count)
        if only_local and skipped_count:
            print "WARNING: the [only_local] option is enabled,"\
                    " and %i files have been skipped!" % skipped_count

