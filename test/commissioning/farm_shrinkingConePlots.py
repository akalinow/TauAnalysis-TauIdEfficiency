import time
import subprocess

from shrinkingConePlots import distributions
from shrinkingConePlots import efficiencies
from shrinkingConePlots import efficiency_versus


MAX_JOBS = 8
TO_SUBMIT = {}
CURRENT_JOBS = {}
FINISHED_JOBS = {}

desired_distributions = ['jetPt', 'jetEta', 'decayMode', 'decayMode_pt40']
desired_efficiencies = ['TaNCOne', 'TaNCHalf', 'TaNCQuarter']
desired_efficiencies_versus = ['Pt', 'Eta']

for dist in distributions.keys():
    if dist not in desired_distributions:
        continue
    TO_SUBMIT[dist] = {
        'args' : ['python', 'shrinkingConePlots.py', '-p', dist],
        'ret' : -1
    }

for eff in efficiencies.keys():
    if eff not in desired_efficiencies:
        continue
    for eff_vs in efficiency_versus.keys():
        if eff_vs not in desired_efficiencies_versus:
            continue
        TO_SUBMIT["%s_%s" % (eff, eff_vs)] = {
            'args' : ['python', 'shrinkingConePlots.py',  '-p', eff, '-v', eff_vs ],
            'ret' : -1
        }
    
def wait_for_jobs(current_jobs, finished_jobs, max_jobs):
    while len(current_jobs.keys()) >= max_jobs:
        time.sleep(1.0)
        for job in current_jobs.keys():
            job_info = current_jobs[job]
            ret = job_info['proc'].poll()
            job_info['ret'] = ret
            if ret is not None:
                print job, "is done with", ret
                finished_jobs[job] = current_jobs[job]
                del current_jobs[job]

def submit_jobs(to_submit, current_jobs, finished_jobs, max_jobs = MAX_JOBS):
    # Run until all are submitted
    for plot_name in to_submit.keys():
        plot_info = to_submit[plot_name]
        print "Submitting", plot_name
        plot_info['proc'] = subprocess.Popen(
            plot_info['args'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        current_jobs[plot_name] = plot_info
        del to_submit[plot_name]
        wait_for_jobs(current_jobs, finished_jobs, max_jobs)

if __name__ == "__main__":
    submit_jobs(TO_SUBMIT, CURRENT_JOBS, FINISHED_JOBS)


