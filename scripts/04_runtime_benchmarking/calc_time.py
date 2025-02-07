import sys

import time
import subprocess
import numpy as np
import pandas as pd
from memory_profiler import memory_usage

def run_bac2feature(test_seq_path, method) -> None:
    cmd = []
    if method == 'homology':
        # Homology based prediction
        cmd = ['bac_to_feature.py',
               '-s', test_seq_path,
               '-o', '../../data/time_calculation/b2f_output_homology.tsv',
               '-m', 'homology',
               '--threads', '1']
    elif method == 'phylogeny':
        # Phylogeny based prediction
        cmd = ['bac_to_feature.py',
               '-s', test_seq_path,
               '-o', '../../data/time_calculation/b2f_output_phylogeny.tsv',
               '-m', 'phylogeny',
               '--threads', '1']
    elif method == 'taxonomy':
        # Taxonomy based prediction
        cmd = ['bac_to_feature.py',
               '-s', test_seq_path,
               '-o', '../../data/time_calculation/b2f_output_taxonomy.tsv',
               '-m', 'taxonomy',
               '--threads', '1']
    subprocess.run(cmd)
    return

nseqs = np.arange(start=1000, stop=10000 + 1000, step=1000)

n_threads = 1
args = sys.argv
method = args[1]

# Run process
l_wall_time, l_cpu_time, l_mem_usage = [], [], []
# for n in nseqs:
for n in nseqs:

    test_seq_path = f"../../data/time_calculation/test/test_seqs_{n}.fasta"
    test_qza_path = f"../../data/time_calculation/test/test_seqs_{n}.qza"

    # Start
    wall_start_time = time.perf_counter()
    cpu_start_time = time.process_time()

    # Heavy process
    mem_usage = memory_usage((run_bac2feature, (test_seq_path, method), {}), interval=1, timeout=100)

    # End
    wall_end_time = time.perf_counter()
    cpu_end_time = time.process_time()

    # Time
    wall_time = wall_end_time - wall_start_time
    cpu_time = cpu_end_time - cpu_start_time

    l_wall_time += [wall_time]
    l_cpu_time += [cpu_time]
    l_mem_usage += [np.max(mem_usage)]

result_df = pd.DataFrame(data={'nseqs': nseqs, 'wall_time': l_wall_time, 'cpu_time': l_cpu_time, 'memory': l_mem_usage})
result_df.to_csv(f"../../data/time_calculation/time_{method}.tsv", sep='\t', index=False)
