# PhyKIT performance baseline

Baseline source is `cf8c76525b634983995860f6594dff30b429f687` (PhyKIT
2.9.0). Run the commands from the repository root with the dependencies
installed in an active virtual environment.

- Ubuntu 24.04.2 LTS, Linux 6.8.0-55-generic, x86-64
- AMD Ryzen Threadripper PRO 7995WX, 192 logical CPUs
- CPython 3.12.14, Biopython 1.88, NumPy 2.5.5, SciPy 1.18.1,
  Matplotlib 3.11.2, scikit-learn 1.9.1, Numba 0.68.1, umap-learn 0.6.0,
  tqdm 4.70.1
- PhyKIT, OpenBLAS, OpenMP, MKL, and Numba each limited to one thread

## Bundled workloads

```sh
env PHYKIT_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  MKL_NUM_THREADS=1 NUMBA_NUM_THREADS=1 \
  python benchmarks/profile_commands.py --repeat 7 --warmup 1
```

Each sample starts a fresh Python process, and includes CLI startup, data
loading, computation, and any output written by the command. The warmup
primes the operating system's file cache. The `version` row estimates
Python and CLI startup under the same conditions. Standard output is
discarded but the command still formats and writes its output.

| Command | Median (s) | Min (s) | Max (s) |
| --- | ---: | ---: | ---: |
| `phylo_gwas_multigroup` | 0.990819 | 0.978012 | 1.030014 |
| `covarying_evolutionary_rates` | 0.152266 | 0.149890 | 0.158245 |
| `long_branch_score` | 0.042541 | 0.041815 | 0.044040 |
| `long_branch_score_verbose` | 0.042832 | 0.041992 | 0.044421 |
| `dvmc` | 0.035256 | 0.033515 | 0.036795 |
| `tip_labels` | 0.034449 | 0.033966 | 0.035898 |
| `version` | 0.030495 | 0.029143 | 0.030874 |
| `monophyly_check_true` | 0.027653 | 0.027300 | 0.028770 |
| `variable_sites` | 0.027533 | 0.026095 | 0.028416 |
| `pairwise_identity` | 0.026975 | 0.026213 | 0.028495 |
| `monophyly_check_false` | 0.027082 | 0.026384 | 0.027622 |
| `pairwise_identity_verbose` | 0.026929 | 0.026304 | 0.027954 |
| `gc_content_verbose` | 0.027508 | 0.025918 | 0.027950 |
| `compositional_bias_per_site` | 0.026411 | 0.026196 | 0.028307 |
| `evolutionary_rate_per_site` | 0.026952 | 0.025962 | 0.027514 |
| `parsimony_informative_sites` | 0.026271 | 0.026245 | 0.027973 |
| `alignment_entropy` | 0.026342 | 0.026086 | 0.028020 |
| `gc_content` | 0.026253 | 0.025995 | 0.027910 |
| `alignment_length_no_gaps` | 0.026103 | 0.025400 | 0.027592 |
| `treeness` | 0.026114 | 0.024403 | 0.026766 |
| `terminal_branch_stats` | 0.025351 | 0.024757 | 0.026887 |
| `patristic_distances` | 0.025055 | 0.024410 | 0.026699 |
| `composition_per_taxon` | 0.025530 | 0.023964 | 0.025917 |
| `column_score` | 0.025073 | 0.024132 | 0.026011 |
| `internal_branch_stats` | 0.024590 | 0.024261 | 0.026249 |
| `sum_of_pairs_score` | 0.024688 | 0.024426 | 0.025589 |
| `alignment_length` | 0.024534 | 0.024364 | 0.025199 |
| `bipartition_support_stats` | 0.024594 | 0.024508 | 0.024996 |
| `total_tree_length` | 0.023947 | 0.023786 | 0.025790 |
| `evolutionary_rate` | 0.023889 | 0.023764 | 0.024744 |

## Scaled workloads

```sh
python benchmarks/profile_scaled_workloads.py --repeat 5 \
  --gnu-time /usr/bin/time --profile-dir /tmp/phykit-scaled-prof
```

The synthetic workloads use 32 or 1,024 taxa. Pairwise identity takes a
300- or 1,000-site alignment respectively. Patristic distances and
phylogenetic signal use balanced Newick trees, with 100 permutations for
the signal test. Inputs are generated outside the measured process. Each
row uses five fresh CLI processes after an unmeasured warmup. GNU time
reports maximum resident set size for the CLI, which does not represent
the combined concurrent memory of a process pool. The one- and four-thread
settings are caps passed through `--threads`. A workload may use fewer
threads. The optional cProfile captures are separate from reported timings.

| Command | Taxa | Thread cap | Median (s) | Min (s) | Max (s) | Median peak RSS (MiB) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `pairwise_identity` | 32 | 1 | 0.1001 | 0.0994 | 0.1035 | 28.5 |
| `pairwise_identity` | 32 | 4 | 0.1053 | 0.1011 | 0.1080 | 30.0 |
| `patristic_distances` | 32 | 1 | 0.1056 | 0.1045 | 0.1068 | 30.0 |
| `patristic_distances` | 32 | 4 | 0.1030 | 0.1000 | 0.1046 | 30.0 |
| `phylogenetic_signal` | 32 | 1 | 0.3019 | 0.3012 | 0.3029 | 57.0 |
| `phylogenetic_signal` | 32 | 4 | 0.3028 | 0.3026 | 0.3052 | 57.0 |
| `pairwise_identity` | 1024 | 1 | 0.4781 | 0.4726 | 0.4807 | 66.0 |
| `pairwise_identity` | 1024 | 4 | 0.4745 | 0.4719 | 0.4807 | 66.0 |
| `patristic_distances` | 1024 | 1 | 0.3684 | 0.3622 | 0.3826 | 58.1 |
| `patristic_distances` | 1024 | 4 | 0.3741 | 0.3652 | 0.3789 | 58.1 |
| `phylogenetic_signal` | 1024 | 1 | 0.5006 | 0.4976 | 0.5029 | 108.8 |
| `phylogenetic_signal` | 1024 | 4 | 0.4330 | 0.4305 | 0.4339 | 107.4 |
