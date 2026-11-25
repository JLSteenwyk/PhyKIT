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
