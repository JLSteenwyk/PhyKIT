# Performance Benchmarks

Benchmarks are local microbenchmarks for optimized hot paths. They use fixed
random seeds and synthetic alignments to make before/after comparisons stable
enough for regression tracking; absolute timings will vary by machine.

## 2026-06-22

Environment:

- Python 3.11
- NumPy from the active local environment
- macOS development workstation

Results:

| Function | Dataset | Baseline mean | Optimized mean | Speedup |
| --- | --- | ---: | ---: | ---: |
| `calculate_summary_statistics_from_arr` | 1M floating-point values | 1.1086s | 0.0448s | 24.8x |
| `calculate_summary_statistics_from_arr` combined quantiles | 1M floating-point values, optimized helper baseline | 0.0327s | 0.0187s | 1.75x |
| `calculate_summary_statistics_from_arr` variance/std reuse | 1M floating-point values, optimized helper baseline | 0.005185s | 0.004021s | 1.29x |
| `calculate_summary_statistics_from_arr` constant-array shortcut | 1M identical floating-point values, side-by-side previous percentile/variance path | 0.049040s | 0.000533s | 92.01x |
| `calculate_summary_statistics_from_arr` ndarray mean reduction | 10 / 100 / 1000 / 100k / 1M floating-point values, side-by-side previous `np.mean` wrapper | 0.000008843s / 0.000008838s / 0.000007437s / 0.000078838s / 0.000916197s | 0.000006742s / 0.000005305s / 0.000006767s / 0.000061452s / 0.000516858s | 1.31x / 1.67x / 1.10x / 1.28x / 1.77x |
| `calculate_summary_statistics_from_arr` ndarray extrema/variance reductions | 5 / 50 / 1000 / 100k floating-point values, side-by-side previous `np.min`/`np.max`/`np.var` wrappers | 2.558751s / 1.764250s / 0.233506s / 0.628135s | 1.982603s / 0.715148s / 0.220886s / 0.295009s | 1.29x / 2.47x / 1.06x / 2.13x |
| `calculate_summary_statistics_from_dict` | 1M floating-point dictionary values | 1.1465s | 0.0539s | 21.3x |
| `calculate_summary_statistics_from_dict` fromiter setup | 1M floating-point dictionary values, optimized helper baseline | 0.0395s | 0.0276s | 1.4x |
| `calculate_summary_statistics_from_dict` combined percentiles | 1M floating-point dictionary values, optimized helper baseline | 0.0373s | 0.0356s | 1.05x |
| `stats_summary` module import without eager NumPy | cold subprocess import of shared summary-statistics helper | 0.118467s | 0.061148s | 1.94x |
| `print_summary_statistics` batched output | 100k captured summary reports, identical stdout text | 0.437117s | 0.342506s | 1.28x |
| `print_summary_statistics` template formatting | 100k captured summary reports, identical stdout text, side-by-side previous f-string join formatter | 1.132736s | 0.931583s | 1.22x |
| `stats_summary._print_no_values_message` batched output | 100k captured no-values diagnostics, identical stdout text | 0.063022s | 0.030512s | 2.07x |
| `stats_summary` sized no-value early return | 20k empty/single array and dictionary summary requests, identical diagnostic text | 0.023061s / 0.039072s / 0.039999s / 0.073569s | 0.008748s / 0.007520s / 0.005851s / 0.007079s | 2.64x / 5.20x / 6.84x / 10.39x |
| `print_json` builtin payload serialization | 100k row dictionaries with nested builtin lists | 0.268362s | 0.062000s | 4.3x |
| `print_json` mixed builtin/NumPy default hook | 100k builtin row dictionaries plus NumPy scalar summary values, identical serialized JSON | 0.541561s | 0.109570s | 4.94x |
| `to_builtin_json_types` mixed builtin/NumPy payload | 300k builtin row dictionaries plus one NumPy scalar summary, identical converted payload | 1.133849s | 0.801193s | 1.42x |
| `to_builtin_json_types` builtin scalar fast path | 300k builtin row dictionaries plus one NumPy scalar summary, side-by-side previous scalar NumPy provenance check | 0.526141s | 0.376409s | 1.40x |
| `to_builtin_json_types` copy-on-write containers | 300k builtin row dictionaries plus one NumPy scalar summary, side-by-side previous eager container allocation | 0.402334s | 0.362144s | 1.11x |
| `json_output` module import without eager NumPy | cold subprocess import of shared JSON output helper | 0.088862s | 0.058492s | 1.52x |
| `json_output` module import without eager stdlib JSON | cold subprocess import of shared JSON output helper, serialization still imports JSON on demand | 0.019455s | 0.018319s | 1.06x |
| `helpers.caching` module import without eager serialization/cache setup | median cold subprocess import after lazy pickle/json and lazy global caches | 0.012563s | 0.004362s | 2.88x |
| `helpers.caching` module import without `typing` startup | median cold subprocess import after replacing annotation-only `Any`/`Callable` aliases with built-in annotations | 0.025169s | 0.023913s | 1.05x |
| `helpers.caching` module import without eager `hashlib` | cold subprocess import of shared caching helper, cache-key generation still imports hashing on demand | 0.021160s | 0.018515s | 1.14x |
| `ResultCache._get_cache_key` cached md5 helper | 200k repeated primitive-argument cache-key generations with kwargs, identical keys | 1.540026s | 1.143370s | 1.35x |
| `ResultCache.clear` scandir cleanup loop | 1M fake cache directory entries, half `.pkl`, identical removed paths while isolating Python loop overhead | 0.663190s | 0.444320s | 1.49x |
| `trait_parsing` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in annotations | 0.022946s | 0.020371s | 1.13x |
| `trait_parsing.parse_multi_trait_file` single-pass parser | 300k-row multi-trait TSV, 3 numeric trait columns, 100k shared taxa | 0.727592s | 0.634490s | 1.15x |
| `trait_parsing.parse_multi_trait_file` valid-row float conversion | 300k-row multi-trait TSV, 3 numeric trait columns, 100k shared taxa | 0.634490s | 0.537817s | 1.18x |
| `trait_parsing.parse_multi_trait_file` all-shared parser fast path | 300k-row multi-trait TSV, 3 numeric trait columns, all taxa shared | 0.423491s | 0.277321s | 1.53x |
| `trait_parsing.parse_multi_trait_file` narrow-row float conversion | 300k-row multi-trait TSV, 3 numeric trait columns, all taxa shared, identical parsed traits and nonnumeric fallback messages | 0.367998s | 0.331910s | 1.11x |
| `trait_parsing.trait_column_from_rows` direct column vector | 120k taxa x 12 parsed trait columns, one selected trait column | 0.013424s | 0.007369s | 1.82x |
| `trait_parsing.trait_matrix_from_rows` direct row matrix | 120k taxa x 12 numeric traits, ordered trait rows to NumPy matrix | 0.137167s | 0.038031s | 3.61x |
| `trait_parsing.response_predictor_arrays` selected-column design matrix | 180k taxa x 10 parsed trait columns, one response plus four predictors | 0.147410s | 0.060780s | 2.43x |
| `trait_parsing.response_predictor_arrays` narrow design matrix | 180k taxa x 20 parsed trait columns, one response plus one predictor | 0.063154s | 0.024249s | 2.60x |
| Regression-style trait-name index resolution | 40k parsed trait columns, one response plus 800 predictors, first duplicate index preserved | 0.461402s | 0.002792s | 165.26x |
| `trait_parsing.subset_traits_to_ordered_shared_taxa` ordered shared comparison | 300k ordered taxa, all shared / 150k shared subset, avoiding set allocation | 0.026659s / 0.034705s | 0.000216s / 0.015925s | 123.5x / 2.18x |
| `phykit.phykit` CLI startup without registry `typing` imports | median cold subprocess import after converting registry/factory annotation-only aliases to built-in annotations | 0.035697s | 0.032062s | 1.11x |
| `Phykit.__init__` normal-command dispatch without top-level help parser | mocked `alignment_length` command dispatch, parser construction bypassed | 0.000442517s | 0.000000397s | 1114.7x |
| `Phykit._dispatch_command` single normal-command lookup | 2M mocked `alignment_length` dispatches, side-by-side previous `hasattr` plus `getattr` path | 1.701180s | 1.476757s | 1.15x |
| `phykit.phykit` CLI startup without eager plot config helper | cold subprocess import after wrapping `add_plot_arguments`; plot commands still import on parser construction | 0.029375s | 0.024250s | 1.21x |
| `phykit.phykit` CLI startup without eager logging setup | cold subprocess import after removing unused module logger/handler construction | 0.024059s | 0.020968s | 1.15x |
| `phykit.phykit` CLI startup without eager `argparse` | cold subprocess import after deferring parser imports until `_new_parser()` and lazy-loading boolean parse errors | 0.020950s | 0.019993s | 1.05x |
| `boolean_argument_parsing` import without eager `argparse` | cold subprocess import of `str2bool`; invalid values still raise `argparse.ArgumentTypeError` | 0.019988s | 0.018408s | 1.09x |
| `geological_timescale.get_timescale_for_range` epoch intervals | 300k auto epoch-range calls, side-by-side previous overwritten epoch-list build | 1.001253s | 0.425282s | 2.35x |
| `geological_timescale.get_timescale_for_range` period early stop | 160k auto period-range calls, identical included period bands | 0.191914s | 0.170999s | 1.12x |
| `geological_timescale.get_timescale_for_range` epoch early stop | 300k mixed auto epoch-range calls, identical included epoch bands, side-by-side previous full epoch-table scan | 1.779906s | 0.629607s | 2.83x |
| `phykit.phykit` CLI startup without eager `textwrap` | cold subprocess import after routing parser-description dedent calls through a lazy `_dedent()` helper | 0.020041s | 0.019360s | 1.04x |
| `Phykit.run_alias` invalid-command cached banner | 3k repeated invalid-alias banner format calls, identical banner text | 0.544822s | 0.000186s | 2931.12x |
| `Phykit.version` cached banner | 3k repeated version banner format calls, identical version text | 0.540932s | 0.000232s | 2327.00x |
| `IdentityMatrix._compute_identity_matrix` | 180 taxa x 1200 sites, alphabet `ACGT-?NX*` | 1.7554s | 0.1377s | 12.7x |
| `IdentityMatrix._compute_identity_matrix` matrix products | 700 taxa x 900 sites, alphabet `ACGT-?NX*` | 1.7909s | 0.8584s | 2.1x |
| `IdentityMatrix._compute_identity_matrix` BLAS-backed float products | 700 taxa x 900 sites, alphabet `ACGT-?NX*` | 0.829083s | 0.069862s | 11.9x |
| `IdentityMatrix` clustered heatmap linkage setup | 1200 taxa, precomputed identity matrix and `--sort cluster` | 0.0311s | 0.0164s | 1.9x |
| `IdentityMatrix._determine_order` tree terminal-name extraction | parsed balanced 65536-tip tree, `--sort tree` terminal order | 0.1344s | 0.0126s | 10.7x |
| `IdentityMatrix._compute_partition_identities` | 220 taxa x 5000 sites split into 20 partitions, alphabet `ACGT-?NX*` | 2.9923s | 0.5447s | 5.5x |
| `IdentityMatrix._compute_partition_identities` precomputed valid float mask | 260 taxa x 6000 sites split into 30 partitions, alphabet `ACGT-?NX*`, identical partition identity array | 0.139126s | 0.128939s | 1.08x |
| `IdentityMatrix._compute_partition_identities` condensed partition extraction | 700 taxa x 30 partitions, representative compared/matches matrices, side-by-side previous reused `np.triu_indices` extraction | 0.839308s | 0.487465s | 1.72x |
| `IdentityMatrix._compute_partition_identities` identical-sequence shortcut | 700 taxa x 900 identical DNA sites split into 30 partitions with gaps/ambiguous symbols, side-by-side previous matrix-product path | 1.533377s | 0.010975s | 139.72x |
| `IdentityMatrix._compute_identity_matrix` clean ASCII direct comparison | 260 taxa x 6000 sites, alphabet `ACGT`, side-by-side previous validity-mask matrix-product path | 0.178639s | 0.122918s | 1.45x |
| `IdentityMatrix._compute_identity_matrix` clean ASCII direct threshold | 150 x 250 / 150 x 500 / 200 x 1000 / 300 x 250 clean DNA matrices, side-by-side previous validity-mask matrix-product path under the old length-only threshold | 0.127482s / 0.280151s / 0.633633s / 0.337570s | 0.002105s / 0.029026s / 0.022576s / 0.040665s | 60.55x / 9.65x / 28.07x / 8.30x |
| `IdentityMatrix._compute_identity_matrix` identical-sequence shortcut | 700 taxa x 900 identical DNA sites with gaps/ambiguous symbols, side-by-side previous matrix-product path | 0.048009s | 0.000289s | 166.26x |
| `IdentityMatrix._all_sequences_identical` no-slice taxon scan | 1M taxon-name sequence dictionary, identical / early-different / late-different cases, side-by-side previous `taxa_names[1:]` shortcut predicate | 1.024048s / 0.011813s / 0.659557s | 0.581272s / 0.000000167s / 0.630424s | 1.76x / 70736.53x / 1.05x |
| `IdentityMatrix._identical_sequence_identity_value` late valid character | 20 repeated 3.5M-site identical DNA sequences with final valid site, side-by-side previous per-character membership scan | 18.405654s | 0.083401s | 220.69x |
| `IdentityMatrix._identical_sequence_identity_value` all invalid characters | 20 repeated 3.5M-site identical DNA sequences with only invalid symbols, side-by-side previous per-character membership scan | 20.856280s | 0.428930s | 48.62x |
| `IdentityMatrix` pairwise fallback boolean counts | 100k-site Unicode fallback pair comparison, side-by-side previous boolean `np.sum` valid/match counts | 0.000246s | 0.000206s | 1.19x |
| `IdentityMatrix._partition_identity_strip` | 350 taxa x 24 partitions, 61075 pairwise partition identities | 2.2786s | 0.0229s | 99.7x |
| `IdentityMatrix._partition_identity_strip` condensed row scan | 1200 taxa x 32 partitions, 719400 pairwise partition identities, side-by-side previous `np.triu_indices` plus `np.add.at` accumulation | 2.144036s | 0.117468s | 18.25x |
| `IdentityMatrix._parse_partitions` partition parser | 500k RAxML-style partition rows with comments/blanks, comma/no-comma forms, and whitespace tolerance | 0.644653s | 0.564478s | 1.14x |
| `IdentityMatrix._summarize_identity_matrix` condensed upper triangle | 3500 x 3500 symmetric identity matrix, identical mean/min/max values and min/max taxon pairs | 0.120000s | 0.019253s | 6.23x |
| `IdentityMatrix._summarize_identity_matrix` condensed-vector reductions | 700 / 1200 / 2500 / 3500-taxon symmetric identity matrices, side-by-side previous generic `np.mean`/`np.min`/`np.max`/`np.argmin`/`np.argmax` reductions after `squareform` | 0.001703s / 0.006129s / 0.020625s / 0.048776s | 0.000981s / 0.005080s / 0.016065s / 0.036241s | 1.74x / 1.21x / 1.28x / 1.35x |
| `IdentityMatrix._print_text_output` batched summary | 100k captured identity-matrix text summaries, identical stdout text | 0.267546s | 0.186169s | 1.44x |
| `IdentityMatrix._print_text_output` template formatting | 100k captured identity-matrix text summaries, identical stdout text, side-by-side previous f-string join formatter | 0.186500s | 0.149763s | 1.25x |
| `identity_matrix` module import without eager SciPy clustering | cold process import for identity-matrix command module | 0.459178s | 0.249947s | 1.8x |
| `IdentityMatrix` cached SciPy clustering wrappers | 2k repeated 16-taxon squareform/linkage/leaves-list cluster ordering calls, SciPy already warm, side-by-side previous import-on-call wrappers | 0.200257s | 0.169233s | 1.18x |
| `helpers.files` module import without eager Bio.AlignIO | cold subprocess import of shared alignment file helper | 0.159080s | 0.060833s | 2.62x |
| `helpers.files` module import without `hashlib` startup | median cold subprocess import after localizing cache-key hashing to `_get_file_hash` | 0.003301s | 0.000753s | 4.38x |
| `helpers.files._get_file_hash` raw stat cache key | 100k cache-key generations for one alignment file | 0.000002628s | 0.000001821s | 1.44x |
| `helpers.files.is_protein_alignment` ASCII nucleotide detection | five scans of 50k records x 112 bp lowercase nucleotide alignment | 0.171807s | 0.049896s | 3.44x |
| `helpers.files.read_single_column_file_to_list` bulk splitlines trim | 1M single-column rows with leading/trailing spaces, identical stripped list | 0.072114s | 0.052382s | 1.38x |
| `helpers.files.read_single_column_file_to_list` streaming trim | 2M single-column rows with leading/trailing spaces, identical stripped list from real file I/O | 0.577303s | 0.491627s | 1.17x |
| `helpers.files._detect_format_by_content` PHYLIP header split | 1M mixed first-line headers, identical detected formats | 2.294802s | 1.625386s | 1.41x |
| `alignment.base` module import without eager Bio.AlignIO | cold subprocess import after lazy shared alignment reader | 0.154313s | 0.113593s | 1.36x |
| `alignment.base` module import without eager NumPy lookup tables | cold subprocess import after lazy RCV NumPy lookup construction | 0.076903s | 0.022717s | 3.39x |
| `alignment.base` module import without eager file helper | median cold subprocess import after lazy shared alignment-reader wrapper | 0.042439s | 0.032811s | 1.29x |
| `identity_matrix` module import without eager Bio.SeqIO/Phylo/AlignIO | cold subprocess import after lazy Biopython alignment/tree readers | 0.197827s | 0.119863s | 1.65x |
| `identity_matrix` module import without eager NumPy/json/plot config | cold subprocess import after lazy NumPy proxy plus JSON and plot config helper deferral | 0.090216s | 0.025531s | 3.53x |
| `identity_matrix` module import without tree-base/typing startup | median cold subprocess import after localizing tree-sort helper import and replacing annotation-only typing aliases | 0.005714s | 0.001001s | 5.71x |
| `PairwiseIdentity._process_pair_batch` with `exclude_gaps=True` | 160 taxa x 1000 sites, alphabet `ACGT-?NX*` | 0.5102s | 0.0710s | 7.2x |
| `PairwiseIdentity._process_pair_batch` `count_nonzero` identity kernel | 1600 fallback pair comparisons, 80 taxa x 3000 sites, alphabet `ACGT-?NX*`, `exclude_gaps=True`, side-by-side previous `np.sum` kernel | 0.015628s | 0.012908s | 1.21x |
| `PairwiseIdentity.calculate_pairwise_identities` equal-length exclude-gaps matrix path | 160 taxa x 1000 sites, alphabet `ACGT-?NX*`, multiprocessing disabled | 0.064975s | 0.051714s | 1.26x |
| `PairwiseIdentity.calculate_pairwise_identities` clean ASCII exclude-gaps shortcut | 400 taxa x 1000 sites, alphabet `ACGT`, multiprocessing disabled | 0.314763s | 0.130653s | 2.41x |
| `PairwiseIdentity.calculate_pairwise_identities` byte sequence arrays | 180 taxa x 1500 sites, alphabet `ACGT`, multiprocessing disabled | 0.0772s | 0.0655s | 1.2x |
| `PairwiseIdentity.calculate_pairwise_identities` default matrix path | 180 taxa x 1500 sites, alphabet `ACGT-?NX*`, multiprocessing disabled | 0.093963s | 0.027521s | 3.41x |
| `PairwiseIdentity.calculate_pairwise_identities` matrix output condensed summaries | synthetic identity-count matrix output assembly for 100 / 250 / 500 / 1000 / 1500 taxa, side-by-side previous matrix indexing plus dict summary | 0.013223s / 0.158875s / 0.497973s / 1.937952s / 5.230754s | 0.003542s / 0.019972s / 0.121561s / 0.636426s / 1.944651s | 3.73x / 7.95x / 4.10x / 3.05x / 2.69x |
| `PairwiseIdentity.run` summary-only matrix stats | 180 taxa x 1500 sites, alphabet `ACGT-?NX*`, non-verbose/no-plot output | 0.024218s | 0.008288s | 2.92x |
| `PairwiseIdentity.calculate_pairwise_identity_stats` clean ASCII exclude-gaps shortcut | 400 taxa x 1000 sites, alphabet `ACGT`, summary-only stats path | 0.263493s | 0.031507s | 8.36x |
| `PairwiseIdentity.calculate_pairwise_identity_stats` matrix block index reuse | 700 taxa x 900 sites, alphabet `ACGT-?NX*`, side-by-side previous per-block column index allocation | 0.131023s | 0.100191s | 1.31x |
| `PairwiseIdentity.calculate_pairwise_identity_stats` gappy matrix condensed extraction | 1800 taxa x 1200 sites synthetic identity-count matrix, side-by-side previous `np.triu_indices` stats extraction | 0.141980s | 0.058894s | 2.41x |
| `PairwiseIdentity.calculate_pairwise_identities` identical-sequence shortcut | 400 taxa x 1000 identical DNA sites with gaps/ambiguous symbols, `exclude_gaps=True`, lowercase/uppercase variants, side-by-side previous matrix path | 0.299976s | 0.014772s | 20.31x |
| `PairwiseIdentity.calculate_pairwise_identity_stats` identical-sequence shortcut | 400 taxa x 1000 identical DNA sites with gaps/ambiguous symbols, `exclude_gaps=True`, summary-only stats path | 0.297837s | 0.000242s | 1231.58x |
| `PairwiseIdentity.calculate_pairwise_identities` single-record early return | 4.5M-site single-record DNA alignment, side-by-side previous fallback sequence-array setup before no-pair result | 0.002617416s | 0.000002709s | 966.28x |
| `PairwiseIdentity.calculate_pairwise_identity_stats` single-record early return | 4.5M-site single-record DNA alignment, side-by-side previous full-result fallback before no-values stats | 0.004708167s | 0.000002458s | 1915.27x |
| `PairwiseIdentity._identity_for_identical_sequence` gap count | 20 repeated 4.5M-site identical DNA sequences, `exclude_gaps=True`, side-by-side previous per-character membership loop | 9.847190s | 0.367208s | 26.82x |
| `AlignmentLengthNoGaps`/`PairwiseIdentity` identical-row no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.259744s / 0.018993s / 0.082823s | 0.126212s / 0.000004s / 0.038411s | 2.06x / 5301.12x / 2.16x |
| `PairwiseIdentity.calculate_pairwise_identities` sequential fallback pair streaming | 500 mixed-length records, multiprocessing disabled, worker and stats helper held constant | 1.010443s | 0.595901s | 1.70x |
| `PairwiseIdentity.calculate_pairwise_identities` multiprocessing fallback streaming chunks | 2500 fallback records, 3,123,750 index-pair chunk setup, side-by-side previous full pair-list slicing | 0.720987s | 0.338193s | 2.13x |
| `PairwiseIdentity.run` verbose batched text output | 100k pair rows, mocked alignment/read and identical stdout text | 0.054454s | 0.040773s | 1.34x |
| `PairwiseIdentity.run` summary-only taxa-list elision | 1M mocked records, scoring and summary output mocked, non-verbose/no-plot output | 0.116345s | 0.000001s | 90068.29x |
| `PairwiseIdentity._plot_pairwise_identity_heatmap` canonical matrix fill | 1200 taxa, 719400 canonical pair identities, identical symmetric heatmap matrix with arbitrary-order fallback retained | 0.119255s | 0.055930s | 2.13x |
| `PairwiseIdentity._pairwise_identity_matrix_from_pairs` squareform canonical fill | 2500 taxa, 3123750 canonical pair identities, identical symmetric float32 heatmap matrix | 0.404774s | 0.057102s | 7.09x |
| `pairwise_identity` module import without eager SciPy clustering | cold process import for non-plot pairwise-identity command module | 0.423027s | 0.273540s | 1.5x |
| `pairwise_identity` module import without eager Bio.Align/tqdm | cold subprocess import after lazy annotation/progress imports | 0.208326s | 0.134796s | 1.55x |
| `pairwise_identity` module import without eager NumPy lookup tables | cold subprocess import after lazy class-level gap lookup construction | 0.088633s | 0.039021s | 2.27x |
| `pairwise_identity` module import without eager multiprocessing/stats/json/plot helpers | cold subprocess import after lazy helper wrappers and localized `partial` import | 0.039676s | 0.025821s | 1.54x |
| `pairwise_identity` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.036227s | 0.032552s | 1.11x |
| `VariableSites.calculate_variable_sites` | 220 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.1652s | 0.0621s | 2.7x |
| `VariableSites.calculate_variable_sites` byte matrix setup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.1697s | 0.0491s | 3.5x |
| `VariableSites.calculate_variable_sites` gap lookup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.048267s | 0.018903s | 2.55x |
| `VariableSites.calculate_variable_sites` ASCII protein min/max columns | 1200 taxa x 6000 sites, alphabet `ACDEFGHIKLMNPQRSTVWY-X?*` | 0.073575s | 0.032941s | 2.23x |
| `VariableSites.calculate_variable_sites` min/max valid-column sentinel | 1200 taxa x 6000 sites, alphabet `ACGT-?NX*`, side-by-side previous valid-count column pass | 0.021288s | 0.015788s | 1.35x |
| `VariableSites.calculate_variable_sites` gap-code mask construction | 1200 taxa x 12000 sites, alphabet `ACGT-?NX*`, side-by-side previous lookup-mask gather | 0.092339s | 0.070243s | 1.31x |
| `VariableSites.calculate_variable_sites` protein no-gap direct min/max | 1200 taxa x 10000 sites, 20 amino-acid symbols, side-by-side previous masked min/max path | 0.027451s | 0.015857s | 1.73x |
| `VariableSites.calculate_variable_sites` DNA no-gap direct min/max | 1200 taxa x 12000 clean DNA sites, side-by-side previous masked min/max path with identical variable-site result | 0.046271s | 0.017388s | 2.66x |
| `VariableSites.calculate_variable_sites` clean ASCII gap-byte precheck | 80 taxa x 1M clean DNA sites, side-by-side previous full invalid-mask setup | 2.584151s | 1.359966s | 1.90x |
| `VariableSites.calculate_variable_sites` identical-sequence shortcut | 1200 taxa x 12000 identical ASCII DNA sites, lowercase/uppercase variants, side-by-side previous matrix path | 0.014955s | 0.005998s | 2.49x |
| `VariableSites.calculate_variable_sites` raw-identical normalization scan | 300k raw-identical DNA rows, side-by-side previous eager uppercase sequence setup with zero variable sites | 0.182161s | 0.112015s | 1.63x |
| `VariableSites.calculate_variable_sites` identical-sequence no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.323234s / 0.005263s / 0.141766s | 0.034174s / 0.000003s / 0.038053s | 9.46x / 1884.87x / 3.73x |
| `VariableSites.calculate_variable_sites` single-record early return | 4.5M-site single-record DNA alignment, side-by-side previous sequence materialization before zero return | 0.004295166s | 0.000001084s | 3962.11x |
| `VariableSites.calculate_variable_sites` Unicode final variable-column count | 1M-site fallback valid-symbol count vector, side-by-side previous boolean `np.sum` final count | 0.000235s | 0.000063s | 3.75x |
| `variable_sites` module import without eager NumPy/Bio.Align | cold subprocess import after lazy NumPy lookup construction and annotation-only Bio.Align import | 0.119398s | 0.023496s | 5.08x |
| `variable_sites` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006383s | 0.005138s | 1.24x |
| `variable_sites` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.003050s | 0.000955s | 3.19x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` | 220 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.1827s | 0.0590s | 3.1x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` byte matrix setup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.1776s | 0.0575s | 3.1x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` gap lookup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.058787s | 0.028491s | 2.06x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` ASCII protein column histograms | 1200 taxa x 6000 sites, alphabet `ACDEFGHIKLMNPQRSTVWY-X?*` | 0.151177s | 0.107077s | 1.41x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` gap-code mask construction | 1200 taxa x 12000 sites, alphabet `ACGT-?NX*`, side-by-side previous lookup-mask gather | 0.269284s | 0.237955s | 1.13x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` clean ASCII block counting | 1200 taxa x 12000 clean DNA sites, side-by-side previous always-masked block counter with identical PI-site result | 0.157442s | 0.142948s | 1.10x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` clean ASCII gap-byte precheck | 80 taxa x 1M clean DNA sites, side-by-side previous full invalid-mask setup | 15.602653s | 14.708943s | 1.06x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` identical-sequence shortcut | 1200 taxa x 12000 identical ASCII DNA sites, lowercase/uppercase variants, side-by-side previous block-count path | 0.142054s | 0.005836s | 24.34x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` identical-row no-slice scan | 1M identical ASCII DNA rows, side-by-side previous `sequences[1:]` equality scan | 0.390887s | 0.279042s | 1.40x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` combined sequence-build identity scan | 1M mixed-symbol DNA rows, identical / late-different setup cases, side-by-side previous sequence-list build plus identity pass | 0.230124s / 0.285500s | 0.174743s / 0.217506s | 1.32x / 1.31x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` single-record early return | 4.5M-site single-record DNA alignment, side-by-side previous sequence materialization before zero return | 0.005767500s | 0.000001250s | 4614.61x |
| `ParsimonyInformative.calculate_parsimony_informative_sites` Unicode final PI-site count | 1M-site fallback recurrent-state count vector, side-by-side previous boolean `np.sum` final count | 0.000235s | 0.000063s | 3.73x |
| `ParsimonyInformative.get_number_of_occurrences_per_character` record-wise direct count loop | 200 sampled columns from 5000 taxa x 2000 sites, alphabet `ACGT-?NX*`, side-by-side previous column slicing path with identical `Counter` output | 0.834731s | 0.524695s | 1.59x |
| `ParsimonyInformative.is_parsimony_informative` early recurrent-state exit | 20k repeated checks over 1000 recurrent and 1000 singleton states, identical truth value | 0.930823s | 0.002989s | 311.41x |
| `parsimony_informative_sites` module import without eager NumPy/Bio.Align | cold subprocess import after lazy NumPy lookup construction and annotation-only Bio.Align import | 0.101964s | 0.023527s | 4.33x |
| `parsimony_informative_sites` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006218s | 0.005062s | 1.23x |
| `parsimony_informative_sites` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.002538s | 0.000957s | 2.65x |
| `AlignmentEntropy.calculate_site_entropies` | 220 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.1706s | 0.0549s | 3.1x |
| `AlignmentEntropy.calculate_site_entropies` byte matrix setup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.1702s | 0.0564s | 3.0x |
| `AlignmentEntropy.calculate_site_entropies` gap lookup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.064573s | 0.031327s | 2.06x |
| `AlignmentEntropy.calculate_site_entropies` protein ASCII counts | 1000 taxa x 5000 sites, 20 amino-acid symbols | 0.076837s | 0.051144s | 1.50x |
| `AlignmentEntropy.calculate_site_entropies` protein ASCII column-major indexing | 1200 taxa x 6000 sites, alphabet `ACDEFGHIKLMNPQRSTVWY-X?*` | 0.156230s | 0.128932s | 1.21x |
| `AlignmentEntropy._entropy_from_ascii_codes` all-valid protein blocks | 1200 taxa x 12000 sites, 20 amino-acid symbols, side-by-side previous boolean-index path | 0.100458s | 0.090770s | 1.11x |
| `AlignmentEntropy.calculate_site_entropies` clean protein gap-mask elision | 1200 taxa x 12000 clean protein sites, side-by-side previous full valid-mask setup with identical entropy values | 0.431117s | 0.378430s | 1.14x |
| `AlignmentEntropy._entropy_from_ascii_codes` array-to-list conversion | 12 taxa x 800k clean DNA sites, full helper output conversion with identical Python float list | 0.750628s | 0.476712s | 1.57x |
| `AlignmentEntropy._entropy_from_counts` masked log terms | 20 x 200k sparse protein count matrix, side-by-side previous `np.where(probs > 0, probs * log2(probs), 0)` terms | 0.040059s | 0.027744s | 1.44x |
| `AlignmentEntropy`/`MaskAlignment` protein entropy column reductions | count matrices shaped 20x5000 / 64x20000, side-by-side previous in-place probability/log-probability product plus `np.sum(..., axis=0)` | 2.243765s / 3.478205s | 1.991224s / 2.424359s | 1.13x / 1.43x |
| Shared boolean mask `any`/`all` reductions | empty / 10 / 1000 / 1M / 1000x1000 boolean masks, side-by-side previous top-level `np.any`/`np.all` dispatch | `any`: 4.323206s / 3.550042s / 1.254590s / 0.002394s / 0.001399s; `all`: 3.013299s / 2.108792s / 1.183967s / 0.001637s / 0.001208s | `any`: 0.750272s / 0.762969s / 0.429710s / 0.000900s / 0.000662s; `all`: 0.505534s / 1.391566s / 0.265627s / 0.001174s / 0.000669s | `any`: 5.76x / 4.65x / 2.92x / 2.66x / 2.11x; `all`: 5.96x / 1.52x / 4.46x / 1.39x / 1.80x |
| Alignment plotting/GWAS no-axis mask guards | sparse/dense boolean masks sized 100 / 1000 / 100k / 1M, side-by-side previous top-level `np.any(mask)` dispatch used before plotting or contingency tests | sparse: 2.580420s / 1.544515s / 0.016977s / 0.000714s; dense: 3.705966s / 1.403889s / 0.018599s / 0.001354s | sparse: 1.436175s / 0.655550s / 0.005228s / 0.000345s; dense: 1.056168s / 0.697029s / 0.026244s / 0.000341s | sparse: 1.80x / 2.36x / 3.25x / 2.07x; dense: 3.51x / 2.01x / 0.71x / 3.97x |
| `AlignmentEntropy.calculate_site_entropies` ASCII DNA entropy counts | 3000 taxa x 8000 sites, four observed DNA symbols, side-by-side previous boolean `np.sum(..., axis=0)` counts | 0.306331s | 0.138240s | 2.22x |
| `AlignmentEntropy.calculate_site_entropies` gap-code mask construction | 1200 taxa x 12000 sites, protein alphabet plus gaps/ambiguous symbols, side-by-side previous lookup-mask gather | 0.392355s | 0.335875s | 1.17x |
| `AlignmentEntropy.calculate_site_entropies` single valid-symbol shortcut | 1200 taxa x 12000 sites, conserved ASCII DNA alignment, side-by-side previous count/probability path | 0.088211s | 0.072725s | 1.21x |
| `AlignmentEntropy.calculate_site_entropies` identical-sequence shortcut | 1200 taxa x 12000 sites, identical ASCII DNA alignment, side-by-side previous byte-matrix path | 0.092795s | 0.005340s | 17.38x |
| `AlignmentEntropy.calculate_site_entropies` identical-row no-slice scan | 1M identical ASCII DNA rows, side-by-side previous `sequences[1:]` equality scan | 0.470263s | 0.391903s | 1.20x |
| `AlignmentEntropy.calculate_site_entropies` deferred Unicode gap-char set | 400 taxa x 2000 ASCII DNA sites, side-by-side previous unconditional `get_gap_chars` set construction | 0.389624s | 0.279371s | 1.39x |
| `AlignmentEntropy.calculate_site_entropies` entropy column dot | 2 / 4 / 8 / 20 states by 12000 / 12000 / 12000 / 5000 sites, side-by-side previous small-alphabet multiply plus `np.sum(axis=0)` path | 0.000073515s / 0.000106187s / 0.000208391s / 0.000220500s | 0.000051922s / 0.000071867s / 0.000164666s / 0.000112586s | 1.42x / 1.48x / 1.27x / 1.96x |
| `AlignmentEntropy.run` verbose batched text output | 100k site rows, mocked alignment/read and identical stdout text | 0.067358s | 0.055155s | 1.22x |
| `AlignmentEntropy.run` verbose text direct rows | 100k site rows, mocked alignment/read/calculation, captured stdout and identical text, side-by-side previous row-dict formatter | 0.057309s | 0.047691s | 1.20x |
| `AlignmentEntropy.run` nonverbose summary output | 100k site entropies, mocked alignment/read/calculation | 0.044581s | 0.000471s | 94.7x |
| `AlignmentEntropy.run` nonverbose plot-only series preparation | 1M site entropies, identical rounded plotted site/value arrays without temporary row dictionaries | 0.404570s | 0.200756s | 2.02x |
| `alignment_entropy` module import without eager NumPy lookup tables | cold subprocess import after lazy NumPy lookup construction | 0.078578s | 0.029387s | 2.67x |
| `alignment_entropy` module import without eager JSON/plot config helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.011941s | 0.005361s | 2.23x |
| `alignment_entropy` module import without shared `hashlib`/`typing` startup | median cold subprocess import after localizing helper hashing and replacing annotation-only typing aliases | 0.006227s | 0.001422s | 4.38x |
| `SumOfPairsScore.determine_number_of_matches_and_total_pairs` | 180 taxa x 1200 sites, complete equal-length pair set, alphabet `ACGT-` | 0.1360s | 0.0025s | 54.4x |
| `SumOfPairsScore.run` complete equal-length setup | 1200 taxa x 300 sites, complete equal-length pair set from reference IDs | 0.447407s | 0.000721s | 620.5x |
| `SumOfPairsScore._calculate_equal_length_complete_records` identical records | 1200 taxa x 1000 sites, complete equal-length query/reference records, side-by-side previous matrix stack path | 0.001236s | 0.000309s | 4.00x |
| `SumOfPairsScore._calculate_equal_length_complete_records` matching taxa count | 4000 taxa x 2000 sites, 2% changed residues, side-by-side previous boolean `np.sum(..., axis=0)` | 0.008258s | 0.005935s | 1.39x |
| `SumOfPairsScore._calculate_equal_length_complete_records` match-vector sum | 100k-site per-site matching-pair vector, side-by-side previous `np.sum` total | 0.000053360s | 0.000038879s | 1.37x |
| `SumOfPairsScore._calculate_equal_length_complete_pairs` ordered pair-set check | 1400 taxa x 300 sites, complete ordered `itertools.combinations` pair list, side-by-side previous full expected-set comparison | 0.528270s | 0.059062s | 8.94x |
| `SumOfPairsScore._has_complete_pair_set` streamed ordered validation | 2400 taxa, 2878800 ordered pair IDs, side-by-side previous nested slice validation | 0.355139s | 0.165752s | 2.14x |
| `sum_of_pairs_score` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.212225s | 0.119312s | 1.78x |
| `sum_of_pairs_score` module import without eager NumPy | cold subprocess import after lazy NumPy proxy | 0.085850s | 0.029637s | 2.90x |
| `sum_of_pairs_score` module import without eager multiprocessing | cold subprocess import after lazy multiprocessing proxy and localized `partial` import | 0.011085s | 0.006683s | 1.66x |
| `sum_of_pairs_score` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006175s | 0.005132s | 1.20x |
| `sum_of_pairs_score` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only typing aliases to built-in annotations | 0.003356s | 0.001074s | 3.12x |
| `SumOfPairsScore._process_pair_batch` | 3570 mixed-length incomplete pair-set comparisons, 120 taxa x ~1600 sites | 0.5732s | 0.0141s | 40.7x |
| `SumOfPairsScore._process_pair_batch` tuple cache and `count_nonzero` masks | 3570 mixed-length incomplete pair-set comparisons, 120 taxa x ~1600 sites, side-by-side previous nested-dict worker | 0.014967s | 0.007266s | 2.06x |
| `SumOfPairsScore.determine_number_of_matches_and_total_pairs` sequential equal-length cache | 49 incomplete pair comparisons, 14 taxa x 50k sites, side-by-side previous per-pair array conversion | 0.001350s | 0.000487s | 2.77x |
| `SumOfPairsScore.determine_number_of_matches_and_total_pairs` sequential mixed-length path | 8 taxa x 50k reference sites x 49k query sites, 28 pairs | 0.141108s | 0.000271s | 520.8x |
| `SumOfPairsScore.determine_number_of_matches_and_total_pairs` unchanged incomplete pairs | 49 incomplete mixed-length pair comparisons, 14 taxa x ~50k sites, side-by-side previous sequential array path | 0.000806s | 0.000244s | 3.30x |
| `SumOfPairsScore._read_fasta` | 50k FASTA records, mixed-case 120 bp each | 0.1076s | 0.0342s | 3.1x |
| `SumOfPairsScore._read_fasta` shared unique first-token parser | 50k FASTA records, mixed-case 120 bp each, legacy `SimpleFastaParser` baseline | 0.037212s | 0.033738s | 1.10x |
| `CompositionPerTaxon.calculate_composition_per_taxon` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.1262s | 0.0548s | 2.3x |
| `CompositionPerTaxon.calculate_composition_per_taxon` matrix counts | 800 taxa x 4000 sites, alphabet `ACGT-?NX*` | 0.0844s | 0.0534s | 1.6x |
| `CompositionPerTaxon.calculate_composition_per_taxon` uint8 row bincounts | 800 taxa x 4000 sites, alphabet `ACGT-?NX*` | 0.0947s | 0.0406s | 2.3x |
| `CompositionPerTaxon.calculate_composition_per_taxon` invalid lookup | 2000 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.143750s | 0.122099s | 1.18x |
| `CompositionPerTaxon.calculate_composition_per_taxon` small-alphabet counts | 2000 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.115404s | 0.083495s | 1.38x |
| `CompositionPerTaxon.calculate_composition_per_taxon` short protein count strategy | 50k taxa x 50 sites, alphabet `ACDEFGHIKLMNPQRSTVWY-X?` | 0.120957s | 0.049710s | 2.43x |
| `CompositionPerTaxon.calculate_composition_per_taxon` no-gap DNA valid lengths | 2000 taxa x 5000 sites, alphabet `ACGT`, side-by-side previous full valid-mask path | 0.054157s | 0.039753s | 1.36x |
| `CompositionPerTaxon.calculate_composition_per_taxon` no-gap protein valid lengths | 2000 taxa x 5000 sites, 20 amino-acid symbols, side-by-side previous full valid-mask path | 0.069665s | 0.047561s | 1.46x |
| `CompositionPerTaxon.calculate_composition_per_taxon` single valid-symbol shortcut | 3000 taxa x 10000 sites, conserved ASCII DNA alignment, side-by-side previous count/frequency path | 0.109004s | 0.099306s | 1.10x |
| `CompositionPerTaxon.calculate_composition_per_taxon` identical multi-symbol shortcut | 1200 taxa x 12000 identical protein sites, lowercase/uppercase variants, side-by-side previous matrix path | 0.055857s | 0.006509s | 8.58x |
| `CompositionPerTaxon.calculate_composition_per_taxon` identical frequency sum | 20-symbol identical-sequence count vector, side-by-side previous `np.sum` normalization | 0.000007918s | 0.000003764s | 2.10x |
| `CompositionPerTaxon.calculate_composition_per_taxon` identical-row no-slice scan | 300k conserved 20-symbol protein records, side-by-side previous `sequences[1:]` equality scan | 0.614861s | 0.460860s | 1.33x |
| `CompositionPerTaxon.calculate_composition_per_taxon` identical-row direct record scan | 300k conserved 20-symbol protein records, side-by-side previous sequence-list setup with identical symbols and frequency rows | 0.452163s | 0.366624s | 1.23x |
| `CompositionPerTaxon.calculate_composition_per_taxon` raw-identical normalization scan | 300k conserved raw-identical 20-symbol protein records, side-by-side previous eager uppercase record setup | 0.766295s | 0.508587s | 1.51x |
| `CompositionPerTaxon.run` text output formatting | 100k taxon rows x 4 composition symbols, mocked alignment/read and identical stdout text | 0.853874s | 0.194774s | 4.38x |
| `CompositionPerTaxon.run` JSON payload formatting | 100k taxon rows x 4 composition symbols, mocked calculation rows, side-by-side previous index lookup loop | 0.147151s | 0.137843s | 1.07x |
| `composition_per_taxon` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.085520s | 0.023586s | 3.63x |
| `composition_per_taxon` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006286s | 0.004998s | 1.26x |
| `composition_per_taxon` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.002980s | 0.001067s | 2.79x |
| `alignment_outlier_taxa` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.003062s | 0.001096s | 2.79x |
| `GCContent.calculate_gc_per_sequence_data` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.0698s | 0.0106s | 6.6x |
| `GCContent.calculate_gc_total_value` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.0704s | 0.0101s | 7.0x |
| `GCContent.calculate_gc_per_sequence_data` byte lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?*XN` | 0.1035s | 0.0461s | 2.2x |
| `GCContent.calculate_gc_total_value` byte lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?*XN` | 0.1031s | 0.0465s | 2.2x |
| `GCContent.calculate_gc_per_sequence_data` case-insensitive byte lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?*XN` | 0.047277s | 0.044028s | 1.07x |
| `GCContent.calculate_gc_total_value` case-insensitive byte lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?*XN` | 0.047642s | 0.043335s | 1.10x |
| `GCContent.calculate_gc_per_sequence_data` ASCII matrix counts | 2000 taxa x 5000 sites, alphabet `ACGTN-?X*` | 0.037366s | 0.025886s | 1.44x |
| `GCContent.calculate_gc_per_sequence_data` no-gap ASCII valid lengths | 2000 taxa x 5000 sites, alphabet `ACGT`, side-by-side previous valid lookup count path | 0.057049s | 0.033166s | 1.72x |
| `GCContent.calculate_gc_per_sequence_data` variable-length ASCII fallback | 50k records x 90-106 bp, alphabet `ACGTN-?X*` | 0.273960s | 0.195534s | 1.40x |
| `GCContent.calculate_gc_per_sequence_data` identical-sequence shortcut | 1200 taxa x 12000 identical DNA sites with ambiguous/gap symbols, lowercase/uppercase variants, side-by-side previous matrix path | 0.062592s | 0.005746s | 10.89x |
| `GCContent._gc_counts_from_ascii_matrix` raw-identical normalized shortcut | 500k identical DNA records with ambiguous/gap symbols, side-by-side previous per-row uppercase equality scan | 0.454749s | 0.375646s | 1.21x |
| `GCContent.calculate_gc_total_value` ASCII flat counts | 2000 taxa x 5000 sites, alphabet `ACGTN-?X*` | 0.037335s | 0.022379s | 1.67x |
| `GCContent.calculate_gc_total_value` no-gap ASCII valid length | 2000 taxa x 12000 sites, alphabet `ACGT`, side-by-side previous valid lookup count path | 0.066701s | 0.045941s | 1.45x |
| `GCContent.calculate_gc_total_value` mixed Unicode/ASCII fallback | 50k ASCII records x 120 bp plus one Unicode record, alphabet `ACGTN-?X*` | 0.307134s | 0.263488s | 1.17x |
| `GCContent.calculate_gc_total_value` mixed Unicode batched ASCII fallback | 50k ASCII records x 120 bp plus one Unicode record, alphabet `ACGTN-?X*` | 0.133722s | 0.034256s | 3.90x |
| `GCContent.calculate_gc_total_value` identical-sequence shortcut | 1200 taxa x 12000 identical DNA sites with ambiguous/gap symbols, lowercase/uppercase variants, side-by-side previous flat byte path | 0.059239s | 0.007356s | 8.05x |
| `GCContent._gc_counts_from_upper_sequence` identical-row count helper | 4M clean uppercase DNA chars / 4M gappy uppercase DNA chars / 1.05M uppercase Unicode DNA chars, side-by-side previous helper re-uppercase and repeated `str.count` output | 0.105214s / 0.238117s / 0.052446s | 0.037932s / 0.051851s / 0.020863s | 2.77x / 4.59x / 2.51x |
| `GCContent._gc_total_from_ascii` raw-identical normalized shortcut | 500k identical DNA records with ambiguous/gap symbols, side-by-side previous per-row uppercase equality scan | 0.384826s | 0.095858s | 4.01x |
| `GCContent.calculate_gc_per_sequence` batched text output | 50k sequence rows, mocked per-sequence data and identical stdout text | 0.025783s | 0.019660s | 1.31x |
| `gc_content` module import without eager NumPy/Bio.Align | cold subprocess import after lazy NumPy lookup construction and annotation-only Bio.Align import | 0.111378s | 0.023406s | 4.76x |
| `gc_content` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006135s | 0.004972s | 1.23x |
| `gc_content` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.002604s | 0.000947s | 2.75x |
| `gc_content` module import without eager file helper | median cold subprocess import after lazy alignment-reader wrappers in service and alignment base | 0.048802s | 0.030473s | 1.60x |
| `AlignmentLength.run` FASTA parser fast path | 50k FASTA records x 304 bp, cold alignment cache each baseline run | 0.605968s | 0.094259s | 6.43x |
| `AlignmentLength._get_fasta_alignment_length` direct scanner | 50k wrapped FASTA records x 304 bp, side-by-side previous `SimpleFastaParser` length helper | 0.068333s | 0.046054s | 1.48x |
| `AlignmentLength._get_fasta_alignment_length` binary direct scanner | 50k wrapped ASCII FASTA records x 304 bp, side-by-side previous text direct scanner | 0.026663s | 0.024407s | 1.09x |
| `alignment_length` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.219706s | 0.117138s | 1.88x |
| `alignment_length` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006248s | 0.005071s | 1.23x |
| `alignment_length` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only typing aliases to built-in annotations | 0.002915s | 0.001239s | 2.35x |
| `AlignmentLengthNoGaps.get_sites_no_gaps_count` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.1058s | 0.0075s | 14.1x |
| `AlignmentLengthNoGaps.get_sites_no_gaps_count` gap lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?NX*` | 0.086420s | 0.023064s | 3.75x |
| `AlignmentLengthNoGaps.get_sites_no_gaps_count` gap-code column reduction | 1200 taxa x 12000 sites, alphabet `ACGT-?NX*`, side-by-side previous lookup-mask gather | 0.040526s | 0.014469s | 2.80x |
| `AlignmentLengthNoGaps.get_sites_no_gaps_count` no-gap byte short-circuit | 2000 taxa x 12000 sites DNA `ACGT` / 1200 taxa x 12000 sites 20 amino-acid symbols, side-by-side previous column-reduction path | 0.019192s / 0.008347s | 0.011522s / 0.005520s | 1.67x / 1.51x |
| `AlignmentLengthNoGaps.get_sites_no_gaps_count` no-gap constant byte precheck | 1200 taxa x 12000 no-gap DNA sites, side-by-side previous NumPy gap-code precheck | 0.012226s | 0.011693s | 1.05x |
| `AlignmentLengthNoGaps.get_sites_no_gaps_count` no-gap byte membership precheck | 1200 taxa x 12000 non-identical no-gap DNA sites, side-by-side previous per-code one-byte allocation precheck | 0.006010s | 0.005257s | 1.14x |
| `AlignmentLengthNoGaps.get_sites_no_gaps_count` identical-sequence shortcut | 2000 taxa x 12000 DNA sites, identical no-gap/gappy sequences; 1200 taxa x 12000 identical protein sites, side-by-side previous alignment-wide byte path | 0.032510s / 0.063387s / 0.008216s | 0.000157s / 0.000284s / 0.000088s | 207.13x / 223.03x / 93.32x |
| `AlignmentLengthNoGaps` identical Unicode fallback counts | 100k-site mixed-case Unicode identical sequence helper, DNA / protein, side-by-side previous `upper()` plus uppercase gap counts | 0.000486s / 0.000401s | 0.000379s / 0.000251s | 1.28x / 1.60x |
| `AlignmentLengthNoGaps` identical-row no-slice scan | see shared `PairwiseIdentity` row above for the common helper benchmark | 0.259744s / 0.018993s / 0.082823s | 0.126212s / 0.000004s / 0.038411s | 2.06x / 5301.12x / 2.16x |
| `alignment_length_no_gaps` module import without eager NumPy/Bio.Align | cold subprocess import after lazy NumPy lookup construction and annotation-only Bio.Align import | 0.116950s | 0.024571s | 4.76x |
| `alignment_length_no_gaps` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.007086s | 0.006015s | 1.18x |
| `alignment_length_no_gaps` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.004736s | 0.002167s | 2.18x |
| `MaskAlignment.calculate_keep_mask` with entropy threshold | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.1572s | 0.0315s | 5.0x |
| `MaskAlignment.calculate_keep_mask` gap lookup with entropy threshold | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.070534s | 0.035500s | 1.99x |
| `MaskAlignment.calculate_keep_mask` protein entropy block counts | 500 taxa x 8000 sites, protein alphabet plus gaps/ambiguous symbols, entropy threshold enabled | 0.114204s | 0.071811s | 1.59x |
| `MaskAlignment.calculate_keep_mask` protein entropy column-major indexing | 1200 taxa x 6000 sites, alphabet `ACDEFGHIKLMNPQRSTVWY-X?*`, entropy threshold enabled | 0.151865s | 0.125777s | 1.21x |
| `MaskAlignment.calculate_keep_mask` clean protein entropy mask elision | 1200 taxa x 6000 clean protein sites, entropy threshold enabled, side-by-side previous occupancy and valid-mask entropy path | 0.072934s | 0.052166s | 1.40x |
| `MaskAlignment._column_entropies_from_ascii_codes` masked log terms | 1000 taxa x 12000 clean protein sites with two symbols per site, side-by-side previous boolean-indexed `log2` terms | 0.141908s | 0.124990s | 1.14x |
| `MaskAlignment.calculate_keep_mask` ASCII DNA entropy counts | 3000 taxa x 8000 sites, four observed DNA symbols, side-by-side previous boolean `np.sum(..., axis=0)` counts | 0.125266s | 0.109594s | 1.14x |
| `MaskAlignment.calculate_keep_mask` single-symbol entropy shortcut | 2000 taxa x 12000 conserved clean DNA sites, entropy threshold enabled, side-by-side previous count/probability path | 0.104621s | 0.080733s | 1.30x |
| `MaskAlignment.calculate_keep_mask` all-pass threshold shortcut | DNA 2000 taxa x 12000 sites, alphabet `ACGT-?NX*` / protein 1200 taxa x 12000 sites, protein alphabet plus gaps/ambiguous symbols, `max_gap=1`, `min_occupancy=0`, no entropy | 0.065244s / 0.035698s | 0.013126s / 0.006714s | 4.97x / 5.32x |
| `MaskAlignment.calculate_keep_mask` all-pass pre-materialization return | 2000 taxa x 12000 DNA sites, all-pass thresholds, side-by-side previous sequence materialization before all-true mask | 0.070370167s | 0.000251167s | 280.17x |
| `MaskAlignment.calculate_keep_mask` clean no-entropy shortcut | 2000 taxa x 12000 clean DNA sites, `max_gap=0.2`, `min_occupancy=0.8`, no entropy threshold, side-by-side previous occupancy mask path | 0.060156s | 0.017106s | 3.52x |
| `MaskAlignment.calculate_keep_mask` identical-sequence shortcut | 1200 taxa x 12000 identical DNA sites with gaps/ambiguous symbols, entropy threshold enabled, lowercase/uppercase variants, side-by-side previous matrix path | 0.100434s | 0.005945s | 16.89x |
| `MaskAlignment.calculate_keep_mask` identical-sequence no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.196770s / 0.006034s / 0.084992s | 0.040581s / 0.000003s / 0.042880s | 4.85x / 1765.72x / 1.98x |
| `MaskAlignment.calculate_keep_mask` identical-sequence byte mask | 8 taxa x 1,000,000 identical DNA sites with gaps/ambiguous symbols, entropy threshold enabled, side-by-side previous Python generator mask path | 0.416832s | 0.009614s | 43.35x |
| `MaskAlignment.apply_mask` | 1200 taxa x 3000 sites, 2000 kept sites | 0.1349s | 0.0044s | 30.7x |
| `MaskAlignment.apply_mask` all-keep shortcut | 2000 taxa x 5000 sites, lowercase ASCII sequences and all-true keep mask, side-by-side previous matrix-slice path | 0.032670s | 0.006820s | 4.79x |
| `MaskAlignment.apply_mask` all-trim shortcut | 2000 taxa x 5000 sites, lowercase ASCII sequences and all-false keep mask, side-by-side previous matrix-slice path | 0.008348s | 0.001225s | 6.81x |
| `MaskAlignment.run` batched FASTA text output | 100k masked FASTA records, mocked alignment/read and identical stdout text | 0.021122s | 0.011388s | 1.85x |
| `MaskAlignment.run` JSON kept-site count | 2M-site boolean keep mask, side-by-side previous boolean `np.sum` count | 0.000407500s | 0.000067708s | 6.02x |
| `mask_alignment` module import without eager NumPy lookup tables | cold subprocess import after lazy NumPy lookup construction and postponed annotations | 0.081088s | 0.023500s | 3.45x |
| `mask_alignment` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006133s | 0.005002s | 1.23x |
| `mask_alignment` module import without `typing` startup | median cold subprocess import after removing annotation-only typing aliases under postponed annotations | 0.002835s | 0.000953s | 2.98x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.1942s | 0.0347s | 5.6x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` gap lookup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.066770s | 0.036214s | 1.84x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` protein block counts | 500 taxa x 8000 sites, protein alphabet plus gaps/ambiguous symbols | 0.093095s | 0.062521s | 1.49x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` protein column-major block counts | 500 taxa x 8000 sites, protein alphabet plus gaps/ambiguous symbols | 0.060433s | 0.045616s | 1.32x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` DNA observed-symbol validity | 500 taxa x 8000 sites, alphabet `ACGT-?NX*`, side-by-side previous full valid-mask path | 0.074886s | 0.044362s | 1.69x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` DNA no-gap observed-symbol validity | 1000 taxa x 12000 sites, alphabet `ACGT`, side-by-side previous full valid-mask path | 0.142018s | 0.118764s | 1.20x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` protein no-gap mask elision | 1000 taxa x 5000 sites, 20 amino-acid symbols, side-by-side previous full valid-mask path | 0.099297s | 0.069440s | 1.43x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` single valid-symbol shortcut | 1200 taxa x 12000 sites, conserved ASCII DNA alignment, side-by-side previous count/frequency path | 0.058379s | 0.044176s | 1.32x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` identical-sequence shortcut | 1200 taxa x 12000 identical ASCII DNA sites, lowercase/uppercase variants, side-by-side previous matrix path | 0.076192s | 0.006297s | 12.10x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` raw-identical normalization scan | 300k raw-identical DNA rows, side-by-side previous eager uppercase sequence setup with zero PIC values | 0.158051s | 0.099522s | 1.59x |
| `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` identical-sequence no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.223614s / 0.006108s / 0.073953s | 0.104555s / 0.000004s / 0.049539s | 2.14x / 1610.87x / 1.49x |
| `EvolutionaryRatePerSite`/`CompositionalBiasPerSite` column count sum-of-squares | count matrices shaped 4x12000 / 8x12000 / 20x5000 / 64x20000, side-by-side previous `np.sum(counts * counts, axis=0)` | 0.426963s / 0.656234s / 0.739237s / 0.859598s | 0.317205s / 0.516822s / 0.375821s / 0.573535s | 1.35x / 1.27x / 1.97x / 1.50x |
| `EvolutionaryRatePerSite.remove_gap_characters` cached translate deletion | 2M-character mixed-case sequence, gap/ambiguous symbols `-?*XxNn`, identical uppercase filtered output | 0.307932s | 0.012577s | 24.48x |
| `EvolutionaryRatePerSite.get_number_of_occurrences_per_character` record-wise direct count loop | 200 sampled columns from 5000 taxa x 2000 sites, alphabet `ACGT-?NX*`, side-by-side previous column slicing path with identical `Counter` output | 1.473608s | 0.948019s | 1.55x |
| `EvolutionaryRatePerSite.get_number_of_occurrences_per_character` multi-character gap fallback | 100k records, non-ASCII column, gap tokens `--` and `?`, identical `Counter` output | 0.020868s | 0.016616s | 1.26x |
| `EvolutionaryRatePerSite.run` batched text output | 100k site rows, mocked alignment/read and identical stdout text | 0.069242s | 0.055913s | 1.24x |
| `EvolutionaryRatePerSite.run` direct terminal text output | 100k site values, mocked alignment/read and identical stdout text | 0.083225s | 0.064262s | 1.30x |
| `EvolutionaryRatePerSite.run` plot-only series preparation | 1M site rates, identical rounded plotted site/value arrays without temporary row dictionaries | 0.376116s | 0.196265s | 1.92x |
| `evolutionary_rate_per_site` module import without eager NumPy/Bio.Align | cold subprocess import after lazy NumPy lookup construction and annotation-only Bio.Align import | 0.114353s | 0.029664s | 3.86x |
| `evolutionary_rate_per_site` module import without eager JSON/plot config helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.013219s | 0.005953s | 2.22x |
| `evolutionary_rate_per_site` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.002573s | 0.000905s | 2.84x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.6591s | 0.0349s | 18.9x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` gap lookup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.086270s | 0.038467s | 2.24x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` protein ASCII counts | 1000 taxa x 5000 sites, 20 amino-acid symbols | 0.078547s | 0.050084s | 1.57x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` protein column-major block counts | 1000 taxa x 5000 sites, protein alphabet plus gaps/ambiguous symbols | 0.078616s | 0.059675s | 1.32x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` DNA observed-symbol validity | 500 taxa x 8000 sites, alphabet `ACGT-?NX*`, side-by-side previous full valid-mask path | 0.056998s | 0.033467s | 1.70x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` DNA no-gap observed-symbol validity | 1000 taxa x 12000 sites, alphabet `ACGT`, side-by-side previous full valid-mask path | 0.137986s | 0.083428s | 1.65x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` protein no-gap mask elision | 1000 taxa x 5000 sites, 20 amino-acid symbols, side-by-side previous full valid-mask path | 0.052347s | 0.043303s | 1.21x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` single valid-symbol shortcut | 1200 taxa x 12000 sites, conserved ASCII DNA alignment, side-by-side previous count/statistic path | 0.066130s | 0.047264s | 1.40x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` identical-sequence shortcut | 1200 taxa x 12000 identical ASCII DNA sites, lowercase/uppercase variants, side-by-side previous matrix/statistic path | 0.087052s | 0.010425s | 8.35x |
| `CompositionalBiasPerSite.calculate_compositional_bias_per_site` identical-row no-slice scan | 1M identical ASCII DNA rows x 32 sites, side-by-side previous `sequences[1:]` equality scan | 0.486066s | 0.234605s | 2.07x |
| `CompositionalBiasPerSite.get_number_of_occurrences_per_character` record-wise direct count loop | 200 sampled columns from 5000 taxa x 2000 sites, alphabet `ACGT-?NX*`, side-by-side previous column slicing path with identical first-seen count order | 1.300166s | 0.822473s | 1.58x |
| `CompositionalBiasPerSite._erfc_array` vectorized square roots | 500k chi-square half-statistics, side-by-side previous scalar `sqrt` + `erfc` loop without SciPy imports | 0.048987s | 0.040776s | 1.20x |
| `CompositionalBiasPerSite` corrected p-value reconstruction | 100k sites, 9,916 valid p-values interleaved with `"nan"` slots | 2.063417s | 0.004771s | 432.50x |
| `CompositionalBiasPerSite` corrected p-value all-valid restore | 1M corrected p-values with no `"nan"` slots | 0.194807s | 0.000001334s | 146019.11x |
| `CompositionalBiasPerSite` all-valid p-value extraction | 3 repeated 1M-site p-value arrays with no NaN slots, side-by-side previous boolean-index extraction | 0.084985s | 0.077307s | 1.10x |
| `CompositionalBiasPerSite` statistic result construction | 1M per-site statistic/p-value pairs converted to `Power_divergenceResult` rows | 2.121716s | 1.699578s | 1.25x |
| `CompositionalBiasPerSite._false_discovery_control` small-list path without NumPy startup | cold subprocess, 4 p-values through Benjamini-Hochberg helper | 0.066106s | 0.024334s | 2.72x |
| `CompositionalBiasPerSite.run` batched text output | 100k site rows, mocked alignment/read and identical stdout text | 0.072149s | 0.061021s | 1.18x |
| `CompositionalBiasPerSite.run` direct terminal text output | 100k site rows, mocked alignment/read and identical stdout text | 0.200112s | 0.168921s | 1.18x |
| `CompositionalBiasPerSite.run` plot-only series preparation | 1M corrected p-values with interleaved `"nan"` slots, identical rounded plotted site/value arrays without temporary row dictionaries | 0.935595s | 0.241023s | 3.88x |
| `compositional_bias_per_site` module import | cold subprocess import, avoid eager `scipy.stats` import | 0.610058s | 0.141978s | 4.3x |
| `compositional_bias_per_site` module import without eager NumPy/Bio.Align | cold subprocess import after lazy NumPy lookup construction and annotation-only Bio.Align import | 0.117581s | 0.029701s | 3.96x |
| `compositional_bias_per_site` module import without eager JSON/plot config helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.012085s | 0.005828s | 2.07x |
| `compositional_bias_per_site` module import without typing startup | median cold subprocess import after converting annotation-only typing names to built-in postponed annotations | 0.006420s | 0.004329s | 1.48x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.0756s | 0.0286s | 2.6x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` invalid lookup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.058200s | 0.025501s | 2.28x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` protein byte bincounts | 2000 taxa x 5000 sites, protein alphabet plus gaps/ambiguous symbols | 0.102812s | 0.095350s | 1.08x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` observed-symbol validity DNA | 2000 taxa x 5000 sites, alphabet `ACGT-?NX*`, side-by-side previous full valid-mask path | 0.100830s | 0.066724s | 1.51x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` observed-symbol validity protein | 2000 taxa x 5000 sites, protein alphabet plus gaps/ambiguous symbols, side-by-side previous full valid-mask path | 0.134780s | 0.085772s | 1.57x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` no-gap observed-symbol validity | 2000 taxa x 5000 sites, DNA `ACGT` / 20 amino-acid symbols, side-by-side previous full valid-mask path | 0.078720s / 0.091700s | 0.052968s / 0.046669s | 1.49x / 1.96x |
| RCV/RCVT/AlignmentOutlierTaxa valid-length mask counts | 5000 taxa x 4000-site boolean valid mask, side-by-side previous `np.sum(mask, axis=1).astype(float)` | 0.008596s | 0.004423s | 1.94x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` many-short protein count table | 50k taxa x 50 no-gap protein sites, side-by-side previous one `np.bincount` per taxon row | 1.300922s | 0.921645s | 1.41x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` count-matrix column totals | count matrices shaped 260x4 / 1200x20 / 2000x20 / 50000x20, side-by-side previous `np.sum(..., axis=0)` wrapper | 3.497999s / 3.163409s / 4.144007s / 4.768939s | 1.929073s / 2.924800s / 3.469921s / 3.756192s | 1.81x / 1.08x / 1.19x / 1.27x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` identical-sequence shortcut | 1200 taxa x 12000 identical DNA sites, lowercase/uppercase variants, side-by-side previous matrix path | 0.069846s | 0.010896s | 6.41x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` raw-identical normalization scan | 300k raw-identical DNA rows, side-by-side previous eager uppercase sequence setup with identical zero rows | 0.298070s | 0.116078s | 2.57x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` identical-sequence no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.058848s / 0.007070s / 0.211013s | 0.042901s / 0.000004s / 0.037338s | 1.37x / 1844.67x / 5.65x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` single-record early return | 5 repeated 4.5M-site single-record RCVT calls, side-by-side previous sequence materialization before zero row | 0.001927875s | 0.000000416s | 4630.96x |
| `RelativeCompositionVariabilityTaxon.calculate_rows` zip-based row assembly | 500k mocked taxon records and computed RCVT values, side-by-side previous index lookups | 0.317740s | 0.284361s | 1.12x |
| `RelativeCompositionVariabilityTaxon.run` batched text output | 50k taxon rows, mocked alignment/read and identical stdout text | 0.017426s | 0.011493s | 1.52x |
| `RelativeCompositionVariabilityTaxon._plot_rcvt` plot series preparation | 300k taxon rows, repeated RCVT values, identical stable descending taxon order and plotted values | 0.068288s | 0.042674s | 1.60x |
| `rcvt` module import without eager NumPy lookup tables | cold subprocess import after lazy NumPy lookup construction | 0.079652s | 0.029906s | 2.66x |
| `rcvt` module import without eager JSON/plot config helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.011765s | 0.005049s | 2.33x |
| `AlignmentOutlierTaxa.calculate_outliers` | 180 taxa x 1200 sites, alphabet `ACGT-?NX*` | 0.734118s | 0.008124s | 90.36x |
| `AlignmentOutlierTaxa.calculate_outliers` ASCII numeric core uint8 lookup | 400 taxa x 1200 sites, alphabet `ACGT-?NX*` | 0.546744s | 0.176515s | 3.10x |
| `AlignmentOutlierTaxa.calculate_outliers` long-branch matrix products | 400 taxa x 1200 sites, alphabet `ACGT-?NX*` | 0.165889s | 0.035340s | 4.69x |
| `AlignmentOutlierTaxa._blocked_long_branch_proxy` large ambiguous DNA fallback | 3000 taxa x 300 sites, alphabet `ACGTN-`, side-by-side previous per-taxon row loop above full-matrix cutoff | 5.048598s | 2.537640s | 1.99x |
| `AlignmentOutlierTaxa.calculate_outliers` comparable-pair row counts | 2000 x 2000 comparable-pair boolean matrix, side-by-side previous `np.sum(..., axis=1)` row counts | 0.001769s | 0.001283s | 1.38x |
| `AlignmentOutlierTaxa.calculate_outliers` constant-composition shortcut | 1000 taxa x 5000 sites, conserved ASCII DNA alignment, side-by-side previous full feature pipeline | 0.429936s | 0.042515s | 10.11x |
| `AlignmentOutlierTaxa.calculate_outliers` identical multi-symbol shortcut | 1000 taxa x 5000 mixed-symbol DNA sites, lowercase/uppercase variants, side-by-side previous full feature pipeline | 8.073867s | 0.014392s | 561.01x |
| `AlignmentOutlierTaxa` identical Unicode valid length | 100k-site uppercase Unicode identical-sequence helper, DNA / protein, side-by-side previous Python character-membership loop | 0.006644s / 0.005970s | 0.000494s / 0.000388s | 13.45x / 15.38x |
| `AlignmentOutlierTaxa.calculate_outliers` identical-sequence no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.081116s / 0.008015s / 0.090683s | 0.046170s / 0.000001s / 0.059985s | 1.76x / 12036.36x / 1.51x |
| `AlignmentOutlierTaxa.calculate_outliers` single-record early return | 4.5M-site all-invalid single-record DNA alignment, side-by-side previous full feature pipeline | 0.195257542s | 0.020914958s | 9.34x |
| `AlignmentOutlierTaxa.calculate_outliers` all-valid ASCII shortcut | 400 taxa x 1200 sites, variable-composition ASCII DNA alignment without invalid symbols | 0.143346s | 0.103475s | 1.39x |
| `AlignmentOutlierTaxa` ASCII symbol count setup | five repeated 1000-taxon x 5000-site protein row/site count builds, side-by-side previous per-symbol equality reductions | 3.406926s | 0.871172s | 3.91x |
| `AlignmentOutlierTaxa._symbol_counts_by_row` large-short ASCII counts | 20000 taxa x 128 sites / 30000 taxa x 96 sites / 50000 taxa x 64 sites / 80000 taxa x 128 sites, 20 valid symbols, side-by-side previous per-row `bincount` loop | 4.549077s / 6.862011s / 7.782192s / 6.304361s | 1.803506s / 1.586349s / 1.926483s / 2.222659s | 2.52x / 4.33x / 4.04x / 2.84x |
| `AlignmentOutlierTaxa._symbol_counts_by_site` expanded ASCII histogram threshold | 80000 taxa x 128 sites / 1200 taxa x 12000 sites / 3000 taxa x 8000 sites, 20 valid symbols, side-by-side previous repeated equality scans above the old 8M-cell cutoff | 5.065210s / 2.208961s / 4.031251s | 0.814091s / 0.832620s / 1.356734s | 6.22x / 2.65x / 2.97x |
| `AlignmentOutlierTaxa.calculate_outliers` all-valid protein long-branch formula | two repeated 220-taxon x 2000-site protein analyses, side-by-side previous all-valid pairwise matrix-product long-branch path | 9.241269s | 0.075789s | 121.93x |
| `AlignmentOutlierTaxa.calculate_outliers` composition-distance row norms | 400 taxa x 1200 DNA sites and 1000 taxa x 5000 protein sites, side-by-side previous `np.linalg.norm(..., axis=1)` with identical rows | 0.039359s / 0.229932s | 0.016265s / 0.188438s | 2.42x / 1.22x |
| `AlignmentOutlierTaxa.calculate_outliers` entropy column dot | site probability/log-probability matrices shaped 4x12000 / 8x12000 / 20x5000 / 64x20000, side-by-side previous `np.sum(site_probs * log_probs, axis=0)` | 0.420023s / 0.620534s / 0.777484s / 1.564784s | 0.380694s / 0.425398s / 0.542617s / 1.111065s | 1.10x / 1.46x / 1.43x / 1.41x |
| `AlignmentOutlierTaxa.calculate_outliers` zipped row assembly | 100k synthetic taxa with six feature arrays and nested reason rows | 0.437748s | 0.370843s | 1.18x |
| `AlignmentOutlierTaxa.run` batched text output | 100k outlier rows, mocked alignment/read and identical stdout text | 0.103591s | 0.090365s | 1.15x |
| `alignment_outlier_taxa` module import without eager NumPy/json helpers | cold subprocess import after lazy NumPy proxy and JSON helper wrapper | 0.069810s | 0.022303s | 3.13x |
| `PlotAlignmentQC` composition-distance scatter panel | 5000 taxa, 250 flagged taxa, Matplotlib Agg setup only | 9.036485s | 0.015365s | 588.14x |
| `PlotAlignmentQC.run` shared plot-array preparation | 200k synthetic taxa, six QC feature arrays, occupancy/gap ordering, composition panel inputs, and heatmap z-score setup with identical arrays | 0.171046s | 0.097432s | 1.76x |
| `PlotAlignmentQC` heatmap equal-feature sigma shortcut | 200k / 500k identical finite feature values, side-by-side previous median/MAD plus `np.std` fallback | 0.001939s / 0.005776s | 0.000040728s / 0.000188017s | 47.61x / 30.72x |
| `PlotAlignmentQC._prepare_plot_arrays` one-pass flagged mask | 500k synthetic taxa, six QC feature arrays, and 13.5k flagged taxa, side-by-side previous `np.fromiter` mask pass | 0.363184s | 0.292066s | 1.24x |
| `PlotAlignmentQC._prepare_plot_arrays` empty-outlier mask shortcut | 500k synthetic taxa, six QC feature arrays, and no flagged taxa, side-by-side previous per-row set membership path | 1.277690s | 1.108565s | 1.15x |
| `PlotAlignmentQC._flag_colors` vectorized mask mapping | 1M ordered taxa flags, identical normal/flagged color sequence | 0.036504s | 0.014694s | 2.48x |
| `plot_alignment_qc` module import without eager NumPy/outlier/json/plot helpers | cold subprocess import after lazy NumPy proxy and forwarding helper wrappers | 0.075926s | 0.026252s | 2.89x |
| `plot_alignment_qc` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.036055s | 0.034670s | 1.04x |
| `ColumnScore.get_columns_from_alignments` + column matching | 260 taxa x 5000 sites query/reference alignments, alphabet `ACGT-?NX*` | 0.4111s | 0.0232s | 17.7x |
| `ColumnScore.run` direct byte-column matching | 260 taxa x 5000 sites query/reference FASTA alignments, alphabet `ACGT-?NX*` | 0.122301s | 0.053935s | 2.27x |
| `ColumnScore._calculate_matches_between_alignments_direct` same-object shortcut | 260 taxa x 5000 sites, alphabet `ACGT-?NX*`, side-by-side previous double-unique/intersection path | 0.008894s | 0.002266s | 3.93x |
| `ColumnScore._calculate_matches_between_alignments_direct` identical-sequence shortcut | 260 taxa x 5000 identical query/reference sites, lowercase/uppercase variants, side-by-side previous double-unique/intersection path | 0.006354s | 0.003171s | 2.00x |
| `ColumnScore._calculate_matches_between_alignments_direct` repeated-row same-object shortcut | 1200 taxa x 12000 sites, every row the same mixed-symbol ASCII sequence, side-by-side previous byte-matrix unique-column path | 0.094739s | 0.011452s | 8.27x |
| `ColumnScore._calculate_matches_between_alignments_direct` repeated-row separate-alignments shortcut | 1200 taxa x 12000 sites, reference/query each have repeated mixed-symbol ASCII rows with partially overlapping symbols, side-by-side previous double-unique/intersection path | 0.382126s | 0.046693s | 8.18x |
| `ColumnScore._calculate_matches_between_alignments_direct` taxon-count mismatch shortcut | 260-reference-taxon x 5000-site alignment against one-query-taxon x 5000-site alignment, side-by-side previous sequence materialization before zero-match result | 0.001336417s | 0.000003041s | 439.43x |
| `ColumnScore._repeated_sequence_symbols_ascii` no-slice row scan | 2M repeated ASCII rows, side-by-side previous `sequences[1:]` equality scan | 0.238253s | 0.085275s | 2.79x |
| `column_score` module import without eager NumPy/Bio.AlignIO | cold subprocess import after lazy NumPy, AlignIO, annotation, and JSON helpers | 0.124730s | 0.022308s | 5.59x |
| `column_score` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.002598s | 0.000871s | 2.98x |
| `DNAThreader.normalize_n_seq` | 96k amino acids with gaps/stops/unknowns, 288k nucleotide output | 0.0467s | 0.0133s | 3.5x |
| `DNAThreader._thread_sequence` | 5k protein/nucleotide pairs, 240 amino acids each with gaps/stops/unknowns | 2.0228s | 1.0611s | 1.9x |
| `DNAThreader._thread_sequence` full-length mask iteration | 5k protein/nucleotide pairs, 240 amino acids each, full-length mixed mask | 0.4463s | 0.4073s | 1.10x |
| `DNAThreader._thread_sequence` no-stop masked output | 5k protein/nucleotide pairs, 240 amino acids each, mixed mask without terminal stop restoration | 1.291884s | 1.173806s | 1.10x |
| `DNAThreader._thread_sequence` grouped mixed mask plan | 5k protein/nucleotide pairs, 240 amino acids each with gaps/unknowns, repeated ClipKIT-style mask | 0.981049s | 0.762770s | 1.29x |
| `DNAThreader._thread_sequence` no-stop direct chunk threading | 5k protein/nucleotide pairs, 240 amino acids each with gaps/unknowns, mixed mask without terminal stop restoration | 0.349655s | 0.338762s | 1.03x |
| `DNAThreader._thread_sequence` terminal-stop direct threading | 5k protein/nucleotide pairs, 240 amino acids each with gaps/unknowns, mixed mask and terminal stop restoration | 0.452564s | 0.349399s | 1.3x |
| `DNAThreader._thread_sequence` precomputed emit plan | 5k protein/nucleotide pairs, 240 amino acids each with gaps/unknowns, repeated mixed mask without terminal stop restoration | 0.754961s | 0.332402s | 2.27x |
| `DNAThreader._thread_sequence` terminal-stop emit plan | 5k protein/nucleotide pairs, 240 amino acids each with gaps/unknowns, repeated mixed mask and terminal stop restoration | 0.757118s | 0.330857s | 2.29x |
| `DNAThreader._thread_sequence` uniform emit-plan run grouping | 5k protein/nucleotide pairs, 240 amino acids each with gaps/unknowns, repeated mixed mask without terminal stop restoration | 0.334219s | 0.225659s | 1.48x |
| `DNAThreader._thread_sequence` terminal-stop uniform emit-plan run grouping | 5k protein/nucleotide pairs, 240 amino acids each with gaps/unknowns, repeated mixed mask and terminal stop restoration | 0.324190s | 0.224287s | 1.45x |
| `DNAThreader._thread_sequence` partial-tail uniform prefix emit plan | 5k calls, 2200 amino acids with gaps/unknowns, repeated mixed mask, no terminal stop, final partial emitter group | 1.009760s | 0.701732s | 1.44x |
| `DNAThreader._thread_sequence` no-gap masked prefix | 5k protein/nucleotide pairs, 240 amino acids each, mixed mask without gap or stop characters | 0.405581s | 0.135780s | 3.0x |
| `DNAThreader._thread_sequence` no-gap masked prefix compression | 5k protein/nucleotide pairs, 240 amino acids each, mixed mask without gap or stop characters | 0.127036s | 0.047419s | 2.7x |
| `DNAThreader._thread_sequence` no-gap prefix gap-symbol scan | 5k calls, 240-aa no-gap protein, repeated mixed mask, optimized helper baseline | 0.125954s | 0.093268s | 1.35x |
| `DNAThreader._thread_sequence` prefix materialization | 200k amino acids + stop, triplet-aligned DNA, ClipKIT-like mask | 0.1382s | 0.0844s | 1.6x |
| `DNAThreader._thread_sequence` all-keep mask | 200k amino acids + stop, triplet-aligned DNA, default all-true mask | 0.0908s | 0.0331s | 2.7x |
| `DNAThreader._thread_sequence` no-gap all-keep mask | 300k amino acids, triplet-aligned DNA, default all-true mask | 0.0465s | 0.0085s | 5.5x |
| `DNAThreader._thread_sequence` precomputed all-keep mask hint | 5k protein/nucleotide pairs, 240 amino acids each, no ClipKIT mask | 0.022675s | 0.013099s | 1.73x |
| `DNAThreader._thread_sequence` all-keep gappy fallback | 5k calls, 2k amino acids with gaps/stops/unknowns, default all-true mask | 3.146158s | 0.877293s | 3.59x |
| `DNAThreader._thread_all_sites_kept` run grouping | 5k calls, 2k amino acids with clustered gaps/stops/unknowns, default all-true mask | 2.065739s | 0.669193s | 3.09x |
| `DNAThreader._create_thread_mask_plan` full/empty group emitters | 60k amino acids, alternating full-keep and full-trim 9-site groups | 0.022028s | 0.012502s | 1.76x |
| `DNAThreader._thread_sequence` full-trim plan groups | 60k amino acids, repeated full-trim/full-keep/gap groups with shared mask plan | 0.006454s | 0.005004s | 1.29x |
| `DNAThreader.create_mask` ClipKIT log setup | 200k ClipKIT rows expanded to nucleotide mask | 0.1014s | 0.0649s | 1.6x |
| `DNAThreader.create_mask` streaming ClipKIT status parse | 200k ClipKIT rows expanded to nucleotide mask | 0.091483s | 0.049973s | 1.8x |
| `DNAThreader.create_mask` binary ClipKIT status parse | 200k ClipKIT rows expanded to nucleotide mask | 0.051766s | 0.049098s | 1.05x |
| `DNAThreader.create_mask` constant triplet expansion | 200k ClipKIT rows expanded to nucleotide mask | 0.066043s | 0.057353s | 1.15x |
| `DNAThreader.thread` ClipKIT all-keep detection | 500k all-keep ClipKIT rows expanded to nucleotide mask and all-true flag | 0.128622s | 0.123556s | 1.04x |
| `DNAThreader.create_mask` bounded byte status comparison | 200k all-keep ClipKIT rows expanded to nucleotide mask, identical single-space status parsing | 0.128046s | 0.104242s | 1.23x |
| `DNAThreader.clipkit_log_data` streaming row split | 300k ClipKIT rows, identical parsed row lists | 0.058648s | 0.054238s | 1.08x |
| `DNAThreader.run` batched FASTA text output | 100k threaded FASTA records, mocked parser/threading and identical stdout text | 0.030984s | 0.011490s | 2.70x |
| `dna_threader` module import without eager Bio.SeqIO/Seq | cold subprocess import after lazy Biopython sequence imports | 0.187897s | 0.084045s | 2.24x |
| `dna_threader` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006269s | 0.004938s | 1.27x |
| `dna_threader` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.032400s | 0.031734s | 1.02x |
| `MemoryEfficientAlignmentProcessor.calculate_column_stats_streaming` one-pass column stats | 300 FASTA records x 1500 sites, alphabet `ACGT-N` | 1.183170s | 0.026735s | 44.26x |
| `StreamingFastaReader.get_sequence_count` mmap header scan | 300k FASTA records, 20 bp each, identical sequence count | 0.123812s | 0.066868s | 1.85x |
| `streaming` module import without eager Bio.SeqIO | cold subprocess import after localizing FASTA parser imports | 0.099270s | 0.002314s | 42.90x |
| `streaming` module import without `typing` startup | median cold subprocess import after replacing annotation-only typing aliases with built-in annotations | 0.022866s | 0.021220s | 1.08x |
| `parallel` module import without eager multiprocessing/NumPy/executors | cold subprocess import after lazy worker, array, and executor proxies | 0.070840s | 0.001929s | 36.72x |
| `parallel` module import without `typing` startup | median cold subprocess import after converting annotation-only collection/callable/optional aliases to built-in postponed annotations | 0.024022s | 0.020614s | 1.17x |
| `ParallelProcessor.get_optimal_workers` cached CPU count | 500k repeated worker-selection calls with varying data sizes, identical worker counts | 2.401188s | 0.699817s | 3.43x |
| `NumpyParallel.parallel_pairwise_operation` explicit sequential symmetric path | 1200 items, symmetric pairwise absolute-difference matrix, `num_workers=1`, identical matrix output | 0.328311s | 0.157943s | 2.08x |
| `ParallelProcessor.parallel_reduce` no-initial iterator reduction | 1M-item sequential identity map plus additive reduce, side-by-side previous `results[1:]` slice | 0.209717s | 0.151342s | 1.39x |
| `NumpyParallel.parallel_apply_along_axis` lazy index range | 1 x 1M array, column-wise sequential apply with identical output, side-by-side previous eager `list(range(...))` index stream | 1.033323s | 0.788665s | 1.31x |
| `CreateConcatenationMatrix._process_alignment_file` | 2500-record FASTA, 2750 requested taxa, 120 sites | 0.0566s | 0.0017s | 33.3x |
| `CreateConcatenationMatrix._process_alignment_file` all-present FASTA return | 8000-record FASTA, all requested taxa present, 120 sites | 0.008992s | 0.006960s | 1.29x |
| `CreateConcatenationMatrix._process_alignment_file` lightweight FASTA parser | 8000-record FASTA, all requested taxa present, 20 sites split across two lines | 0.005270s | 0.004516s | 1.17x |
| `CreateConcatenationMatrix.get_list_of_taxa_and_records` | 50k FASTA records, mixed-case 120 bp each | 0.2495s | 0.0659s | 3.8x |
| `CreateConcatenationMatrix.get_list_of_taxa_and_records` slots records | 50k FASTA records with descriptions, 120 bp each | 0.052496s | 0.041054s | 1.3x |
| `CreateConcatenationMatrix.get_list_of_taxa_and_records` shared first-token records parser | 50k FASTA records with descriptions, mixed-case 120 bp split across two lines, legacy `SimpleFastaParser` baseline | 0.080893s | 0.069486s | 1.16x |
| `CreateConcatenationMatrix._compute_effective_occupancy` | 900 taxa x 80 genes x 300 sites, threshold filtering invalid symbols | 0.6085s | 0.1985s | 3.1x |
| `CreateConcatenationMatrix._compute_effective_occupancy` joined counts | 800 taxa x 300 genes x 180 sites, threshold filtering invalid symbols | 0.4558s | 0.2637s | 1.7x |
| `CreateConcatenationMatrix._compute_effective_occupancy` byte deletion | 800 taxa x 300 genes x 180 sites, threshold filtering invalid symbols | 0.2376s | 0.1540s | 1.5x |
| `CreateConcatenationMatrix._build_occupancy_state_matrix` | 800 taxa x 220 genes x 180 sites, occupancy plot state matrix | 2.8920s | 0.2351s | 12.3x |
| `CreateConcatenationMatrix._build_occupancy_state_matrix` batched gene lookup | 800 taxa x 220 genes x 180 sites, 80% present, occupancy plot state matrix | 0.249490s | 0.061994s | 4.0x |
| `CreateConcatenationMatrix._plot_concatenation_occupancy` gene-boundary rendering | 4096 gene boundary lines, real Matplotlib Agg render | 0.796482s | 0.055704s | 14.30x |
| `CreateConcatenationMatrix._plot_concatenation_occupancy` represented-row counts | 6000 taxa x 4000 concatenated-position state matrix, side-by-side previous boolean `np.sum(..., axis=1)` | 0.010258s | 0.007614s | 1.35x |
| `CreateConcatenationMatrix.add_to_occupancy_info` cached taxa | 800 occupancy rows over 6000 taxa, sorted missing-taxa lists | 0.2575s | 0.1625s | 1.6x |
| `CreateConcatenationMatrix.process_taxa_sequences` cached taxa set | 600 sequential alignments, 8000 taxa, 80% present per alignment | 0.7268s | 0.6294s | 1.2x |
| `CreateConcatenationMatrix.process_taxa_sequences` combined present scan | 600 sequential alignments, 8000 taxa, 80% present per alignment | 0.7095s | 0.5664s | 1.3x |
| `CreateConcatenationMatrix.process_taxa_sequences` ordered unique-taxa missing scan | 600 sequential alignments, 8000 taxa, 80% present per alignment | 0.585364s | 0.548891s | 1.07x |
| `CreateConcatenationMatrix.process_taxa_sequences` parsed present-taxa reuse | 600 sequential alignments, 8000 taxa, 80% present per alignment | 1.668222s | 1.438374s | 1.2x |
| `CreateConcatenationMatrix.process_taxa_sequences` all-present fast return | 600 sequential alignments, 8000 taxa, all taxa present per alignment | 1.575316s | 1.338536s | 1.2x |
| `CreateConcatenationMatrix.process_taxa_sequences` slots records | 600 sequential alignments, 8000 taxa, 80% present per alignment | 0.487060s | 0.377785s | 1.3x |
| `CreateConcatenationMatrix.process_taxa_sequences` all-present slots records | 600 sequential alignments, 8000 taxa, all taxa present per alignment | 0.442465s | 0.302022s | 1.5x |
| `CreateConcatenationMatrix.process_taxa_sequences` parsed-record append fast path | 600 sequential alignments, 8000 taxa, parser-created records, all taxa present | 0.298402s | 0.198637s | 1.5x |
| `CreateConcatenationMatrix.process_taxa_sequences` parsed-record append fast path with missing taxa | 600 sequential alignments, 8000 taxa, parser-created records, 80% present per alignment | 0.346466s | 0.261893s | 1.3x |
| `CreateConcatenationMatrix` shared missing-taxa list | 600 sequential alignments, 8000 taxa, parser-created records, 80% present per alignment plus occupancy rows | 0.487739s | 0.397005s | 1.23x |
| `CreateConcatenationMatrix` ordered append cached taxon lists | 600 ordered alignment dictionaries, 8000 taxa, 80% represented per alignment | 0.221014s | 0.195735s | 1.13x |
| `CreateConcatenationMatrix.fasta_file_write` batched FASTA rows | 80k taxa x 80 sequence chunks, identical FASTA text via generated `writelines` rows | 0.194861s | 0.177925s | 1.10x |
| `CreateConcatenationMatrix.fasta_file_write` chunked FASTA row writes | 80k taxa x 80 sequence chunks, identical FASTA text, side-by-side previous generator `writelines` writer | 0.231785s | 0.178997s | 1.29x |
| `CreateConcatenationMatrix` threshold exclusion text output | 100k excluded taxa rows, captured stdout and identical text | 0.040259s | 0.028425s | 1.42x |
| `CreateConcatenationMatrix._get_taxa_from_alignment` | 50k FASTA records, 12 bp each | 0.0441s | 0.0238s | 1.9x |
| `CreateConcatenationMatrix._get_taxa_from_alignment` header-only parser | 50k FASTA records, mixed-case 120 bp each, legacy `SimpleFastaParser` baseline | 0.060259s | 0.038296s | 1.57x |
| `create_concatenation_matrix` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.206190s | 0.135351s | 1.52x |
| `create_concatenation_matrix` module import without eager NumPy/concurrency helpers | cold subprocess import after lazy NumPy occupancy lookup, concurrency, JSON, and plot config helpers | 0.096497s | 0.025028s | 3.86x |
| `create_concatenation_matrix` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.035103s | 0.033782s | 1.04x |
| `create_concatenation_matrix` module import without eager file helper | median cold subprocess import, interleaved lazy import vs eager-equivalent `phykit.helpers.files` preload | 0.028115s | 0.028138s | 1.00x |
| `create_concatenation_matrix` module import without eager `textwrap` | median cold subprocess import after lazy verbose-message dedent helper | 0.047605s | 0.041895s | 1.14x |
| `TaxonGroups._extract_taxa` FASTA mode | 50k FASTA records, 12 bp each | 0.0808s | 0.0182s | 4.4x |
| `TaxonGroups._extract_taxa` FASTA header-only parser | 50k FASTA records, mixed-case 120 bp each, legacy `SimpleFastaParser` baseline | 0.039226s | 0.022543s | 1.74x |
| `TaxonGroups._extract_taxa` tree terminal-name extraction | parsed balanced 65536-tip Newick tree, collect terminal names | 0.1344s | 0.0129s | 10.4x |
| `TaxonGroups._read_file_list` streaming parse | 500k relative paths with comments/blanks | 1.617507s | 1.484327s | 1.09x |
| `TaxonGroups._read_file_list` Path-compatible string resolver | 500k relative/absolute paths with comments/blanks and redundant separators, old `Path`-per-row parser baseline | 1.587105s | 0.320406s | 4.95x |
| `TaxonGroups._read_file_list` simple-path normalization fast path | 500k relative/absolute paths with comments/blanks and redundant separators, side-by-side previous string resolver comparison | 0.313885s | 0.178136s | 1.76x |
| `TaxonGroups._read_file_list` separator-guarded normalization checks | 500k simple relative paths with extensions, side-by-side previous string resolver comparison | 0.171839s | 0.134575s | 1.28x |
| `TaxonGroups._read_file_list` no-separator fast path | 500k simple relative paths with extensions and comments/blanks, side-by-side previous separator-guarded parser comparison | 0.139625s | 0.104450s | 1.34x |
| `TaxonGroups.run` direct taxon-set grouping | 50k mocked file paths, 100 recurring taxon sets, grouping work isolated | 0.047093s | 0.026018s | 1.81x |
| `TaxonGroups.run` group-size sort key | 200k mocked taxon groups with varied group sizes, side-by-side previous negated-length sort key | 0.063989s | 0.052508s | 1.22x |
| `TaxonGroups.run` batched text report output | 10k one-file taxon groups, mocked extraction and identical stdout text | 0.014605s | 0.010450s | 1.40x |
| `TaxonGroups._extract_taxa` existence guard | 50k existing small file paths, taxa parser mocked | 0.257713s | 0.115714s | 2.23x |
| `taxon_groups` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006413s | 0.004990s | 1.29x |
| `taxon_groups` module import without tree-base startup | median cold subprocess import after localizing tree-only helpers to tree extraction | 0.004458s | 0.000601s | 7.42x |
| `OccupancyFilter._extract_taxa` FASTA mode | 50k FASTA records, 12 bp each | 0.0744s | 0.0195s | 3.8x |
| `OccupancyFilter._extract_taxa` FASTA header-only parser | 50k FASTA records, mixed-case 120 bp each, legacy `SimpleFastaParser` baseline | 0.039226s | 0.027756s | 1.41x |
| `OccupancyFilter._extract_taxa` tree terminal-name extraction | parsed balanced 65536-tip Newick tree, collect terminal names | 0.1349s | 0.0128s | 10.5x |
| `OccupancyFilter._extract_taxa` existence guard | 50k existing small file paths, taxa parser mocked | 0.257713s | 0.115714s | 2.23x |
| `OccupancyFilter._read_file_list` streaming parse | 500k relative paths with comments/blanks | 1.608151s | 1.482792s | 1.08x |
| `OccupancyFilter._read_file_list` Path-compatible string resolver | 500k relative/absolute paths with comments/blanks and redundant separators, old `Path`-per-row parser baseline | 1.587105s | 0.312125s | 5.08x |
| `OccupancyFilter._read_file_list` simple-path normalization fast path | 500k relative/absolute paths with comments/blanks and redundant separators, side-by-side previous string resolver comparison | 0.315469s | 0.177995s | 1.77x |
| `OccupancyFilter._read_file_list` separator-guarded normalization checks | 500k simple relative paths with extensions, side-by-side previous string resolver comparison | 0.171840s | 0.135726s | 1.27x |
| `OccupancyFilter._read_file_list` no-separator fast path | 500k simple relative paths with extensions and comments/blanks, side-by-side previous separator-guarded parser comparison | 0.135912s | 0.088443s | 1.54x |
| `OccupancyFilter._filter_fasta` | 30k FASTA records x 120 bp, 15k retained | 0.0507s | 0.0281s | 1.8x |
| `OccupancyFilter._filter_fasta` streaming retained-sequence parser | 50k FASTA records x 120 bp, 10k retained, legacy `SimpleFastaParser` baseline | 0.050271s | 0.036378s | 1.38x |
| `OccupancyFilter._write_wrapped_fasta_sequence` | 25k retained sequences x 1200 bp, 60-char wrapping | 0.0686s | 0.0553s | 1.2x |
| `OccupancyFilter._write_wrapped_fasta_sequence` list-backed chunks | 25k retained sequences x 1200 bp, 60-char wrapping | 0.0799s | 0.0553s | 1.4x |
| `OccupancyFilter.run` threshold classification without dead sort | 500k shuffled taxa occupancy counts, classify kept/removed sets | 0.186998s | 0.057841s | 3.23x |
| `OccupancyFilter.run` batched occupancy Counter updates | 200 input files x 20k taxa sampled from 500k taxa, identical occupancy counts | 0.843763s | 0.414411s | 2.04x |
| `OccupancyFilter._filter_tree` | balanced 2048-tip Newick tree, prune 1024 tips and write filtered tree | 0.5579s | 0.4044s | 1.4x |
| `OccupancyFilter._filter_tree` batch standard-tree pruning | balanced 8192-tip Newick tree, prune 4096 tips and write filtered tree | 3.1946s | 0.0453s | 70.5x |
| `OccupancyFilter._filter_tree` terminal scan child push | balanced 131072-tip tree, collect terminals and 65536 prune targets, side-by-side previous `reversed(children)` scan | 0.072382s | 0.061099s | 1.18x |
| `OccupancyFilter._print_text` batched report output | 50k taxa, 1k output files, mocked occupancy data and identical stdout text | 0.024754s | 0.016128s | 1.53x |
| `OccupancyFilter._print_text` direct occupancy-key sorting | 200k taxa x 12 files report, side-by-side previous `sorted(occupancy.keys())` loop | 0.494521s | 0.454472s | 1.09x |
| `occupancy_filter` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.184065s | 0.065460s | 2.81x |
| `occupancy_filter` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006182s | 0.004883s | 1.27x |
| `occupancy_filter` module import without tree-base/typing startup | median cold subprocess import after localizing tree-only helpers and replacing annotation-only typing aliases | 0.005686s | 0.000560s | 10.16x |
| `RenameFastaEntries.replace_ids_in_file_and_write` | 30k FASTA records x 120 bp, 15k renamed | 0.0710s | 0.0376s | 1.9x |
| `RenameFastaEntries.replace_ids_in_file_and_write` combined record write | 30k FASTA records x 120 bp, 15k renamed | 0.060598s | 0.052972s | 1.14x |
| `RenameFastaEntries.replace_ids_in_file_and_write` short-sequence no-wrap path | 120k FASTA records x 50 bp, half renamed, identical headers and sequence output | 0.087410s | 0.072161s | 1.21x |
| `RenameFastaEntries.load_idmap` explicit split loop | 500k two-column ID-map rows, identical whitespace parsing and duplicate-key behavior | 1.323420s | 1.219563s | 1.09x |
| `RenameFastaEntries._write_wrapped_fasta_sequence` | 25k renamed sequences x 1200 bp, 60-char wrapping | 0.0692s | 0.0597s | 1.2x |
| `RenameFastaEntries._write_wrapped_fasta_sequence` list-backed chunks | 25k renamed sequences x 1200 bp, 60-char wrapping | 0.0598s | 0.0541s | 1.1x |
| `rename_fasta_entries` module import without eager Bio.SeqIO | cold subprocess import after lazy Bio.SeqIO/FastaIO imports | 0.186867s | 0.117185s | 1.59x |
| `rename_fasta_entries` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006248s | 0.005029s | 1.24x |
| `rename_fasta_entries` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.002786s | 0.000976s | 2.86x |
| `Faidx.run` multi-entry text output | 50k indexed FASTA entries x 120 bp, in-memory text stream | 0.015406s | 0.009832s | 1.57x |
| `Faidx.run` streaming FASTA fetch | 50k FASTA records x 120 bp, three requested entries | 0.333717s | 0.145065s | 2.30x |
| `Faidx._fetch_entries` selected-entry parser | 50k FASTA records x 120 bp, three requested entries, legacy `SimpleFastaParser` baseline | 0.043153s | 0.026709s | 1.62x |
| `Faidx.run` direct sequence mapping | 100k requested FASTA entries x 120 bp, side-by-side previous `_FastaEntry` wrapper allocation/output path | 0.204819s | 0.068155s | 3.01x |
| shared `_fasta._clean_sequence` single-line fast path | 80k FASTA records x 120 bp, single-line / wrapped two-line sequences, identical first-token parser output | 0.159723s / 0.209468s | 0.116201s / 0.163365s | 1.37x / 1.28x |
| `faidx` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.181946s | 0.079052s | 2.30x |
| `faidx` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006110s | 0.004992s | 1.22x |
| `faidx` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only typing aliases to built-in annotations | 0.002548s | 0.000926s | 2.75x |
| `AlignmentSubsample._read_alignment` | 50k FASTA records, 12 bp each | 0.0508s | 0.0244s | 2.1x |
| `AlignmentSubsample._read_alignment` shared case-preserving parser | 50k FASTA records, mixed-case 120 bp each, legacy `SimpleFastaParser` baseline | 0.047439s | 0.036622s | 1.30x |
| `AlignmentSubsample._run_sites` site selection | 500 taxa x 10k sites, 5k sampled sites without replacement | 0.0911s | 0.0324s | 2.8x |
| `AlignmentSubsample._run_sites` full-site non-bootstrap shortcut | 500 taxa x 10000 sites, selected site count equals alignment length, mocked FASTA write and summary output | 0.076459s | 0.000031s | 2459.84x |
| `AlignmentSubsample._assemble_partition_subsample` | 900 taxa x 800 partitions x 80 sites, 600 selected partitions with duplicates | 0.1243s | 0.0502s | 2.5x |
| `AlignmentSubsample._print_summary` batched text output | seven captured summaries with 100k output-file rows each, identical stdout text | 0.118185s | 0.050466s | 2.34x |
| `AlignmentSubsample._write_fasta` chunked output | 1M FASTA records x 40 bp, identical output file text, previous per-record write baseline | 0.273654s | 0.193250s | 1.42x |
| `AlignmentSubsample._write_partition_file` chunked output | 1M RAxML-style partition rows, identical output file text, previous per-row write baseline | 0.276847s | 0.200761s | 1.38x |
| `AlignmentSubsample._parse_partition_file` split parser | 100k RAxML-style partition rows, preserving comments/invalid-row skips and trailing text tolerance | 0.264859s | 0.155592s | 1.70x |
| `alignment_subsample` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.188059s | 0.110484s | 1.70x |
| `alignment_subsample` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.008061s | 0.006792s | 1.19x |
| `alignment_subsample` module import without eager random | median cold subprocess import after localizing command-run sampling import | 0.007269s | 0.005914s | 1.23x |
| `alignment_subsample` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.035616s | 0.028305s | 1.26x |
| `AlignmentSubsample._run_genes` selected path writer | 100k selected gene paths, identical output file text | 0.007651s | 0.001694s | 4.52x |
| `AlignmentRecoding.recode_alignment` | 500 taxa x 10k sites, RY nucleotide recoding with lowercase/gap symbols | 0.5274s | 0.0180s | 29.3x |
| `AlignmentRecoding.run` string-backed recoded FASTA assembly | 500 taxa x 10k sites, RY nucleotide recoding with lowercase/gap symbols and mocked stdout | 0.039963s | 0.004749s | 8.42x |
| `AlignmentRecoding.run` batched FASTA text output | 100k recoded FASTA records, mocked alignment/read and identical stdout text | 0.066343s | 0.055510s | 1.20x |
| `AlignmentRecoding.read_recoding_table` bounded split | 1M custom recoding rows with ignored trailing columns, identical recoding dictionary | 0.213408s | 0.163705s | 1.30x |
| `AlignmentRecoding._build_translation_table` direct case inserts | 200k translation-table builds from uppercase/lowercase recoding symbols plus DNA gap symbols, side-by-side previous per-symbol two-item set builder | 6.007475s | 5.653971s | 1.06x |
| `alignment_recoding` module import without eager Bio.Align/json helpers | cold subprocess import after postponed annotations and lazy JSON helper wrapper | 0.114186s | 0.026325s | 4.34x |
| `alignment_recoding` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.002459s | 0.000894s | 2.75x |
| `PhyloGwas._read_fasta` | 50k FASTA records, lowercase 12 bp each | 0.0539s | 0.0244s | 2.2x |
| `PhyloGwas._read_fasta` shared first-token parser | 50k FASTA records, lowercase 120 bp each, legacy `SimpleFastaParser` baseline | 0.051468s | 0.040648s | 1.27x |
| `PhyloGwas._read_phenotype` bounded row split | 1M phenotype TSV rows with comments/blanks and three ignored trailing columns, side-by-side previous full-row split | 0.959486s | 0.867587s | 1.11x |
| `PhyloGwas._read_phenotype` binary bounded row split | 1M phenotype TSV rows with comments/blanks and three ignored trailing columns, identical first-two-column/duplicate output | 0.402401s | 0.382109s | 1.05x |
| `PhyloGwas` allele extraction | 600 taxa x 5k ASCII sites, no ambiguous columns | 0.2494s | 0.0264s | 9.4x |
| `PhyloGwas` ASCII valid-column prefilter | 600 taxa x 50k ASCII sites, 1% ambiguous character rate | 0.116165s | 0.033381s | 3.48x |
| `PhyloGwas._benjamini_hochberg` | 1M synthetic p-values | 0.3779s | 0.1207s | 3.1x |
| `PhyloGwas.run` zip-based FDR result annotation | 1M site result rows with adjusted p-value array, side-by-side previous index lookups | 0.263524s | 0.181904s | 1.45x |
| `PhyloGwas` partition gene lookup | 10k result positions across 10k sorted partitions | 1.3198s | 0.0028s | 472.5x |
| `PhyloGwas._parse_partition_file` cached regex matcher | 300k RAxML-style partition rows with comments/blanks, whitespace tolerance, and trailing text tolerance | 0.404678s | 0.344136s | 1.18x |
| `PhyloGwas._test_site_categorical` | 5k biallelic sites x 600 taxa x 6 phenotype groups | 1.8773s | 0.3523s | 5.3x |
| `PhyloGwas._test_site_categorical` ASCII allele counting | 5k biallelic ASCII sites x 600 taxa x 6 phenotype groups | 0.207769s | 0.177623s | 1.17x |
| `PhyloGwas._test_site_categorical` ASCII byte-column path | 5k biallelic ASCII sites x 600 taxa x 6 phenotype groups | 0.096297s | 0.063694s | 1.5x |
| `PhyloGwas._test_site_categorical` byte histogram allele counting | 5k biallelic ASCII sites x 600 taxa x 6 phenotype groups | 0.165865s | 0.133463s | 1.24x |
| `PhyloGwas._test_site_categorical` grouped byte histogram | 5k biallelic ASCII sites x 600 taxa x 6 phenotype groups | 0.131545s | 0.123580s | 1.06x |
| `PhyloGwas._test_site_categorical` two-group Fisher exact p-value | 5k biallelic ASCII sites x 600 taxa x 2 phenotype groups | 1.793960s | 0.442959s | 4.05x |
| `PhyloGwas._test_site_categorical` cached two-group Fisher exact p-values | 5k biallelic ASCII sites x 600 taxa x 2 phenotype groups, repeated row totals | 0.452611s | 0.315611s | 1.43x |
| `PhyloGwas._test_site_categorical` two-group byte offsets | 5k biallelic ASCII sites x 600 taxa x 2 phenotype groups, repeated row totals | 0.293545s | 0.273126s | 1.07x |
| `PhyloGwas._test_site_categorical` cached Fisher log-combinations | 5k biallelic ASCII sites x 600 taxa x 2 phenotype groups, repeated row totals | 0.092839s | 0.082627s | 1.12x |
| `PhyloGwas._test_site_categorical` prepared two-group byte counts | 5k biallelic ASCII sites x 600 taxa x 2 phenotype groups, repeated row totals | 0.095170s | 0.060784s | 1.57x |
| `PhyloGwas` byte major/minor ndarray minmax | 20 / 100 / 600 / 5k / 10k / 100k byte alleles, shared categorical/continuous helper, side-by-side previous `np.min`/`np.max` wrappers | 5.198329s / 2.533913s / 1.377197s / 0.196970s / 0.126234s / 0.012328s | 2.893385s / 1.602335s / 0.941029s / 0.052535s / 0.054851s / 0.009868s | 1.80x / 1.58x / 1.46x / 3.75x / 2.30x / 1.25x |
| `PhyloGwas._test_site_continuous` | 5k biallelic sites x 600 taxa, continuous phenotype | 1.1932s | 0.1200s | 9.9x |
| `PhyloGwas._test_site_continuous` ASCII byte-column path | 5k biallelic ASCII sites x 600 taxa, continuous phenotype | 0.122781s | 0.094730s | 1.3x |
| `PhyloGwas._test_site_continuous` byte-count formula | 5k biallelic ASCII sites x 600 taxa, continuous phenotype | 0.094610s | 0.056158s | 1.7x |
| `PhyloGwas._test_site_continuous` byte min/max allele detection | 5k biallelic ASCII sites x 600 taxa, continuous phenotype | 0.069616s | 0.060067s | 1.16x |
| `PhyloGwas._test_site_continuous` scalar math p-value path | 5k biallelic ASCII sites x 600 taxa, continuous phenotype | 0.084417s | 0.075367s | 1.12x |
| `PhyloGwas._test_site_continuous` byte phenotype dot product | 6k biallelic ASCII sites x 10k taxa, continuous phenotype | 0.351814s | 0.182585s | 1.93x |
| `PhyloGwas._prune_tree_to_taxa` | balanced 1024-tip tree, prune 512 tips before phylogenetic classification | 0.8222s | 0.1153s | 7.1x |
| `PhyloGwas._prune_tree_to_taxa` no-op shared-taxa setup | balanced 65536-tip tree, all taxa retained before phylogenetic classification | 0.176819s | 0.040796s | 4.33x |
| `PhyloGwas._terminal_clades` order-preserving child push | balanced 131072-tip tree, terminal helper for prune setup, side-by-side previous `reversed(children)` helper | 0.038431s | 0.033913s | 1.13x |
| `PhyloGwas._classify_phylo_pattern` | balanced 1024-tip tree, 5000 significant-site derived-taxon classifications | 7.7410s | 0.0062s | 1245.4x |
| `PhyloGwas._build_phylo_pattern_index` direct postorder | balanced 32768-tip tree, descendant taxon-set index for monophyly classification | 0.276455s | 0.187328s | 1.48x |
| `PhyloGwas._build_phylo_pattern_index` reverse-preorder helper | balanced 32768-tip tree, descendant taxon-set index for monophyly classification | 0.068807s | 0.055296s | 1.24x |
| `PhyloGwas._build_phylo_pattern_index` binary child union | balanced 65536-tip tree, descendant taxon-set index for monophyly classification, side-by-side previous generic child-set union | 0.189455s | 0.179548s | 1.06x |
| `PhyloGwas._benjamini_hochberg` in-place adjustment allocation | 1M synthetic p-values | 0.108093s | 0.106394s | 1.02x |
| `PhyloGwas._create_manhattan_plot` partition span rendering | 4096 partitions / 2048 alternating gray spans, real Matplotlib Agg render | 0.570093s | 0.089715s | 6.35x |
| `PhyloGwas._create_manhattan_plot` point series preparation | 1M GWAS result rows with mixed FDR-significant phylogenetic patterns, identical positions, colors, transformed p-values, and FDR threshold | 0.155766s | 0.135860s | 1.15x |
| `PhyloGwas._print_text_output` top-hit selection and batched report output | 50k categorical text reports with 1k shuffled significant sites, captured stdout and identical text | 5.129769s | 3.727684s | 1.38x |
| `PhyloGwas._print_text_output` cached p-value key | 50k categorical text reports with 1k shuffled significant sites, captured stdout and identical text | 3.838282s | 3.414573s | 1.12x |
| `PhyloGwas` JSON result payload assembly | 200k categorical GWAS result rows, identical `json.dumps(sort_keys=True)` output | 0.575563s | 0.402213s | 1.43x |
| `PhyloGwas` significant-result summary | 1M GWAS result rows, two-thirds FDR-significant with mixed phylogenetic patterns | 0.082049s | 0.066148s | 1.24x |
| `PhyloGwas._write_csv` tuple-row writer | 100k categorical GWAS rows with ignored extra fields and quoted gene names | 0.222516s | 0.170786s | 1.30x |
| `PhyloGwas._write_csv` required-field fast path | 100k categorical GWAS rows with ignored extra fields and quoted gene names, side-by-side previous `get()` row writer with incomplete-row fallback preserved | 0.698594s | 0.646321s | 1.08x |
| `PhyloGwas.run` ASCII biallelic column prefilter | 400 taxa x 50,000 sites, continuous phenotype scan, 90% invariant / 5% multiallelic / 5% biallelic ASCII columns | 0.352346s | 0.119668s | 2.94x |
| `PhyloGwas.run` phylo-pattern minor-taxa scan | 400 taxa x 50,000 sites with 5000 significant tree-classified sites, side-by-side previous temporary allele-list filter | 2.555218s | 1.682205s | 1.52x |
| `PhyloGwas.run` shared-taxon setup | 1M alignment taxa x 1M phenotype taxa with 750k overlap, identical sorted shared taxa | 3.113728s | 2.794556s | 1.11x |
| `PhyloGwas._test_site_categorical` Unicode multiallelic skip | 300k non-ASCII alleles, first three alleles multiallelic, categorical phenotype | 0.036496s | 0.000002125s | 17174.6x |
| `PhyloGwas._test_site_continuous` Unicode multiallelic skip | 300k non-ASCII alleles, first three alleles multiallelic, continuous phenotype | 0.012269s | 0.000001333s | 9204.1x |
| `phylo_gwas` module import without eager Bio.Phylo/FASTA parser | cold subprocess import after lazy Biopython parser imports | 0.225119s | 0.116552s | 1.93x |
| `phylo_gwas` module import without eager NumPy/json/plot config | cold subprocess import after lazy NumPy proxy plus JSON and plot config helper deferral | 0.075720s | 0.024089s | 3.14x |
| `phylo_gwas` module import without shared typing startup | median cold subprocess import after removing annotation-only `typing` imports from shared error/alignment helpers | 0.007420s | 0.006305s | 1.18x |
| `Alignment.calculate_rcv` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.0766s | 0.0292s | 2.6x |
| `Alignment.calculate_rcv` invalid lookup | 500 taxa x 8000 sites, alphabet `ACGT-?NX*` | 0.063539s | 0.028435s | 2.23x |
| `Alignment.calculate_rcv` protein byte bincounts | 1200 taxa x 12000 sites, protein alphabet plus gaps/ambiguous symbols, side-by-side previous per-symbol scan | 0.157945s | 0.130001s | 1.21x |
| `Alignment.calculate_rcv` observed-symbol validity DNA | 2000 taxa x 5000 sites, alphabet `ACGT-?NX*`, side-by-side previous full valid-mask path | 0.080328s | 0.047165s | 1.70x |
| `Alignment.calculate_rcv` observed-symbol validity protein | 2000 taxa x 5000 sites, protein alphabet plus gaps/ambiguous symbols, side-by-side previous full valid-mask path | 0.099660s | 0.078707s | 1.27x |
| `Alignment.calculate_rcv` no-gap observed-symbol validity | 2000 taxa x 5000 sites, DNA `ACGT` / 20 amino-acid symbols, side-by-side previous full valid-mask path | 0.051012s / 0.063966s | 0.034910s / 0.042129s | 1.46x / 1.52x |
| `Alignment.calculate_rcv` identical-sequence shortcut | 1200 taxa x 12000 sites, lowercase/uppercase identical DNA records, side-by-side previous matrix path | 0.046277s | 0.006825s | 6.78x |
| `Alignment.calculate_rcv` identical-sequence no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.163491s / 0.008818s / 0.075810s | 0.050280s / 0.000004s / 0.031403s | 3.25x / 2377.98x / 2.41x |
| `Alignment.calculate_rcv` final RCV total | median per-call reduction of 260 / 1200 / 2000 per-taxon RCV values, side-by-side previous `np.sum` wrapper | 0.000003516s / 0.000005445s / 0.000004957s | 0.000002338s / 0.000002423s / 0.000003238s | 1.50x / 2.25x / 1.53x |
| `Alignment.calculate_rcv` clean large-short ASCII count matrix | 10000 taxa x 128 sites / 50000 taxa x 64 sites / 200000 taxa x 32 sites, 20 valid symbols, side-by-side previous per-row `bincount` loop | 5.113480s / 5.967428s / 6.979634s | 1.557019s / 1.341158s / 1.625248s | 3.28x / 4.45x / 4.29x |
| `Alignment.calculate_rcv` single-record early return | 5 repeated 4.5M-site single-record RCV calls, side-by-side previous sequence materialization before zero return | 0.014365s | 0.000000417s | 34428.75x |
| `rcv` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006077s | 0.004840s | 1.26x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` | 260 taxa x 5000 sites, alphabet `ACGT-?NX*` | 0.0483s | 0.0037s | 13.1x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` byte lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?*XN` | 0.0753s | 0.0264s | 2.9x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` case-insensitive byte lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?*XN` | 0.025249s | 0.021711s | 1.16x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` ASCII matrix counts | 2000 taxa x 5000 sites, alphabet `ACGTN-?X*` | 0.020122s | 0.014212s | 1.42x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` variable-length ASCII fallback | 50k records x 90-106 bp, alphabet `ACGTN-?X*` | 0.202738s | 0.177142s | 1.14x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` variable-length byte translate | 50k records x 90-106 bp, alphabet `ACGTN-?X*` | 0.073642s | 0.033598s | 2.19x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` no-gap DNA ASCII shortcut | 2000 taxa x 5000 sites, alphabet `ACGT`, side-by-side previous full matrix lookup path | 0.025015s | 0.004553s | 5.49x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` no-gap protein ASCII shortcut | 2000 taxa x 5000 sites, 20 amino-acid symbols, side-by-side previous full matrix lookup path | 0.027621s | 0.003884s | 7.11x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` identical gappy ASCII shortcut | 1200 taxa x 12000 sites, identical mixed-symbol DNA records with gaps/ambiguous symbols, side-by-side previous byte-matrix lookup path | 0.126362s | 0.000168s | 750.48x |
| `OccupancyPerTaxon.calculate_occupancy_per_taxon` single-record direct count | 4.5M-site single-record DNA alignment, side-by-side previous record-data and matrix-helper setup | 0.022511792s | 0.010062833s | 2.24x |
| `OccupancyPerTaxon._occupancy_from_ascii_matrix` identical-row no-slice scan | 1M identical mixed-symbol DNA records, side-by-side previous `sequences[1:]` equality scan | 0.569269s | 0.392282s | 1.45x |
| `OccupancyPerTaxon._occupancy_from_ascii_matrix` combined length/identity scan | 1M mixed-symbol DNA records, identical / late-different / late variable-length cases, side-by-side previous sequence-list length pass plus identity pass | 0.288323s / 0.292319s / 0.085046s | 0.138665s / 0.282464s / 0.055508s | 2.08x / 1.03x / 1.53x |
| `OccupancyPerTaxon.run` batched text output | 50k taxon rows, mocked alignment/read and identical stdout text | 0.025741s | 0.019400s | 1.33x |
| `occupancy_per_taxon` module import without eager NumPy/json helpers | cold subprocess import after lazy NumPy proxy, lookup construction, and JSON helper wrapper | 0.081361s | 0.024702s | 3.29x |
| `occupancy_per_taxon` module import without `typing` startup | median cold subprocess import after removing annotation-only typing aliases under postponed annotations | 0.002602s | 0.000946s | 2.75x |
| `Dstatistic._run_alignment_mode` | 400k sites, ABBA/BABA/invariant/ambiguous synthetic alignment, block size 1000 | 0.4105s | 0.0089s | 46.1x |
| `Dstatistic._count_site_patterns` all-valid ASCII shortcut | 2M sites, clean ASCII ABBA/BABA/invariant synthetic alignment, block size 1000 | 0.017517s | 0.007029s | 2.49x |
| `Dstatistic._count_site_patterns` all-invariant shortcut | 2M sites, P1/P2/P3 identical to outgroup, block size 1000, identical zero totals and block arrays | 0.052473s | 0.000005s | 10494.6x |
| `Dstatistic._count_site_patterns` identical-ingroup shortcut | 2M sites, P1/P2/P3 identical and outgroup fixed different, block size 1000, identical zero totals and block arrays | 0.008103s | 0.000002s | 4051.5x |
| `Dstatistic._count_site_patterns` impossible sister-pair shortcut | 5M sites, P1/P2 identical with P3/outgroup different, block size 1000, side-by-side previous byte-array scan path | 0.013529s | 0.000002958s | 4573.80x |
| `Dstatistic._count_site_patterns` single-pattern mask shortcut | 5M sites, P1/outgroup identical and P2/P3 identical, ABBA possible and BABA impossible, side-by-side previous two-mask vector path | 0.014601s | 0.009187s | 1.59x |
| `Dstatistic._count_site_patterns_scalar` direct skip checks | 1M Unicode-containing scalar fallback sites, block size 1000, identical ABBA/BABA totals and block arrays | 0.967857s | 0.443901s | 2.18x |
| `Dstatistic._normal_two_tailed_p_value` | cold process, alignment-mode jackknife z-score p-value | 0.556567s | 0.000003125s | 178101.4x |
| `Dstatistic._jackknife_d_values` | 300k ABBA/BABA jackknife blocks | 0.0972s | 0.0022s | 43.3x |
| `Dstatistic._jackknife_d_values` block total reductions | 300k ABBA/BABA jackknife blocks, side-by-side previous `np.sum` totals | 0.000460071s | 0.000216358s | 2.13x |
| `Dstatistic` jackknife mean reduction | 10 / 20 / 50 / 100 / 1000 / 4000 / 10000 / 300k jackknife D values, side-by-side previous `np.mean` wrapper | 3.081938167s / 4.086735041s / 2.810416792s / 2.284828209s / 0.727887375s / 0.418584084s / 0.231934833s / 0.050983125s | 2.317377792s / 2.801813542s / 3.321993875s / 2.227213833s / 0.376665792s / 0.149053916s / 0.121398167s / 0.047169209s | 1.33x / 1.46x / 0.85x / 1.03x / 1.93x / 2.81x / 1.91x / 1.08x |
| `Dstatistic` jackknife standard-error sum of squares | 10 / 20 / 50 / 100 / 1000 / 4000 / 10000 jackknife D values, side-by-side previous `np.sum((x - mean) ** 2)` | 0.000008163s / 0.000007648s / 0.000007759s / 0.000007694s / 0.000008656s / 0.000017993s / 0.000020941s | 0.000002708s / 0.000002864s / 0.000002712s / 0.000003855s / 0.000003504s / 0.000004506s / 0.000017397s | 3.01x / 2.67x / 2.86x / 2.00x / 2.47x / 3.99x / 1.20x |
| `Dstatistic._read_fasta` / `Dfoil._read_fasta` | 50k FASTA records, lowercase 120 bp each | 0.0638s | 0.0360s | 1.8x |
| `Dstatistic._read_fasta` / `Dfoil._read_fasta` shared first-token parser | 50k FASTA records, lowercase 120 bp each, legacy `SimpleFastaParser` baseline | 0.047871s | 0.040750s | 1.17x |
| `Dfoil._count_site_patterns` | 475k sites, DFOIL informative/invariant/ambiguous/non-biallelic synthetic alignment | 0.3978s | 0.0073s | 54.5x |
| `Dfoil._count_site_patterns` pairwise derived-allele predicate | 1M valid ASCII sites, synthetic DFOIL informative/invariant/non-biallelic alignment | 0.008275s | 0.007467s | 1.11x |
| `Dfoil._count_site_patterns` all-valid ASCII shortcut | 2M sites, clean ASCII DFOIL informative/uninformative synthetic alignment | 0.027733s | 0.021082s | 1.32x |
| `Dfoil._count_site_patterns` small skip-code lookup mask | 1k / 5k / 10k / 100k ASCII sites with 2% skip codes, side-by-side previous six-code validity loop | 2.403774s / 1.269534s / 0.692042s / 0.224607s | 1.063043s / 0.773009s / 0.649486s / 0.209380s | 2.26x / 1.64x / 1.07x / 1.07x |
| `Dfoil._count_site_patterns` all-invariant shortcut | 2M sites, P1/P2/P3/P4 identical to outgroup, identical all-zero pattern dictionary | 0.028345s | 0.000002s | 14172.5x |
| `Dfoil._count_site_patterns_scalar` pattern-code fallback | 1M Unicode-containing scalar fallback sites, identical 16-pattern count dictionary | 2.091503s | 0.917063s | 2.28x |
| `Dstatistic._print_alignment_text_output` batched report output | 50k alignment-mode D-statistic text reports, captured stdout and identical text | 0.249669s | 0.141300s | 1.77x |
| `Dstatistic._print_gene_tree_text_output` batched report output | 50k gene-tree-mode D-statistic text reports, captured stdout and identical text | 0.243808s | 0.134732s | 1.81x |
| `Dstatistic._print_alignment_text_output` single-report formatting | 50k alignment-mode D-statistic text reports, captured stdout and identical text, side-by-side previous line-list comparison | 0.095364s | 0.074071s | 1.29x |
| `Dstatistic._print_gene_tree_text_output` single-report formatting | 50k gene-tree-mode D-statistic text reports, captured stdout and identical text, side-by-side previous line-list comparison | 0.093211s | 0.076612s | 1.22x |
| `Dstatistic._collect_clade_taxa_and_nonterminals_direct` reverse-preorder helper | balanced 32768-tip gene tree, descendant taxon sets plus nonterminal preorder | 0.056586s | 0.046768s | 1.21x |
| `Dstatistic._collect_clade_taxa_and_nonterminals_direct` binary descendant aggregation | balanced 32768-tip gene tree, descendant taxon sets plus nonterminal preorder | 0.046967s | 0.037408s | 1.26x |
| `dstatistic` module import without eager Bio.Phylo/FASTA parser | cold subprocess import after lazy Biopython parser imports | 0.219520s | 0.114720s | 1.91x |
| `dstatistic` module import without eager NumPy/json helpers | cold subprocess import after lazy NumPy proxy and JSON helper wrapper | 0.078477s | 0.025417s | 3.09x |
| `dstatistic` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.035487s | 0.031252s | 1.14x |
| `dfoil` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.207642s | 0.119953s | 1.73x |
| `dfoil` module import without eager NumPy/json helpers | cold subprocess import after lazy NumPy proxy and JSON helper wrapper | 0.074428s | 0.025339s | 2.94x |
| `dfoil` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.003288s | 0.001646s | 2.00x |
| `Dfoil._chi2_sf_df1` | cold process, four DFOIL chi-square p-values | 0.542926s | 0.000004833s | 112337.3x |
| `Dfoil._print_text_output` batched report output | 50k DFOIL text reports, captured stdout and identical text | 0.439878s | 0.310966s | 1.41x |
| `Dfoil._print_text_output` cached pattern groups and single report formatting | 50k DFOIL text reports, captured stdout and identical text, side-by-side previous pattern-list report comparison | 0.264166s | 0.232821s | 1.13x |
| `Dstatistic._get_quartet_topology` | 80 balanced gene trees x 128 taxa, quartet topology classification | 0.1713s | 0.0648s | 2.6x |
| `Dstatistic._get_quartet_topology` combined direct traversal | 80 balanced gene trees x 128 taxa, quartet topology classification | 0.059130s | 0.014816s | 3.99x |
| `Dstatistic._get_quartet_topology` quartet-only split check | balanced 32768-tip gene tree, side-by-side previous full complement set per internal bipartition | 0.174334s | 0.082270s | 2.12x |
| `Dstatistic._chi2_sf_df1` | cold process, gene-tree mode chi-square p-value | 0.551261s | 0.000002458s | 224272.2x |
| `TipToTipDistance.calculate_all_pairwise_distances` | balanced tree with 220 tips | 4.3503s | 0.0269s | 161.7x |
| `TipToTipDistance.calculate_all_pairwise_distances` tip-name setup | balanced 65536-tip tree, all-pairs terminal-name extraction | 0.129071s | 0.017473s | 7.39x |
| `Tree.calculate_pairwise_tip_distances_fast` deep-tree LCA index | pectinate 1200-tip tree, 719400 all-pairs combo/distance rows, copied old path baseline | 6.667369s | 0.862153s | 7.73x |
| `TipToTipDistance._build_distance_matrix` sorted all-pairs heatmap fill | 2500 taxa, 3,123,750 sorted upper-triangle all-pairs rows, side-by-side previous taxon-index dictionary fill | 6.848783s | 4.044856s | 1.69x |
| `TipToTipDistance.run` all-pairs text output | 200k pairwise distance rows, mocked tree/read and identical stdout text | 0.082166s | 0.057857s | 1.42x |
| `TipToTipDistance.calculate_tip_to_tip_distance` | balanced 32768-tip tree, opposite terminal tips | 0.1331s | 0.0189s | 7.1x |
| `TipToTipDistance.calculate_tip_to_tip_distance` child-list terminal check | balanced 65536-tip tree, opposite terminal tips, optimized helper baseline | 0.061879s | 0.052576s | 1.18x |
| `TipToTipDistance.calculate_tip_to_tip_distance` same-tip shortcut | pectinate 20000-tip tree, deepest terminal compared to itself after validation | 1.666261s | 0.856954s | 1.94x |
| `TipToTipDistance.calculate_tip_to_tip_distance` same-tip validation scan | pectinate 20000-tip tree, deepest terminal compared to itself, side-by-side previous same-tip shortcut | 1.326063s | 0.515092s | 2.57x |
| `TipToTipDistance.run` cached read-only tree setup | balanced 32768-tip cached tree, single-pair calculation and output mocked | 0.376594s | 0.000098s | 3823.29x |
| `tip_to_tip_distance` module import without eager SciPy clustering | cold process import for non-plot tip-distance command module | 0.427668s | 0.234465s | 1.8x |
| `tip_to_tip_distance` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy and TreeMixin imports | 0.135689s | 0.031125s | 4.36x |
| `tip_to_tip_distance` module import without eager JSON/plot config helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.011718s | 0.004924s | 2.38x |
| `tip_to_tip_distance` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005712s | 0.003511s | 1.63x |
| `TipToTipNodeDistance.calculate_tip_to_tip_node_distance` | balanced 32768-tip tree, opposite terminal tips | 0.1359s | 0.0130s | 10.5x |
| `TipToTipNodeDistance.calculate_tip_to_tip_node_distance` child-list terminal check | balanced 65536-tip tree, opposite terminal tips, side-by-side previous `is_terminal()` traversal | 0.025455s | 0.020033s | 1.27x |
| `TipToTipNodeDistance.calculate_tip_to_tip_node_distance` ancestor-distance path calculation | balanced 262144-tip tree, opposite terminal tips, side-by-side previous root-path-list helper | 0.275265s | 0.240866s | 1.14x |
| `TipToTipNodeDistance.calculate_tip_to_tip_node_distance` same-tip validation scan | pectinate 20000-tip tree, deepest terminal compared to itself, side-by-side previous ancestor-distance path | 0.355335s | 0.124900s | 2.84x |
| `TipToTipNodeDistance.run` cached read-only tree setup | balanced 32768-tip cached tree, node-distance calculation and output mocked | 0.355596s | 0.000087s | 4108.99x |
| `tip_to_tip_node_distance` module import without eager Bio.Phylo | cold subprocess import after lazy TreeMixin proxy and postponed annotations | 0.129296s | 0.026274s | 4.92x |
| `tip_to_tip_node_distance` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006248s | 0.004857s | 1.29x |
| `tip_to_tip_node_distance` module import without `typing` startup | median cold subprocess import after converting the annotation-only typing alias to a built-in postponed annotation | 0.008182s | 0.004153s | 1.97x |
| `PatristicDistances.calculate_distance_between_pairs` | balanced tree with 220 tips | 4.3503s | 0.0211s | 206.2x |
| `PatristicDistances.calculate_patristic_distance_stats` stats-only path | balanced tree with 768 tips | 0.095535s | 0.078120s | 1.22x |
| `PatristicDistances.calculate_pairwise_tip_distance_values_fast` child push | balanced 768-tip tree, stats-only pair distances, side-by-side previous `reversed(children)` setup | 0.505635s | 0.377701s | 1.34x |
| `PatristicDistances.calculate_pairwise_tip_distance_values_fast` deep-tree LCA index | pectinate 1200-tip tree, 719400 stats-only pair distances, copied old path baseline | 6.610513s | 0.820027s | 8.06x |
| `PatristicDistances.calculate_distance_values_between_pairs` stats-only fallback streaming pairs | 2000 nonstandard-tree tips, 1,999,000 pair values, side-by-side previous returned-combo path setup | 1.020654s | 0.610121s | 1.67x |
| `PatristicDistances.run` verbose text output | 200k pairwise distance rows, mocked tree/read and identical stdout text | 0.134971s | 0.111984s | 1.21x |
| `patristic_distances` module import without eager Bio.Phylo/tqdm | cold subprocess import of patristic-distances command module | 0.199123s | 0.121178s | 1.64x |
| `patristic_distances` module import without eager stats NumPy | cold subprocess import after lazy shared summary helper | 0.121178s | 0.073601s | 1.65x |
| `patristic_distances` module import without eager multiprocessing/pickle | cold subprocess import after lazy `mp` and `pickle` proxies plus localized `partial` import | 0.036541s | 0.030312s | 1.21x |
| `patristic_distances` module import without eager stats/json helpers | cold subprocess import after lazy forwarding wrappers for summary and JSON helpers | 0.028098s | 0.022213s | 1.26x |
| `patristic_distances` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005310s | 0.003553s | 1.49x |
| `PatristicDistances.run` cached read-only tree setup | balanced 32768-tip cached tree, non-verbose stats/output mocked | 0.388442s | 0.000097s | 3994.26x |
| `SpuriousSequence.get_branch_lengths_and_their_names` | balanced 65536-tip tree, terminal branch-length collection | 0.1404s | 0.0184s | 7.6x |
| `SpuriousSequence._iter_terminal_clades` order-preserving child push | balanced 131072-tip tree, terminal clade list with identical tip order, side-by-side previous `reversed(children)` helper | 0.045395s | 0.040772s | 1.11x |
| `SpuriousSequence.get_branch_lengths_and_their_names` one-pass collection | balanced 131072-tip standard tree, terminal branch lengths and name map, side-by-side previous terminal-list materialization | 0.095943s | 0.075439s | 1.27x |
| `SpuriousSequence.identify_spurious_sequence` large median | balanced 131072-tip tree, varied terminal branch lengths and identical threshold | 0.075865s | 0.057356s | 1.32x |
| `SpuriousSequence.run` batched text output | 50k flagged terminal rows, mocked tree/read and identical stdout text | 0.087150s | 0.031513s | 2.77x |
| `SpuriousSequence.run` cached read-only tree path | balanced 32768-tip cached tree with 1% long terminal branches, output mocked | 0.122502s | 0.005828s | 21.02x |
| `spurious_sequence` module import without eager Bio.Phylo | cold subprocess import of spurious-sequence command module | 0.157770s | 0.070032s | 2.25x |
| `spurious_sequence` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.010024s | 0.009378s | 1.07x |
| `spurious_sequence` module import without eager statistics | median cold subprocess import after local median helper | 0.010221s | 0.005496s | 1.86x |
| `spurious_sequence` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.007591s | 0.003433s | 2.21x |
| `DensityMap._terminal_clades` layout setup | balanced 65536-tip tree, terminal y-order mapping | 0.1362s | 0.0157s | 8.7x |
| `DensityMap._iter_preorder` order-preserving child push | balanced 131072-tip tree, preorder helper, side-by-side previous `reversed(children)` helper | 0.038833s | 0.026046s | 1.49x |
| `DensityMap._terminal_clades` order-preserving child push | balanced 131072-tip tree, terminal helper, side-by-side previous `reversed(children)` helper | 0.041157s | 0.037719s | 1.09x |
| `DensityMap._plot_density_map` phylogram layout setup | balanced 32768-tip tree, precomputed mappings and parent map | 0.8647s | 0.1049s | 8.2x |
| `DensityMap._plot_density_map` cladogram layout setup | balanced 32768-tip tree, precomputed mappings and parent map | 0.9129s | 0.1336s | 6.8x |
| `DensityMap._plot_density_map` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.051980s | 0.041415s | 1.26x |
| `DensityMap._plot_density_map` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.080175s | 0.066198s | 1.21x |
| `DensityMap._plot_density_map` rectangular batched posterior segments | balanced 256-tip tree, 50 posterior segments per branch, real Matplotlib Agg render | 6.629881s | 0.889630s | 7.45x |
| `DensityMap._plot_density_map` circular batched posterior segments | balanced 256-tip tree, 50 posterior segments per branch plus internal arcs, real Matplotlib Agg render | 6.455050s | 0.687876s | 9.38x |
| `DensityMap._compute_branch_posteriors` one-pass history scan | 5000 stochastic mappings, 80 posterior segments, 4 states, 8 transitions per branch history | 0.700506s | 0.221886s | 3.16x |
| `DensityMap.run` all-tip-state read-only setup | balanced 32768-tip cached tree, trait data for every tip, fitting/simulation/plot/output mocked | 0.258584s | 0.048538s | 5.33x |
| `DensityMap._print_text_output` batched summary | 100k captured density-map text summaries, identical stdout text | 0.104672s | 0.063514s | 1.65x |
| `density_map` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.084363s | 0.032070s | 2.63x |
| `density_map` module import without eager JSON/plot/color helpers | median cold subprocess import after localizing output, PlotConfig, circular-layout, and color-annotation helpers | 0.014782s | 0.005735s | 2.58x |
| `density_map` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.007946s | 0.005545s | 1.43x |
| `RobinsonFouldsDistance.run` rooting-tip setup | balanced 65536-tip tree, first terminal lookup before rooting | 0.1232s | 0.000005s | 27122.6x |
| `RobinsonFouldsDistance._first_terminal_name` direct leftmost descent | balanced depth-17 tree, first terminal lookup, side-by-side previous direct stack helper | 0.000002650s | 0.000001039s | 2.55x |
| `RobinsonFouldsDistance.calculate_robinson_foulds_distance` | balanced 256-tip tree pair, rooted descendant-set RF semantics | 0.2184s | 0.0022s | 101.6x |
| `RobinsonFouldsDistance.calculate_robinson_foulds_distance` direct terminal count | balanced 4096-tip identical tree pair, rooted descendant-set RF semantics | 0.068054s | 0.055135s | 1.23x |
| `RobinsonFouldsDistance._terminal_count_direct` localized stack operations | balanced 262144-tip tree, terminal count for RF normalization, side-by-side previous direct helper | 0.027284s | 0.026463s | 1.03x |
| `RobinsonFouldsDistance.get_all_bipartitions` direct postorder | balanced 8192-tip tree, rooted descendant-set RF split extraction | 0.045320s | 0.015135s | 2.99x |
| `RobinsonFouldsDistance._get_all_bipartitions_direct` binary child union | balanced 32768-tip tree, rooted descendant-set RF split extraction, side-by-side previous direct helper | 0.049324s | 0.038265s | 1.29x |
| `RobinsonFouldsDistance._get_all_bipartition_id_sets_direct` binary child union | balanced 32768-tip tree, compact rooted descendant split ids, side-by-side previous direct helper | 0.048783s | 0.039282s | 1.24x |
| `RobinsonFouldsDistance._get_all_bipartitions_direct` reverse-preorder postorder | balanced 32768-tip tree, rooted descendant-set RF split extraction, side-by-side binary-union helper | 0.039855s | 0.034580s | 1.15x |
| `RobinsonFouldsDistance._get_all_bipartition_id_sets_direct` reverse-preorder postorder | balanced 32768-tip tree, compact rooted descendant split ids, side-by-side binary-union helper | 0.037405s | 0.033784s | 1.11x |
| `RobinsonFouldsDistance.calculate_robinson_foulds_distance` same-object shortcut | balanced 32768-tip tree compared to itself, side-by-side previous compact split-id path and unchanged RF normalization denominator | 1.189548s | 0.004597s | 258.79x |
| `rf_distance` module import without annotation-only Bio.Phylo | cold subprocess import after postponed annotations and removing `Newick` import | 0.138859s | 0.043381s | 3.20x |
| `rf_distance` module import without eager concurrent futures | cold subprocess import after lazy `ProcessPoolExecutor` proxy | 0.043023s | 0.025542s | 1.68x |
| `rf_distance` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.007791s | 0.006287s | 1.24x |
| `rf_distance` module import without eager pickle | median cold subprocess import after lazy pickle proxy | 0.006557s | 0.004961s | 1.32x |
| `rf_distance` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005337s | 0.003524s | 1.51x |
| `KuhnerFelsensteinDistance.run` rooting-tip setup | balanced 65536-tip tree, first terminal lookup before rooting | 0.1238s | 0.000004s | 27506.1x |
| `KuhnerFelsensteinDistance.calculate_kf_distance` | balanced 1024-tip tree pair, branch score split map | 0.0472s | 0.0105s | 4.5x |
| `KuhnerFelsensteinDistance._get_splits_with_lengths` direct postorder | balanced 32768-tip tree, identical branch-score split map | 0.299108s | 0.200414s | 1.49x |
| `KuhnerFelsensteinDistance._get_splits_with_lengths` binary child union | balanced 32768-tip tree, branch-score split map, side-by-side previous direct split helper | 0.049567s | 0.035310s | 1.40x |
| `KuhnerFelsensteinDistance.calculate_kf_distance` same-object shortcut | balanced 8192-tip tree compared to itself, side-by-side previous double split-map path | 0.508294s | 0.000000s | >1e6x |
| `KuhnerFelsensteinDistance.calculate_kf_distance` split accumulation | balanced 16384-tip tree pair, same topology with perturbed branch lengths, side-by-side previous split-map union accumulation | 0.388520s | 0.204482s | 1.90x |
| `kf_distance` module import without eager Bio.Phylo | cold subprocess import of KF-distance command module | 0.164549s | 0.066013s | 2.49x |
| `kf_distance` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006883s | 0.005195s | 1.32x |
| `kf_distance` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005886s | 0.004585s | 1.28x |
| `TreeSpace._get_shared_taxa` terminal-name intersections | 80 balanced trees x 512 shared taxa | 0.0770s | 0.0095s | 8.1x |
| `TreeSpace`/`SpectralDiscordance._get_shared_taxa` no-slice gene-tree scan | 200k cached-tip gene-tree objects, identical 64 shared taxa, side-by-side previous `gene_trees[1:]` loop | 2.905650s | 2.567368s | 1.13x |
| `TreeSpace._build_distance_matrix` no-prune setup checks | 80 balanced trees x 512 shared taxa | 0.0771s | 0.0087s | 8.9x |
| `TreeSpace._build_distance_matrix` copied-tree batch pruning setup | balanced 8192-tip tree, prune 4096 copied tips before split extraction | 3.2511s | 0.0324s | 100.5x |
| `TreeSpace._prune_to_taxa` batch standard-tree pruning | balanced 4096-tip tree, prune to 2048 retained tips during shared-taxa normalization | 6.366311s | 0.102619s | 62.04x |
| `TreeSpace._build_distance_matrix` RF | 80 balanced trees x 128 shared taxa, rooted split matrix | 0.1168s | 0.0764s | 1.5x |
| `TreeSpace._build_distance_matrix` KF | 80 balanced trees x 128 shared taxa, branch-score matrix | 0.2151s | 0.0737s | 2.9x |
| `TreeSpace._extract_splits` direct postorder | 20 balanced 512-tip trees, RF splits filtered to 256 shared taxa | 0.056101s | 0.027321s | 2.1x |
| `TreeSpace._extract_splits` reverse-preorder postorder helper | 20 balanced 512-tip trees, RF splits filtered to 256 shared taxa | 0.030492s | 0.026352s | 1.16x |
| `TreeSpace._extract_splits` binary child union | balanced 1024-tip tree, RF splits filtered to 512 shared taxa, side-by-side previous reverse-preorder helper | 0.002574s | 0.002395s | 1.07x |
| `TreeSpace._extract_splits_with_lengths` binary child union | balanced 1024-tip tree, KF split lengths filtered to 512 shared taxa, side-by-side previous reverse-preorder helper | 0.002544s | 0.002405s | 1.06x |
| `TreeSpace._canonical_split` equal-size tiebreak | 3k equal-size 600-vs-600 bipartitions over 1200 taxa, identical sorted-lexicographic canonical side | 0.688258s | 0.257655s | 2.67x |
| `TreeSpace._build_distance_matrix` direct split extraction | 20 balanced 512-tip trees, prune to 256 shared taxa before RF splits | 0.101763s | 0.085802s | 1.2x |
| `TreeSpace._build_distance_matrix` RF shared-taxa direct extraction | 20 balanced 512-tip trees, 256 shared taxa before RF split matrix | 0.066800s | 0.021859s | 3.06x |
| `TreeSpace._auto_detect_k` normalized Laplacian scaling | 900-item synthetic affinity matrix, diagonal scaling step | 0.029895s | 0.002207s | 13.5x |
| `TreeSpace._auto_detect_k` normalized Laplacian row sums | 900 / 2500 / 5000-item synthetic affinity matrices, side-by-side previous `np.sum(..., axis=1)` degree-vector reduction | 0.000159105s / 0.003991750s / 0.015878400s | 0.000153917s / 0.001476250s / 0.008788600s | 1.03x / 2.70x / 1.81x |
| `TreeSpace._spectral_cluster` condensed distance setup | 2500 x 2500 symmetric distance matrix, side-by-side previous triangular-index upper-distance gather for bandwidth median | 0.340237s | 0.207055s | 1.64x |
| `TreeSpace._print_text_output` batched cluster rows | 100k cluster rows, captured stdout and identical text | 0.026882s | 0.016951s | 1.59x |
| `TreeSpace._write_distance_matrix` row-slice CSV formatting | 800 x 800 NumPy distance matrix, identical six-decimal CSV text | 0.199160s | 0.169827s | 1.17x |
| `TreeSpace._parse_trees_from_source` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.090451s | 0.067223s | 1.35x |
| `TreeSpace._parse_trees_from_source` path-list resolver | 50k existing relative tree paths, tree parsing mocked | 0.720531s | 0.507855s | 1.42x |
| `TreeSpace._parse_trees_from_source` streamed row parsing | 200k rows with comments/blanks, old cleaned-list parser vs streaming parser with tree parsing mocked | 0.237905s | 0.220375s | 1.08x |
| `TreeSpace._parse_trees_from_source` inline-Newick streamed row parsing | 200k inline Newick rows with comments/blanks, old cleaned-list parser vs streaming parser with tree parsing mocked | 0.316756s | 0.214420s | 1.48x |
| `TreeSpace._plot_heatmap` condensed distance vector setup | 4000 x 4000 symmetric distance matrix, identical upper-triangle condensed vector | 0.808634s | 0.005716s | 141.46x |
| `tree_space` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy proxy and lazy Phylo reader | 0.112687s | 0.031951s | 3.53x |
| `tree_space` module import without eager JSON/plot helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.015061s | 0.007018s | 2.15x |
| `tree_space` module import without eager pickle | median cold subprocess import after lazy pickle proxy | 0.006942s | 0.005695s | 1.22x |
| `tree_space` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.008027s | 0.003640s | 2.21x |
| `SpectralDiscordance._get_shared_taxa` terminal-name intersections | 80 balanced gene trees x 512 taxa, no species tree | 0.0787s | 0.0091s | 8.6x |
| `SpectralDiscordance._get_shared_taxa` with species tree | 80 balanced gene trees x 512 taxa plus species tree | 0.0808s | 0.0092s | 8.7x |
| `SpectralDiscordance._build_bipartition_matrix` NRF | 80 balanced trees x 128 shared taxa, bipartition presence matrix | 0.1008s | 0.0752s | 1.3x |
| `SpectralDiscordance._build_bipartition_matrix` WRF | 80 balanced trees x 128 shared taxa, branch-length weighted matrix | 0.1087s | 0.0806s | 1.3x |
| `SpectralDiscordance._extract_splits` direct postorder | balanced 2048-tip tree, identical canonical bipartition set | 0.055506s | 0.050864s | 1.09x |
| `SpectralDiscordance._extract_splits_with_lengths` direct postorder | balanced 2048-tip tree, identical canonical branch-length map | 0.056510s | 0.050592s | 1.12x |
| `SpectralDiscordance._canonical_split` equal-size tiebreak | 3k equal-size 600-vs-600 bipartitions over 1200 taxa, identical sorted-lexicographic canonical side | 0.487977s | 0.186022s | 2.62x |
| `SpectralDiscordance._canonical_split` size-first complement avoidance | 9k mixed 20-vs-1180, 1180-vs-20, and 600-vs-600 bipartitions over 1200 taxa | 0.269033s | 0.236251s | 1.14x |
| `SpectralDiscordance._spectral_cluster` normalized Laplacian scaling | 900-item synthetic affinity matrix, diagonal scaling step | 0.029895s | 0.002207s | 13.5x |
| `SpectralDiscordance._spectral_cluster` normalized Laplacian row sums | 900 / 2500 / 5000-item synthetic affinity matrices, side-by-side previous `np.sum(..., axis=1)` degree-vector reduction | 0.000159105s / 0.003991750s / 0.015878400s | 0.000153917s / 0.001476250s / 0.008788600s | 1.03x / 2.70x / 1.81x |
| `SpectralDiscordance._spectral_cluster` condensed distance reuse | 2500 gene-tree PCA score rows x 8 dimensions, side-by-side previous squareform plus triangular-index distance gather for bandwidth median | 0.303434s | 0.091461s | 3.32x |
| `SpectralDiscordance._spectral_cluster` eigenvector row normalization | 200x4 / 1000x8 / 5000x12 / 10000x20 eigenvector matrices, side-by-side previous `np.linalg.norm(..., axis=1, keepdims=True)` | 0.000003750s / 0.000011625s / 0.000133792s / 0.000351958s | 0.000002667s / 0.000005625s / 0.000054666s / 0.000056291s | 1.41x / 2.07x / 2.45x / 6.25x |
| `SpectralDiscordance._kmeans` vectorized distance and center updates | 80k rows x 8 dimensions, 12 clusters, fixed RandomState seed and identical labels | 12.232717s | 0.789952s | 15.49x |
| `SpectralDiscordance._kmeans` k-means++ initialization distances | 80k rows x 8 dimensions, 12 seeded initial centers, identical centers and closest-distance vector | 0.047207s | 0.033778s | 1.40x |
| `SpectralDiscordance._get_top_loadings` partial top-N selection | 20 PCs x 300k bipartitions, top 10 loadings per PC, identical reported entries | 0.442687s | 0.050885s | 8.70x |
| `SpectralDiscordance._run_pca` singular-value total variance | 2 / 3 / 4 / 8 / 16 / 32 / 128 / 1024 singular values, side-by-side previous `np.sum(S ** 2)` | 0.000005887s / 0.000006054s / 0.000005545s / 0.000005994s / 0.000006021s / 0.000006877s / 0.000004525s / 0.000006631s | 0.000001281s / 0.000001132s / 0.000001423s / 0.000001156s / 0.000001312s / 0.000001259s / 0.000000879s / 0.000001432s | 4.60x / 5.35x / 3.90x / 5.19x / 4.59x / 5.46x / 5.15x / 4.63x |
| `SpectralDiscordance._parse_gene_trees` path-list existence guard | 50k existing absolute tree paths, tree parsing mocked | 0.263681s | 0.121357s | 2.17x |
| `spectral_discordance` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy proxy, postponed annotations, and localized gene-tree parser import | 0.142026s | 0.031427s | 4.52x |
| `spectral_discordance` module import without eager JSON/plot helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.015141s | 0.007166s | 2.11x |
| `spectral_discordance` module import without eager pickle | median cold subprocess import after lazy pickle proxy | 0.006739s | 0.004936s | 1.37x |
| `spectral_discordance` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005708s | 0.003931s | 1.45x |
| `SpectralDiscordance._print_text` PC-score rows | 100k gene-tree score rows x 5 PCs, identical stdout text | 0.240505s | 0.197126s | 1.22x |
| `SpectralDiscordance._format_json` row-slice score payload | 100k gene-tree score rows x 5 PCs, identical nested payload | 0.157021s | 0.103268s | 1.52x |
| `SpectralDiscordance.run` cached read-only species tree setup | balanced 32768-tip cached species tree, gene-tree parsing, bipartition matrix, PCA, clustering, and output mocked | 0.100765s | 0.000038s | 2651.7x |
| `QuartetNetwork._normalize_taxa` identical taxa setup | 80 balanced trees x 512 shared taxa | 0.0784s | 0.0095s | 8.3x |
| consensus/quartet `_prune_to_taxa` direct terminal objects | balanced 8192-tip tree, prune to 4096 retained tips | 0.029969s | 0.013771s | 2.18x |
| `QuartetNetwork._compute_quartet_cfs` | 50 balanced trees x 24 taxa, all quartet topology counts | 1.0438s | 0.2805s | 3.7x |
| `QuartetNetwork._compute_quartet_cfs` precomputed quartet masks | 50 balanced trees x 24 taxa, all quartet topology counts | 0.363644s | 0.286839s | 1.27x |
| `QuartetNetwork._classify_quartet` chi-square p-values | cold process, one hybrid-like quartet count vector | 0.494692s | 0.000012084s | 40937.8x |
| `QuartetNetwork._compute_nanuq_distance` topology index helpers | 42 taxa, 111930 quartet classifications, mixed tree/hybrid/unresolved results | 0.111250s | 0.096283s | 1.16x |
| `QuartetNetwork._extract_bipartitions` direct ordered traversal | balanced 4096-tip tree, public frozenset bipartition helper | 0.389116s | 0.366308s | 1.06x |
| `QuartetNetwork._extract_bipartition_masks` | balanced 32768-tip tree, bitmask split extraction setup | 0.2149s | 0.1446s | 1.5x |
| `QuartetNetwork._extract_bipartition_masks` direct postorder | balanced 32768-tip tree, bitmask split extraction setup | 0.175142s | 0.081568s | 2.1x |
| `QuartetNetwork._extract_bipartition_masks` reverse-preorder postorder helper | balanced 32768-tip tree, bitmask split extraction setup | 0.104120s | 0.067357s | 1.55x |
| `QuartetNetwork._extract_bipartition_masks` binary child masks | balanced 32768-tip tree, bitmask split extraction setup, side-by-side reverse-preorder helper | 0.053298s | 0.049024s | 1.09x |
| `QuartetNetwork._compute_circular_split_weights` equal-size tiebreak | 80 / 160 / 220-taxon circular orderings, identical positive split weights, side-by-side previous sorted-half tiebreak | 0.075409s / 0.734834s / 1.975865s | 0.042086s / 0.551053s / 1.674072s | 1.79x / 1.33x / 1.18x |
| `QuartetNetwork._compute_circular_split_weights` rolling arc sets | 40 / 80 / 160 / 220-taxon circular orderings, identical split order and weights, side-by-side previous per-arc set reconstruction | 0.030619s / 0.182034s / 0.589527s / 1.769014s | 0.005200s / 0.020853s / 0.555364s / 0.464726s | 5.89x / 8.73x / 1.06x / 3.81x |
| `QuartetNetwork._compute_split_directions` cached one-pass split centers | 500 circular splits over 2000 taxa with cached gap positions, side-by-side previous per-split trigonometric center sums | 0.098739s | 0.041051s | 2.41x |
| `QuartetNetwork._build_splits_graph` | 80 taxa, 20 circular splits, 33-node quartet graph | 6.1411s | 0.0016s | 3910.4x |
| `QuartetNetwork._draw_quartet_network` edge rendering | 80 taxa, 20 circular splits, real Matplotlib Agg internal and pendant edge render | 0.027545s | 0.004864s | 5.66x |
| `QuartetNetwork.run` batched text quartet output | 100k quartet rows, captured stdout and identical text | 0.180601s | 0.166519s | 1.08x |
| `QuartetNetwork._parse_trees_from_source` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.091090s | 0.067767s | 1.34x |
| `QuartetNetwork._parse_trees_from_source` path-list resolver | 50k existing relative tree paths, tree parsing mocked | 0.720531s | 0.507855s | 1.42x |
| `quartet_network` module import without eager Bio.Phylo | cold subprocess import after lazy `Phylo.read` proxy | 0.134003s | 0.032005s | 4.19x |
| `quartet_network` module import without eager JSON/plot helpers | median cold subprocess import after localizing `PlotConfig` and lazy JSON wrapper | 0.017726s | 0.006822s | 2.60x |
| `quartet_network` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in annotations | 0.006016s | 0.004340s | 1.39x |
| `quartet_utils.compute_gcf_per_node` | 80 balanced gene trees x 128 taxa plus species tree, gCF/gDF counts | 0.2140s | 0.1058s | 2.0x |
| `quartet_utils.compute_gcf_per_node` direct traversal | 80 balanced gene trees x 128 taxa plus species tree, cached gCF/gDF counts | 0.1542s | 0.0944s | 1.6x |
| `quartet_utils.compute_gcf_per_node` split-frequency lookups | 80 balanced gene trees x 128 taxa plus species tree, cached gCF/gDF counts | 0.103750s | 0.057895s | 1.8x |
| `quartet_utils.canonical_split` tuple container | 50k repeated canonical split constructions over 64 taxa, side-by-side previous list container | 0.074187s | 0.067238s | 1.10x |
| `quartet_utils._collect_clade_tip_sets_direct` reverse-preorder helper | balanced 32768-tip tree, descendant tip sets for all taxa | 0.049942s | 0.044048s | 1.13x |
| `quartet_utils._collect_clade_tip_sets_direct` reverse-preorder shared-taxa helper | balanced 32768-tip tree, descendant tip sets filtered to half the taxa | 0.038905s | 0.035033s | 1.11x |
| `quartet_utils._collect_clade_tip_sets_direct` binary child merge | balanced 65536-tip tree, descendant tip sets for all taxa, side-by-side previous set-update helper | 0.179555s | 0.118090s | 1.52x |
| `quartet_utils._collect_clade_tip_sets_direct` filtered binary child merge | balanced 65536-tip tree, descendant tip sets filtered to half the taxa, side-by-side previous set-update helper | 0.136503s | 0.104268s | 1.31x |
| `quartet_utils._preorder_clades_direct` order-preserving child push | balanced 131072-tip tree, direct preorder helper, optimized helper baseline | 0.033150s | 0.028746s | 1.15x |
| `quartet_utils.parse_astral_annotations` direct preorder | balanced 65536-tip tree, q1/q2/q3 labels on every internal node | 0.264422s | 0.097986s | 2.70x |
| `quartet_utils.parse_astral_branch_info` direct preorder | balanced 65536-tip tree, f1/pp1 labels on every internal node | 0.274189s | 0.107192s | 2.56x |
| `quartet_utils` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in annotations | 0.022979s | 0.020954s | 1.10x |
| `ConcordanceAsr._compute_gcf_per_node` | 80 balanced gene trees x 128 taxa plus species tree, weighted-ASR gCF/gDF proportions | 0.1748s | 0.0843s | 2.1x |
| `ConcordanceAsr._compute_gcf_per_node` direct preorder iterators | 80 balanced gene trees x 512 taxa plus species tree, weighted-ASR gCF/gDF proportions, side-by-side previous `get_nonterminals`/`find_clades` iterators | 0.500742s | 0.417360s | 1.20x |
| `ConcordanceAsr._canonical_split` equal-size tiebreak | 3k equal-size 600-vs-600 bipartitions over 1200 taxa, identical sorted-lexicographic canonical side | 0.822011s | 0.289217s | 2.84x |
| `ConcordanceAsr._canonical_split` size-first complement avoidance | 9k mixed 20-vs-1180, 1180-vs-20, and 600-vs-600 bipartitions over 1200 taxa | 0.262131s | 0.237326s | 1.10x |
| `ConcordanceAsr._normalize_taxa` setup | 80 balanced gene trees x 512 taxa plus species tree | 0.0819s | 0.0101s | 8.1x |
| `ConcordanceAsr._parse_gene_trees` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.094455s | 0.071398s | 1.32x |
| `ConcordanceAsr._parse_gene_trees` path-list resolver | 50k relative tree path rows, tree parsing mocked | 0.081335s | 0.014871s | 5.47x |
| `ConcordanceAsr._parse_gene_trees` streaming rows | 500k mixed relative path and inline-Newick rows with comments/blanks, tree parsing mocked, side-by-side previous cleaned-list parser | 1.160041s | 0.920783s | 1.26x |
| `ConcordanceAsr._run_asr_on_tree` no-prune setup | balanced 32768-tip tree, all tips shared with trait data | 0.1338s | 0.0078s | 17.2x |
| `ConcordanceAsr._run_asr_on_tree` no-prune copy skip | balanced 32768-tip tree, all tips shared with trait data, ASR internals mocked | 0.640786s | 0.385320s | 1.66x |
| `ConcordanceAsr._run_asr_on_tree` prune-needed setup | balanced 32768-tip tree, half tips shared with trait data | 0.1367s | 0.0145s | 9.4x |
| `ConcordanceAsr.run` species-copy prune setup | balanced 32768-tip species tree, all tips shared with trait data | 0.0693s | 0.0076s | 9.2x |
| `ConcordanceAsr.run` all-shared species-copy skip | balanced 32768-tip cached species tree, all species tips have trait values, gene/ASR/output mocked | 0.652220s | 0.238411s | 2.74x |
| `ConcordanceAsr._run_distribution` result assembly | balanced 2048-tip species tree, ASR/gCF stubs, 20 gene-tree result dictionaries | 0.1196s | 0.0494s | 2.4x |
| `ConcordanceAsr._run_distribution` scalar estimate summaries | 2047 species-tree nodes x 20 gene-tree estimates, identical mean/population variance and sigma2 average | 0.078126s | 0.003796s | 20.58x |
| `ConcordanceAsr._law_of_total_variance` scalar small-list path | three weighted topology estimates, identical total/within/between variances | 0.000039046s | 0.000004145s | 9.42x |
| `ConcordanceAsr._collect_uncertainty_node_data` | balanced 2048-tip species tree, uncertainty entries for every internal node | 0.1004s | 0.0149s | 6.7x |
| `ConcordanceAsr._plot_concordance_contmap` rectangular setup | balanced 32768-tip tree, precomputed ancestral estimates and gCF values | 1.3490s | 0.1657s | 8.1x |
| `ConcordanceAsr._plot_concordance_contmap` circular setup | balanced 32768-tip tree, precomputed ancestral estimates and gCF values | 0.9902s | 0.1985s | 5.0x |
| `ConcordanceAsr._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.034202s | 0.023919s | 1.43x |
| `ConcordanceAsr._collect_clade_tip_sets` reverse-preorder binary merge | balanced 32768-tip tree, descendant tip-set cache used by gCF, ASR re-keying, and uncertainty output | 0.343929s | 0.172670s | 1.99x |
| `ConcordanceAsr._plot_concordance_contmap` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.051516s | 0.039994s | 1.29x |
| `ConcordanceAsr._plot_concordance_contmap` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.055215s | 0.044103s | 1.25x |
| `ConcordanceAsr._plot_concordance_contmap` rectangular batched base branches | balanced 512-tip tree, gCF markers present, real Matplotlib Agg render | 2.872509s | 1.955351s | 1.47x |
| `ConcordanceAsr._plot_concordance_contmap` batched gCF markers | 4096 gCF markers with variable sizes/colors, real Matplotlib Agg scatter render | 9.035145s | 0.063602s | 142.06x |
| `ConcordanceAsr._plot_uncertainty` mean markers | 2047 uncertainty rows x 20 estimates, identical marker x-coordinates | 0.035963s | 0.000745s | 48.26x |
| `ConcordanceAsr._print_text_output` batched estimate output | 100k ancestral-estimate rows, mixed CI/non-CI rows, captured stdout and identical text | 0.167087s | 0.151078s | 1.11x |
| `concordance_asr` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy proxy and lazy Phylo reader | 0.146693s | 0.033499s | 4.38x |
| `concordance_asr` module import without eager ASR/plot helpers | cold subprocess import after localizing ASR helper, pickle, and plotting-helper imports | 0.033717s | 0.026711s | 1.26x |
| `concordance_asr` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006562s | 0.005254s | 1.25x |
| `concordance_asr` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006545s | 0.004335s | 1.51x |
| `EvoTempoMap._classify_gene_trees` | 40 balanced 256-tip gene trees plus species tree, branch concordance classification | 0.2659s | 0.1611s | 1.7x |
| `EvoTempoMap._parse_and_validate_gene_trees` branch-length validation | 80 balanced 2048-tip gene trees, identical validated branch count | 0.494464s | 0.048419s | 10.21x |
| `EvoTempoMap._parse_gene_trees` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.097501s | 0.072305s | 1.35x |
| `EvoTempoMap._parse_gene_trees` path-list resolver | 50k relative tree path rows, tree parsing mocked | 0.081335s | 0.014871s | 5.47x |
| `EvoTempoMap._classify_gene_trees` shared-taxa/no-prune setup | 80 balanced gene trees x 512 taxa plus species tree | 0.1580s | 0.0088s | 17.9x |
| `EvoTempoMap._classify_gene_trees` single-pass gene splits | 40 balanced 256-tip gene trees plus species tree, branch concordance classification | 0.0817s | 0.0621s | 1.3x |
| `EvoTempoMap._canonical_split` equal-size tiebreak | 3k equal-size 600-vs-600 bipartitions over 1200 taxa, identical sorted-lexicographic canonical side | 0.871998s | 0.283520s | 3.08x |
| `EvoTempoMap._compute_global_treeness` | 40 balanced 512-tip gene trees plus species tree, global split concordance and treeness | 0.5462s | 0.2757s | 2.0x |
| `EvoTempoMap._compute_global_treeness` direct split traversal | 40 balanced 512-tip gene trees plus species tree, global split concordance and treeness | 0.217001s | 0.127218s | 1.71x |
| `EvoTempoMap._compute_global_treeness` treeness summaries | 10k mean/median summaries of 40 treeness values, identical summary values | 0.365470s | 0.016978s | 21.53x |
| `EvoTempoMap._collect_clade_taxa` reverse-preorder postorder helper | balanced 8192-tip tree, descendant taxon sets for branch classification | 0.028838s | 0.019994s | 1.44x |
| `EvoTempoMap._build_parent_map` direct map traversal | balanced 65536-tip tree, branch-classification parent map, optimized helper baseline | 0.026751s | 0.017737s | 1.51x |
| `EvoTempoMap._iter_preorder_clades` binary-child fast path | balanced 131072-tip tree, branch-classification preorder helper, side-by-side previous `reversed(children)` helper | 0.039015s | 0.034670s | 1.13x |
| `EvoTempoMap._compute_treeness` batch | 40 balanced 4096-tip gene trees, helper-only treeness values | 0.9191s | 0.0330s | 27.8x |
| `EvoTempoMap._test_branch` insufficient-data summaries | 10k singleton concordant vs singleton discordant length summaries, identical early-return stats | 0.676468s | 0.022739s | 29.75x |
| `EvoTempoMap._fdr` | 1M synthetic p-values | 0.647786s | 0.122101s | 5.3x |
| `EvoTempoMap._fdr` small-list path without NumPy startup | cold subprocess, 7 p-values through Benjamini-Hochberg helper | 0.071440s | 0.023255s | 3.07x |
| `evo_tempo_map` module import without eager `scipy.stats` | cold process import for evo-tempo-map command module | 0.690765s | 0.209734s | 3.3x |
| `evo_tempo_map` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy proxy and lazy Phylo reader | 0.129755s | 0.031645s | 4.10x |
| `evo_tempo_map` module import without eager `pathlib` | median cold subprocess import after lazy gene-tree-list `Path` wrapper | 0.040799s | 0.039905s | 1.02x |
| `evo_tempo_map` module import without eager JSON/plot config helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.011676s | 0.004933s | 2.37x |
| `evo_tempo_map` module import without `typing` startup | median cold subprocess import after replacing annotation-only typing aliases with built-in postponed annotations | 0.035356s | 0.034328s | 1.03x |
| `EvoTempoMap.run` cached read-only species tree setup | balanced 32768-tip cached species tree, gene-tree parsing, branch classification, global treeness, and text output mocked | 0.128927s | 0.000028s | 4604.5x |
| `EvoTempoMap._output_text` batched verbose branch output | 100k branch rows plus verbose concordant/discordant length details, captured stdout and identical text | 0.294768s | 0.247926s | 1.19x |
| `EvoTempoMap._output_text` single-pass verbose formatting | 100k branch rows plus verbose concordant/discordant length details, captured stdout and identical text, side-by-side previous second verbose pass comparison | 0.217947s | 0.198972s | 1.10x |
| `EvoTempoMap._plot` strip-point rendering | 1024 branch groups, 4 concordant and 4 discordant strip points per branch, real Matplotlib Agg scatter render | 3.675574s | 0.042251s | 86.99x |
| `EvoTempoMap._plot` significant star rendering | 4096 significant branch markers, real Matplotlib Agg star render | 0.656416s | 0.018294s | 35.88x |
| `Hybridization._count_topologies` | 40 balanced 256-tip gene trees plus species tree, NNI topology counts | 0.0838s | 0.0782s | 1.1x |
| `Hybridization`/`DiscordanceAsymmetry._canonical_split` equal-size tiebreak | 3k equal-size 600-vs-600 bipartitions over 1200 taxa, identical sorted-lexicographic canonical side | 0.601252s | 0.234306s | 2.57x |
| `Hybridization._canonical_split` size-first complement avoidance | 9k mixed 20-vs-1180, 1180-vs-20, and 600-vs-600 bipartitions over 1200 taxa | 0.260723s | 0.230776s | 1.13x |
| `DiscordanceAsymmetry._canonical_split` size-first complement avoidance | 9k mixed 20-vs-1180, 1180-vs-20, and 600-vs-600 bipartitions over 1200 taxa | 0.269212s | 0.215084s | 1.25x |
| `Hybridization`/`DiscordanceAsymmetry._count_split_matches` single-pass split counts | 200k gene-tree split sets x 3 target topologies, identical concordant/alt1/alt2 counts | 0.089838s | 0.041921s | 2.14x |
| `Hybridization._count_topologies` species-taxon setup | balanced 32768-tip species tree | 0.0673s | 0.0077s | 8.8x |
| `Hybridization._parse_gene_trees` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.093972s | 0.070400s | 1.33x |
| `Hybridization._parse_gene_trees` path-list resolver | 50k relative tree path rows, tree parsing mocked | 0.081335s | 0.014871s | 5.47x |
| `Hybridization._fdr` | 1M synthetic p-values | 0.738633s | 0.164136s | 4.5x |
| `Hybridization._fdr` small-list path without NumPy startup | cold subprocess, 7 p-values through Benjamini-Hochberg helper | 0.077783s | 0.023435s | 3.32x |
| `hybridization` module import without `scipy.stats` | cold process import for reticulation command module | 0.629033s | 0.183021s | 3.4x |
| `hybridization` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy proxy and localized gene-tree parser import | 0.147853s | 0.032011s | 4.62x |
| `hybridization` module import without eager JSON/plot/circular/color helpers | median cold subprocess import after localizing PlotConfig, plotting helpers, circular/color helpers, and lazy JSON wrapper | 0.012903s | 0.005302s | 2.43x |
| `hybridization` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.034640s | 0.030366s | 1.14x |
| `Hybridization.run` cached read-only species tree setup | balanced 32768-tip cached species tree, gene-tree parsing, topology counting, and text output mocked | 0.097596s | 0.000028s | 3485.6x |
| `Hybridization._output_text` batched significant-branch output | three captured text reports with 100k significant branch rows each, identical stdout text | 0.480525s | 0.435568s | 1.10x |
| `Hybridization._output_text` single-pass row formatting | 100k significant branch rows, captured stdout and identical text, side-by-side previous list-filter/f-string row formatter comparison | 0.147229s | 0.134598s | 1.09x |
| `Hybridization`/`DiscordanceAsymmetry._build_parent_map` unordered child push | balanced 65536-tip tree, mirrored object parent-map helper, optimized helper baseline | 0.024100s | 0.018367s | 1.31x |
| `Hybridization`/`DiscordanceAsymmetry._get_terminal_clades` order-preserving binary push | balanced 131072-tip tree, mirrored terminal-clade helper with identical tip order, optimized helper baseline | 0.021994s | 0.017857s | 1.23x |
| `Hybridization._collect_clade_taxa` reverse-preorder binary merge | balanced 65536-tip tree, mirrored clade-taxon map helper, side-by-side previous visited-stack set-update helper | 0.150542s | 0.098837s | 1.52x |
| `DiscordanceAsymmetry._collect_clade_taxa` reverse-preorder binary merge | balanced 65536-tip tree, mirrored clade-taxon map helper, side-by-side previous visited-stack set-update helper | 0.150542s | 0.083811s | 1.80x |
| `Hybridization._plot` traversal reuse | balanced 32768-tip tree, rectangular coordinate/result/branch/star traversal setup | 0.712564s | 0.160012s | 4.45x |
| `Hybridization._plot` rectangular batched branch rendering | balanced 2048-tip tree, per-internal-node hybrid scores, real Matplotlib Agg branch render | 4.688986s | 2.918381s | 1.61x |
| `Hybridization._plot` circular batched branch rendering | balanced 2048-tip tree, per-internal-node hybrid scores, real Matplotlib Agg branch/arc render | 1.927446s | 0.729285s | 2.64x |
| `Hybridization._plot` significant star markers | 4096 significant branch markers, real Matplotlib Agg scatter render | 6.537405s | 0.013889s | 470.69x |
| `Hybridization._plot` repeated hybrid-score color cache | balanced 2048-tip tree, rectangular Agg plot with labels/title/legend disabled and repeated per-internal-node hybrid score | 0.485382s | 0.425973s | 1.14x |
| `DiscordanceAsymmetry._count_topologies` | 40 balanced 256-tip gene trees plus species tree, NNI topology counts | 0.0830s | 0.0789s | 1.1x |
| `DiscordanceAsymmetry._count_topologies` species-taxon setup | balanced 32768-tip species tree | 0.0673s | 0.0070s | 9.6x |
| `DiscordanceAsymmetry._parse_gene_trees` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.093972s | 0.070400s | 1.33x |
| `DiscordanceAsymmetry._parse_gene_trees` path-list resolver | 50k relative tree path rows, tree parsing mocked | 0.081335s | 0.014871s | 5.47x |
| `DiscordanceAsymmetry._fdr` | 1M synthetic p-values | 0.712219s | 0.164136s | 4.3x |
| `DiscordanceAsymmetry._fdr` small-list path without NumPy startup | cold subprocess, 7 p-values through Benjamini-Hochberg helper | 0.067908s | 0.023407s | 2.90x |
| `discordance_asymmetry` module import without `scipy.stats` | cold process import for NNI-asymmetry command module | 0.638335s | 0.187172s | 3.4x |
| `discordance_asymmetry` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy proxy and localized gene-tree parser import | 0.153823s | 0.032043s | 4.80x |
| `discordance_asymmetry` module import without eager JSON/plot/circular/color helpers | median cold subprocess import after localizing PlotConfig, plotting helpers, circular/color helpers, and lazy JSON wrapper | 0.012810s | 0.005459s | 2.35x |
| `discordance_asymmetry` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.007310s | 0.004311s | 1.70x |
| `DiscordanceAsymmetry.run` cached read-only species tree setup | balanced 32768-tip cached species tree, gene-tree parsing, topology counting, and text output mocked | 0.094878s | 0.000032s | 2964.9x |
| `DiscordanceAsymmetry._output_text` batched verbose branch output | 100k branch rows plus verbose per-branch details, captured stdout and identical text | 0.230438s | 0.197002s | 1.17x |
| `DiscordanceAsymmetry._output_text` single-pass verbose formatting | 100k branch rows plus verbose per-branch details, captured stdout and identical text, side-by-side previous second verbose pass comparison | 0.194901s | 0.170617s | 1.14x |
| `DiscordanceAsymmetry._plot` traversal reuse | balanced 32768-tip tree, rectangular coordinate/result/branch/star traversal setup | 0.647252s | 0.153974s | 4.20x |
| `DiscordanceAsymmetry._plot` rectangular batched branch rendering | balanced 2048-tip tree, per-internal-node asymmetry ratios, real Matplotlib Agg branch render | 4.725613s | 2.922123s | 1.62x |
| `DiscordanceAsymmetry._plot` circular batched branch rendering | balanced 2048-tip tree, per-internal-node asymmetry ratios, real Matplotlib Agg branch/arc render | 1.912087s | 0.724776s | 2.64x |
| `DiscordanceAsymmetry._plot` significant star markers | 4096 significant branch markers, real Matplotlib Agg scatter render | 6.537405s | 0.013889s | 470.69x |
| `DiscordanceAsymmetry._plot` repeated asymmetry-ratio color cache | balanced 2048-tip tree, rectangular Agg plot with labels/title/legend disabled and repeated per-internal-node asymmetry ratio | 0.482074s | 0.427223s | 1.13x |
| `Hybridization`/`DiscordanceAsymmetry` plot setup | balanced 4096-tip species tree, rectangular branch-result lookup setup | 0.1725s | 0.0696s | 2.5x |
| `Hybridization`/`DiscordanceAsymmetry` child-y coordinate mean | balanced 32768-tip species tree, postorder node-y coordinate setup with identical positions | 0.466845s | 0.046268s | 10.09x |
| `PolytomyTest._evaluate_tree_triplets_fast` | balanced 1024-tip tree, 8k group triplets | 5.0266s | 0.0481s | 104.4x |
| `PolytomyTest.examine_all_triplets_and_sister_pairing` large-dispatch triplet setup | 90x90x90 group product with tree read/prep/fast evaluator mocked | 0.137157s | 0.000028416s | 4826.7x |
| `PolytomyTest._evaluate_tree_triplets_fast` triplet/group-set reuse | balanced 1024-tip tree, 8k group triplets, side-by-side previous per-triplet group conversion | 0.037052s | 0.025503s | 1.45x |
| `PolytomyTest._common_ancestor_from_path_cache` triplet-node check | 2M zipped triplet path node rows, same / second-different / third-different cases, side-by-side previous `nodes[1:]` generator predicate | 2.809178s / 1.873774s / 1.963923s | 0.185879s / 0.076666s / 0.159605s | 15.11x / 24.44x / 12.30x |
| `PolytomyTest.chisquare_tests` p-value helper | cold process, two three-category equal-frequency chi-square tests | 0.487223s | 0.000007208s | 67594.8x |
| `PolytomyTest.determine_sisters_and_add_to_counter` polytomy check | 5000 resolved three-tip counter evaluations | 0.250329s | 0.201572s | 1.24x |
| `PolytomyTest._has_exactly_three_terminals` triplet check | 5000 resolved three-tip trees before legacy sister counting | 0.027710s | 0.002097s | 13.21x |
| `PolytomyTest._has_exactly_three_terminals` fallback terminal scan | nonstandard tree yielding 500k terminals, false case | 0.043168s | 0.000001167s | 36990.6x |
| `PolytomyTest.check_if_triplet_is_a_polytomy` direct internal count | 5000 resolved plus 5000 polytomy three-tip trees | 0.043043s | 0.004448s | 9.68x |
| `PolytomyTest.check_if_triplet_is_a_polytomy` fallback nonterminal scan | nonstandard tree yielding 500k nonterminals, false case | 0.059486s | 0.000000958s | 62093.9x |
| `PolytomyTest.count_number_of_groups_in_triplet` precomputed group-set reuse | 5k triplet group-representation checks over three 1k-taxon groups | 1.415379s | 0.617556s | 2.29x |
| `PolytomyTest.sister_relationship_counter` direct nested increment | 3M repeated sister-count increments over one tree summary | 1.674614s | 0.906683s | 1.85x |
| `PolytomyTest.set_branch_lengths_in_tree_to_one` direct traversal | balanced 65536-tip tree, reset every branch length to one | 0.202605s | 0.015667s | 12.93x |
| `PolytomyTest.print_gene_support_freq_res` batched text output | 100k captured gene-support frequency summaries, identical stdout text | 0.234641s | 0.155517s | 1.51x |
| `PolytomyTest.get_triplet_tree` prune-list setup | 200 mocked legacy triplet-tree preparations over 8195 tips, tree read/root/prune mocked | 0.194432s | 0.121392s | 1.60x |
| `PolytomyTest` triplet group first-identifier lookup | 500k group identifiers, per-lookup median from batched repetitions | 0.033081725s | 0.000000108s | 306771.8x |
| `polytomy_test` module import without eager Bio.Phylo | cold subprocess import with lazy `Phylo.read` proxy | 0.190013s | 0.095192s | 2.00x |
| `polytomy_test` module import without eager unittest.mock | cold subprocess import after lightweight mock-object detection | 0.056568s | 0.035368s | 1.60x |
| `polytomy_test` module import without eager concurrency/copy helpers | cold subprocess import after lazy multiprocessing/thread-pool proxies and deferred tree-copy imports | 0.035875s | 0.025992s | 1.38x |
| `polytomy_test` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006955s | 0.005686s | 1.22x |
| `polytomy_test` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.007281s | 0.004536s | 1.60x |
| `polytomy_test` module import without eager file helper | median cold subprocess import after lazy tree-list reader wrapper | 0.031790s | 0.021889s | 1.45x |
| `TransferAnnotations.run` taxa validation setup | two balanced 32768-tip trees | 0.1371s | 0.0148s | 9.3x |
| `TransferAnnotations.run` cached read-only source setup | two balanced 32768-tip trees, taxa validation, annotation extraction, transfer, write, and text output mocked | 0.227289s | 0.148054s | 1.54x |
| `TransferAnnotations._extract_annotations` + `_transfer` | 1024-tip balanced source/target trees, annotations on all internal source nodes | 0.0585s | 0.0386s | 1.5x |
| `TransferAnnotations._extract_annotations` + `_transfer` single postorder pass | 1024-tip balanced source/target trees, annotations on all internal source nodes | 0.0260s | 0.0207s | 1.3x |
| `TransferAnnotations._extract_annotations` + `_transfer` direct postorder and lazy complement | 4096-tip balanced source/target trees, annotations on all internal source nodes | 0.295462s | 0.013271s | 22.26x |
| `TransferAnnotations._collect_clade_taxa` unary child fast path | 65536-tip mixed unary/binary tree, identical clade-to-taxa map | 0.090294s | 0.073181s | 1.23x |
| `TransferAnnotations._print_text_output` batched summary | 100k captured transfer-annotation text summaries, identical stdout text | 0.065141s | 0.044902s | 1.45x |
| `TransferAnnotations._write_annotated_tree` single output write | 8192-tip balanced annotated tree, identical normalized Newick output | 0.099131s | 0.031941s | 3.10x |
| `transfer_annotations` module import without eager Bio.Phylo | cold subprocess import after lazy `Phylo.read`/`Phylo.write` proxy | 0.113289s | 0.025829s | 4.39x |
| `transfer_annotations` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006235s | 0.004812s | 1.30x |
| `transfer_annotations` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only typing aliases to built-in annotations | 0.011073s | 0.007502s | 1.48x |
| `color_annotations.get_clade_tip_ids` | balanced 65536-tip clade, collect terminal object ids for highlighted ranges/clades | 0.1280s | 0.0200s | 6.4x |
| `color_annotations.get_clade_branch_ids` | balanced 32768-tip clade, collect highlighted clade descendant node ids | 0.098573s | 0.013878s | 7.10x |
| `color_annotations.get_clade_branch_ids` unordered child push | balanced 32768-tip clade, descendant id set for highlighted clades, optimized helper baseline | 0.009671s | 0.007221s | 1.34x |
| `color_annotations.resolve_mrca` early-stop taxa validation | balanced 32768-tip tree, MRCA lookup for first two tips | 0.014703s | 0.0000155s | 946.92x |
| `color_annotations._terminal_clades` order-preserving child push | balanced 131072-tip clade, terminal list for highlighted ranges/clades, optimized helper baseline | 0.038837s | 0.034347s | 1.13x |
| `color_annotations._valid_mrca_taxa` order-preserving child push | balanced 131072-tip tree, validate first and last requested taxa, optimized helper baseline | 0.040392s | 0.037106s | 1.09x |
| `color_annotations._range_wedge_angle_bounds` | 100 repeated 100k-tip circular range angle-bound calculations | 1.148038s | 0.468314s | 2.45x |
| `color_annotations._range_rect_tip_bounds` | 500k highlighted tips, identical rectangular x/y range bounds without temporary coordinate lists | 0.328038s | 0.295813s | 1.11x |
| `color_annotations.apply_label_colors` | 4096 label colors matched against 4096 real Matplotlib text artists | 1.220676s | 0.006927s | 176.21x |
| `color_annotations.parse_color_file` streaming parser | 500k mixed label/range/clade TSV rows with comments/blanks | 0.679407s | 0.642493s | 1.06x |
| `color_annotations.parse_color_file` lowercase fast path | 500k mixed label/range/clade TSV rows with comments/blanks, side-by-side original parser comparison | 0.735540s | 0.712200s | 1.03x |
| `color_annotations.parse_color_file` consolidated type dispatch | 500k mixed exact/mixed-case label/range/clade TSV rows with comments/blanks | 1.298780s | 1.181437s | 1.10x |
| tree color-file label recoloring via `apply_label_colors` | 28 converted plotter label-color blocks, each matching 1024 label colors against 1024 real Matplotlib text artists | 2.116338s | 0.046892s | 45.13x |
| `circular_layout.compute_circular_coords` | balanced 32768-tip tree, circular coordinates for all nodes | 0.4025s | 0.0746s | 5.4x |
| `circular_layout.compute_circular_coords` reverse-preorder tip ranges | balanced 32768-tip tree, circular coordinates for all nodes | 0.094112s | 0.069021s | 1.36x |
| `circular_layout._terminal_clades_direct` order-preserving child push | balanced 131072-tip tree, direct terminal helper, optimized helper baseline | 0.022347s | 0.017749s | 1.26x |
| `circular_layout._preorder_clades_direct` order-preserving child push | balanced 131072-tip tree, direct preorder helper, optimized helper baseline | 0.024254s | 0.020663s | 1.17x |
| `circular_layout.draw_circular_branches` | balanced 32768-tip tree, no-op axis with precomputed circular coordinates | 0.7247s | 0.4274s | 1.7x |
| `circular_layout.draw_circular_branches` cached arc fractions | balanced 32768-tip tree, no-op axis with precomputed circular coordinates | 0.393488s | 0.322528s | 1.2x |
| `circular_layout.draw_circular_branches` vectorized arc points | balanced 32768-tip tree, no-op axis with precomputed circular coordinates | 0.316143s | 0.157609s | 2.0x |
| `circular_layout.draw_circular_branches` direct single-pass drawing | balanced 32768-tip tree, no-op axis with precomputed circular coordinates | 0.178953s | 0.151419s | 1.18x |
| `circular_layout.draw_circular_branches` fallback child push | balanced 8192-tip tree, no-op axis with precomputed circular coordinates, identical plot-call count | 0.045360s | 0.043200s | 1.05x |
| `circular_layout.draw_circular_branches` batched LineCollections | balanced 2048-tip tree, precomputed circular coordinates, real Matplotlib Agg render | 1.072863s | 0.068903s | 15.57x |
| `circular_layout.circular_branch_points` trig hoist | 100k radial branches x 30 interpolated points | 0.635819s | 0.259429s | 2.45x |
| `circular_layout.draw_circular_gradient_branch` batched LineCollections | 1024 gradient radial branches x 30 color segments, real Matplotlib Agg render | 5.785908s | 1.013460s | 5.71x |
| `circular_layout.draw_circular_gradient_branches` whole-tree LineCollection | 1024 gradient radial branches x 30 color segments, real Matplotlib Agg render | 0.973249s | 0.373876s | 2.60x |
| `circular_layout.draw_circular_colored_arcs` whole-tree LineCollection | 1024 colored internal arcs x 61 polyline points, real Matplotlib Agg render | 0.159383s | 0.029367s | 5.43x |
| `circular_layout.draw_circular_multi_segment_branch` batched LineCollections | 1024 radial branches x 4 discrete color segments, real Matplotlib Agg render | 0.680789s | 0.657507s | 1.04x |
| `circular_layout.draw_circular_tip_labels` | balanced 32768-tip tree, no-op text axis for circular label placement | 0.0892s | 0.0253s | 3.5x |
| `circular_layout.draw_circular_tip_labels` localized lookups | balanced 32768-tip tree, no-op text axis for circular label placement, identical text-call count | 0.023353s | 0.021039s | 1.11x |
| `circular_layout` module import without `typing` startup | median cold subprocess import after replacing local annotation-only typing aliases with built-in generics | 0.025111s | 0.022998s | 1.09x |
| `plot_config` module import without `typing` startup | median cold subprocess import after replacing annotation-only typing aliases with built-in generics | 0.031106s | 0.028606s | 1.09x |
| `plot_config.build_parent_map` | balanced 65536-tip tree, parent map for shared tree plotting helpers | 0.1946s | 0.0266s | 7.3x |
| `plot_config.build_parent_map` unordered child push | balanced 131072-tip tree, parent map for shared tree plotting helpers, optimized helper baseline | 0.058822s | 0.047896s | 1.23x |
| `plot_config.compute_node_positions` phylogram layout | balanced 32768-tip tree, rectangular coordinates for all nodes | 0.3831s | 0.0494s | 7.8x |
| `plot_config.compute_node_positions` cladogram layout | balanced 32768-tip tree, rectangular coordinates for all nodes | 0.4833s | 0.0602s | 8.0x |
| `plot_config.compute_node_positions` phylogram reverse-preorder pass | balanced 32768-tip tree, rectangular coordinates for all nodes | 0.056965s | 0.043426s | 1.31x |
| `plot_config.compute_node_positions` cladogram reverse-preorder pass | balanced 32768-tip tree, rectangular coordinates for all nodes | 0.074095s | 0.065698s | 1.13x |
| `plot_config.compute_node_positions` phylogram inline child-y mean | balanced 65536-tip tree, rectangular coordinates with identical internal midpoint y-coordinates | 0.200336s | 0.151505s | 1.32x |
| `plot_config.compute_node_positions` cladogram inline child-y mean | balanced 65536-tip tree, rectangular coordinates with identical internal midpoint y-coordinates | 0.175507s | 0.137023s | 1.28x |
| `plot_config._terminal_clades_direct` order-preserving child push | balanced 131072-tip tree, direct terminal helper, optimized helper baseline | 0.039957s | 0.035492s | 1.13x |
| `plot_config._preorder_clades_direct` order-preserving child push | balanced 131072-tip tree, direct preorder helper, optimized helper baseline | 0.045027s | 0.040956s | 1.10x |
| `plot_config.draw_tree_branches` | balanced 32768-tip tree, no-op axis with precomputed rectangular coordinates | 0.1674s | 0.0645s | 2.6x |
| `plot_config.draw_tree_branches` batched LineCollections | balanced 2048-tip tree, precomputed rectangular coordinates, real Matplotlib Agg render | 1.400047s | 0.052574s | 26.63x |
| `plot_config.draw_tree_branches` lookup hoist | balanced 2048-tip tree, precomputed rectangular coordinates, real Matplotlib Agg render, side-by-side previous per-branch lookups | 0.346691s | 0.242124s | 1.43x |
| `plot_config.draw_tip_labels` | balanced 32768-tip tree, no-op text axis for rectangular label placement | 0.0760s | 0.0205s | 3.7x |
| `BipartitionSupportStats.get_bipartition_support_vals` | 4096-tip balanced tree, support on every internal node | 0.0932s | 0.0290s | 3.2x |
| `BipartitionSupportStats.get_bipartition_support_vals` direct output-order traversal | 4096-tip balanced tree, support on every internal node | 0.023487s | 0.003529s | 6.65x |
| `BipartitionSupportStats._get_bipartition_support_vals_direct` reverse-preorder names | balanced 4096-tip tree, support on every internal node, side-by-side previous visited-stack helper | 0.003063s | 0.002692s | 1.14x |
| `BipartitionSupportStats.get_bipartition_support_values` direct non-verbose extraction | 32768-tip balanced tree, support on every internal node | 0.081002s | 0.007867s | 10.30x |
| `BipartitionSupportStats._get_bipartition_support_values_array_direct` order-preserving child push | balanced 131072-tip tree, side-by-side previous `reversed(children)` array helper | 0.024698s | 0.021426s | 1.15x |
| `BipartitionSupportStats.calculate_threshold_stats` multi-threshold counts | 1M support values, 17 thresholds | 1.107047s | 0.146290s | 7.57x |
| `BipartitionSupportStats.calculate_threshold_stats` NumPy support values | 1M NumPy support values, 17 thresholds | 0.312667s | 0.037389s | 8.36x |
| `BipartitionSupportStats.calculate_threshold_stats` single NumPy threshold | 2M NumPy support values, one threshold cutoff | 0.175425s | 0.000868s | 202.1x |
| `BipartitionSupportStats.calculate_threshold_stats` large single-threshold Python list | 1M support values, one threshold cutoff, side-by-side previous Python generator count | 0.034144s | 0.028985s | 1.18x |
| `BipartitionSupportStats.calculate_threshold_stats` many-threshold Python list | 500k support values, 1000 thresholds | 0.133931s | 0.036207s | 3.70x |
| `BipartitionSupportStats.calculate_threshold_stats` empty thresholds | 2M support values, default no-threshold path, side-by-side previous unnecessary sort | 2.736035s | 0.000004s | 691083.37x |
| `BipartitionSupportStats._to_builtin` already-builtin verbose payload | 300k bipartition JSON rows, recursive NumPy conversion helper with identical values | 0.969412s | 0.563180s | 1.72x |
| `BipartitionSupportStats._to_builtin` builtin scalar fast path | 300k already-builtin verbose JSON rows, side-by-side previous scalar NumPy provenance check | 0.284150s | 0.251574s | 1.13x |
| `BipartitionSupportStats.build_json_output` verbose row construction | 300k bipartition JSON rows, side-by-side previous index-based `dict(...)` builder | 0.033431s | 0.028094s | 1.19x |
| `BipartitionSupportStats.run` verbose text output | 200k bipartition rows, mocked tree/read and identical stdout text | 0.065119s | 0.038501s | 1.69x |
| `BipartitionSupportStats._print_threshold_stats` batched threshold output | 100k threshold summary rows, captured stdout and identical text | 0.108712s | 0.098320s | 1.11x |
| `BipartitionSupportStats.run` cached read-only tree path | balanced 32768-tip cached tree with support on every internal node and two thresholds, non-verbose output mocked | 0.118359s | 0.004085s | 28.97x |
| `bipartition_support_stats` module import without eager Bio.Phylo/NumPy | cold subprocess import of bipartition-support-stats command module | 0.171310s | 0.069630s | 2.46x |
| `bipartition_support_stats` module import without eager stdlib JSON | median cold subprocess import after lazy `json.dumps` wrapper | 0.010258s | 0.009065s | 1.13x |
| `bipartition_support_stats` module import without eager stats helper | median cold subprocess import after lazy forwarding wrappers for summary helpers and JSON dumps | 0.028963s | 0.025592s | 1.13x |
| `bipartition_support_stats` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006199s | 0.003971s | 1.56x |
| `InternalBranchStats.get_internal_branch_lengths` | 4096-tip balanced tree, branch length on every internal node | 0.0922s | 0.0251s | 3.7x |
| `InternalBranchStats.calculate_internal_branch_stats` non-verbose | 4096-tip balanced tree, branch length on every internal node | 0.0259s | 0.0103s | 2.5x |
| `InternalBranchStats.calculate_internal_branch_stats` direct array non-verbose | 32768-tip balanced tree, branch length on every internal node | 0.010820s | 0.006748s | 1.60x |
| `InternalBranchStats._get_internal_branch_lengths_array_direct` order-preserving child push | balanced 131072-tip tree, side-by-side previous `reversed(children)` array helper | 0.024481s | 0.020200s | 1.21x |
| `InternalBranchStats.get_internal_branch_lengths` non-verbose order-preserving child push | balanced 131072-tip tree, `include_names=False`, side-by-side previous `reversed(children)` helper | 0.029400s | 0.024940s | 1.18x |
| `InternalBranchStats.get_internal_branch_lengths` direct non-verbose traversal | balanced 65536-tip tree, branch length on every internal node | 0.1472s | 0.0180s | 8.2x |
| `InternalBranchStats.get_internal_branch_lengths` direct verbose traversal | balanced 4096-tip tree, branch length and descendant tips for every internal node | 0.0258s | 0.0050s | 5.2x |
| `InternalBranchStats.get_internal_branch_lengths` reverse-preorder verbose helper | balanced 4096-tip tree, branch length and descendant tips for every internal node | 0.004030s | 0.003545s | 1.14x |
| `InternalBranchStats._get_internal_branch_lengths_direct` verbose child push | balanced 16384-tip tree, branch length and descendant tips for every internal node, side-by-side previous `reversed(children)` helper | 0.031354s | 0.015960s | 1.96x |
| `InternalBranchStats.run` verbose batched text output | 100k internal-branch rows, mocked tree/read and identical stdout text | 0.056698s | 0.043318s | 1.31x |
| `InternalBranchStats.run` cached read-only tree path | balanced 32768-tip cached tree with varied internal lengths, non-verbose summary output mocked | 0.111553s | 0.003668s | 30.41x |
| `internal_branch_stats` module import without eager Bio.Phylo | cold subprocess import of internal-branch-stats command module | 0.172597s | 0.113385s | 1.52x |
| `internal_branch_stats` module import without eager stats NumPy | cold subprocess import after lazy shared summary helper | 0.113385s | 0.069735s | 1.63x |
| `internal_branch_stats` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.011026s | 0.008896s | 1.24x |
| `internal_branch_stats` module import without eager stats helper | median cold subprocess import after lazy forwarding wrappers for summary and JSON helpers | 0.028662s | 0.023435s | 1.22x |
| `internal_branch_stats` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.012579s | 0.006258s | 2.01x |
| `TerminalBranchStats.calculate_terminal_branch_stats` non-verbose | 32768-tip balanced tree, branch length on every terminal branch | 0.0895s | 0.0868s | 1.03x |
| `TerminalBranchStats.calculate_terminal_branch_stats` direct array non-verbose | 32768-tip balanced tree, branch length on every terminal branch | 0.011780s | 0.007040s | 1.67x |
| `TerminalBranchStats._get_terminal_branch_lengths_array_direct` order-preserving child push | balanced 131072-tip tree, side-by-side previous `reversed(children)` array helper | 0.024257s | 0.020686s | 1.17x |
| `TerminalBranchStats.get_terminal_branch_lengths` direct non-verbose traversal | balanced 65536-tip tree, branch length on every terminal branch | 0.1379s | 0.0187s | 7.4x |
| `TerminalBranchStats.get_terminal_branch_lengths` direct verbose traversal | balanced 65536-tip tree, branch length and taxon name for every terminal branch | 0.1697s | 0.0480s | 3.5x |
| `TerminalBranchStats._get_terminal_branch_lengths_direct` order-preserving child push | balanced 65536-tip tree, branch length and taxon name for every terminal branch, side-by-side previous `reversed(children)` helper | 0.016234s | 0.012539s | 1.29x |
| `TerminalBranchStats.run` verbose batched text output | 100k terminal rows, mocked tree/read and identical stdout text | 0.052963s | 0.038226s | 1.39x |
| `TerminalBranchStats.run` non-verbose skipped length-list return | 1M terminal branch lengths, stats helper mocked, identical summary payload | 0.068538458s | 0.000000317s | 216448.85x |
| `TerminalBranchStats.run` cached read-only tree path | balanced 32768-tip cached tree with varied terminal lengths, non-verbose summary output mocked | 0.108551s | 0.003696s | 29.37x |
| `terminal_branch_stats` module import without eager Bio.Phylo | cold subprocess import of terminal-branch-stats command module | 0.169211s | 0.118903s | 1.42x |
| `terminal_branch_stats` module import without eager stats NumPy | cold subprocess import after lazy shared summary helper | 0.118903s | 0.070102s | 1.70x |
| `terminal_branch_stats` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.010540s | 0.009028s | 1.17x |
| `terminal_branch_stats` module import without eager stats helper | median cold subprocess import after lazy forwarding wrappers for summary and JSON helpers | 0.027848s | 0.024292s | 1.15x |
| `terminal_branch_stats` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006072s | 0.003834s | 1.58x |
| `Tree.calculate_total_branch_length_fast` | balanced 65536-tip tree, total branch length | 0.2511s | 0.0106s | 23.7x |
| `Tree.calculate_total_branch_length_fast` localized stack operations | balanced 65536-tip tree, total branch length, optimized helper baseline | 0.014608s | 0.012725s | 1.15x |
| `Tree.calculate_total_branch_length_and_terminal_count_fast` | balanced 65536-tip tree, total branch length plus terminal count | 0.3955s | 0.0119s | 33.2x |
| `Tree.calculate_total_branch_length_and_terminal_count_fast` localized count traversal | balanced 131072-tip tree, total branch length plus terminal count, optimized helper baseline | 0.041796s | 0.036840s | 1.13x |
| `Tree.calculate_terminal_count_fast` | balanced 65536-tip tree, terminal count only | 0.136221s | 0.010934s | 12.46x |
| `Tree.calculate_terminal_count_fast` clade root | balanced 32768-tip clade, terminal count only | 0.065122s | 0.003916s | 16.63x |
| `Tree.calculate_terminal_count_fast` localized stack operations | balanced 262144-tip tree, terminal count only, optimized helper baseline | 0.060390s | 0.055346s | 1.09x |
| `Tree.calculate_terminal_root_distances_fast` binary child push | balanced 262144-tip tree, root-to-terminal distances, side-by-side previous `reversed(children)` helper | 0.061517s | 0.056853s | 1.08x |
| `Tree.calculate_treeness` | balanced 65536-tip tree, internal branch length / total branch length | 0.3920s | 0.0127s | 31.0x |
| `Treeness.run` cached read-only tree path | balanced 32768-tip cached tree with varied branch lengths, output mocked | 0.122251s | 0.004065s | 30.08x |
| `treeness` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006082s | 0.004897s | 1.24x |
| `treeness` module import without `typing` startup | median cold subprocess import after postponing annotations and converting the annotation-only typing alias to a built-in annotation | 0.005941s | 0.003868s | 1.54x |
| `Tree.calculate_first_terminal_name_fast` direct leftmost descent | balanced depth-17 tree, first terminal lookup, side-by-side previous direct stack helper | 0.000002652s | 0.000001049s | 2.53x |
| `Tree.get_tip_names_from_tree` | balanced 65536-tip tree, terminal name extraction | 0.1356s | 0.0122s | 11.1x |
| `Tree.calculate_terminal_names_fast` child push loop | balanced 65536-tip tree, terminal name extraction with identical order | 0.018008s | 0.010621s | 1.70x |
| `Tree.calculate_terminal_names_fast` binary child push | balanced 131072-tip tree, terminal name extraction with identical order, side-by-side previous `reversed(children)` helper | 0.021721s | 0.018765s | 1.16x |
| `Tree.calculate_terminal_clades_fast` child push loop | balanced 65536-tip tree, terminal clade extraction with identical objects/order | 0.017121s | 0.010629s | 1.61x |
| `Tree.calculate_terminal_clades_fast` binary child push | balanced 131072-tip tree, terminal clade extraction with identical objects/order, side-by-side previous `reversed(children)` helper | 0.021770s | 0.018708s | 1.16x |
| `Tree.shared_tips` single intersection | 500k names per side with 250k shared names, side-by-side previous duplicate `set.intersection` call | 0.374889s | 0.119761s | 3.13x |
| `Tree` base module import without eager Bio.Phylo/NumPy | cold subprocess import of `phykit.services.tree.base` | 0.160683s | 0.064937s | 2.47x |
| `Tree` base module import without typing startup | median cold subprocess import after converting annotation-only typing names to built-in postponed annotations | 0.006148s | 0.003180s | 1.93x |
| `Tree` base module import without `hashlib` startup | median cold subprocess import after localizing cache-key hashing to `_get_file_hash` | 0.003903s | 0.000873s | 4.47x |
| `Tree._get_file_hash` raw stat cache key | 100k cache-key generations for one tree file | 0.000002628s | 0.000001821s | 1.44x |
| `TipLabels.run` tip-name extraction | balanced 65536-tip tree, terminal labels before output | 0.135883s | 0.017736s | 7.66x |
| `TipLabels.run` cached read-only tree path | balanced 32768-tip cached tree, output mocked | 0.112631s | 0.004871s | 23.12x |
| `TipLabels.run` JSON row literals | 500k terminal labels, JSON row payload construction only, side-by-side previous `dict(taxon=...)` rows | 0.092065s | 0.063262s | 1.45x |
| `tip_labels` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006193s | 0.004810s | 1.29x |
| `tip_labels` module import without `typing` startup | median cold subprocess import after postponing annotations and converting the annotation-only typing alias to a built-in annotation | 0.011704s | 0.007020s | 1.67x |
| `NeighborNet._compute_distances_from_alignment` | 180 taxa x 1200 sites, default p-distance, alphabet `ACGT-?NX` | 1.0017s | 0.0792s | 12.6x |
| `NeighborNet._compute_distance_matrix_from_equal_length_sequences` float count products | 700 taxa x 1200 sites, alphabet `ACGT-?NX` | 0.964211s | 0.047624s | 20.25x |
| `NeighborNet._read_alignment` | 50k FASTA records, mixed-case 120 bp each | 0.1525s | 0.0471s | 3.2x |
| `NeighborNet._read_alignment` shared ordered first-token parser | 50k FASTA records, mixed-case 120 bp each, legacy `SimpleFastaParser` baseline | 0.049768s | 0.044420s | 1.12x |
| `neighbor_net` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.228679s | 0.120682s | 1.89x |
| `neighbor_net` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.076428s | 0.027498s | 2.78x |
| `neighbor_net` module import without eager JSON/plot helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.009942s | 0.002834s | 3.51x |
| `neighbor_net` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006431s | 0.002320s | 2.77x |
| `NeighborNet._nj_circular_ordering` terminal extraction | balanced 65536-tip NJ-like tree, leaf ordering extraction after NJ construction mocked/completed | 0.236324s | 0.028959s | 8.16x |
| `NeighborNet._enumerate_circular_splits` all-taxa reuse | 160-taxon circular ordering, enumerate canonical non-trivial splits | 0.242429s | 0.203312s | 1.19x |
| `NeighborNet._enumerate_circular_splits` rolling side set | 160-taxon circular ordering, enumerate canonical non-trivial splits | 0.357076s | 0.105253s | 3.39x |
| `NeighborNet._enumerate_circular_splits` direct complement arcs | 160-taxon circular ordering, enumerate canonical non-trivial splits with identical split order | 0.077278s | 0.063128s | 1.22x |
| `NeighborNet._enumerate_circular_splits` complement reuse | 220-taxon circular ordering, enumerate canonical non-trivial splits with identical split order | 0.294846s | 0.200605s | 1.47x |
| `NeighborNet._enumerate_circular_splits` equal-size tiebreak | 220-taxon circular ordering, enumerate canonical non-trivial splits with identical split order, side-by-side previous sorted-half tiebreak | 0.136059s | 0.121730s | 1.12x |
| `NeighborNet._canonical_split` size-first min tiebreak | 9k mixed 20-vs-1180, 1180-vs-20, and 600-vs-600 bipartitions over 1200 taxa, side-by-side previous sorted-half tiebreak | 0.728545s | 0.230656s | 3.16x |
| `NeighborNet._is_circular_split` modulo-free boundary scan | 4096-taxon circular ordering, 610 contiguous/non-contiguous split checks with identical results | 0.230614s | 0.131844s | 1.75x |
| `NeighborNet._is_circular_split` third-boundary early rejection | 4096-taxon circular ordering, mixed 746 circular/non-circular split checks with identical accepted count | 0.141770s | 0.076098s | 1.86x |
| `NeighborNet._compute_split_directions` modulo-free gap scan | 4096-taxon circular ordering, 316 circular split direction vectors with identical output | 0.142790s | 0.079448s | 1.80x |
| `NeighborNet._compute_split_directions` cached taxon coordinates | 50 taxa/200 splits, 200 taxa/1000 splits, 500 taxa/3000 splits, identical direction dictionaries with cached gap positions | 0.001240s / 0.042574s / 0.371985s | 0.000936s / 0.025069s / 0.056481s | 1.33x / 1.70x / 6.59x |
| `NeighborNet._compute_split_directions` one-pass split centers | 500 circular splits over 2000 taxa with cached gap positions, side-by-side previous two generator sums | 0.053674s | 0.044812s | 1.20x |
| `NeighborNet` circular split filtering plus direction setup | 8192-taxon circular ordering, mixed 1490 circular/non-circular split checks with identical 745 accepted direction vectors | 0.631656s | 0.345524s | 1.83x |
| `NeighborNet._compute_distance_matrix_from_equal_length_sequences` clean ASCII direct comparison | 260 taxa x 6000 sites, alphabet `ACGT`, p-distance matrix | 0.161045s | 0.120611s | 1.34x |
| `NeighborNet._read_distance_matrix` row-slice CSV fill | 800 x 800 labeled CSV distance matrix, identical taxa and matrix values | 0.124120s | 0.094420s | 1.31x |
| `NeighborNet._read_distance_matrix` simple CSV `fromstring` parser | 1200 x 1200 labeled CSV distance matrix, identical taxa and matrix values, quoted-field fallback preserved | 0.215591s | 0.126250s | 1.71x |
| `NeighborNet._nj_circular_ordering` lower-triangle row conversion | 1500 x 1500 symmetric distance matrix, identical BioPython lower-triangle payload | 0.186708s | 0.035887s | 5.20x |
| `NeighborNet._estimate_split_weights` | 80 taxa, 3080 circular splits, dense NNLS design matrix | 3.5731s | 2.9731s | 1.2x |
| `NeighborNet._build_splits_graph` | 80 taxa, first 20 circular splits, 21-node Buneman graph | 3.3135s | 0.0010s | 3262.5x |
| `NeighborNet` / `ConsensusNetwork` / `QuartetNetwork` split-graph extent helper | 500k node positions, identical maximum x/y extent without temporary coordinate lists | 0.070550s | 0.029240s | 2.41x |
| `NeighborNet` / `ConsensusNetwork` network edge rendering | 80 taxa, 20 circular splits, real Matplotlib Agg internal and pendant edge render | 0.025885s | 0.012323s | 2.10x |
| `NeighborNet` / `ConsensusNetwork` unlabeled fallback point rendering | 4096 taxa without accepted splits, real Matplotlib Agg point render | 0.756119s | 0.026598s | 28.43x |
| `NeighborNet.run` batched text split output | 100k positive split rows, captured stdout and identical text | 0.055873s | 0.043923s | 1.27x |
| `ParsimonyScore._parse_alignment` | 50k FASTA records, mixed-case 120 bp each | 0.1198s | 0.0382s | 3.1x |
| `ParsimonyScore._parse_alignment` shared first-token parser | 50k FASTA records, mixed-case 120 bp each, legacy `SimpleFastaParser` baseline | 0.051323s | 0.037320s | 1.38x |
| `parsimony_score` module import without eager FASTA parser | cold subprocess import after lazy Bio.SeqIO.FastaIO import | 0.227281s | 0.124062s | 1.83x |
| `parsimony_score` module import without eager NumPy | cold subprocess import after lazy NumPy proxy | 0.077951s | 0.025578s | 3.05x |
| `RelativeRateTest._run_single` | 140 ingroup taxa x 1500 sites plus outgroup, alphabet `ACGT-?NX` | 1.4043s | 0.0275s | 51.1x |
| `RelativeRateTest._tajima_test` long ASCII helper | five 1M-site Tajima tests, alphabet `ACGTn` | 1.111929s | 0.012380s | 89.82x |
| `RelativeRateTest._tajima_test` clean ASCII shortcut | 2M-site Tajima triplet, alphabet `ACGT`, side-by-side previous validity-mask vector path | 0.008916s | 0.007849s | 1.14x |
| `RelativeRateTest._identify_outgroup` singleton scan | rooted tree with 32768-tip ingroup plus singleton outgroup | 0.063596s | 0.000001s | 47668.9x |
| `RelativeRateTest._chi2_sf` df=1 scalar p-value | cold process, Tajima relative-rate chi-square survival probability | 0.538501s | 0.000002958s | 182049.0x |
| `RelativeRateTest._run_pairwise_tests_vectorized` df=1 vector p-values | 1M synthetic chi-square statistics, large informative-pair branch | 0.151878s | 0.010071s | 15.08x |
| `RelativeRateTest._run_pairwise_tests_vectorized` zip-based result rows | 1,124,250 pairwise result rows from 1500 taxa, side-by-side previous index lookups | 0.698618s | 0.558125s | 1.25x |
| `RelativeRateTest._run_pairwise_tests_vectorized` row-wise pairwise results | 1,124,250 pairwise result rows from 1500 taxa, side-by-side previous `np.triu_indices` gather path | 3.519097s | 2.335869s | 1.51x |
| `RelativeRateTest._ingroup_ascii_matrix` joined byte-buffer setup | 1000 repeated 140-ingroup-taxon x 1500-site ASCII matrix builds, side-by-side previous per-taxon `np.frombuffer` plus `np.vstack` setup | 0.845237s | 0.296323s | 2.85x |
| `RelativeRateTest._bonferroni` | 1M synthetic p-values | 0.096721s | 0.033725s | 2.9x |
| `RelativeRateTest._fdr` | 1M synthetic p-values | 0.698518s | 0.164136s | 4.3x |
| `RelativeRateTest._add_multiple_testing_corrections` zip-based correction assignment | 1M result rows with Bonferroni/FDR arrays, side-by-side previous index lookups | 0.293454s | 0.222238s | 1.32x |
| `RelativeRateTest` small multiple-testing corrections without NumPy startup | cold subprocess, 7 p-values through Bonferroni and FDR helpers | 0.087782s | 0.027210s | 3.23x |
| `RelativeRateTest._output_single` batched text output | 100k pairwise relative-rate rows, captured stdout and identical text | 0.124079s | 0.111355s | 1.11x |
| `RelativeRateTest._output_batch` batched text output | 100k pair summary rows, captured stdout and identical text | 0.121325s | 0.109700s | 1.11x |
| `RelativeRateTest._output_batch` gene-result summary | 100k taxon-pair result groups x 8 gene results, identical reject counts and median chi-square values | 0.116062s | 0.107860s | 1.08x |
| `RelativeRateTest.run` batch list cleanup | 500k alignment-list rows with comments/blanks, cleanup before per-alignment analysis | 0.088629s | 0.067571s | 1.31x |
| `RelativeRateTest.run` batch list path resolver | 100k relative alignment-list rows, alignment analysis mocked | 0.461438s | 0.029321s | 15.74x |
| `RelativeRateTest._plot_heatmap` matrix/significance setup | 80 taxa, 3160 pairwise FDR results | 0.9333s | 0.0015s | 605.8x |
| `RelativeRateTest._plot_heatmap` significance marker rendering | 4096 significant heatmap cells, real Matplotlib Agg star render | 0.650818s | 0.022971s | 28.33x |
| `RelativeRateTest._plot_heatmap` no-significance marker lookup | 2500 x 2500 FDR matrix with diagonal NaNs and no significant cells, side-by-side previous coordinate extraction | 0.021465s | 0.001585s | 13.54x |
| `RelativeRateTest._significant_heatmap_cells` sparse marker coordinates | 2500 x 2500 FDR matrix with diagonal NaNs and sparse symmetric significant cells, side-by-side previous `np.where` coordinate extraction | 0.030663s | 0.008484s | 3.61x |
| `relative_rate_test` module import without eager NumPy | cold subprocess import after lazy NumPy proxy | 0.088372s | 0.031916s | 2.77x |
| `relative_rate_test` module import without eager JSON/plot helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.014537s | 0.005749s | 2.53x |
| `relative_rate_test` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in annotations | 0.006280s | 0.004398s | 1.43x |
| `relative_rate_test` module import without eager file helper | median cold subprocess import after lazy alignment-reader wrapper | 0.048166s | 0.029581s | 1.63x |
| `relative_rate_test` module import without eager `pathlib` | median cold subprocess import after lazy batch-list `Path` wrapper | 0.044962s | 0.038405s | 1.17x |
| `RelativeRateTest.run` cached read-only tree setup | balanced 32768-tip rooted cached tree, outgroup identification with alignment analysis and output mocked | 0.393320s | 0.000271s | 1451.37x |
| `RelativeRateTest._identify_outgroup` generic fallback singleton scan | nonstandard root with 300k-tip ingroup, singleton outgroup, and second 300k-tip ingroup | 0.317954s | 0.000007459s | 42627.3x |
| `CovaryingEvolutionaryRates.correct_branch_lengths` | exact-matching balanced 256-tip gene trees plus reference | 0.3197s | 0.0044s | 72.4x |
| `CovaryingEvolutionaryRates.correct_branch_lengths` same-tree branch-length reuse | balanced 8192-tip tree passed as both gene trees and reference, side-by-side previous duplicate branch-length extraction | 4.022530s | 2.395759s | 1.68x |
| `CovaryingEvolutionaryRates.correct_branch_lengths` same-tree/reference shortcut | balanced 32768-tip tree passed as both gene trees and reference, side-by-side previous same-gene branch-map reuse path | 4.946255s | 2.663574s | 1.86x |
| `CovaryingEvolutionaryRates._correct_branch_lengths_from_exact_splits` direct order traversal | exact-matching balanced 4096-tip gene trees plus reference | 0.105679s | 0.090155s | 1.17x |
| `CovaryingEvolutionaryRates._branch_lengths_by_tipset` direct postorder | balanced 65536-tip tree, descendant-tip branch-length map | 0.276784s | 0.140246s | 1.97x |
| `CovaryingEvolutionaryRates._branch_lengths_by_tipset` reverse-preorder helper | balanced 65536-tip tree, descendant-tip branch-length map | 0.405383s | 0.360719s | 1.12x |
| `CovaryingEvolutionaryRates._branch_lengths_by_tipset` binary frozenset union | balanced 32768-tip tree, descendant-tip branch-length map | 0.144733s | 0.126943s | 1.14x |
| `CovaryingEvolutionaryRates._correct_branch_lengths_from_exact_splits` reference direct postorder | exact-matching balanced 16384-tip gene trees plus reference | 0.188390s | 0.153965s | 1.22x |
| `CovaryingEvolutionaryRates._reference_tipsets_and_order` reverse-preorder helper | exact-matching balanced 16384-tip reference tree metadata | 0.064371s | 0.053343s | 1.21x |
| `CovaryingEvolutionaryRates._reference_tipsets_and_order` order-preserving child push | balanced 32768-tip reference tree metadata, optimized helper baseline | 0.105540s | 0.092782s | 1.14x |
| `CovaryingEvolutionaryRates._zscore` centered sum-of-squares path | 260 / 1000 / 5000 / 50000 branch-length values, side-by-side previous separate `np.mean` and `np.std` reductions | 0.000019512s / 0.000015718s / 0.000019248s / 0.000085509s | 0.000008080s / 0.000006586s / 0.000009209s / 0.000049810s | 2.41x / 2.39x / 2.09x / 1.72x |
| `CovaryingEvolutionaryRates._pearsonr` dot-product correlation | 10k corrected branch-length pairs, SciPy already warm, side-by-side previous normalized-vector helper | 0.000086246s | 0.000045948s | 1.88x |
| `CovaryingEvolutionaryRates._tips_to_prune_for_shared` ordered prune-list scan | 400k tree-zero tips, tree-one tips, reference tips, and 300k shared tips | 0.277122s | 0.069708s | 3.98x |
| `CovaryingEvolutionaryRates.get_indices_of_outlier_branch_lengths` flat outlier indices | 3M corrected branch lengths with sparse `abs(x) > 5` and `NaN` outliers, side-by-side previous `np.where(...)[0]` extraction | 0.012458s | 0.005730s | 2.17x |
| `CovaryingEvolutionaryRates.get_indices_of_outlier_branch_lengths` empty-prior direct return | 3M corrected branch lengths with sparse `abs(x) > 5` and `NaN` outliers, no existing outlier indices | 0.105751s | 0.063676s | 1.66x |
| `CovaryingEvolutionaryRates.remove_outliers_based_on_indices` direct array filter | 1M corrected branch lengths with 10 outlier indices, side-by-side previous NumPy mask plus Python enumerate filter | 0.128672s | 0.031693s | 4.06x |
| `CovaryingEvolutionaryRates.remove_outliers_based_on_indices` sparse list slices | 1M branch-length list with 10 outlier indices, preserving order and negative-index handling | 0.047255s | 0.007984s | 5.92x |
| `covarying_evolutionary_rates` module import | cold subprocess import, avoid eager `scipy.stats` import | 0.780821s | 0.191193s | 4.1x |
| `covarying_evolutionary_rates` module import without eager NumPy | cold subprocess import after lazy NumPy proxy | 0.108575s | 0.048410s | 2.24x |
| `covarying_evolutionary_rates` module import without eager concurrent futures | cold subprocess import after lazy executor and completion proxies | 0.048199s | 0.031029s | 1.55x |
| `covarying_evolutionary_rates` module import without eager pickle/plot config | cold subprocess import after lazy pickle proxy and localized `PlotConfig` import | 0.035961s | 0.025588s | 1.41x |
| `covarying_evolutionary_rates` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.007304s | 0.005600s | 1.30x |
| `CovaryingEvolutionaryRates.run` verbose text output | 100k verbose branch rows, captured stdout and identical text | 0.095201s | 0.082940s | 1.15x |
| `LastCommonAncestorSubtree.run` | balanced 4096-tip tree, write 1024-tip MRCA subtree | 0.2959s | 0.2702s | 1.1x |
| `LastCommonAncestorSubtree.run` LCA lookup | balanced 4096-tip tree, 1024-tip MRCA subtree lookup with output stubbed | 0.2210s | 0.0320s | 6.9x |
| `LastCommonAncestorSubtree._find_lca_subtree` one-pass direct lookup | balanced 4096-tip tree, 1024-tip MRCA subtree lookup | 0.0207s | 0.0019s | 10.9x |
| `LastCommonAncestorSubtree._find_lca_subtree` parent-depth lookup | 200 repeated balanced 4096-tip tree, 1024-tip MRCA subtree lookups | 0.712718s | 0.358302s | 1.99x |
| `LastCommonAncestorSubtree._find_lca_subtree` selected terminal lookup | balanced 65536-tip tree, 1024 spread taxa before MRCA lookup | 0.070156s | 0.047469s | 1.48x |
| `LastCommonAncestorSubtree._find_lca_subtree` early stop after selected terminals | balanced 32768-tip tree, first 1024 taxa before MRCA lookup | 0.031224s | 0.001963s | 15.91x |
| `LastCommonAncestorSubtree._find_lca_subtree` order-preserving stack push | balanced 65536-tip tree, 1024 spread taxa, optimized helper baseline | 0.047902s | 0.039079s | 1.23x |
| `LastCommonAncestorSubtree._find_lca_subtree` parent-depth target scan | 20 repeated balanced 32768-tip tree, 2048 spread-taxa MRCA lookups, side-by-side previous `targets[1:]` loop | 1.827804s | 1.504977s | 1.21x |
| `LastCommonAncestorSubtree._find_lca_subtree` single-taxon direct lookup | balanced 32768-tip tree, rightmost single terminal lookup, side-by-side previous parent-depth path | 1.893951s | 0.781276s | 2.42x |
| `LastCommonAncestorSubtree.run` cached read-only tree setup | balanced 4096-tip cached tree, 1024-tip MRCA subtree lookup with output stubbed | 0.039566s | 0.008036s | 4.92x |
| `LastCommonAncestorSubtree.run` no-copy large-tree LCA setup | balanced 32768-tip cached tree, opposite terminal MRCA lookup, side-by-side previous pickle-copy `common_ancestor` core | 0.182315s | 0.013709s | 13.30x |
| `last_common_ancestor_subtree` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006235s | 0.005155s | 1.21x |
| `last_common_ancestor_subtree` module import without `typing` startup | median cold subprocess import after postponing annotations and converting the annotation-only typing alias to a built-in annotation | 0.006208s | 0.004222s | 1.47x |
| `last_common_ancestor_subtree` module import without eager file helper | median cold subprocess import after lazy taxa-list reader wrapper | 0.044738s | 0.022908s | 1.95x |
| `MonophylyCheck.run` exact-clade path | balanced 4096-tip tree, 1024-taxon exact clade, output stubbed | 3.4868s | 0.0199s | 174.9x |
| `MonophylyCheck.run` cached read-only tree setup | balanced 32768-tip cached tree, taxa read, clade resolution, bootstrap stats, and output mocked | 0.095587s | 0.001858s | 51.4x |
| `MonophylyCheck._find_exact_clade_by_taxa` count-based lookup | balanced 32768-tip tree, all 32768 taxa exactly matching root clade | 0.182382s | 0.126754s | 1.44x |
| `MonophylyCheck._find_exact_clade_by_taxa` direct count traversal | balanced 65536-tip tree, all taxa exactly matching root clade, side-by-side previous generic postorder count path | 0.364151s | 0.098613s | 3.69x |
| `MonophylyCheck._find_exact_clade_by_taxa_direct` stack-frame counts | balanced 65536-tip tree, first-quarter exact clade, side-by-side previous dictionary-backed direct count path | 0.032201s | 0.015213s | 2.12x |
| `MonophylyCheck.get_bootstrap_statistics` direct support traversal | balanced 32768-tip clade, support on every internal node | 0.066024s | 0.006910s | 9.56x |
| `MonophylyCheck._collect_bootstrap_values_direct` unordered support scan | balanced 131072-tip clade, support on every internal node, optimized helper baseline | 0.028672s | 0.021346s | 1.34x |
| `MonophylyCheck.print_results` batched text output | 200k mixed status/support rows, captured stdout and identical text | 0.301773s | 0.271484s | 1.11x |
| `monophyly_check` module import without annotation-only Bio.Phylo | cold subprocess import after postponed annotations and removing `Newick` import | 0.137385s | 0.030312s | 4.53x |
| `monophyly_check` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.010405s | 0.009056s | 1.15x |
| `monophyly_check` module import without eager stats helper | median cold subprocess import after lazy forwarding wrapper for bootstrap summary helper | 0.030237s | 0.024028s | 1.26x |
| `monophyly_check` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006314s | 0.004114s | 1.53x |
| `monophyly_check` module import without eager file helper | median cold subprocess import after lazy taxa-list reader wrapper | 0.070020s | 0.061728s | 1.13x |
| `HiddenParalogyCheck.run` sequential exact-clade path | balanced 4096-tip tree, eight 512-taxon exact clades, output stubbed | 24.3482s | 0.0442s | 551.0x |
| `HiddenParalogyCheck.run` cached read-only master tree setup | balanced 4096-tip cached tree, eight 512-taxon exact clades, output stubbed | 0.043373s | 0.026190s | 1.66x |
| `HiddenParalogyCheck._process_clade_batch` exact-clade path | balanced 2048-tip tree, four 256-taxon exact clades in one batch | 3.0635s | 0.0180s | 170.6x |
| `HiddenParalogyCheck._process_clade_batch` shared exact-clade index | balanced 2048-tip tree, four 256-taxon exact clades in one batch | 0.0169s | 0.000055s | 307.3x |
| `HiddenParalogyCheck._build_exact_clade_index` direct postorder | balanced 32768-tip tree, exact descendant-taxon index for every clade | 0.262937s | 0.131806s | 2.00x |
| `HiddenParalogyCheck._build_exact_clade_index` binary child-set union | balanced 32768-tip tree, exact descendant-taxon index for every clade | 0.131806s | 0.121338s | 1.09x |
| `HiddenParalogyCheck._terminal_names_direct` | balanced 32768-tip clade, collect terminal names in batch fallback | 0.065931s | 0.008021s | 8.22x |
| `HiddenParalogyCheck._terminal_names_direct` set traversal | balanced 65536-tip clade, collect identical terminal-name set | 0.019487s | 0.012567s | 1.55x |
| `HiddenParalogyCheck` requested-clade shared-taxa intersection | 500k three-taxon requested clades against 100k master tree tips, duplicate/off-tree semantics preserved | 0.247522s | 0.190846s | 1.30x |
| `HiddenParalogyCheck.print_results` text output | 200k clade status rows, identical stdout text | 0.028718s | 0.004742s | 6.06x |
| `HiddenParalogyCheck.read_clades_file` bulk read split | 500k three-taxon clade rows with blank-line compatibility | 0.572789s | 0.435850s | 1.31x |
| `HiddenParalogyCheck.read_clades_file` direct line iteration | 500k three-taxon clade rows with blank-line compatibility, side-by-side previous `read().splitlines()` parser | 0.094809s | 0.092886s | 1.02x |
| `hidden_paralogy_check` module import without eager Bio.Phylo | cold subprocess import after lazy Phylo reader proxy | 0.135706s | 0.029766s | 4.56x |
| `hidden_paralogy_check` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.012157s | 0.010456s | 1.16x |
| `hidden_paralogy_check` module import without eager multiprocessing | median cold subprocess import after lazy multiprocessing proxy and localized `partial` import | 0.030027s | 0.023741s | 1.26x |
| `hidden_paralogy_check` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.035490s | 0.032260s | 1.10x |
| `BranchLengthMultiplier.run` | balanced 8192-tip tree, multiply every branch by 2 | 0.1126s | 0.1015s | 1.1x |
| `BranchLengthMultiplier.multiply_branch_lengths_by_factor` traversal | balanced 65536-tip tree, multiply every branch by 2 | 0.1801s | 0.0117s | 15.4x |
| `BranchLengthMultiplier.multiply_branch_lengths_by_factor` inline scaling | balanced 65536-tip tree, multiply every branch by 2.5 during direct stack traversal | 0.039696s | 0.032010s | 1.24x |
| `BranchLengthMultiplier._multiply_standard_tree_branch_lengths` localized stack operations | balanced 131072-tip tree, multiply every branch by 2.5, side-by-side previous direct helper median | 0.068762s | 0.055513s | 1.24x |
| `BranchLengthMultiplier.run` factor-one read-only cached tree path | balanced 32768-tip cached tree, factor `1.0`, output stubbed | 0.212082s | 0.000031s | 6841.35x |
| `BranchLengthMultiplier._count_standard_tree_branch_lengths` localized stack operations | balanced 131072-tip tree, factor-one JSON branch count, side-by-side previous direct helper median | 0.015458s | 0.014903s | 1.04x |
| `branch_length_multiplier` module import without eager Bio.Phylo/NumPy | cold subprocess import of branch-length multiplier command module | 0.169961s | 0.065100s | 2.61x |
| `branch_length_multiplier` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006216s | 0.004949s | 1.26x |
| `branch_length_multiplier` module import without `typing` startup | median cold subprocess import after converting the remaining annotation-only typing alias to a built-in postponed annotation | 0.006049s | 0.003979s | 1.52x |
| `RenameTreeTips.run` | balanced 32768-tip tree, rename 16384 tips from id map | 0.4165s | 0.2528s | 1.6x |
| `RenameTreeTips.run` non-JSON match detection | balanced 32768-tip cached tree, rename 16384 tips, output stubbed | 0.895882s | 0.821195s | 1.09x |
| `RenameTreeTips.replace_tip_names` terminal traversal | balanced 65536-tip tree, rename half the tips from id map | 0.1275s | 0.0124s | 10.3x |
| `RenameTreeTips.replace_tip_names` single-pass direct rename | balanced 65536-tip tree, rename half the tips from id map | 0.042480s | 0.037260s | 1.14x |
| `RenameTreeTips._replace_standard_tree_tip_names` unordered child push | balanced 65536-tip tree, rename half the tips from id map, side-by-side previous order-preserving helper | 0.031663s | 0.024518s | 1.29x |
| `RenameTreeTips.run` no-match read-only cached tree path | balanced 32768-tip cached tree, id map with no matching tips, output stubbed | 0.209356s | 0.008288s | 25.26x |
| `RenameTreeTips.count_matching_tip_names` unordered read-only scan | balanced 131072-tip tree, every 16th tip in id map, optimized helper baseline | 0.023922s | 0.016833s | 1.42x |
| `RenameTreeTips.has_matching_tip_name` unordered read-only scan | balanced 131072-tip tree, no matching tips, optimized helper baseline | 0.022242s | 0.015908s | 1.40x |
| `RenameTreeTips` empty id-map short-circuit | balanced 32768-tip tree, empty id map for has/count/replace helper paths | 0.011177500s / 0.006682375s / 0.013954125s | 0.000005042s / 0.000002667s / 0.000003791s | 2216.80x / 2505.73x / 3680.91x |
| `rename_tree_tips` module import without eager Bio.Phylo | cold subprocess import of rename-tree-tips command module | 0.162538s | 0.065407s | 2.49x |
| `rename_tree_tips` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006052s | 0.004723s | 1.28x |
| `RootTree.run` | balanced 32768-tip tree, root with one-tip outgroup | 0.3064s | 0.1872s | 1.6x |
| `root_tree` module import without eager Bio.Phylo | cold subprocess import after lazy root-with-outgroup proxy | 0.120389s | 0.025771s | 4.67x |
| `root_tree` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006306s | 0.005138s | 1.23x |
| `root_tree` module import without eager file helper | median cold subprocess import after lazy outgroup-list reader wrapper | 0.040920s | 0.038210s | 1.07x |
| `PruneTree.run` | balanced 32768-tip tree, prune one tip | 0.3089s | 0.1792s | 1.7x |
| `PruneTree.run` keep-mode complement | 8192-tip tree, keep 4096 named tips, pruning/output stubbed | 0.1859s | 0.0026s | 70.4x |
| `PruneTree.run` keep-mode terminal-name selection | balanced 32768-tip tree, keep half of tips, pruning/output stubbed | 0.0671s | 0.0072s | 9.4x |
| `PruneTree.run --ignore-branch-labels` terminal-name selection | balanced 32768-tip tree, keep half of tips with every third tip labeled | 0.0726s | 0.0122s | 5.9x |
| `PruneTree._strip_branch_label` plain-name regex guard | 32768 tip names x 40 scans, every third tip labeled | 0.358330s | 0.246274s | 1.45x |
| `PruneTree.run` keep-all read-only cached tree path | balanced 32768-tip cached tree, keep all tips, output stubbed | 0.220430s | 0.014711s | 14.98x |
| `prune_tree` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006864s | 0.005426s | 1.27x |
| `prune_tree` module import without `typing` startup | median cold subprocess import after postponing annotations and converting the annotation-only typing alias to a built-in annotation | 0.006264s | 0.004162s | 1.50x |
| `prune_tree` module import without eager file helper | median cold subprocess import, interleaved lazy import vs eager-equivalent `phykit.helpers.files` preload | 0.026866s | 0.026566s | 1.01x |
| `RobinsonFouldsDistance.calculate_robinson_foulds_distance` compact split ids | identical balanced 16384-tip trees, rooted descendant split comparison | 0.069501s | 0.052131s | 1.33x |
| `NearestNeighborInterchange._build_parent_map` direct traversal | balanced 32768-tip tree, parent map for NNI generation | 0.101185s | 0.010405s | 9.72x |
| `NearestNeighborInterchange._build_parent_map` unordered child push | balanced 65536-tip tree, parent map for NNI generation, optimized helper baseline | 0.020223s | 0.014583s | 1.39x |
| `NearestNeighborInterchange.get_neighbors` direct level-order nonterminals | balanced 8192-tip tree, NNI generation with tree-copy creation stubbed, side-by-side previous `get_nonterminals(order="level")` loop | 0.021400s | 0.011773s | 1.82x |
| `nearest_neighbor_interchange` module import without eager Bio.Phylo | cold subprocess import after lazy `Phylo.write` proxy and postponed annotations | 0.123042s | 0.025484s | 4.83x |
| `nearest_neighbor_interchange` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.007619s | 0.006372s | 1.20x |
| `nearest_neighbor_interchange` module import without eager pickle | median cold subprocess import after lazy pickle proxy | 0.006408s | 0.005053s | 1.27x |
| `nearest_neighbor_interchange` module import without `typing` startup | median cold subprocess import after converting the runtime branch-spec alias and annotations to built-in generics | 0.028164s | 0.023944s | 1.18x |
| `DVMC.determine_dvmc` terminal count reuse | balanced 32768-tip tree with unit branch lengths | 0.1459s | 0.0779s | 1.9x |
| `DVMC.determine_dvmc` one-pass terminal distances | balanced 65536-tip tree with unit branch lengths | 0.1630s | 0.0171s | 9.5x |
| `DVMC.determine_dvmc` scalar terminal-distance stats | balanced 65536-tip tree with varied terminal branch lengths | 0.014995s | 0.013724s | 1.09x |
| `DVMC.determine_dvmc` NumPy-free fallback stats | cold subprocess fallback tree with 200k terminal distances and identical scalar result | 0.129657s | 0.051419s | 2.52x |
| `dvmc` module import without eager Bio.Phylo/NumPy | cold subprocess import of DVMC command module | 0.166131s | 0.064778s | 2.56x |
| `dvmc` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006003s | 0.004850s | 1.24x |
| `dvmc` module import without `typing` startup | median cold subprocess import after converting the annotation-only typing alias to a built-in postponed annotation | 0.007384s | 0.005897s | 1.25x |
| `DVMC.run` cached read-only tree path | balanced 32768-tip cached tree with varied root-to-tip lengths, output mocked | 0.109749s | 0.004151s | 26.44x |
| `TotalTreeLength.run` cached read-only tree path | balanced 32768-tip cached tree, output mocked | 0.106133s | 0.003879s | 27.36x |
| `EvolutionaryRate.run` cached read-only tree path | balanced 32768-tip cached tree, output mocked | 0.111459s | 0.004104s | 27.16x |
| `evolutionary_rate` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006227s | 0.004876s | 1.28x |
| `evolutionary_rate` module import without `typing` startup | median cold subprocess import after postponing annotations and converting the annotation-only typing alias to a built-in annotation | 0.011586s | 0.008064s | 1.44x |
| `total_tree_length` module import without eager Bio.Phylo/NumPy | cold subprocess import of total-tree-length command module | 0.153291s | 0.065309s | 2.35x |
| `total_tree_length` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006164s | 0.004823s | 1.28x |
| `total_tree_length` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.011637s | 0.007487s | 1.55x |
| `LTT` shared terminal/depth setup | balanced 32768-tip tree, compute gamma and LTT data | 0.3442s | 0.2626s | 1.3x |
| `LTT._terminal_clades` setup | balanced 65536-tip tree, terminal clade list for gamma/LTT helpers | 0.1272s | 0.0172s | 7.4x |
| `LTT._terminal_clades` order-preserving binary push | balanced 131072-tip tree, terminal clade list with identical tip order, optimized helper baseline | 0.022452s | 0.017732s | 1.27x |
| `LTT._depths_from_root` direct stack traversal | balanced 32768-tip tree, root-depth map setup | 0.014824s | 0.010569s | 1.40x |
| `LTT._compute_gamma` + `_compute_ltt` internal-depth scans | balanced 32768-tip tree, shared terminal list and depth map | 0.1966s | 0.0398s | 4.9x |
| `LTT._compute_gamma` cumulative-stat loop | balanced 65536-tip tree, shared terminal list and depth map, identical gamma, p-value, branching times, and intervals | 0.039136s | 0.034603s | 1.13x |
| `LTT._compute_gamma` combined ST/stat accumulation | 65536 internode intervals, side-by-side previous separate `ST` generator pass plus cumulative-stat loop | 0.021182s | 0.006768s | 3.13x |
| `LTT._compute_gamma` p-value calculation | cold process, two-tailed standard-normal p-value for gamma statistic | 0.557722s | 0.000003708s | 150410.5x |
| `LTT._compute_gamma` / `_compute_ltt` streaming tip-height max | 1M terminal depths, identical root-relative tree height from shared depth map | 0.070200s | 0.044918s | 1.56x |
| `LTT.run` cached read-only tree path | balanced 16384-tip cached tree, text output mocked | 0.171729s | 0.034802s | 4.93x |
| `LTT._output_text` batched verbose rows | 100k branching-time rows plus 100k LTT rows, captured stdout and identical text | 0.073615s | 0.052625s | 1.40x |
| `LTT._compute_gamma_and_ltt` shared internal-depth calculation | balanced 32768-tip tree, shared terminal list and depth map, identical gamma/LTT output | 0.028798s | 0.018729s | 1.54x |
| `LTT._compute_gamma_and_ltt` internal-depth child push | balanced 131072-tip tree, shared terminal list and depth map, identical gamma/LTT output | 0.092145s | 0.078426s | 1.17x |
| `LTT._compute_gamma_and_ltt` LTT row construction | 1M sorted internal depths, identical LTT rows, side-by-side previous `internal_depths[1:]` loop | 0.363777s | 0.274531s | 1.33x |
| `ltt` module import without eager JSON/plot helpers | median cold subprocess import after localizing PlotConfig and lazy JSON wrapper | 0.013677s | 0.006328s | 2.16x |
| `ltt` module import without `typing` startup | median cold subprocess import after converting annotation-only typing names to built-in postponed annotations | 0.006019s | 0.004281s | 1.41x |
| `TreenessOverRCV.run` read-only treeness setup | balanced 32768-tip cached tree, RCV calculation and output mocked, side-by-side previous copied tree read | 0.241418s | 0.008151s | 29.62x |
| `treeness_over_rcv` module import | cold subprocess import, defer alignment RCV helper import until run | 0.165134s | 0.150374s | 1.10x |
| `treeness_over_rcv` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006141s | 0.004849s | 1.27x |
| `treeness_over_rcv` module import without `typing` startup | median cold subprocess import after postponing annotations and converting the annotation-only typing alias to a built-in annotation | 0.006743s | 0.003997s | 1.69x |
| `InternodeLabeler.run` | balanced 32768-tip tree, label all internal nodes | 0.4292s | 0.2612s | 1.6x |
| `InternodeLabeler.add_labels_to_tree` | balanced 65536-tip tree, preorder labels on every internal node | 0.1264s | 0.0168s | 7.5x |
| `InternodeLabeler._add_labels_to_standard_tree` order-preserving child push | balanced 131072-tip tree, direct preorder label assignment, optimized helper baseline | 0.027385s | 0.021706s | 1.26x |
| `internode_labeler` module import without eager Bio.Phylo/NumPy | cold subprocess import of internode-labeler command module | 0.156960s | 0.065330s | 2.40x |
| `internode_labeler` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006050s | 0.005131s | 1.18x |
| `internode_labeler` module import without `typing` startup | median cold subprocess import after converting the annotation-only typing alias to a built-in postponed annotation | 0.011143s | 0.007092s | 1.57x |
| `CollapseBranches.run` | supported balanced 32768-tip tree, no branches below threshold | 0.4744s | 0.3769s | 1.3x |
| `CollapseBranches.run` non-JSON count skipping | supported balanced 32768-tip tree, no branches below threshold, output stubbed | 0.2158s | 0.0866s | 2.5x |
| `CollapseBranches.run` no-op collapse pre-scan | supported balanced 32768-tip tree, no branches below threshold, output stubbed | 0.0996s | 0.0133s | 7.5x |
| `CollapseBranches.run` read-only no-op cached tree path | supported balanced 32768-tip cached tree, no branches below threshold, output stubbed | 0.197230s | 0.007733s | 25.50x |
| `CollapseBranches.count_internal_nodes` | balanced 65536-tip tree, JSON collapsed-branch count setup | 0.1224s | 0.0099s | 12.4x |
| `CollapseBranches.run` JSON no-op setup scan | supported balanced 65536-tip tree, no branches below threshold | 0.074100s | 0.046069s | 1.6x |
| `CollapseBranches.run` non-JSON weak-branch pre-scan | weak-root balanced 262144-tip tree, pre-collapse setup only | 0.054011s | 0.000009708s | 5563.49x |
| `collapse_branches` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006097s | 0.004794s | 1.27x |
| `collapse_branches` module import without `typing` startup | median cold subprocess import after converting the annotation-only typing alias to a built-in postponed annotation | 0.007718s | 0.004077s | 1.89x |
| `Spr.run` single-taxon terminal validation and subtree lookup | balanced 65536-tip tree, single terminal subtree | 0.2741s | 0.0171s | 16.0x |
| `Spr._generate_spr_trees` | balanced 64-tip tree, two-tip subtree, 122 generated SPR rearrangements | 0.1695s | 0.0773s | 2.2x |
| `Spr._generate_spr_trees` direct copy traversal | balanced 64-tip tree, two-tip subtree, 122 generated SPR rearrangements | 0.077625s | 0.026085s | 3.0x |
| `Spr._clades_and_parent_map` order-preserving binary push | balanced 65536-tip tree, preorder clade list plus parent map, optimized helper baseline | 0.031078s | 0.020642s | 1.51x |
| `Spr._iter_preorder` binary-child fast path | balanced 131072-tip tree, subtree preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.060582s | 0.038110s | 1.59x |
| `Spr._iter_postorder` reverse-preorder helper | balanced 32768-tip tree, direct postorder helper, side-by-side previous visited-tuple helper | 0.011335s | 0.004941s | 2.29x |
| `Spr._collect_clade_taxa` reverse-preorder and binary union | balanced 32768-tip tree, descendant taxon cache for SPR regraft descriptions | 0.070923s | 0.054266s | 1.31x |
| `Spr._print_output_summary` batched summary output | 100k captured SPR output-file summaries, identical stdout text | 0.086666s | 0.051927s | 1.67x |
| `Spr._print_spr_trees` batched stdout output | 100k generated tree strings, fake writer, captured stdout and identical text | 0.039031s | 0.028162s | 1.39x |
| `spr` module import without eager Bio.Phylo | cold subprocess import after lazy `Phylo.write` and `Clade` proxies | 0.122829s | 0.025707s | 4.78x |
| `spr` module import without eager JSON/pickle helpers | median cold subprocess import after lazy JSON wrapper, removing unused `json`, and lazy pickle proxy | 0.007684s | 0.004994s | 1.54x |
| `spr` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only typing aliases to built-in annotations | 0.007167s | 0.004261s | 1.68x |
| `NearestNeighborInterchange._build_parent_map` | balanced 2048-tip tree, parent map for NNI generation | 1.3251s | 0.0056s | 234.9x |
| `NearestNeighborInterchange._generate_targeted_nnis` terminal-name setup | balanced 65536-tip tree, targeted branch taxon validation | 0.1371s | 0.0143s | 9.6x |
| `IndependentContrasts.run` | balanced 32768-tip tree with one trait value per tip, JSON output captured | 0.9886s | 0.8625s | 1.1x |
| `IndependentContrasts.run` tip-name setup | balanced 32768-tip tree, all tips shared with trait data | 0.0692s | 0.0091s | 7.6x |
| `IndependentContrasts._print_json` direct serialization | 32767 contrasts with descendant-tip lists, stdout captured | 0.160804s | 0.051197s | 3.1x |
| `IndependentContrasts._compute_pic` scalar/label merge traversal | balanced tree with 32768 tips, full descendant-tip labels | 0.094674s | 0.072704s | 1.3x |
| `IndependentContrasts._resolve_polytomies` direct no-op traversal | balanced 32768-tip binary tree, no polytomies to resolve | 0.106671s | 0.009941s | 10.7x |
| `IndependentContrasts.run` direct polytomy traversal | balanced 32768-tip tree with one trait value per tip, JSON output captured | 0.645442s | 0.502745s | 1.3x |
| `IndependentContrasts.run` no-list binary polytomy pass | balanced 32768-tip tree with one trait value per tip, JSON output captured | 0.167298s | 0.154622s | 1.08x |
| `IndependentContrasts._resolve_polytomies` unordered no-op scan | balanced 131072-tip binary tree, no-op dichotomy setup scan, optimized helper baseline | 0.030895s | 0.015197s | 2.03x |
| `IndependentContrasts._print_text` batched output | 200k contrast rows, identical stdout text | 0.160924s | 0.135159s | 1.19x |
| `IndependentContrasts` contrast summary reductions | 100 / 1000 / 10000 / 100000 / 200000 contrast values, side-by-side previous list-based `np.mean(np.abs(...))` and `np.var(..., ddof=1)` wrappers | 0.000019139s / 0.000086641s / 0.000621055s / 0.006865s / 0.013458s | 0.000015530s / 0.000043070s / 0.000297021s / 0.002686s / 0.005619s | 1.23x / 2.01x / 2.09x / 2.56x / 2.39x |
| `IndependentContrasts` text-output compact descendant summaries | pectinate 12000-tip tree, compute PIC and captured text output with identical displayed tip previews/counts | 0.520021s | 0.030507s | 17.05x |
| `IndependentContrasts._postorder_clades_fast` reverse-preorder helper | balanced 32768-tip tree, direct PIC postorder traversal | 0.010762s | 0.003350s | 3.21x |
| `IndependentContrasts.run` trait parsing/prune setup | balanced 32768-tip tree, all tips shared with trait data | 0.030738s | 0.018286s | 1.68x |
| `IndependentContrasts._parse_trait_data` exact-order all-shared fast path | 500k two-column trait rows with comments/blanks, file order matching tree-tip order, identical dict order | 0.599368s | 0.372073s | 1.61x |
| `IndependentContrasts.run` all-shared read-only setup | balanced 32768-tip cached binary tree, one trait value for every tip, PIC/output mocked | 0.281312s | 0.050213s | 5.60x |
| `IndependentContrasts`/`FitDiscrete._needs_default_branch_lengths` unordered scan | balanced 131072-tip tree, complete branch lengths, optimized helper baseline | 0.023833s | 0.016582s | 1.44x |
| `independent_contrasts` module import without eager NumPy | cold subprocess import after lazy NumPy proxy | 0.081544s | 0.025828s | 3.16x |
| `independent_contrasts` module import without eager stdlib JSON | median cold subprocess import after lazy `json.dumps` wrapper | 0.006503s | 0.005449s | 1.19x |
| `independent_contrasts` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only aliases to built-in annotations | 0.006418s | 0.004098s | 1.57x |
| `ParsimonyScore.run` | balanced 32768-tip tree with 4-site FASTA, JSON output captured | 0.7835s | 0.5557s | 1.4x |
| `ParsimonyScore.run` taxon-set setup | balanced 32768-tip tree, all tips shared with alignment | 0.0695s | 0.0091s | 7.6x |
| `ParsimonyScore.run` verbose JSON serialization | 200k per-site scores, stdout captured | 0.054060s | 0.008258s | 6.5x |
| `ParsimonyScore.run` verbose text output | 200k per-site scores, mocked tree/alignment and identical stdout text | 0.050216s | 0.027030s | 1.86x |
| `ParsimonyScore._resolve_polytomies` | balanced 65536-tip binary tree, no-op dichotomy setup scan | 0.1759s | 0.0219s | 8.0x |
| `ParsimonyScore._postorder_clades_fast` reverse-preorder helper | balanced 32768-tip tree, direct parsimony postorder traversal | 0.009456s | 0.003729s | 2.54x |
| `ParsimonyScore._resolve_polytomies` reverse-preorder helper | balanced 32768-tip binary tree, no-op dichotomy setup scan | 0.010886s | 0.005248s | 2.07x |
| `ParsimonyScore._resolve_polytomies` unordered no-op scan | balanced 131072-tip binary tree, no-op dichotomy setup scan, optimized helper baseline | 0.028605s | 0.015738s | 1.82x |
| `ParsimonyScore.run` no-list binary polytomy pass | balanced 32768-tip tree with 4-site FASTA, JSON output captured | 0.088947s | 0.083471s | 1.07x |
| `ParsimonyScore.run` all-shared read-only setup | balanced 32768-tip cached binary tree, four-site FASTA for every tip, Fitch scoring mocked | 0.286847s | 0.057090s | 5.02x |
| `IndependentContrasts`/`ParsimonyScore._has_polytomies` unordered scan | balanced 131072-tip binary tree, no polytomies, optimized helper baseline | 0.023511s | 0.017187s | 1.37x |
| `parsimony_score` module import without eager stdlib JSON | median cold subprocess import after lazy `json.dumps` wrapper | 0.005933s | 0.004855s | 1.22x |
| `parsimony_score` module import without `typing` startup | median cold subprocess import after postponing annotations and replacing annotation-only typing aliases | 0.031268s | 0.028986s | 1.08x |
| `parsimony_utils.build_parent_map` direct traversal | balanced 65536-tip tree, parent map for every non-root clade | 0.190389s | 0.021707s | 8.77x |
| `parsimony_utils.build_parent_map` unordered child push | balanced 131072-tip tree, parent map for every non-root clade, optimized helper baseline | 0.046072s | 0.038417s | 1.20x |
| `parsimony_utils.resolve_polytomies` direct no-op scan | balanced 131072-tip binary tree, no-op dichotomy setup scan, optimized helper baseline | 0.022617s | 0.014873s | 1.52x |
| `parsimony_utils.fitch_downpass` direct postorder | balanced 32768-tip tree, eight discrete characters | 0.367538s | 0.295333s | 1.24x |
| `parsimony_utils.fitch_downpass` binary-node specialization | balanced 32768-tip tree, eight discrete characters | 0.647346s | 0.448144s | 1.44x |
| `parsimony_utils.fitch_downpass` bitmask resolved-tree path | balanced 32768-tip tree, eight binary discrete characters | 0.373862s | 0.162835s | 2.30x |
| `parsimony_utils.fitch_downpass` reverse-preorder postorder helper | balanced 32768-tip tree, eight binary discrete characters | 0.130832s | 0.110949s | 1.18x |
| `parsimony_utils._postorder_clades_direct` localized stack operations | balanced 131072-tip tree, direct postorder helper materialized as a list, side-by-side previous reverse-preorder helper | 0.071677s | 0.052405s | 1.37x |
| `parsimony_utils.fitch_downpass` repeated terminal mask cache | balanced 32768-tip tree, eight binary discrete characters with repeated terminal state vectors | 0.155554s | 0.124489s | 1.25x |
| `parsimony_utils.fitch_uppass_acctran` direct preorder | balanced 32768-tip tree, eight-character state sets after downpass | 0.162745s | 0.077135s | 2.11x |
| `parsimony_utils._preorder_clades_direct` binary-child fast path | balanced 32768-tip tree, direct preorder helper for parsimony uppass/change traversal | 0.013750s | 0.006956s | 1.98x |
| `parsimony_utils.classify_changes` plain-dict transition counts | balanced 32768-tip tree, 24-character synthetic branch-change map with 1.26M changes, identical classifications | 1.033947s | 0.943712s | 1.10x |
| `parsimony_utils.retention_index` large-alphabet ASCII state counts | 1000 characters x 500 taxa / 6000 characters x 1200 taxa / 10000 characters x 256 taxa with 22-46 observed states, side-by-side repeated equality scans | 0.790840s / 1.323010s / 1.573656s | 0.333097s / 0.462331s / 0.455574s | 2.37x / 2.86x / 3.45x |
| `parsimony_utils` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in annotations | 0.022870s | 0.020804s | 1.10x |
| `PhyloLogistic._build_logistic_vcv` | balanced 512-tip tree, OU-transformed VCV at alpha=0.05 | 0.1049s | 0.0970s | 1.1x |
| `discrete_models.felsenstein_pruning` transition matrix cache | balanced 4096-tip tree, two-state ER model, repeated unit branch lengths | 0.0725s | 0.0244s | 3.0x |
| `discrete_models.felsenstein_pruning` direct postorder | balanced 4096-tip tree, two-state ER model, repeated unit branch lengths | 0.023353s | 0.014722s | 1.59x |
| `discrete_models._prepare_felsenstein_context` direct postorder | balanced 4096-tip tree, two-state ER model, prepared pruning metadata | 0.017437s | 0.008280s | 2.1x |
| `discrete_models.felsenstein_pruning` reverse-preorder postorder helper | balanced 4096-tip tree, two-state ER model, repeated unit branch lengths | 0.017259s | 0.015753s | 1.10x |
| `discrete_models._prepare_felsenstein_context` reverse-preorder postorder helper | balanced 4096-tip tree, two-state ER model, prepared pruning metadata | 0.006522s | 0.004503s | 1.45x |
| `discrete_models._felsenstein_loglik_prepared` two-state ARD scalar transitions | balanced 512-tip tree, binary unequal-rate Q matrix, prepared pruning metadata | 0.001701s | 0.000725s | 2.35x |
| `discrete_models.build_q_matrix` direct ER diagonal | two-state ER Q matrix, repeated optimizer-style calls | 0.00000250s | 0.00000186s | 1.34x |
| `discrete_models.build_q_matrix` cached SYM off-diagonal indices | 32-state SYM Q matrix, repeated optimizer-style calls | 0.00014780s | 0.00000732s | 20.19x |
| `discrete_models.build_q_matrix` cached ARD off-diagonal mask | 32-state ARD Q matrix, repeated optimizer-style calls | 0.00016262s | 0.00000408s | 39.83x |
| `discrete_models.fit_q_matrix` prepared pruning context | balanced 512-tip tree, two-state ER model, full Q fit | 1.674459s | 0.802447s | 2.1x |
| `discrete_models.fit_q_matrix` cached internal pruning rows | balanced 512-tip tree, two-state ER model, full Q fit | 0.811803s | 0.740218s | 1.1x |
| `discrete_models.fit_q_matrix` two-state ER analytic transition matrices | balanced 512-tip tree, two-state ER model, full Q fit | 0.824497s | 0.689012s | 1.20x |
| `discrete_models.fit_q_matrix` scalar two-state ER prepared pruning | balanced 512-tip tree, two-state ER model, full Q fit | 0.707479s | 0.324087s | 2.18x |
| `discrete_models.fit_q_matrix` scalar two-state ER rate objective | balanced 512-tip tree, two-state ER model, full Q fit | 0.353477s | 0.321415s | 1.10x |
| `discrete_models.fit_q_matrix` two-state ER scalar optimizer | balanced 512-tip tree, two-state ER model, full Q fit with equal log-likelihood | 0.296045s | 0.022109s | 13.39x |
| `discrete_models` two-state scalar exp/log primitives | one scalar transition decay and one scalar root log-likelihood operation, side-by-side previous NumPy ufunc dispatch | 0.000000403s | 0.000000145s | 2.79x |
| `discrete_models` generic root likelihood total | 2 / 3 / 4 / 8 / 16 / 64-state prior-weighted likelihood vectors, side-by-side previous `np.sum(pi * root_lik)` | 0.000006310s / 0.000006381s / 0.000006961s / 0.000006135s / 0.000005135s / 0.000005092s | 0.000001144s / 0.000001084s / 0.000001175s / 0.000001519s / 0.000001047s / 0.000001155s | 5.52x / 5.89x / 5.93x / 4.04x / 4.90x / 4.41x |
| `FitDiscrete.run` shared pruning context | balanced 8192-tip tree, ER/SYM/ARD setup context reuse | 0.045272s | 0.014442s | 3.13x |
| `FitDiscrete.run` all-shared read-only setup | balanced 32768-tip cached tree, trait state for every tip, model fitting/output mocked | 0.319352s | 0.102667s | 3.11x |
| `FitDiscrete._print_text` batched model table | captured model comparison table with 100k synthetic rows, identical stdout text | 0.182404s | 0.168677s | 1.08x |
| `fit_discrete` module import without eager NumPy/discrete helper | cold subprocess import after lazy helper wrappers and local model constant | 0.090685s | 0.026126s | 3.47x |
| `fit_discrete` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006547s | 0.005270s | 1.24x |
| `fit_discrete` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only aliases to built-in annotations | 0.006143s | 0.003948s | 1.56x |
| `discrete_models` module import without eager SciPy linalg/optimize | cold process import for shared discrete-model utilities | 0.427970s | 0.121082s | 3.5x |
| `discrete_models` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.121082s | 0.019515s | 6.21x |
| `discrete_models` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.022798s | 0.020829s | 1.09x |
| `discrete_models.parse_discrete_traits` two-column streaming parser | 300k-row two-column discrete trait TSV, 100k shared taxa | 0.306287s | 0.274550s | 1.12x |
| `discrete_models.parse_discrete_traits` multi-column streaming parser | 300k-row multi-column discrete trait TSV, 300k tree tips | 0.210558s | 0.201712s | 1.04x |
| `discrete_models.parse_discrete_traits` two-column all-shared parser fast path | 300k-row two-column discrete trait TSV, all taxa shared | 0.229224s | 0.129126s | 1.78x |
| `discrete_models.parse_discrete_traits` two-column partition parser | 300k-row two-column discrete trait TSV with comments/blanks, all taxa shared, side-by-side previous split parser comparison | 0.131254s | 0.127151s | 1.03x |
| `discrete_models.parse_discrete_traits` multi-column all-shared parser fast path | 300k-row multi-column discrete trait TSV, all taxa shared | 0.244142s | 0.135837s | 1.80x |
| `discrete_models.parse_discrete_traits` multi-column bounded selected-column split | 300k-row x 24-trait discrete TSV, all taxa shared, selected early trait column, side-by-side previous full row split | 0.284023s | 0.214614s | 1.32x |
| `discrete_models.matrix_exp` cached SciPy expm wrapper | 20k four-state ARD transition matrices, SciPy already warm | 0.136597s | 0.125134s | 1.09x |
| `discrete_models.matrix_exp` two-state ER analytic path | 20k two-state ER transition matrices | 0.105901s | 0.034737s | 3.05x |
| `discrete_models.matrix_exp` two-state ARD analytic path | 20k binary unequal-rate transition matrices, side-by-side previous SciPy `expm` path | 0.125131s | 0.042341s | 2.96x |
| `Tree.validate_tree` required branch lengths | balanced 65536-tip tree, min-tip and branch-length validation | 0.2985s | 0.0184s | 16.2x |
| `Tree.validate_tree` default branch lengths | balanced 65536-tip tree, fill half of terminal branch lengths | 0.2983s | 0.0149s | 20.0x |
| `Tree._validate_standard_tree` unordered child push, required branch lengths | balanced 262144-tip tree, min-tip and branch-length validation, optimized helper baseline | 0.053449s | 0.031135s | 1.72x |
| `Tree._validate_standard_tree` unordered child push, default branch lengths | balanced 65536-tip tree, fill half of terminal branch lengths, optimized helper baseline | 0.016180s | 0.010987s | 1.47x |
| `Tree.validate_tree` generic fallback min-tip check | nonstandard tree yielding 500k terminals, `min_tips=3` | 0.036222s | 0.000001292s | 28035.6x |
| `Tree.prune_tree_using_taxa_list` | balanced 2048-tip tree, prune 1024 named tips | 0.2327s | 0.2070s | 1.1x |
| `Tree.prune_tree_using_taxa_list` batch standard-tree pruning | balanced 2048-tip tree, prune 1024 resolved terminal targets | 0.2030s | 0.0022s | 91.9x |
| `Tree._prune_terminal_objects_batch_standard_tree` reverse-preorder pass | balanced 2048-tip tree, prune 1024 resolved terminal targets | 0.001459s | 0.001125s | 1.30x |
| `Tree.prune_tree_using_taxa_list` terminal target setup | balanced 8192-tip tree, resolve 512 named prune targets | 0.0283s | 0.0040s | 7.0x |
| `Tree.prune_tree_using_taxa_list` selected target setup | balanced 65536-tip tree, resolve 512 named prune targets | 0.024874s | 0.018461s | 1.35x |
| `Tree._terminal_targets_and_count_fast` unordered child push | balanced 131072-tip tree, resolve 512 named prune targets, optimized helper baseline | 0.024863s | 0.014603s | 1.70x |
| `Tree._terminal_by_name_fast` unordered child push | balanced 131072-tip tree, full terminal-name index for prune fallback, optimized helper baseline | 0.032400s | 0.021422s | 1.51x |
| `ThresholdModel._prune_tree_to_taxa` | balanced 2048-tip tree, prune 1024 tips before VCV/MCMC setup | 3.1806s | 0.4040s | 7.9x |
| `ThresholdModel._prune_tree_to_taxa` batch keep-pruning | balanced 2048-tip tree, keep 1024 tips before VCV/MCMC setup | 0.417008s | 0.002112s | 197.5x |
| `TreeSpace._build_distance_matrix` prune step | 20 balanced 512-tip trees, prune to 256 shared taxa before RF splits | 0.7711s | 0.6363s | 1.2x |
| `ConsensusNetwork._prune_to_taxa` | balanced 2048-tip tree, prune 1024 tips to shared taxa | 0.5093s | 0.4056s | 1.3x |
| `ConsensusNetwork._prune_to_taxa` batch standard-tree pruning | balanced 2048-tip tree, prune 1024 tips to shared taxa | 0.2286s | 0.0071s | 32.3x |
| `QuartetNetwork._prune_to_taxa` | balanced 2048-tip tree, prune 1024 tips to shared taxa | 0.5070s | 0.4073s | 1.2x |
| `build_vcv_matrix` descendant-index accumulation | balanced 2048-tip tree, all terminal names ordered | 0.056182s | 0.048285s | 1.16x |
| `vcv_utils.build_vcv_matrix` broadcast block indexing | balanced 2048-tip tree, descendant-index VCV branch accumulation | 0.041264s | 0.037537s | 1.10x |
| `vcv_utils.build_vcv_matrix` single-tip diagonal updates | balanced 2048-tip tree, descendant-index VCV branch accumulation | 0.151061s | 0.054459s | 2.77x |
| `build_discordance_vcv` shared-taxa setup | 80 balanced gene trees x 512 taxa plus species tree | 0.0810s | 0.0093s | 8.7x |
| `build_discordance_vcv` prune step | 10 balanced 1024-tip gene trees, prune to 512 shared taxa before averaging VCVs | 2.7773s | 2.1955s | 1.3x |
| `build_discordance_vcv` copied gene-tree batch pruning | 10 balanced 1024-tip gene trees, prune copies to 512 shared taxa | 1.0501s | 0.0440s | 23.9x |
| `build_discordance_vcv` direct pruned-subset VCV | 10 balanced 1024-tip gene trees, 512 shared taxa before averaging VCVs | 0.118993s | 0.086442s | 1.38x |
| `build_discordance_vcv` pruned-subset single-tip diagonal updates | balanced 4096-tip gene tree, 2048 retained taxa | 0.069937s | 0.051761s | 1.35x |
| `build_discordance_vcv` pruned-subset broadcast block indexing | 10 balanced 1024-tip gene trees, 512 shared taxa before averaging VCVs | 0.059726s | 0.052204s | 1.14x |
| `build_discordance_vcv` gene-tree branch-length validation | 20 balanced 4096-tip gene trees | 0.2277s | 0.0156s | 14.6x |
| `build_discordance_vcv` combined gene-tree tip scan and branch validation | 20 balanced 4096-tip gene trees | 0.035306s | 0.018838s | 1.87x |
| `build_discordance_vcv` PSD fast path | 10 balanced 1024-tip gene trees, all taxa shared before VCV averaging | 0.376386s | 0.301010s | 1.25x |
| `build_discordance_vcv` all-shared gene-tree copy skip | 10 balanced 512-tip gene trees, all taxa shared before VCV averaging | 0.124170s | 0.096836s | 1.3x |
| `vcv_utils._copy_prune_gene_tree_to_shared_taxa` all-shared copy skip | balanced 32768-tip gene tree, all taxa shared | 0.425988s | 0.008369s | 50.9x |
| `vcv_utils`/`SpectralDiscordance` all-shared prune preflight scan | balanced 131072-tip standard tree, all taxa shared, optimized direct-scan baseline | 0.034115s | 0.021787s | 1.57x |
| `vcv_utils.parse_gene_trees` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.096195s | 0.072211s | 1.33x |
| `vcv_utils.parse_gene_trees` path-list resolver | 50k relative tree path rows, tree parsing mocked | 0.081335s | 0.014871s | 5.47x |
| `vcv_utils._gene_tree_has_missing_branch_lengths` unordered scan | balanced 131072-tip gene tree, complete branch lengths, optimized helper baseline | 0.023896s | 0.016475s | 1.45x |
| `vcv_utils._get_tip_names_and_missing_branch_lengths` order-preserving scan | balanced 131072-tip gene tree, terminal names plus branch-length validation, optimized helper baseline | 0.037289s | 0.029137s | 1.28x |
| `vcv_utils._nearest_psd` smallest-eigenvalue check | 1024 x 1024 already-PSD averaged VCV matrix | 0.106880s | 0.061465s | 1.74x |
| `vcv_utils._nearest_psd` cached SciPy eigvalsh wrapper | 2k repeated 8 x 8 already-PSD matrices, SciPy already warm | 0.027812s | 0.026592s | 1.05x |
| `vcv_utils._nearest_psd` correction reconstruction | 900 x 900 eigenvector matrix with clipped eigenvalues | 0.029000s | 0.015881s | 1.8x |
| `vcv_utils._nearest_psd_from_eigendecomposition` eigenvalue minimum | 8 / 64 / 1024 / 100k / 1M eigenvalues, side-by-side previous `np.min` wrapper | 2.048513s / 0.970033s / 0.288490s / 0.105668s / 0.029085s | 1.375394s / 0.674878s / 0.170990s / 0.095815s / 0.025196s | 1.49x / 1.44x / 1.69x / 1.10x / 1.15x |
| `vcv_utils` module import without eager SciPy linalg | cold process import for shared VCV helpers | 0.323420s | 0.183398s | 1.8x |
| `vcv_utils` module import without eager Bio.Phylo | cold subprocess import with lazy `Phylo.read` proxy | 0.174747s | 0.112800s | 1.55x |
| `vcv_utils` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.112800s | 0.024335s | 4.64x |
| `vcv_utils` module import without eager pickle | median cold subprocess import after lazy pickle proxy | 0.006726s | 0.005047s | 1.33x |
| `vcv_utils` module import without typing/path parser startup | median cold subprocess import after converting annotation-only typing names and localizing `Path`/`StringIO` to gene-tree parsing | 0.013944s | 0.004052s | 3.44x |
| `SpectralDiscordance._copy_prune_if_needed` | balanced 2048-tip tree, prune 1024 tips before bipartition extraction | 0.5258s | 0.4127s | 1.3x |
| `SpectralDiscordance._copy_prune_if_needed` direct terminal scans | balanced 2048-tip tree, prune 1024 tips before bipartition extraction | 0.0142s | 0.0076s | 1.9x |
| `SpectralDiscordance._copy_prune_if_needed` no-prune direct scan | balanced 2048-tip tree, all taxa shared | 0.00365s | 0.00047s | 7.8x |
| `SpectralDiscordance._copy_prune_if_needed` copied-tree target scan child push | balanced 262144-tip tree, collect 131072 prune targets, side-by-side previous `reversed(children)` scan | 0.082522s | 0.071608s | 1.15x |
| `PrintTree.run` with `--remove --json` | balanced 8192-tip tree, branch lengths removed before Newick formatting | 0.1028s | 0.0963s | 1.1x |
| `PrintTree.run` branch-length removal traversal | balanced 65536-tip tree, clear every branch length | 0.1727s | 0.0077s | 22.5x |
| `PrintTree._remove_standard_tree_branch_lengths` one-pass traversal | balanced 65536-tip tree, clear every branch length | 0.014211s | 0.012182s | 1.17x |
| `PrintTree._remove_standard_tree_branch_lengths` localized stack operations | balanced 131072-tip tree, clear every branch length, side-by-side previous direct helper median | 0.031891s | 0.026640s | 1.20x |
| `print_tree` module import without eager Bio.Phylo | cold subprocess import after lazy `Phylo.draw_ascii` proxy | 0.130837s | 0.025049s | 5.22x |
| `PrintTree.run` no-remove cached read-only tree setup | balanced 32768-tip cached tree, ASCII rendering mocked | 0.389771s | 0.000138s | 2824.43x |
| `print_tree` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006284s | 0.005009s | 1.25x |
| `print_tree` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only aliases to built-in annotations | 0.006847s | 0.004296s | 1.59x |
| `ConsensusTree._prune_to_taxa` | balanced 1024-tip tree, prune 512 tips to shared taxa | 0.1857s | 0.1485s | 1.3x |
| `ConsensusTree._tips` / `ConsensusNetwork._tips` | 40 balanced 4096-tip trees, taxon-set extraction | 0.3272s | 0.0356s | 9.2x |
| `ConsensusTree` / `ConsensusNetwork` / `QuartetNetwork` `_normalize_taxa` identical taxon sets | 80 precomputed taxon sets x 1024 shared taxa, defer shared intersection until needed | 2.673776s | 0.475789s | 5.62x |
| `ConsensusTree` / `ConsensusNetwork` / `QuartetNetwork` `_normalize_taxa` identical tip-set no-slice scan | 1M precomputed tip sets, identical / early-different / late-different cases, side-by-side previous `tip_sets[1:]` shortcut predicate | 0.267546s / 0.006743s / 0.105294s | 0.226317s / 0.000000208s / 0.077291s | 1.18x / 32392.41x / 1.36x |
| `ConsensusTree._parse_trees_from_source` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.090417s | 0.067520s | 1.34x |
| `ConsensusTree._parse_trees_from_source` path-list resolver | 50k existing relative tree paths, tree parsing mocked | 0.720531s | 0.507855s | 1.42x |
| `consensus_tree` module import without eager Bio.Phylo | cold subprocess import after lazy `Phylo.read` and `Consensus` proxies | 0.129701s | 0.025590s | 5.07x |
| `consensus_tree` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006000s | 0.004952s | 1.21x |
| `consensus_tree` module import without `typing` startup | median cold subprocess import after postponing annotations and converting annotation-only aliases to built-in annotations | 0.005166s | 0.003532s | 1.46x |
| `ConsensusNetwork._count_splits` | 120 balanced trees x 256 taxa, rooted split counting | 0.5522s | 0.1808s | 3.1x |
| `ConsensusNetwork._count_splits` mask counter | 120 balanced trees x 256 taxa, rooted split counting | 0.208316s | 0.042667s | 4.9x |
| `ConsensusNetwork._count_splits` reverse-preorder mask extractor | 80 balanced trees x 1024 taxa, rooted split counting, side-by-side previous visited-tuple mask extractor | 0.111707s | 0.094082s | 1.19x |
| `ConsensusNetwork._count_splits` batched split-set Counter updates | 600 extracted split sets x 700 splits sampled from 50k possible splits, identical split and possible counters | 0.113693s | 0.037813s | 3.01x |
| `ConsensusNetwork._canonical_split` equal-size tiebreak | 3k equal-size 600-vs-600 bipartitions over 1200 taxa, identical sorted-lexicographic canonical side | 1.113460s | 0.330384s | 3.37x |
| `ConsensusNetwork._canonical_split_mask` equal-size tiebreak | 3k equal-size 600-vs-600 bitmask splits over 1200 sorted taxa, identical lexicographic canonical side | 0.509757s | 0.000929s | 548.91x |
| `ConsensusNetwork._count_splits` allow-mode taxon setup | 80 balanced trees x 1024 taxa with one varying taxon | 0.1450s | 0.0173s | 8.4x |
| `ConsensusNetwork._compute_circular_ordering` terminal-order extraction | parsed balanced 65536-tip consensus-like tree, collect terminal names | 0.139837s | 0.018060s | 7.74x |
| `ConsensusNetwork._build_splits_graph` | 80 taxa, 20 circular splits, 33-node consensus graph | 7.2405s | 0.0117s | 621.0x |
| `ConsensusNetwork._build_splits_graph` integer-mask edge pass | 4096 taxa, 12 independent splits, 4096 nodes / 24576 edges, side-by-side previous tuple-slicing edge pass | 0.078513s | 0.062518s | 1.26x |
| `ConsensusNetwork._compute_split_directions` cached taxon coordinates | 50 taxa/200 splits, 200 taxa/1000 splits, 500 taxa/3000 splits, identical direction dictionaries with cached gap positions | 0.000906s / 0.013686s / 0.127458s | 0.000681s / 0.008597s / 0.102900s | 1.33x / 1.59x / 1.24x |
| `ConsensusNetwork._compute_split_directions` one-pass split centers | 500 circular splits over 2000 taxa with cached gap positions, side-by-side previous two generator sums | 0.059449s | 0.041242s | 1.44x |
| `NeighborNet` / `ConsensusNetwork` network edge rendering | 80 taxa, 20 circular splits, real Matplotlib Agg internal and pendant edge render | 0.025885s | 0.012323s | 2.10x |
| `NeighborNet` / `ConsensusNetwork` unlabeled fallback point rendering | 4096 taxa without accepted splits, real Matplotlib Agg point render | 0.756119s | 0.026598s | 28.43x |
| `ConsensusNetwork.run` batched text split output | 100k filtered split rows, captured stdout and identical text | 0.064258s | 0.050664s | 1.27x |
| `ConsensusNetwork._parse_trees_from_source` source cleanup | 500k path-like rows with comments/blanks, cleanup before tree parsing | 0.091090s | 0.067767s | 1.34x |
| `ConsensusNetwork._parse_trees_from_source` path-list resolver | 50k existing relative tree paths, tree parsing mocked | 0.720531s | 0.507855s | 1.42x |
| `consensus_network` module import without eager Bio.Phylo | cold subprocess import after lazy `Phylo.read` and `Consensus.majority_consensus` proxies | 0.121535s | 0.032077s | 3.79x |
| `consensus_network` module import without eager JSON/plot helpers | median cold subprocess import after localizing PlotConfig and lazy JSON wrapper | 0.023345s | 0.006703s | 3.48x |
| `consensus_network` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006103s | 0.004336s | 1.41x |
| `LBScore.calculate_lb_score` | balanced tree with 160 tips | 5.3248s | 0.0474s | 112.3x |
| `LBScore.calculate_lb_score` shared pairwise cache | balanced 512-tip tree | 0.2313s | 0.1520s | 1.5x |
| `LBScore.calculate_lb_score` direct per-tip distance sums | balanced 512-tip tree, shared pairwise cache | 0.1777s | 0.0707s | 2.5x |
| `LBScore.calculate_lb_score` linear tree distance sums | balanced 512-tip tree | 0.054681s | 0.005159s | 10.6x |
| `LBScore._calculate_lb_components_fast` historical denominator arithmetic | balanced 8192-tip tree, linear component helper, side-by-side previous `len(tip_set - set(tip))` denominator | 1.071748s | 0.015581s | 68.79x |
| `LBScore.calculate_average_distance_of_taxon_to_other_taxa` fallback tip-set reuse | 1000 fallback taxon names, side-by-side previous per-tip `set(tips)` construction while preserving `set(tip)` behavior | 0.024459s | 0.009492s | 2.58x |
| `LBScore.calculate_average_distance_between_tips` fallback streaming pair batches | 2000 fallback tips, 1,999,000 pair batch setup, side-by-side previous full pair-list slicing | 0.484847s | 0.317065s | 1.53x |
| `LBScore._calculate_lb_components_fast` postorder child push | balanced 32768-tip tree, linear component helper, side-by-side previous `reversed(children)` setup | 0.159546s | 0.132146s | 1.21x |
| `LBScore.calculate_lb_score_per_taxa` without NumPy startup | cold subprocess, 32768 average-distance values transformed to LB scores | 0.098413s | 0.029658s | 3.32x |
| `LBScore.run` verbose text output | 200k taxon LB-score rows, mocked tree/read and identical stdout text | 0.119431s | 0.093329s | 1.28x |
| `LBScore.run` cached read-only tree setup | balanced 32768-tip cached tree, LB calculation and output mocked | 0.381093s | 0.000095s | 4002.74x |
| `lb_score` module import without eager tqdm | cold subprocess import of `phykit.services.tree.lb_score` | 0.211205s | 0.182790s | 1.16x |
| `lb_score` module import without annotation-only Bio.Phylo | cold subprocess import after postponed annotations and removing `Newick` import | 0.147666s | 0.046242s | 3.19x |
| `lb_score` module import without eager multiprocessing/concurrent futures | cold subprocess import after lazy executor and cpu-count proxies | 0.045831s | 0.030564s | 1.50x |
| `lb_score` module import without eager pickle/stats/json helpers | cold subprocess import after lazy `pickle` proxy and local forwarding wrappers | 0.030135s | 0.022288s | 1.35x |
| `lb_score` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005367s | 0.003556s | 1.51x |
| `Saturation.loop_through_combos_and_calculate_pds_and_pis` | balanced tree with 150 tips x 800 sites, `exclude_gaps=True` | 1.5329s | 0.1745s | 8.8x |
| `Saturation.loop_through_combos_and_calculate_pds_and_pis` byte sequence arrays | balanced tree with 150 tips x 800 sites, `exclude_gaps=True` | 0.1039s | 0.0957s | 1.1x |
| `Saturation._gap_mask_for_array` byte lookup | 1200 taxa x 12000 sites, alphabet `ACGT-?NX*` | 0.043702s | 0.017843s | 2.45x |
| `Saturation.loop_through_combos_and_calculate_pds_and_pis` uint8 valid-position counts | balanced tree with 150 tips x 800 sites, `exclude_gaps=True` | 0.1417s | 0.0304s | 4.7x |
| `Saturation.loop_through_combos_and_calculate_pds_and_pis` ordered cached distances | balanced tree with 150 tips x 800 sites, `exclude_gaps=True` | 0.0340s | 0.0269s | 1.26x |
| `Saturation.loop_through_combos_and_calculate_pds_and_pis` all-pairs matrix distances | balanced tree with 150 tips x 800 sites, `exclude_gaps=True`, ordered cached patristic distances | 0.022602s | 0.010371s | 2.2x |
| `Saturation._calculate_uncorrected_distances_matrix` no-gap direct comparison | 300 taxa x 1200 sites, 44,850 requested pairs, `exclude_gaps=False` | 0.117667s | 0.031655s | 3.72x |
| `Saturation._calculate_uncorrected_distances_matrix` clean `exclude_gaps` handoff | 300 taxa x 1200 clean DNA sites, 44,850 requested pairs, `exclude_gaps=True`, side-by-side previous gappy valid-matrix path | 8.094869s | 0.888346s | 9.11x |
| `Saturation._calculate_uncorrected_distances_matrix` gappy valid-mask construction | 150 taxa x 800 sites / 300 taxa x 1200 sites / 1000 taxa x 1200 sites / 1000 taxa x 8000 sites, side-by-side previous per-tip inverted-mask list before `vstack` | 1.230092s / 0.782973s / 0.323495s / 0.138809s | 0.636202s / 0.347071s / 0.080768s / 0.080888s | 1.93x / 2.26x / 4.01x / 1.72x |
| `Saturation._standard_upper_triangle_gappy_distances` all-valid row division | synthetic gappy distance matrices for 100 / 250 / 500 / 1000 taxa, side-by-side previous per-row valid-mask scratch arrays | 0.000375s / 0.002955s / 0.011395s / 0.079078s | 0.000170s / 0.001857s / 0.004250s / 0.070790s | 2.20x / 1.59x / 2.68x / 1.12x |
| `Saturation._standard_upper_triangle_gappy_distances` zero-length pairs | synthetic gappy distance matrices for 100 / 250 / 500 / 1000 taxa with 5% zero adjusted lengths, side-by-side previous per-row valid-mask scratch arrays | 0.000634s / 0.017047s / 0.016864s / 0.127760s | 0.000369s / 0.000873s / 0.005999s / 0.022578s | 1.72x / 19.53x / 2.81x / 5.66x |
| `Saturation._calculate_uncorrected_distances_matrix` standard pair-order extraction | 1000 taxa, 499500 requested pairs in upper-triangle order, identical extracted matrix distances | 0.265385s | 0.093287s | 2.84x |
| `Saturation._standard_upper_triangle_values` trusted run-order extraction | 1000 taxa, 499500 internally generated upper-triangle pairs, identical extracted matrix distances | 0.134489s | 0.013554s | 9.92x |
| `Saturation.loop_through_combos_and_calculate_pds_and_pis` identical-sequence shortcut | balanced 300-tip tree x 1200 mixed-symbol sites, 44,850 requested pairs, `exclude_gaps=True`, side-by-side previous matrix distance path | 1.866907s | 0.140145s | 13.32x |
| `Saturation` identical Unicode valid-site count | 100k-site uppercase Unicode identical-sequence helper, DNA / protein gap sets, side-by-side previous Python character-membership loop | 0.006550s / 0.007644s | 0.000488s / 0.000397s | 13.42x / 19.24x |
| `Saturation.loop_through_combos_and_calculate_pds_and_pis` identical-sequence no-slice scan | 1M uppercase sequence strings, identical / early-different / late-different cases, side-by-side previous `sequences[1:]` shortcut predicate | 0.058916s / 0.004838s / 0.051696s | 0.033481s / 0.000000334s / 0.045077s | 1.76x / 14490.77x / 1.15x |
| `Saturation._constant_uncorrected_distance_for_identical_sequences` raw-identical normalization scan | 500k raw-identical requested tips with gappy DNA symbols, side-by-side previous eager uppercase tip dictionary and sequence list | 0.824100s | 0.671471s | 1.23x |
| `Saturation.run` cached read-only tree setup | balanced 32768-tip cached tree, alignment parsing, pairwise calculation, and output mocked | 0.346762s | 0.000101s | 3426.23x |
| `Saturation.print_res` verbose text output | 200k pairwise rows, captured stdout and identical text | 0.205054s | 0.166528s | 1.23x |
| `saturation` module import without eager NumPy/Bio.Phylo | cold subprocess import after lazy NumPy and annotation-only Biopython imports | 0.144697s | 0.035127s | 4.12x |
| `saturation` module import without eager multiprocessing | cold subprocess import after lazy multiprocessing proxy and localized `partial` import | 0.035324s | 0.031814s | 1.11x |
| `saturation` module import without eager plot config | cold subprocess import after localizing `PlotConfig` to argument processing | 0.029371s | 0.023684s | 1.24x |
| `saturation` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006523s | 0.005264s | 1.24x |
| `saturation` module import without `typing` startup | median cold subprocess import after removing runtime `TYPE_CHECKING` and converting annotation-only typing aliases to built-in annotations | 0.006053s | 0.003971s | 1.52x |
| `saturation` module import without eager file helper | median cold subprocess import after lazy alignment-reader wrapper | 0.051539s | 0.025990s | 1.98x |
| `vcv_utils.build_vcv_matrix` | balanced tree with 220 tips | 3.8985s | 0.0186s | 209.3x |
| `Tree.calculate_pairwise_tip_distances_fast` path setup | balanced 1024-tip tree, pairwise distances among 768 named tips | 0.2270s | 0.1362s | 1.7x |
| `Tree.calculate_pairwise_tip_distances_fast` prefix LCA lookup | balanced 1024-tip tree, pairwise distances among 768 named tips | 0.132881s | 0.083547s | 1.59x |
| `Tree.calculate_pairwise_tip_distances_fast` binary child push setup | balanced 65536-tip tree, pairwise distances among 64 named tips, optimized helper baseline | 0.028739s | 0.027155s | 1.06x |
| `vcv_utils.build_vcv_matrix` path setup | balanced 1024-tip tree, standard VCV for all tips | 0.3823s | 0.2230s | 1.7x |
| `vcv_utils.build_transformed_vcv_matrix` path setup | balanced 1024-tip tree, transformed VCV for all tips | 0.3930s | 0.2237s | 1.8x |
| `vcv_utils.build_vcv_matrix` branch accumulation | balanced 1024-tip tree, standard VCV for all tips | 0.220607s | 0.018267s | 12.1x |
| `vcv_utils.build_transformed_vcv_matrix` branch accumulation | balanced 1024-tip tree, transformed VCV for all tips | 0.219507s | 0.018513s | 11.9x |
| `tree_paths.build_object_parent_map` | balanced 65536-tip tree, object parent map for VCV/path helpers | 0.1865s | 0.0177s | 10.6x |
| `tree_paths.build_object_parent_map` unordered child push | balanced 131072-tip tree, object parent map for VCV/path helpers, optimized helper baseline | 0.042188s | 0.030392s | 1.39x |
| `tree_paths.build_root_path_map` | balanced 32768-tip tree, root-to-node path lists for all clades | 0.146555s | 0.050841s | 2.9x |
| `LTT._compute_gamma` + `LTT._compute_ltt` | balanced tree with 500 tips | 0.3363s | 0.0077s | 43.8x |
| `Chronogram._compute_root_to_tip` | balanced tree with 5000 tips | 9.8673s | 0.0041s | 2417.5x |
| `Chronogram._print_json` payload construction | balanced tree with 2048 tips, printing stubbed | 0.0577s | 0.0214s | 2.7x |
| `Chronogram._print_json` direct traversal payload construction | balanced tree with 32768 tips, printing stubbed | 0.227351s | 0.072172s | 3.15x |
| `Chronogram._build_descendant_tip_name_cache` exact direct postorder helper | balanced 4096-tip tree, sorted descendant tip names for JSON node ages | 0.005863s | 0.003530s | 1.66x |
| `Chronogram._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.046541s | 0.033467s | 1.39x |
| `Chronogram._preorder_clades_direct` binary-child fast path | balanced 131072-tip tree, direct preorder list helper, side-by-side previous `reversed(children)` helper | 0.073212s | 0.048706s | 1.50x |
| `Chronogram._plot_rectangular` setup | balanced 32768-tip tree, precomputed parent map and root distances | 0.6582s | 0.1047s | 6.3x |
| `Chronogram._plot_rectangular` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, y-coordinate setup | 0.049343s | 0.039908s | 1.24x |
| `Chronogram._plot_rectangular` batched base branches | balanced 2048-tip tree, real Matplotlib Agg branch/label/timescale render | 1.938913s | 0.539964s | 3.59x |
| `Chronogram._plot_rectangular` HPD interval rendering | 4096 rectangular HPD bars, real Matplotlib Agg bar render | 1.311417s | 0.147596s | 8.89x |
| `Chronogram._plot_rectangular` clade-color overlay rendering | 4096 highlighted rectangular clade branches, real Matplotlib Agg line render | 1.466798s | 0.060810s | 24.12x |
| `Chronogram._plot_circular` setup | balanced 32768-tip tree, precomputed parent map and root distances | 0.6870s | 0.1836s | 3.7x |
| `Chronogram._plot_circular` circular coordinate clade-list reuse | balanced 32768-tip tree, root distances, parent map, preorder list, and tips already available | 0.059658s | 0.046128s | 1.29x |
| `Chronogram._plot_circular` batched base branches/arcs | balanced 2048-tip tree, real Matplotlib Agg branch/arc/label/timescale render | 1.627218s | 0.569001s | 2.86x |
| `Chronogram._plot_circular` HPD interval rendering | 4096 circular HPD radial intervals, real Matplotlib Agg line render | 0.700745s | 0.044718s | 15.67x |
| `Chronogram._parse_hpd_intervals` direct traversal and cached regexes | balanced 32768-tip tree with HPD comments on two-thirds of internal nodes | 0.121682s | 0.034318s | 3.55x |
| `chronogram` module import without eager NumPy | cold subprocess import after lazy NumPy proxy | 0.079933s | 0.032269s | 2.48x |
| `chronogram` module import without eager helper modules | median cold subprocess import after lazy wrappers/localized imports for JSON, plot, timescale, circular, and color helpers | 0.013520s | 0.005307s | 2.55x |
| `chronogram` module import without typing startup | median cold subprocess import after removing annotation-only typing imports from tree base and chronogram module | 0.006369s | 0.003916s | 1.63x |
| `Chronogram.run` cached read-only tree setup | balanced 32768-tip cached tree, validation, parent map, root distances, HPD parsing, and plotting mocked | 0.103109s | 0.000023s | 4483.0x |
| `DVMC.determine_dvmc` | balanced tree with 2500 tips | 1.0085s | 0.0117s | 85.9x |
| `PhyloHeatmap.run` initial tip-name setup | balanced 32768-tip tree, all tips shared with data | 0.0667s | 0.0069s | 9.7x |
| `PhyloHeatmap.run` tip-order setup | balanced 32768-tip tree after pruning/ladderize step | 0.0661s | 0.0060s | 10.9x |
| `PhyloHeatmap.run` all-shared read-only setup | balanced 32768-tip cached tree, two numeric traits for every tip, plotting/output mocked | 0.301831s | 0.063210s | 4.78x |
| `PhyloHeatmap._parse_trait_matrix` streaming parser | 200k taxa x 8 numeric traits, comments/blanks before header, all taxa shared | 0.729967s | 0.701524s | 1.04x |
| `PhyloHeatmap._parse_trait_matrix` all-shared parser fast path | 200k taxa x 8 numeric traits, comments/blanks before header, all taxa shared | 0.359596s | 0.247998s | 1.45x |
| `PhyloHeatmap._parse_trait_matrix` numeric hot-loop cleanup | 200k taxa x 8 numeric traits, comments/blanks, all taxa shared, side-by-side previous parser comparison | 0.285750s | 0.256122s | 1.12x |
| `PhyloHeatmap._build_heatmap_matrix` typed matrix construction | 200k taxa x 8 numeric traits, ordered row lookup, side-by-side previous inline `np.array` construction | 0.173474s | 0.133871s | 1.30x |
| `PhyloHeatmap._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.057904s | 0.027449s | 2.11x |
| `PhyloHeatmap._plot_phylo_heatmap` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram tree-panel coordinate setup | 0.049305s | 0.040452s | 1.22x |
| `PhyloHeatmap._plot_phylo_heatmap_circular` phylogram setup | balanced 32768-tip tree, precomputed heatmap matrix | 0.4398s | 0.1925s | 2.3x |
| `PhyloHeatmap._plot_phylo_heatmap_circular` phylogram circular coordinate clade-list reuse | balanced 32768-tip tree, node x positions, parent map, preorder list, and tips already available | 0.056538s | 0.045025s | 1.26x |
| `PhyloHeatmap._plot_phylo_heatmap_circular` cladogram setup | balanced 32768-tip tree, precomputed heatmap matrix | 0.3637s | 0.1921s | 1.9x |
| `PhyloHeatmap._plot_phylo_heatmap_circular` cladogram node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, cladogram coordinate setup | 0.058446s | 0.047804s | 1.22x |
| `PhyloHeatmap._plot_phylo_heatmap_circular` cladogram circular coordinate clade-list reuse | balanced 32768-tip tree, cladogram node x positions, parent map, preorder list, and tips already available | 0.055849s | 0.044909s | 1.24x |
| `PhyloHeatmap._plot_phylo_heatmap_circular` ring-cell rendering | 1024 tips x 12 traits, real Matplotlib Agg heatmap wedge render | 6.315404s | 1.000884s | 6.31x |
| `PhyloHeatmap` rectangular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 1.226909s | 0.041817s | 29.34x |
| `PhyloHeatmap` circular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 0.609568s | 0.029928s | 20.37x |
| `phylo_heatmap` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.096110s | 0.032152s | 2.99x |
| `phylo_heatmap` module import without eager JSON/plot/color helpers | median cold subprocess import after lazy wrappers/localized imports for JSON, plot config, plot helpers, and color annotations | 0.013094s | 0.005342s | 2.45x |
| `phylo_heatmap` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006222s | 0.003853s | 1.62x |
| `Cophylo._get_tip_order` | balanced tree with 4096 tips | 0.0081s | 0.0008s | 10.1x |
| `Cophylo.run` tip-name setup | two balanced trees with 4096 shared tips | 0.0159s | 0.0018s | 8.8x |
| `Cophylo._validate_tree` direct traversal | balanced 65536-tip tree, count tips and fill missing branch lengths | 0.3349s | 0.0211s | 15.9x |
| `Cophylo._validate_standard_tree` unordered validation scan | balanced 131072-tip tree, count tips and fill missing branch lengths, optimized helper baseline | 0.026761s | 0.017404s | 1.54x |
| `Cophylo._rotate_tree` | balanced tree with 4096 tips, reversed matched target order | 0.3578s | 0.1501s | 2.4x |
| `Cophylo._rotate_tree` aggregate position means | balanced tree with 4096 tips, reversed matched target order | 0.0298s | 0.0293s | 1.02x |
| `Cophylo._rotate_tree` direct postorder | balanced 8192-tip tree, reversed matched target order | 0.067452s | 0.011656s | 5.79x |
| `Cophylo._postorder_clades` reverse-preorder helper | balanced 32768-tip tree, direct postorder traversal list | 0.008283s | 0.003330s | 2.49x |
| `Cophylo._postorder_clades` localized stack operations | balanced 262144-tip tree, direct postorder traversal list with identical clade order | 0.083764s | 0.073222s | 1.14x |
| `Cophylo._rotate_tree` reverse-preorder postorder helper | balanced 16384-tip tree, reversed matched target order | 0.025372s | 0.022264s | 1.14x |
| `Cophylo._rotate_tree` single-pass binary swaps | balanced 32768-tip tree, reversed matched target order, side-by-side previous two-phase mean-position path | 0.185358s | 0.142129s | 1.30x |
| `Cophylo._rotate_tree` direct child aggregate lookups | balanced 32768-tip tree, reversed matched target order, side-by-side single-pass helper with `.get()` lookups | 0.074279s | 0.062701s | 1.18x |
| `Cophylo._build_parent_map` unordered child push | balanced 65536-tip tree, object parent map for plotting setup, optimized helper baseline | 0.023673s | 0.018281s | 1.29x |
| `Cophylo._preorder_clades` order-preserving child push | balanced 131072-tip tree, direct preorder helper with identical clade order, side-by-side previous `reversed(children)` helper | 0.051055s | 0.043574s | 1.17x |
| `Cophylo._terminal_clades` order-preserving child push | balanced 131072-tip tree, direct terminal helper with identical terminal order, side-by-side previous `reversed(children)` helper | 0.044053s | 0.040119s | 1.10x |
| `Cophylo._draw_phylogram` rectangular setup | balanced 8192-tip tree, no-op axes for branch/label drawing | 0.1437s | 0.0556s | 2.6x |
| `Cophylo._draw_phylogram` batched LineCollections | balanced 2048-tip tree, one side phylogram, real Matplotlib Agg branch/label render | 1.928144s | 0.522134s | 3.69x |
| `Cophylo._plot_cophylo_rect` association connectors | 4096 mapped taxa, real Matplotlib Agg middle-panel connector render | 0.693452s | 0.098376s | 7.05x |
| `Cophylo` rectangular clade-color overlay rendering | two balanced 2048-tip trees, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 2.652615s | 0.090716s | 29.24x |
| `Cophylo` circular clade-color overlay rendering | two balanced 2048-tip trees, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 1.287788s | 0.057825s | 22.27x |
| `Cophylo._plot_cophylo_circular` circular coordinate and terminal-list reuse | two balanced 32768-tip trees, node x positions, parent maps, preorder lists, and tips already available | 0.136640s | 0.093813s | 1.46x |
| `Cophylo.run` cached tree2 setup | two balanced 32768-tip trees, validation, tip ordering, plotting, and text output mocked | 0.245234s | 0.209092s | 1.17x |
| `Cophylo._print_text_output` batched summary | 100k captured cophylo text summaries, identical stdout text | 0.099878s | 0.057622s | 1.73x |
| `Cophylo._parse_mapping_file` streaming parser | 1M two-column mapped taxa rows with comments/blanks | 0.393696s | 0.391772s | 1.00x |
| `Cophylo._parse_mapping_file` partition parser | 1M two-column mapped taxa rows with comments/blanks, side-by-side previous split parser comparison | 0.392522s | 0.351642s | 1.12x |
| `cophylo` module import without NumPy | cold subprocess import after replacing internal y-position mean with Python arithmetic | 0.086469s | 0.031961s | 2.71x |
| `cophylo` module import without eager color annotations | cold subprocess import after localizing optional color-file helpers | 0.014231s | 0.013016s | 1.09x |
| `cophylo` module import without eager JSON/plot helpers | median cold subprocess import after localizing PlotConfig, cladogram layout helper, and lazy JSON wrapper | 0.013248s | 0.005757s | 2.30x |
| `cophylo` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006165s | 0.002045s | 3.01x |
| `SimmapSummary._summarize_per_branch` direct accumulation | balanced 1024-tip tree, 80 stochastic mappings, 3 states | 0.3717s | 0.1183s | 3.1x |
| `SimmapSummary._summarize_per_branch` single-segment history fast path | 80 mappings x 2046 branches, 90% single-segment histories, 3 states | 0.107847s | 0.039522s | 2.73x |
| `SimmapSummary._summarize_per_branch` direct branch-clade setup | balanced 8192-tip tree, empty branch histories, 2 states | 0.071274s | 0.048508s | 1.47x |
| `SimmapSummary._summarize_per_branch` single branch-data table | balanced 8192-tip tree, empty branch histories, 2 states | 0.048508s | 0.041788s | 1.16x |
| `SimmapSummary._summarize_node_posteriors` direct traversal | balanced 32768-tip tree, 3 mappings, 3 states | 0.140093s | 0.044254s | 3.17x |
| `SimmapSummary._summarize_node_posteriors` ordered child push and accumulation | balanced 32768-tip tree, 3 mappings, 3 states, optimized helper baseline | 0.053613s | 0.047107s | 1.14x |
| `SimmapSummary._collect_clade_tip_names` direct postorder | balanced 65536-tip tree, branch-label descendant tip-name cache | 0.242102s | 0.092232s | 2.62x |
| `SimmapSummary._collect_clade_tip_names` reverse-preorder helper | balanced 65536-tip tree, branch-label descendant tip-name cache | 0.121940s | 0.080830s | 1.51x |
| `SimmapSummary._print_text` batched table output | balanced 32768-tip synthetic tree, 3 states, branch/node summaries for every clade, captured stdout and identical text | 0.467031s | 0.337105s | 1.39x |
| `SimmapSummary._print_text` direct preorder row preparation | balanced 16384-tip synthetic tree, 3 states, branch/node rows for every clade, identical row text | 0.142745s | 0.097380s | 1.47x |
| `SimmapSummary._print_text` row-wise dense tables | 64-state synthetic SIMMAP summary, 600 branch rows, 600 node rows, identical text | 0.031693s | 0.027633s | 1.15x |
| `SimmapSummary._print_text` percent row formatting | balanced 16384-tip synthetic tree, 3 states, branch/node summaries for every clade, captured stdout and identical text, side-by-side previous f-string row-builder comparison | 0.124867s | 0.081454s | 1.53x |
| `SimmapSummary._print_json_output` direct preorder reuse | balanced 16384-tip synthetic tree, 3 states, branch/node summaries for every clade, identical payload data | 0.191861s | 0.092817s | 2.07x |
| `SimmapSummary._print_json_output` Q-matrix row iteration | 32-state synthetic Q matrix, rounded nested JSON payload | 0.000321s | 0.000275s | 1.17x |
| `SimmapSummary._write_csv` direct preorder reuse | balanced 16384-tip synthetic tree, 3 states, branch/node summaries for every clade, identical CSV text | 0.185094s | 0.086272s | 2.15x |
| `SimmapSummary._write_csv` joined branch proportion rows | 120k branch rows x 32 states, identical branch CSV row text | 1.306681s | 1.114563s | 1.17x |
| `SimmapSummary._write_csv` joined node posterior rows | 120k node posterior rows x 32 states, identical node CSV row text | 1.804939s | 1.663184s | 1.09x |
| `SimmapSummary._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.035303s | 0.023099s | 1.53x |
| `SimmapSummary._plot_posterior_pie` setup | balanced 32768-tip tree, precomputed summaries and parent map | 0.4529s | 0.1683s | 2.7x |
| `SimmapSummary._plot_posterior_pie` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.048659s | 0.037910s | 1.28x |
| `SimmapSummary._plot_posterior_pie` batched base branches | balanced 2048-tip tree, precomputed branch summaries, pie placement disabled, real Matplotlib Agg render | 3.526622s | 1.912688s | 1.84x |
| `SimmapSummary.run` all-tip-state read-only setup | balanced 32768-tip cached tree, trait data for every tip, fitting/simulation/summary/output mocked | 0.264782s | 0.047127s | 5.62x |
| `simmap_summary` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.083424s | 0.032221s | 2.59x |
| `simmap_summary` module import without eager JSON/plot helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig`/plot helper imports | 0.013608s | 0.007155s | 1.90x |
| `simmap_summary` module import without inherited `typing` startup | median cold subprocess import after removing annotation-only typing aliases from SIMMAP summary and its stochastic-map parent | 0.005650s | 0.003692s | 1.53x |
| `StochasticCharacterMap._summarize_simulations` direct accumulation | balanced 1024-tip tree, 80 stochastic mappings, 3 states | 0.3331s | 0.0943s | 3.5x |
| `StochasticCharacterMap._summarize_simulations` direct setup traversal | balanced 32768-tip tree, 3 stochastic mappings, 3 states | 0.2189s | 0.1238s | 1.8x |
| `StochasticCharacterMap._summarize_simulations` single-segment history fast path | 100 mappings x 4094 branches, 90% single-segment histories, 4 states | 0.190871s | 0.103495s | 1.84x |
| `StochasticCharacterMap._sample_ancestral_states` transition cache | balanced 512-tip tree, 80 ancestral-state samples, repeated branch lengths | 0.6845s | 0.4208s | 1.6x |
| `StochasticCharacterMap._sample_ancestral_states` precomputed parent-state probabilities | balanced 512-tip tree, 80 ancestral-state samples, repeated branch lengths | 0.160799s | 0.092451s | 1.7x |
| `StochasticCharacterMap._prepare_ancestral_sampling_nodes` vectorized large-state CDFs | 32768 sampling nodes, 8 states, repeated unit branch length | 0.520323s | 0.246565s | 2.11x |
| `StochasticCharacterMap._prepare_ancestral_sampling_nodes` all-positive large-state CDFs | 32768 sampling nodes, 8 states, repeated unit branch length and positive child likelihoods | 0.393498s | 0.220084s | 1.79x |
| `StochasticCharacterMap._simulate_branch_history` precomputed rates | 8k branch-history samples, 4 states, repeated start/end states | 0.3931s | 0.3005s | 1.3x |
| `StochasticCharacterMap` Q diagonal rate extraction | 2 / 4 / 16 / 64 / 256 / 1024 / 2048 states, side-by-side previous `-np.diag(Q)` rate vector copy | 0.000000875s / 0.000001098s / 0.000001258s / 0.000001612s / 0.000000908s / 0.000001260s / 0.000002136s | 0.000000496s / 0.000000304s / 0.000000300s / 0.000000306s / 0.000000346s / 0.000000669s / 0.000001446s | 1.77x / 3.61x / 4.19x / 5.27x / 2.62x / 1.88x / 1.48x |
| `StochasticCharacterMap._simulate_branch_history` precomputed transition CDFs | 8k branch-history samples, 4 states, repeated start/end states | 0.135483s | 0.094216s | 1.44x |
| `StochasticCharacterMap._simulate_branch_history` deterministic two-state transitions | 8k branch-history samples, 2 states, repeated start/end states | 0.021078s | 0.016228s | 1.30x |
| `StochasticCharacterMap._run_single_simulation` traversal metadata | balanced 512-tip tree, 80 stochastic mappings, cached transition matrices | 1.3863s | 1.0821s | 1.3x |
| `StochasticCharacterMap._run_single_simulation` categorical draw helper | balanced 512-tip tree, 80 stochastic mappings, cached metadata and transition matrices | 1.133712s | 0.772678s | 1.47x |
| `StochasticCharacterMap._run_single_simulation` branch-history context reuse | balanced 256-tip tree, 80 stochastic mappings, 8 states, cached metadata and transition matrices | 0.657115s | 0.584082s | 1.13x |
| `StochasticCharacterMap._run_single_simulation` two-state deterministic branch transitions | balanced 512-tip tree, 80 stochastic mappings, cached metadata and branch-history context | 0.215536s | 0.209447s | 1.03x |
| `StochasticCharacterMap._sample_categorical_cdf` small-state branches | balanced 256-tip tree, 80 stochastic mappings, 2 states, cached metadata and transition matrices | 0.173280s | 0.111058s | 1.56x |
| `StochasticCharacterMap._sample_categorical_cdf` 8-state branch | 500k scalar categorical draws from an 8-state CDF | 0.796993s | 0.363060s | 2.20x |
| `StochasticCharacterMap._sample_categorical_cdf` 8-state binary comparisons | 500k scalar categorical draws from an 8-state CDF | 0.217262s | 0.189751s | 1.14x |
| `StochasticCharacterMap._count_missing_tip_states` unordered scan | balanced 131072-tip tree, half the tips present in trait map, optimized helper baseline | 0.026210s | 0.019871s | 1.32x |
| `StochasticCharacterMap._format_result` row-wise payload dictionaries | 32-state synthetic summary, Q matrix, dwelling vectors, and transition matrix | 0.000415s | 0.000246s | 1.69x |
| `StochasticCharacterMap._build_parent_map` direct traversal | balanced 32768-tip tree, id-to-parent map for SIMMAP setup | 0.1006s | 0.0110s | 9.1x |
| `StochasticCharacterMap._build_parent_map` unordered child push | balanced 65536-tip tree, id-to-parent map for SIMMAP setup, optimized helper baseline | 0.024008s | 0.018110s | 1.33x |
| `StochasticCharacterMap._build_simulation_metadata` direct traversal | balanced 32768-tip tree, SIMMAP simulation and branch-node metadata | 0.1226s | 0.0281s | 4.4x |
| `StochasticCharacterMap._print_text_output` batched summary output | 400-state synthetic SIMMAP summary, dense Q/transition tables, captured stdout and identical text | 0.152400s | 0.134238s | 1.14x |
| `StochasticCharacterMap._print_text_output` row-wise dense tables | 200-state synthetic SIMMAP summary, dense Q/transition tables, identical text | 0.035830s | 0.026349s | 1.36x |
| `StochasticCharacterMap._print_text_output` percent Q-row formatting | 400-state synthetic SIMMAP summary, dense Q/transition tables, captured stdout and identical text, side-by-side previous f-string Q-row formatter comparison | 0.108428s | 0.078860s | 1.37x |
| `StochasticCharacterMap._plot_stochastic_map` phylogram layout setup | balanced 32768-tip tree, precomputed mapping and parent map | 0.5905s | 0.0837s | 7.1x |
| `StochasticCharacterMap._plot_stochastic_map` cladogram layout setup | balanced 32768-tip tree, precomputed mapping and parent map | 0.7017s | 0.0978s | 7.2x |
| `StochasticCharacterMap._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.053782s | 0.023900s | 2.25x |
| `StochasticCharacterMap._plot_stochastic_map` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.047556s | 0.037577s | 1.27x |
| `StochasticCharacterMap._plot_stochastic_map` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.063955s | 0.048064s | 1.33x |
| `StochasticCharacterMap._plot_stochastic_map` rectangular batched history segments | balanced 512-tip tree, two mapped segments per branch, real Matplotlib Agg render | 1.831008s | 1.295113s | 1.41x |
| `StochasticCharacterMap._plot_stochastic_map` circular batched history segments/arcs | balanced 512-tip tree, two mapped segments per branch plus internal arcs, real Matplotlib Agg render | 0.901516s | 0.428369s | 2.10x |
| `stochastic_character_map` module import via lazy discrete SciPy helpers | cold process import for stochastic-character-map command module | 0.482006s | 0.202569s | 2.4x |
| `stochastic_character_map` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.086216s | 0.032006s | 2.69x |
| `stochastic_character_map` module import without eager JSON/plot/discrete helpers | median cold subprocess import after localizing output, plotting, color, circular, and discrete-model helpers | 0.014049s | 0.005457s | 2.57x |
| `stochastic_character_map` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005457s | 0.003808s | 1.43x |
| `Phenogram._fast_anc` plus node labeling | balanced tree with 4096 tips, one continuous trait per tip | 0.1058s | 0.0405s | 2.6x |
| `Phenogram._fast_anc` postorder weighted sums | balanced 32768-tip tree, helper-only postorder downpass | 0.065973s | 0.056026s | 1.18x |
| `Phenogram._fast_anc` cached id lookups | balanced 32768-tip tree, full helper with exact node estimates and sigma2 | 0.212574s | 0.171936s | 1.24x |
| `Phenogram._compute_sigma2_from_contrasts` streaming valid children | balanced 32768-tip tree, postorder contrast scan | 0.039013s | 0.025422s | 1.53x |
| `Phenogram._compute_sigma2_from_contrasts` reverse-preorder fallback helper | balanced 32768-tip tree, postorder contrast scan without supplied clade list | 0.035139s | 0.027408s | 1.28x |
| `Phenogram`/`ContMap._build_parent_map` direct map traversal | balanced 65536-tip tree, shared map-only parent helper, optimized helper baseline | 0.030587s | 0.018606s | 1.64x |
| `Phenogram.run` copied-tree prune setup | balanced 32768-tip tree, all tips shared with trait data | 0.0671s | 0.0072s | 9.4x |
| `Phenogram.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, ancestral reconstruction, plotting, and output mocked; protective prune-copy retained | 1.004997s | 0.450786s | 2.23x |
| `Phenogram.run` all-shared read-only setup | balanced 32768-tip cached tree, one trait for every tip, reconstruction/plot/output mocked | 0.300021s | 0.063182s | 4.75x |
| `Phenogram._parse_single_trait_data` streaming valid-row parser | 500k two-column trait rows with comments/blanks, all taxa shared | 0.467365s | 0.463534s | 1.01x |
| `Phenogram._parse_single_trait_data` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.433240s | 0.240995s | 1.80x |
| `Phenogram._parse_single_trait_data` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.235785s | 0.229632s | 1.03x |
| `Phenogram._plot_phenogram` layout and branch setup | balanced 32768-tip tree, precomputed ancestral estimates and parent map | 0.5104s | 0.0577s | 8.8x |
| `Phenogram._plot_phenogram` estimate value range helper | 1M plotted node/tip estimate values, identical min/max range without temporary list | 0.024999s | 0.020917s | 1.20x |
| `Phenogram._plot_phenogram` batched gradient branches | balanced 512-tip tree, 50 color segments per branch, real Matplotlib Agg branch/label render | 11.650694s | 1.956507s | 5.95x |
| `Phenogram._print_text_output` batched summary | 100k captured phenogram text summaries, identical stdout text | 0.087289s | 0.057836s | 1.51x |
| `phenogram` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.088925s | 0.031880s | 2.79x |
| `phenogram` module import without eager pickle/JSON/plot helpers | median cold subprocess import after localizing pickle, PlotConfig, and lazy JSON wrapper | 0.014977s | 0.008669s | 1.73x |
| `phenogram` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005520s | 0.004003s | 1.38x |
| `ContMap._fast_anc` plus node labeling | balanced tree with 4096 tips, one continuous trait per tip | 0.1056s | 0.0404s | 2.6x |
| `ContMap._fast_anc` postorder weighted sums | balanced 32768-tip tree, helper-only postorder downpass | 0.071900s | 0.038321s | 1.88x |
| `ContMap._fast_anc` cached id lookups | balanced 32768-tip tree, full helper with exact node estimates and sigma2 | 0.209880s | 0.169919s | 1.24x |
| `ContMap._compute_sigma2_from_contrasts` streaming valid children | balanced 32768-tip tree, postorder contrast scan | 0.034165s | 0.018992s | 1.80x |
| `ContMap._compute_sigma2_from_contrasts` reverse-preorder fallback helper | balanced 32768-tip tree, postorder contrast scan without supplied clade list | 0.039045s | 0.031494s | 1.24x |
| `ContMap.run` copied-tree prune setup | balanced 32768-tip tree, all tips shared with trait data | 0.0671s | 0.0071s | 9.4x |
| `ContMap.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, ancestral reconstruction, plotting, and output mocked; protective prune-copy retained | 1.460611s | 0.388562s | 3.76x |
| `ContMap.run` all-shared read-only setup | balanced 32768-tip cached tree, one trait for every tip, ladderize off, reconstruction/plot/output mocked | 0.306910s | 0.059638s | 5.15x |
| `ContMap._parse_single_trait_data` streaming valid-row parser | 500k two-column trait rows with comments/blanks, all taxa shared | 0.465422s | 0.456907s | 1.02x |
| `ContMap._parse_single_trait_data` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.436988s | 0.237011s | 1.84x |
| `ContMap._parse_single_trait_data` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.241585s | 0.225893s | 1.07x |
| `ContMap._plot_contmap` phylogram layout setup | balanced 32768-tip tree, precomputed ancestral estimates and parent map | 0.4223s | 0.0633s | 6.7x |
| `ContMap._plot_contmap` cladogram layout setup | balanced 32768-tip tree, precomputed ancestral estimates and parent map | 0.5226s | 0.1030s | 5.1x |
| `ContMap._prepare_contmap_plot_data` single-pass setup | balanced 32768-tip tree, parent map, preorder, estimates, and tips | 0.038240s | 0.027062s | 1.41x |
| `ContMap._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.089880s | 0.055680s | 1.61x |
| `ContMap._plot_contmap` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.049826s | 0.042940s | 1.16x |
| `ContMap._plot_contmap` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.069976s | 0.052805s | 1.33x |
| `ContMap._plot_contmap` estimate value range helper | 1M plotted node/tip estimate values, identical min/max range without temporary list | 0.036248s | 0.032907s | 1.10x |
| `ContMap._plot_contmap` rectangular batched gradient branches | balanced 512-tip tree, 50 color segments per branch, real Matplotlib Agg branch/label/colorbar render | 11.957742s | 1.949793s | 6.13x |
| `ContMap._print_text_output` batched summary | 100k captured contMap text summaries, identical stdout text | 0.089211s | 0.059444s | 1.50x |
| `cont_map` module import without eager NumPy | cold subprocess import after lazy NumPy proxy plus lazy circular-layout arc cache | 0.085155s | 0.038423s | 2.22x |
| `cont_map` module import without eager pickle/JSON/plot helpers | median cold subprocess import after localizing pickle, PlotConfig/layout helpers, circular/color helpers, and lazy JSON wrapper | 0.014370s | 0.005967s | 2.41x |
| `cont_map` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.015563s | 0.004050s | 3.84x |
| `TraitRateMap` branch-rate computation pipeline | balanced tree with 4096 tips, one continuous trait per tip | 0.0779s | 0.0283s | 2.7x |
| `TraitRateMap._ancestral_reconstruction` direct weighted sums | balanced 32768-tip tree, one continuous trait per tip | 0.047155s | 0.025578s | 1.84x |
| `TraitRateMap._ancestral_reconstruction` reverse-preorder helper | balanced 32768-tip tree, one continuous trait per tip | 0.032665s | 0.024686s | 1.32x |
| `TraitRateMap._compute_branch_rates` cached id lookups | balanced 32768-tip tree, precomputed node values and labels, exact branch-rate entries | 0.058582s | 0.041518s | 1.41x |
| `TraitRateMap.run` copied-tree prune setup | balanced 32768-tip tree, all tips shared with trait data | 0.0671s | 0.0071s | 9.5x |
| `TraitRateMap.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, reconstruction, branch-rate calculation, plotting, and output mocked; protective prune-copy retained | 1.578476s | 0.628135s | 2.51x |
| `TraitRateMap.run` all-shared read-only setup | balanced 32768-tip cached tree, one trait for every tip, ladderize off, reconstruction/rates/plot/output mocked | 0.283459s | 0.061932s | 4.58x |
| `TraitRateMap._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.044759s | 0.035023s | 1.28x |
| `TraitRateMap._parse_single_trait_data` streaming valid-row parser | 500k two-column trait rows with comments/blanks, all taxa shared | 0.467069s | 0.460840s | 1.01x |
| `TraitRateMap._parse_single_trait_data` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.435507s | 0.238663s | 1.82x |
| `TraitRateMap._parse_single_trait_data` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.237344s | 0.227041s | 1.05x |
| `TraitRateMap._plot_rate_map` phylogram layout setup | balanced 32768-tip tree, precomputed rates and parent map | 0.6640s | 0.1272s | 5.2x |
| `TraitRateMap._plot_rate_map` cladogram layout setup | balanced 32768-tip tree, precomputed rates and parent map | 0.7293s | 0.3048s | 2.4x |
| `TraitRateMap._plot_rate_map` precomputed layout traversal | balanced 32768-tip tree, cladogram coordinate setup with existing preorder list | 0.070080s | 0.057533s | 1.22x |
| `TraitRateMap._plot_rate_map` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.072342s | 0.056249s | 1.29x |
| `TraitRateMap._plot_rate_map` rectangular batched rate branches | balanced 2048-tip tree, varied branch rates, real Matplotlib Agg branch/label/colorbar render | 2.158820s | 0.663060s | 3.26x |
| `TraitRateMap._plot_rate_map` circular batched rate branches/arcs | balanced 2048-tip tree, varied branch rates, real Matplotlib Agg branch/arc/colorbar render | 1.606744s | 0.352826s | 4.55x |
| `TraitRateMap._plot_rate_map` repeated branch-rate color cache | balanced 2048-tip tree, rectangular Agg plot with labels/title disabled and repeated branch rate | 0.547496s | 0.427287s | 1.28x |
| `TraitRateMap._prepare_rate_map_plot_data` single-pass setup | balanced 32768-tip tree, rate lookup, preorder, tips, and rate list | 0.017025s | 0.012340s | 1.38x |
| `TraitRateMap._print_text_output` batched summary | 100k captured trait-rate-map text summaries, identical stdout text | 0.174115s | 0.110930s | 1.57x |
| `TraitRateMap._parse_multi_trait_data` selected-column streaming | 200k-row TSV, comments/blanks before header, one selected continuous trait | 0.163704s | 0.155640s | 1.05x |
| `TraitRateMap._parse_multi_trait_data` bounded selected-column split | 300k-row TSV, 8 numeric trait columns, selected early continuous trait, all taxa shared | 0.343296s | 0.281438s | 1.22x |
| `trait_rate_map` module import without NumPy | cold subprocess import after replacing summary mean with Python arithmetic | 0.090290s | 0.031843s | 2.84x |
| `trait_rate_map` module import without eager pickle/JSON/plot helpers | median cold subprocess import after localizing pickle, PlotConfig/layout helpers, circular/color helpers, and lazy JSON wrapper | 0.014715s | 0.006266s | 2.35x |
| `trait_rate_map` module import without `typing` startup | median cold subprocess import after postponing annotations and converting typing aliases to built-in annotations | 0.006567s | 0.004268s | 1.54x |
| `PhyloLogistic._build_logistic_vcv` | balanced tree with 600 tips, `alpha=0.05` | 0.3304s | 0.1492s | 2.2x |
| `PhyloLogistic._build_logistic_vcv` diagonal correction update | 8 / 40 / 260 / 900 / 2000 taxa VCV matrices with vector correction, isolated in-place update and reset, side-by-side previous `np.fill_diagonal(np.diag(vcv) + diag_corr)` path | 0.000006873s / 0.000007194s / 0.000017359s / 0.000073396s / 0.000025940s | 0.000002899s / 0.000001468s / 0.000001254s / 0.000002492s / 0.000011652s | 2.37x / 4.90x / 13.84x / 29.45x / 2.23x |
| `PhyloLogistic._root_tip_distances` | balanced 65536-tip tree, ordered OU diagonal distances | 0.1591s | 0.0242s | 6.6x |
| `PhyloLogistic._neg_pen_loglik` reused root distances and scalar branch transform | balanced 256-tip tree, one Firth-correction likelihood evaluation | 0.002322s | 0.002147s | 1.08x |
| `PhyloLogistic._compute_info_matrix` | 420 taxa SPD VCV x 3-predictor design matrix | 0.0052s | 0.0007s | 8.0x |
| `PhyloLogistic._compute_info_matrix_cholesky` cached SciPy linalg wrappers | 2k repeated 8-taxon SPD VCV x 2-column design matrix calls, SciPy already warm | 0.024984s | 0.022490s | 1.11x |
| `PhyloLogistic._compute_info_matrix_inverse` row scaling | 900 taxa SPD VCV x 8-column design matrix, inverse already available | 0.001697s | 0.000365s | 4.7x |
| `PhyloLogistic` standard-error diagonal solve | 200-coefficient SPD information matrix, side-by-side previous explicit inverse diagonal extraction | 0.002265682s | 0.000320565s | 7.07x |
| Tree model read-only diagonal extraction | 2 / 8 / 64 / 512 / 2048 square matrices, side-by-side previous `np.diag(matrix)` wrapper | 0.000001570s / 0.000001338s / 0.000001205s / 0.000000862s / 0.000000824s | 0.000000122s / 0.000000121s / 0.000000160s / 0.000000093s / 0.000000093s | 12.88x / 11.05x / 7.52x / 9.26x / 8.82x |
| `PhyloLogistic._normal_two_tailed_p_values` vectorized special erfc | 200k z statistics, side-by-side previous scalar Python `math.erfc` loop | 0.024712s | 0.001564s | 15.80x |
| `PhyloLogistic._logistic_starting_values` row-scaled IRLS | 1500 taxa x 8-column design matrix, synthetic binary response | 0.018337s | 0.000254s | 72.3x |
| `PhyloLogistic`/`PhylogeneticGLM` starting-value convergence max | 2 / 4 / 8 / 32 / 128 coefficient deltas, side-by-side previous `np.max(np.abs(beta_new - beta))` loop check | 0.000007025s / 0.000007767s / 0.000007604s / 0.000008465s / 0.000007491s | 0.000003384s / 0.000004525s / 0.000003101s / 0.000004361s / 0.000004134s | 2.08x / 1.72x / 2.45x / 1.94x / 1.81x |
| `PhyloLogistic`/`PhylogeneticGLM` binary fallback class count | 80k validated 0/1 responses, side-by-side previous two equality-mask reductions | 0.000089306s | 0.000012742s | 7.01x |
| `phylo_logistic` module import without `scipy.stats` | cold process import for logistic regression command module | 0.624020s | 0.431614s | 1.45x |
| `phylo_logistic` module import without eager SciPy linalg/optimize | cold process import for logistic regression command module | 0.440761s | 0.166163s | 2.7x |
| `phylo_logistic` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.082503s | 0.026603s | 3.10x |
| `phylo_logistic` module import without eager pickle/JSON/trait parsing helpers | median cold subprocess import after localizing pickle and lazy helper wrappers | 0.010470s | 0.006457s | 1.62x |
| `phylo_logistic` module import without `typing` startup | median cold subprocess import after converting annotation-only typing names to built-in postponed annotations | 0.006210s | 0.004573s | 1.36x |
| `PhyloLogistic.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, fitting, and output mocked | 0.516410s | 0.000033s | 15648.79x |
| `PhyloLogistic._print_text_output` batched coefficient output | 100k coefficient rows, captured stdout and identical text | 0.116695s | 0.101976s | 1.14x |
| `PhyloLogistic._print_text_output` inline significance thresholds | 100k coefficient rows, captured stdout and identical text, side-by-side previous `_signif_code()` loop | 0.105801s | 0.098002s | 1.08x |
| `PhyloLogistic._format_result` zip-based coefficient mapping | 5k coefficient JSON entries, side-by-side previous index lookups | 0.148853s | 0.066884s | 2.23x |
| `FitContinuous._vcv_ou` | 900-tip synthetic VCV matrix, `alpha=0.7` | 0.340860s | 0.009100s | 37.5x |
| `FitContinuous.run`/`OUwie.run` tree-height diagonal maximum | 4 / 16 / 64 / 256 / 1024 / 2048 square VCV matrices, side-by-side previous `np.max(np.diag(vcv))` setup | 0.000002088s / 0.000002103s / 0.000002091s / 0.000002138s / 0.000002588s / 0.000003452s | 0.000001150s / 0.000001157s / 0.000001159s / 0.000001211s / 0.000001572s / 0.000002449s | 1.82x / 1.82x / 1.80x / 1.77x / 1.65x / 1.41x |
| `FitContinuous._build_transformed_vcv` branch accumulation | balanced 1024-tip synthetic root-to-tip paths, delta-style transform | 0.216503s | 0.015251s | 14.2x |
| `FitContinuous._build_transformed_vcv` single-tip diagonal updates | balanced 4096-tip synthetic root-to-tip paths, transformed branch lengths | 0.237096s | 0.177280s | 1.34x |
| `FitContinuous._build_parent_map` direct traversal | balanced 65536-tip tree, parent-id map setup | 0.277243s | 0.028661s | 9.67x |
| `FitContinuous._build_parent_map` unordered child push | balanced 65536-tip tree, parent-id map setup, optimized helper baseline | 0.021944s | 0.018356s | 1.20x |
| `FitContinuous.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, VCV/model fitting, and output mocked | 0.483131s | 0.000328s | 1473.15x |
| `FitContinuous._parse_trait_file` streaming valid-row parser | 500k two-column trait rows with comments/blanks, all taxa shared | 0.471328s | 0.458778s | 1.03x |
| `FitContinuous._parse_trait_file` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.433228s | 0.238331s | 1.82x |
| `FitContinuous._parse_trait_file` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.239696s | 0.223276s | 1.07x |
| `FitContinuous._compute_model_comparison` scalar AIC weights | seven synthetic model-result dictionaries, identical normalized weights | 0.000018015s | 0.000006892s | 2.61x |
| `FitContinuous._print_text_output` batched model table | three captured text reports with 100k model rows each, identical stdout text | 0.776887s | 0.714921s | 1.09x |
| `FitContinuous._print_text_output` row-template formatting | 100k model rows, captured stdout and identical text, side-by-side previous f-string row formatter comparison | 0.235230s | 0.205315s | 1.15x |
| `FitContinuous._print_text_output` percent row formatting | 100k model rows, captured stdout and identical text, side-by-side previous `.format()` row formatter comparison | 0.192132s | 0.140039s | 1.37x |
| `FitContinuous._concentrated_ll_cholesky` cached SciPy linalg wrappers | 120 repeated 420-taxon SPD VCV concentrated likelihood evaluations, SciPy already warm, side-by-side previous import-on-call wrappers | 0.077998s | 0.057574s | 1.35x |
| `PhylogeneticGLM._make_ultrametric` | balanced tree with 2500 tips | 2.2726s | 0.0068s | 336.2x |
| `PhylogeneticGLM._root_tip_distances` | balanced 65536-tip tree, ordered ultrametric correction distances | 0.1598s | 0.0245s | 6.5x |
| `PhylogeneticGLM._make_ultrametric` root-height max reduction | 120 / 1200 / 10000 / 65536 root-to-tip distances, side-by-side previous `np.max(heights)` wrapper | 0.000004475s / 0.000006131s / 0.000005739s / 0.000014409s | 0.000002587s / 0.000001470s / 0.000001774s / 0.000011340s | 1.73x / 4.17x / 3.23x / 1.27x |
| `PhylogeneticGLM._poisson_starting_values` row-scaled IRLS | 1200 taxa x 8-column design matrix, synthetic counts | 0.015975s | 0.000268s | 59.5x |
| `PhylogeneticGLM._poisson_gee_information_and_score` row scaling | 1200 taxa SPD correlation inverse x 8-column design matrix | 0.002394s | 0.000952s | 2.5x |
| `PhylogeneticGLM._poisson_gee_information_and_score` RHS-first score multiply | 1200 taxa SPD correlation inverse x 8-column design matrix | 0.001306s | 0.000615s | 2.1x |
| `PhylogeneticGLM._poisson_gee_information_and_score_cholesky` | 100 repeated 800-taxon SPD correlation x 8-column design matrix Poisson GEE information/score evaluations | 0.059989s | 0.030535s | 1.96x |
| `PhylogeneticGLM._poisson_gee_information_cholesky` final covariance solve | 100 repeated 800-taxon SPD correlation x 8-column design matrix final Poisson GEE information evaluations | 0.018674s | 0.016339s | 1.14x |
| `PhylogeneticGLM._fit_poisson_gee` convergence delta reduction | 2 / 4 / 8 / 32 / 128 coefficient deltas, side-by-side previous `np.sum(np.abs(delta))` loop check | 0.000006064s / 0.000004654s / 0.000006191s / 0.000005158s / 0.000006489s | 0.000002508s / 0.000001882s / 0.000002131s / 0.000002579s / 0.000002675s | 2.42x / 2.47x / 2.90x / 2.00x / 2.43x |
| `PhylogeneticGLM` cached SciPy linalg wrappers | 5k repeated 2x2 Cholesky factor/solve calls, SciPy already warm, side-by-side previous import-on-call wrappers | 0.061566s | 0.052477s | 1.17x |
| `PhylogeneticGLM` standard-error diagonal solve | 200-coefficient SPD information matrix, side-by-side previous explicit inverse diagonal extraction | 0.002194492s | 0.000213396s | 10.28x |
| `PhylogeneticGLM._normal_two_tailed_p_values` vectorized special erfc | 200k z statistics, side-by-side previous scalar Python `math.erfc` loop | 0.024712s | 0.001635s | 15.11x |
| `phylogenetic_glm` module import without `scipy.stats` | cold process import for binomial/Poisson GLM command module | 0.625725s | 0.467757s | 1.34x |
| `phylogenetic_glm` module import without eager SciPy optimize | cold process import for binomial/Poisson GLM command module | 0.459322s | 0.156656s | 2.9x |
| `phylogenetic_glm` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.085569s | 0.026700s | 3.20x |
| `phylogenetic_glm` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006674s | 0.005485s | 1.22x |
| `phylogenetic_glm` module import without trait parsing/typing startup | median cold subprocess import after lazy trait parser wrapper and built-in postponed annotations | 0.006709s | 0.005068s | 1.32x |
| `PhylogeneticGLM.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, two model fits, and output mocked | 0.233055s | 0.000038s | 6133.03x |
| `PhylogeneticGLM.run` gene-tree VCV reuse | 1800 trait taxa with 900 shared gene-tree taxa, mocked discordance VCV build returning shared-taxa matrix and metadata | 0.013563s | 0.005657s | 2.40x |
| `PhylogeneticGLM._print_text_output` batched coefficient output | 100k coefficient rows, captured stdout and identical text | 0.119475s | 0.105142s | 1.14x |
| `PhylogeneticGLM._print_text_output` inline significance thresholds | 100k coefficient rows, captured stdout and identical text, side-by-side previous `_signif_code()` loop | 0.114085s | 0.099997s | 1.14x |
| `PhylogeneticGLM` zip-based fitted-value mapping | 500k fitted taxon values, side-by-side previous index lookups | 0.133401s | 0.115596s | 1.15x |
| `PhylogeneticGLM._format_result` zip-based coefficient mapping | 5k coefficient JSON entries, side-by-side previous index lookups | 0.003066s | 0.002548s | 1.20x |
| `Dtt._compute_dtt` | balanced tree with 300 tips x 2 traits | 0.8787s | 0.0472s | 18.6x |
| `Dtt._compute_dtt` lineage event sweep | balanced tree with 300 tips x 2 traits, optimized DTT baseline | 0.0472s | 0.0110s | 4.3x |
| `Dtt._compute_dtt` lineage list means | 5k synthetic lineage events over 20k cached clade disparities, identical DTT values | 0.153168s | 0.062355s | 2.46x |
| `Dtt._simulate_null_avg_sq_batch` selected clade means | 1000 simulations x 5000 clades, selected-column widths 1 / 2 / 4 / 8, side-by-side previous temporary slice plus `np.mean(..., axis=1)` | 0.000005500s / 0.000007084s / 0.000009541s / 0.000013687s | 0.000000333s / 0.000003417s / 0.000005916s / 0.000010459s | 16.52x / 2.07x / 1.61x / 1.31x |
| `Dtt._prepare_dtt_context` direct standard-tree setup | balanced 8192-tip tree, all tips represented in trait order | 0.302738s | 0.273333s | 1.11x |
| `Dtt._prepare_dtt_context` streaming terminal-depth max | 1M terminal depths, identical root-relative tree height without temporary list | 0.132334s | 0.078064s | 1.70x |
| `Dtt._simulate_null` DTT context reuse | balanced 512-tip tree x 2 traits, 50 simulated DTT curves | 7.4703s | 0.5626s | 13.3x |
| `Dtt._simulate_null` batched avg-squared simulations | balanced 512-tip tree x 2 traits, 50 simulated DTT curves | 0.416623s | 0.048629s | 8.6x |
| `Dtt._simulate_null` batched terminal-point extension | 20k simulated DTT rows x 1024 time points, terminal observation grid | 0.046184s | 0.022259s | 2.07x |
| `Dtt._simulate_null` terminal-column masking | 200k simulated DTT rows x 80 time points, terminal observation grid final column | 0.138169s | 0.096509s | 1.43x |
| `Dtt._simulate_null` vectorized MDI reductions | balanced 512-tip tree x 2 traits, 50 simulated DTT curves | 0.119121s | 0.097534s | 1.22x |
| `Dtt._simulate_null` MDI p-value counts | 1M simulated MDI values, side-by-side previous `np.mean(abs(sim_mdis) >= abs(mdi))` reduction | 0.000991s | 0.000390s | 2.54x |
| `Dtt._compute_disparity` observed avg-squared sum-of-squares | observed trait matrices shaped 32x2 / 300x2 / 512x8 / 2048x4, side-by-side previous `np.sum(data * data)` and `np.sum(sums * sums)` reductions | 7.283553s / 2.413355s / 1.135589s / 3.148002s | 1.361637s / 1.171414s / 0.707937s / 2.592290s | 5.35x / 2.06x / 1.60x / 1.21x |
| `Dtt._compute_disparity` observed avg-Manhattan weighted sum | observed trait matrices shaped 50x1k / 100x10k / 100x50k / 1000x10k, side-by-side previous broadcast multiply plus `np.sum` reduction | 0.000028s / 0.001034s / 0.009492s / 0.019283s | 0.000006562s / 0.000119s / 0.000962s / 0.002053s | 4.28x / 8.66x / 9.87x / 9.39x |
| `Dtt._batch_clade_disparities_avg_sq` postorder subtree aggregation | balanced 2048-tip tree x 2 traits, 50 simulated DTT curves | 0.039455s | 0.026875s | 1.47x |
| `Dtt._simulate_null_avg_sq_batch` row sum-of-squares reductions | simulated trait cubes shaped 50x512x2 / 500x512x2 / 100x2048x4 / 1000x128x8, side-by-side previous `np.sum(values * values, ...)` reductions | 0.946699s / 1.249285s / 0.788768s / 0.593306s | 0.636814s / 0.980172s / 0.678222s / 0.359702s | 1.49x / 1.27x / 1.16x / 1.65x |
| `Dtt._simulate_null` observed-DTT reuse | balanced 512-tip tree x 2 traits, 50 simulated DTT curves | 0.071679s | 0.045011s | 1.59x |
| `Dtt._print_text` batched time-point output | 100k DTT time-point rows, captured stdout and identical text | 0.057375s | 0.046449s | 1.24x |
| `Dtt.run` trait matrix setup | 120k taxa x 12 parsed trait columns, selected column / full matrix setup | 0.013424s / 0.134250s | 0.007369s / 0.037870s | 1.82x / 3.54x |
| `Dtt.run` copied-tree prune setup | balanced 32768-tip tree, all tips shared with traits | 0.0671s | 0.0071s | 9.5x |
| `Dtt.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, DTT calculation, and output mocked; protective prune-copy retained | 1.087033s | 0.289034s | 3.76x |
| `Dtt.run` all-shared read-only setup | balanced 32768-tip cached tree, two traits for every tip, DTT/output mocked | 0.361751s | 0.094806s | 3.82x |
| `dtt` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.087831s | 0.031497s | 2.79x |
| `dtt` module import without eager JSON/plot/trait helpers | median cold subprocess import after localizing `PlotConfig` and lazy helper wrappers | 0.016348s | 0.005264s | 3.10x |
| `dtt` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005063s | 0.001347s | 3.76x |
| `Phylomorphospace._reconstruct_ancestral_scores` | balanced tree with 1200 tips x 2 scores | 1.6265s | 0.0587s | 27.7x |
| `Phylomorphospace._reconstruct_ancestral_scores` direct traversal and weighted sums | balanced tree with 4096 tips x 2 scores | 0.066720s | 0.034426s | 1.94x |
| `Phylomorphospace._reconstruct_ancestral_scores` prune setup | balanced 8192-tip tree, all tips retained for 2D scores | 0.2648s | 0.0021s | 128.3x |
| `Phylomorphospace._reconstruct_ancestral_scores` all-shared copy skip | balanced 32768-tip tree, all tips retained for 2D scores | 0.366612s | 0.107335s | 3.42x |
| `Phylomorphospace._plot_phylomorphospace` direct branch traversal | balanced 32768-tip tree, branch segment/color setup with precomputed node estimates | 0.240922s | 0.153599s | 1.57x |
| `Phylomorphospace._plot_phylomorphospace` root-distance max helper | 1M node-distance values, identical maximum-distance fallback behavior without temporary list | 0.028209s | 0.011876s | 2.38x |
| `Phylomorphospace._preorder_clades_direct` order-preserving child push | balanced 131072-tip tree, plotting preorder helper with identical clade order, side-by-side previous `reversed(children)` helper | 0.051598s | 0.047193s | 1.09x |
| `Phylomorphospace.run` trait matrix setup | 120k taxa x 12 parsed trait columns, selected x/y axes plus full color matrix | 0.194092s | 0.038694s | 5.02x |
| `Phylomorphospace`/`PhylogeneticOrdination._parse_color_by` numeric file values | 200k numeric color values, side-by-side previous list comprehension plus `np.array` conversion | 0.029553s | 0.025725s | 1.15x |
| `Phylomorphospace.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, reconstruction, plotting, and output mocked | 0.277017s | 0.000037s | 7486.95x |
| `phylomorphospace` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.089377s | 0.031532s | 2.83x |
| `phylomorphospace` module import without eager pickle/JSON/plot/trait helpers | median cold subprocess import after localizing pickle/PlotConfig and lazy helper wrappers | 0.014470s | 0.005942s | 2.44x |
| `phylomorphospace` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005568s | 0.001530s | 3.64x |
| `PhylogeneticOrdination._reconstruct_ancestral_scores` | balanced tree with 1200 tips x 2 scores | 1.4303s | 0.0361s | 39.7x |
| `PhylogeneticOrdination._reconstruct_ancestral_scores` direct traversal and weighted sums | balanced tree with 8192 tips x 2 scores | 0.172621s | 0.052040s | 3.32x |
| `PhylogeneticOrdination._reconstruct_ancestral_scores` prune setup | balanced 8192-tip tree, all tips retained for ordination scores | 0.2648s | 0.0018s | 150.9x |
| `PhylogeneticOrdination._reconstruct_ancestral_scores` all-shared copy skip | balanced 32768-tip tree, all tips retained for ordination scores | 0.358900s | 0.107630s | 3.33x |
| `PhylogeneticOrdination._draw_tree_overlay` direct branch traversal | balanced 32768-tip tree, branch segment/color setup with precomputed node estimates | 0.167162s | 0.048821s | 3.42x |
| `PhylogeneticOrdination._plot_tree` root-distance max helper | 1M node-distance values, identical maximum-distance fallback behavior without temporary list | 0.030450s | 0.018935s | 1.61x |
| `PhylogeneticOrdination._multi_trait_log_likelihood` | 420 taxa SPD VCV x 8 traits | 0.0084s | 0.0005s | 15.6x |
| `PhylogeneticOrdination._multi_trait_log_likelihood_inverse` combined RHS multiply | 420 taxa SPD VCV x 700 traits, inverse fallback path | 0.074771s | 0.033401s | 2.2x |
| `PhylogeneticOrdination._multi_trait_log_likelihood_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV x 10-trait likelihood evaluations, SciPy already warm | 0.066843s | 0.053844s | 1.24x |
| `PhylogeneticOrdination._multi_trait_lambda` lambda-matrix diagonal restoration | 8 / 40 / 260 / 900 / 2000 taxa VCV transform with precomputed diagonal, side-by-side previous `np.fill_diagonal` restoration per lambda evaluation | 0.000008163s / 0.000009956s / 0.000119735s / 0.003229405s / 0.009533993s | 0.000004753s / 0.000001878s / 0.000031953s / 0.002818906s / 0.007760870s | 1.72x / 5.30x / 3.75x / 1.15x / 1.23x |
| `PhylogeneticOrdination` PCA GLS centering/covariance setup | 420 taxa SPD VCV x 8 traits | 0.005137s | 0.000454s | 11.3x |
| `PhylogeneticOrdination` PCA corr-mode diagonal scaling | 420 taxa x 700 traits, synthetic covariance and centered scores | 0.023722s | 0.001051s | 22.6x |
| `PhylogeneticOrdination._run_pca` eigenvalue total variance | 2 / 3 / 4 / 8 / 16 / 32 / 128 / 1024 eigenvalues, side-by-side previous `np.sum(eigenvalues)` | 0.000005367s / 0.000005370s / 0.000004776s / 0.000004082s / 0.000004451s / 0.000004098s / 0.000004246s / 0.000005265s | 0.000001938s / 0.000002333s / 0.000001965s / 0.000002455s / 0.000002093s / 0.000002365s / 0.000002326s / 0.000002696s | 2.77x / 2.30x / 2.43x / 1.66x / 2.13x / 1.73x / 1.83x / 1.95x |
| `PhylogeneticOrdination._center_traits_by_vcv_inverse` combined RHS multiply | 420 taxa SPD VCV x 700 traits, weighted inverse fallback path | 0.043051s | 0.012594s | 3.4x |
| `PhylogeneticOrdination._center_traits_by_vcv_cholesky` residual solve reuse | 120 repeated 420-taxon SPD VCV x 10-trait weighted centering calls, SciPy already warm | 0.056996s | 0.050926s | 1.12x |
| `PhylogeneticOrdination` broadcast centering products | 420 taxa x 700 traits, centered traits plus weighted centered traits | 0.000867s | 0.000635s | 1.37x |
| `PhylogeneticOrdination._print_pca_text_output` batched text output | 100k taxa x 4 PC scores plus 8 trait loadings, captured stdout and identical text | 0.183670s | 0.169344s | 1.08x |
| `PhylogeneticOrdination._print_pca_text_output` score/loadings row conversion | 100k taxa x 4 PC scores plus 8 trait loadings, identical captured text while avoiding per-cell NumPy indexing | 0.602249s | 0.511011s | 1.18x |
| `PhylogeneticOrdination._format_pca_result` row-slice JSON payload | 100k taxa x 4 PC scores plus 8 trait loadings, identical nested payload | 0.072574s | 0.058033s | 1.25x |
| `PhylogeneticOrdination._print_dimreduce_text_output` batched text output | 100k taxa x 2 embedding dimensions, captured stdout and identical text | 0.111761s | 0.099769s | 1.12x |
| `PhylogeneticOrdination._print_dimreduce_text_output` embedding row conversion | 100k taxa x 2 embedding dimensions, identical captured text while avoiding per-cell NumPy indexing | 0.692168s | 0.517560s | 1.34x |
| `PhylogeneticOrdination._format_dimreduce_result` row-list JSON payload | 5k taxa x 2 embedding dimensions, identical nested payload | 0.290704s | 0.096889s | 3.00x |
| `PhylogeneticOrdination._resolve_tree_color_trait` single-pass color file | 300k-row external tree-color TSV, all taxa covered, reconstruction stubbed to isolate parsing | 0.934379s | 0.522889s | 1.79x |
| `phylogenetic_ordination` module import without eager SciPy linalg/optimize | cold process import for phylogenetic-ordination command module | 0.440223s | 0.168008s | 2.6x |
| `phylogenetic_ordination` module import without eager NumPy/PGLS helper | cold subprocess import after lazy NumPy proxy, postponed annotations, and localized PGLS import | 0.168008s | 0.032468s | 5.17x |
| `phylogenetic_ordination` module import without eager pickle/JSON/plot/trait helpers | median cold subprocess import after localizing pickle/PlotConfig and lazy helper wrappers | 0.016052s | 0.005263s | 3.05x |
| `phylogenetic_ordination` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005693s | 0.004226s | 1.35x |
| `PhylogeneticOrdination.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, VCV build, centering, and PCA output mocked | 0.241108s | 0.000032s | 7534.62x |
| `AncestralReconstruction._anc_ml` | balanced tree with 180 tips, continuous trait, CIs enabled | 18.6911s | 0.4921s | 38.0x |
| `AncestralReconstruction._anc_ml` root-path setup | balanced tree with 180 tips, continuous trait, CIs enabled | 0.0409s | 0.0296s | 1.4x |
| `AncestralReconstruction._anc_ml` batched inverse products | balanced tree with 512 tips, continuous trait, CIs enabled | 0.300064s | 0.232316s | 1.3x |
| `AncestralReconstruction._anc_ml` residual inverse-product reuse | 250 repeated 2048-tip `C_inv @ residuals` products after `C_inv @ x` and `C_inv @ 1` are available | 0.334895s | 0.000888s | 377.3x |
| `AncestralReconstruction._anc_ml` direct cross-covariance setup | balanced 512-tip tree, internal-node cross-covariance setup only | 0.257493s | 0.007714s | 33.38x |
| `AncestralReconstruction._fast_anc` direct traversal setup | balanced 8192-tip tree, CI calculations, log-likelihood excluded | 0.480207s | 0.191153s | 2.51x |
| `AncestralReconstruction._fast_anc` cached id lookups | balanced 16384-tip tree, CI calculations with exact estimates/CIs and log-likelihood excluded | 0.124563s | 0.095632s | 1.30x |
| `AncestralReconstruction._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.034317s | 0.024390s | 1.41x |
| `AncestralReconstruction._compute_sigma2_from_contrasts` streaming valid children | balanced 32768-tip tree, postorder contrast scan | 0.031982s | 0.018191s | 1.76x |
| `AncestralReconstruction._compute_log_likelihood` Cholesky RHS solve | balanced 512-tip tree, BM VCV log-likelihood, side-by-side previous explicit inverse path | 0.029255s | 0.004545s | 6.44x |
| `AncestralReconstruction.run` continuous copied-tree prune setup | balanced 32768-tip tree, all tips shared with trait data | 0.0678s | 0.0079s | 8.6x |
| `AncestralReconstruction.run` discrete copied-tree prune setup | balanced 32768-tip tree, all tips shared with state data | 0.0681s | 0.0081s | 8.4x |
| `AncestralReconstruction.run` cached read-only tree setup | balanced 32768-tip cached tree, continuous trait parsing, reconstruction, formatting, and output mocked; protective prune-copy retained | 1.200916s | 0.450668s | 2.66x |
| `AncestralReconstruction.run` continuous all-shared read-only setup | balanced 32768-tip cached tree, continuous trait parsing, reconstruction/formatting/output mocked, ladderize off | 0.273998s | 0.047376s | 5.78x |
| `AncestralReconstruction._label_internal_nodes` direct traversal | balanced 32768-tip tree, five internal-node labeling runs | 0.537023s | 0.071732s | 7.49x |
| `AncestralReconstruction._build_parent_map` direct traversal | balanced 32768-tip tree, five parent-map builds | 0.529001s | 0.067176s | 7.87x |
| `AncestralReconstruction._build_parent_map` unordered child push | balanced 65536-tip tree, parent-map setup, optimized helper baseline | 0.023996s | 0.017968s | 1.34x |
| `AncestralReconstruction._get_descendant_tips` direct clade traversal | balanced 32768-tip tree, descendant-tip lists for every internal node | 0.935165s | 0.117913s | 7.93x |
| `AncestralReconstruction._discrete_marginal_posteriors` transition cache | balanced 2048-tip tree, 3 states, repeated unit branch length | 0.310393s | 0.104882s | 2.96x |
| `AncestralReconstruction._discrete_marginal_posteriors` ndarray normalization sums | 3-state posterior vectors, per-node upward/posterior normalization reductions | 0.000006715s | 0.000002635s | 2.55x |
| `AncestralReconstruction._format_discrete_result` Q-matrix row iteration | 32-state synthetic Q matrix, nested JSON payload | 0.000128s | 0.000086s | 1.48x |
| `AncestralReconstruction._print_text_output` batched continuous table output | synthetic tree with 100k internal-node estimate rows, captured stdout and identical text | 0.188848s | 0.175623s | 1.08x |
| `AncestralReconstruction._print_discrete_text_output` descendant counts | balanced 32768-tip tree, precomputed posteriors, stdout stubbed | 0.323206s | 0.131654s | 2.45x |
| `AncestralReconstruction._print_discrete_text_output` batched table output | balanced 65536-tip tree, 3 states, posterior rows for every internal node, captured stdout and identical text | 0.289856s | 0.224001s | 1.29x |
| `AncestralReconstruction._print_discrete_text_output` percent row formatting | balanced 32768-tip tree, 3 states, posterior rows for every internal node, captured stdout and identical text, side-by-side previous f-string row concatenation comparison | 0.106752s | 0.079717s | 1.34x |
| `AncestralReconstruction._format_result` descendant-name cache | balanced 16384-tip tree, all internal node result entries | 0.120305s | 0.084881s | 1.42x |
| `AncestralReconstruction._collect_descendant_tip_names` ordered-range concatenation | balanced 32768-tip tree, sorted tip labels | 0.085289s | 0.029060s | 2.94x |
| `AncestralReconstruction._plot_contmap` rectangular setup | balanced 32768-tip tree, precomputed estimates and trait values | 0.9713s | 0.1385s | 7.0x |
| `AncestralReconstruction._plot_contmap` circular setup | balanced 32768-tip tree, precomputed estimates and trait values | 1.1534s | 0.2223s | 5.2x |
| `AncestralReconstruction._plot_contmap` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.047892s | 0.039938s | 1.20x |
| `AncestralReconstruction._plot_contmap` circular coordinate clade-list reuse | balanced 32768-tip tree, node x positions, parent map, preorder list, and tips already available | 0.057251s | 0.044612s | 1.28x |
| `AncestralReconstruction._plot_contmap` estimate value range helper | 1M plotted node/tip estimate values, identical min/max range without temporary list | 0.039032s | 0.030498s | 1.28x |
| `AncestralReconstruction._plot_contmap` rectangular gradient branch rendering | balanced 512-tip tree, 50 color segments per branch, real Matplotlib Agg branch render | 10.877298s | 1.865517s | 5.83x |
| `AncestralReconstruction._plot_contmap` rectangular CI overlay rendering | 2048 CI nodes, real Matplotlib Agg CI bars and point estimates | 8.205795s | 0.036551s | 224.50x |
| `AncestralReconstruction._plot_contmap` circular CI overlay rendering | 512 CI nodes, real Matplotlib Agg tangential CI bars and point estimates | 0.793536s | 0.008322s | 95.35x |
| `AncestralReconstruction._plot_discrete_asr` rectangular setup | balanced 32768-tip tree, precomputed node posteriors and tip states | 0.8363s | 0.1233s | 6.8x |
| `AncestralReconstruction._plot_discrete_asr` circular setup | balanced 32768-tip tree, precomputed node posteriors and tip states | 0.7626s | 0.1767s | 4.3x |
| `AncestralReconstruction._plot_discrete_asr` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.047892s | 0.039938s | 1.20x |
| `AncestralReconstruction._plot_discrete_asr` circular coordinate clade-list reuse | balanced 32768-tip tree, node x positions, parent map, preorder list, and tips already available | 0.057251s | 0.044612s | 1.28x |
| `AncestralReconstruction._plot_discrete_asr` rectangular base branch rendering | balanced 2048-tip tree, precomputed node posteriors and tip states, real Matplotlib Agg branch render | 1.247995s | 0.051896s | 24.05x |
| `AncestralReconstruction._plot_discrete_asr` rectangular node pie rendering | 2048 internal nodes, 3 posterior states per node, real Matplotlib Agg wedge render | 2.924065s | 0.275539s | 10.61x |
| `AncestralReconstruction._plot_discrete_asr` circular node pie rendering | 2048 internal nodes, 3 posterior states per node, real Matplotlib Agg wedge render | 2.652814s | 0.287181s | 9.24x |
| `AncestralReconstruction._plot_discrete_asr` rectangular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 1.287406s | 0.046800s | 27.51x |
| `AncestralReconstruction._plot_discrete_asr` circular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 0.626749s | 0.028314s | 22.14x |
| `AncestralReconstruction._parse_single_trait_data` streaming parser | 500k two-column continuous trait rows with comments/blanks, all taxa shared | 0.467832s | 0.436189s | 1.07x |
| `AncestralReconstruction._parse_single_trait_data` all-shared parser fast path | 500k two-column continuous trait rows with comments/blanks, all taxa shared | 0.434114s | 0.235606s | 1.84x |
| `AncestralReconstruction._parse_multi_trait_data` single-pass parser | 300k-row multi-trait TSV, 3 numeric trait columns, 100k shared taxa | 0.366382s | 0.299851s | 1.22x |
| `AncestralReconstruction._parse_multi_trait_data` all-shared parser fast path | 300k-row multi-trait TSV, 3 numeric trait columns, all taxa shared | 0.273898s | 0.171096s | 1.60x |
| `AncestralReconstruction._parse_discrete_trait_data_single` streaming parser | 500k two-column categorical trait rows with comments/blanks, all taxa shared | 0.384133s | 0.375407s | 1.02x |
| `AncestralReconstruction._parse_discrete_trait_data_single` all-shared parser fast path | 500k two-column categorical trait rows with comments/blanks, all taxa shared | 0.411670s | 0.214221s | 1.92x |
| `AncestralReconstruction._parse_discrete_trait_data_multi` streaming parser | 300k-row categorical TSV, 5 state columns, all taxa shared | 0.263099s | 0.257110s | 1.02x |
| `AncestralReconstruction._parse_discrete_trait_data_multi` all-shared parser fast path | 300k-row categorical TSV, 5 state columns, all taxa shared | 0.242633s | 0.132431s | 1.83x |
| `ancestral_reconstruction` module import via lazy discrete SciPy helpers | cold process import for ancestral-reconstruction command module | 0.475903s | 0.392662s | 1.2x |
| `ancestral_reconstruction` module import without eager NumPy | cold subprocess import after lazy NumPy proxies in command and discrete helper | 0.081936s | 0.034135s | 2.40x |
| `ancestral_reconstruction` module import without eager runtime helpers | cold subprocess import after lazy pickle, VCV, discrete-model, plot, circular-layout, and color-helper proxies | 0.036479s | 0.027797s | 1.31x |
| `ancestral_reconstruction` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.007307s | 0.005960s | 1.23x |
| `ancestral_reconstruction` module import without `typing` startup | median cold subprocess import after removing annotation-only typing import under postponed annotations | 0.006027s | 0.004341s | 1.39x |
| `ancestral_reconstruction` module import without eager `heapq` | median cold subprocess import after lazy `heapq.merge` proxy for descendant-name merging | 0.005179s | 0.004457s | 1.16x |
| `FaithsPD.calculate_faiths_pd` | balanced tree with 2500 tips, 1500-taxon community, `include_root=False` | 0.8434s | 0.4622s | 1.8x |
| `FaithsPD.calculate_faiths_pd` one-pass traversal | balanced tree with 2500 tips, 1500-taxon community, `include_root=False` | 0.5249s | 0.0213s | 24.6x |
| `FaithsPD.calculate_faiths_pd` all-tip community fast path | balanced tree with 32768 tips, all tips selected, `include_root=False` | 0.060900s | 0.030831s | 2.0x |
| `FaithsPD.calculate_faiths_pd` default-root parent bookkeeping | balanced tree with 32768 tips, 1500-taxon community, `include_root=True` | 0.050422s | 0.039901s | 1.26x |
| `FaithsPD.calculate_faiths_pd` default-root selected path walk | balanced tree with 32768 tips, 1500-taxon sparse community, `include_root=True` | 0.050974s | 0.045142s | 1.13x |
| `FaithsPD.calculate_faiths_pd` default-root inline branch sum | balanced 32768-tip tree, 1500-taxon sparse community, `include_root=True` | 0.031601s | 0.026425s | 1.20x |
| `FaithsPD.calculate_faiths_pd` default-root overlapping path cutoff | pectinate 12000-tip tree, every 12th taxon selected, `include_root=True` | 0.522531s | 0.014042s | 37.21x |
| `FaithsPD.calculate_faiths_pd` selected-tip lookup | balanced 32768-tip tree, every 20th taxon selected, `include_root=True` | 0.030370s | 0.027354s | 1.11x |
| `FaithsPD.calculate_faiths_pd` initial child push | balanced 32768-tip tree, every 20th taxon selected, `include_root=True`, side-by-side previous `reversed(children)` traversal | 0.043449s | 0.032935s | 1.32x |
| `FaithsPD.calculate_faiths_pd` binary selected-count aggregation | balanced 8192-tip tree, every third taxon selected, `include_root=False`, side-by-side previous generator-sum child count | 0.020002s | 0.015337s | 1.30x |
| `FaithsPD.calculate_faiths_pd` single-tip exclude-root skip | balanced 32768-tip tree, one selected taxon, `include_root=False` | 0.047379s | 0.029310s | 1.62x |
| `FaithsPD._load_taxa` order-preserving dedupe | 400k taxa rows plus blanks, 180k unique taxa, full file load path | 0.075677s | 0.066871s | 1.13x |
| `FaithsPD.run` cached read-only tree path | balanced 16384-tip cached tree, 2048-taxon community, taxa/output mocked | 0.138387s | 0.022596s | 6.12x |
| `faiths_pd` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006675s | 0.005228s | 1.28x |
| `faiths_pd` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.007869s | 0.003592s | 2.19x |
| `faiths_pd` module import without eager file helper | median cold subprocess import, interleaved lazy import vs eager-equivalent `phykit.helpers.files` preload | 0.029179s | 0.028322s | 1.03x |
| `ParsimonyScore._fitch_parsimony` | balanced tree with 256 tips x 2000 sites, alphabet `ACGT-?NX` | 1.8456s | 0.0488s | 37.8x |
| `ParsimonyScore._fitch_parsimony` byte lookup masks | balanced tree with 256 tips x 2000 sites, alphabet `ACGT-?NX` | 0.0484s | 0.0157s | 3.1x |
| `ParsimonyScore._fitch_parsimony` small-alignment path | balanced 32768-tip tree x 4 sites, alphabet `ACGT-?N` | 0.124124s | 0.044296s | 2.8x |
| `ParsimonyScore._fitch_parsimony` small repeated-sequence cache | balanced 32768-tip tree x 4 sites, four repeated sequence patterns | 0.081057s | 0.050042s | 1.6x |
| `ParsimonyScore._fitch_parsimony` reverse-preorder postorder helper | balanced 32768-tip tree x 4 sites, four repeated sequence patterns | 0.042540s | 0.035535s | 1.20x |
| `ParsimonyScore._fitch_parsimony_sets` postorder reuse | balanced 128-tip tree x 300 sites, 80 observed Unicode states forcing set fallback | 0.202495s | 0.033781s | 5.99x |
| `ParsimonyScore._fitch_parsimony` identical-sequence shortcut | 2048-tip balanced tree, 10k identical DNA sites, identical zero total and per-site scores | 4.821393s | 0.000458s | 10527.1x |
| `ParsimonyScore._fitch_parsimony` non-verbose score-only path | 4-tip tree x 500k sites, vectorized Fitch path, identical total score with/without returned per-site scores | 0.445905s | 0.248677s | 1.79x |
| `IndependentContrasts._compute_pic` | balanced tree with 2500 tips, continuous trait | 0.0219s | 0.0039s | 5.7x |
| `FitContinuous._concentrated_ll` | 420 taxa SPD VCV, single continuous trait | 0.0064s | 0.0005s | 12.7x |
| `FitContinuous._concentrated_ll_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV concentrated likelihood evaluations, SciPy already warm | 0.057255s | 0.042199s | 1.36x |
| `FitContinuous._fit_lambda` lambda-matrix diagonal restoration | 8 / 40 / 260 / 900 / 2000 taxa VCV transform with precomputed diagonal, side-by-side previous `np.fill_diagonal` restoration per lambda evaluation | 0.000008163s / 0.000009956s / 0.000119735s / 0.003229405s / 0.009533993s | 0.000004753s / 0.000001878s / 0.000031953s / 0.002818906s / 0.007760870s | 1.72x / 5.30x / 3.75x / 1.15x / 1.23x |
| `FitContinuous._build_root_to_tip_paths` | balanced 32768-tip tree, paths for 16384 modeled tips | 2.8907s | 0.0613s | 47.1x |
| `fit_continuous` module import without eager SciPy linalg/optimize | cold process import for continuous-model-fitting command module | 0.497876s | 0.214771s | 2.3x |
| `fit_continuous` module import without eager NumPy/PGLS helper | cold subprocess import after lazy NumPy proxy, postponed annotations, and localized PGLS import | 0.081247s | 0.025834s | 3.15x |
| `fit_continuous` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006266s | 0.005033s | 1.25x |
| `fit_continuous` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005514s | 0.001744s | 3.16x |
| `pgls_utils.max_lambda` | balanced tree with 2500 tips | 0.0108s | 0.0005s | 19.9x |
| `pgls_utils.max_lambda` streaming tip min/max | balanced 262144-tip ultrametric tree, direct standard-tree lambda upper-bound helper, side-by-side previous tip-height list helper | 0.060381s | 0.050485s | 1.20x |
| `pgls_utils.pgls_log_likelihood` | 420 taxa SPD VCV x 3-column design matrix | 0.0064s | 0.0006s | 10.7x |
| `pgls_utils.pgls_log_likelihood` cached SciPy linalg wrappers | 2k repeated 8-taxon SPD VCV x 2-column design matrix calls, SciPy already warm | 0.064135s | 0.061040s | 1.05x |
| `pgls_utils.pgls_log_likelihood` combined RHS solve | 160 repeated 420-taxon SPD VCV x 4-column design matrix likelihood evaluations, SciPy already warm | 0.094840s | 0.075566s | 1.26x |
| Cholesky logdet diagonal reductions | 4 / 32 / 256 / 5k / 100k diagonal values, shared likelihood idiom, side-by-side previous `np.sum(np.log(diag))` wrapper | 3.604498s / 1.149306s / 1.004864s / 0.182699s / 0.198437s | 1.484793s / 0.463736s / 0.355035s / 0.090258s / 0.105189s | 2.43x / 2.48x / 2.83x / 2.02x / 1.89x |
| `pgls_utils.estimate_lambda` | 260 taxa SPD VCV x 3-column design matrix, `max_lam=1.0` | 0.7152s | 0.0676s | 10.6x |
| `pgls_utils.estimate_lambda` lambda-matrix diagonal restoration | 8 / 40 / 260 / 900 / 2000 taxa VCV transform with precomputed diagonal, side-by-side previous `np.fill_diagonal` restoration per lambda evaluation | 0.000008163s / 0.000009956s / 0.000119735s / 0.003229405s / 0.009533993s | 0.000004753s / 0.000001878s / 0.000031953s / 0.002818906s / 0.007760870s | 1.72x / 5.30x / 3.75x / 1.15x / 1.23x |
| `pgls_utils.fit_gls` combined RHS multiply | 900 taxa SPD inverse x 13-column design matrix | 0.001092s | 0.000848s | 1.3x |
| `pgls_utils.fit_gls` preallocated combined RHS | 300 repeated 250-taxon SPD inverse x 4-column design matrix fits, side-by-side previous `np.column_stack((X, y))` RHS assembly | 0.084559s | 0.019053s | 4.44x |
| `pgls_utils.fit_gls` RHS-first coefficient multiply | 900 taxa SPD inverse x 13-column design matrix, coefficient solve step | 0.000015s | 0.000007s | 2.27x |
| `pgls_utils.fit_gls` coefficient covariance identity solve | 900 taxa SPD inverse x 13-column design matrix, side-by-side previous explicit coefficient-matrix inverse | 0.006303s | 0.006034s | 1.04x |
| `pgls_utils` module import without eager SciPy linalg/optimize | cold process import for shared PGLS utilities | 0.431510s | 0.120650s | 3.6x |
| `pgls_utils` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.060729s | 0.002053s | 29.58x |
| `pgls_utils` module import without `typing` startup | median cold subprocess import after replacing annotation-only `Tuple` aliases with built-in postponed annotations | 0.023289s | 0.020680s | 1.13x |
| `PhylogeneticRegression._fit_model` | 420 taxa SPD VCV x 3-predictor design matrix | 0.0081s | 0.0007s | 11.5x |
| `PhylogeneticRegression._fit_model` combined Cholesky RHS solve | 120 repeated 420-taxon SPD VCV x 3-predictor design matrix fits, SciPy already warm | 0.071721s | 0.048587s | 1.48x |
| `PhylogeneticRegression._fit_model` combined normal-equation RHS | 420 taxa SPD VCV x 6-predictor design matrix, side-by-side previous explicit normal-matrix inverse | 0.000440s | 0.000406s | 1.08x |
| `phylogenetic_regression` module import without `scipy.stats` | cold process import for PGLS command module | 0.642346s | 0.445118s | 1.4x |
| `phylogenetic_regression` module import without eager SciPy linalg | cold process import for PGLS command module | 0.317047s | 0.173800s | 1.8x |
| `phylogenetic_regression` module import without eager NumPy/PGLS helper | cold subprocess import after lazy NumPy proxy, postponed annotations, and localized PGLS imports | 0.085897s | 0.025865s | 3.32x |
| `phylogenetic_regression` module import without eager JSON/trait parsing helpers | median cold subprocess import after lazy helper wrappers | 0.008477s | 0.005541s | 1.53x |
| `phylogenetic_regression` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005956s | 0.003949s | 1.51x |
| `phylogenetic_regression` cached lazy SciPy special helpers | 200k chunked Student-t and 200k F survival p-value helper calls, SciPy already warm | 0.291163s | 0.238297s | 1.22x |
| `phylogenetic_regression` single-predictor F survival helper | 200k model F p-values with numerator df=1 | 0.288397s | 0.201494s | 1.43x |
| `PhylogeneticRegression.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, VCV build, model fit, p-values, and output mocked | 0.564282s | 0.000063s | 8956.86x |
| `PhylogeneticRegression._print_text_output` batched coefficient rows | 100k PGLS coefficient rows, captured stdout and identical text | 0.113786s | 0.103553s | 1.10x |
| `PhylogeneticRegression._print_text_output` zip-based coefficient rows | 100k PGLS coefficient rows, side-by-side previous index lookups with identical lines | 0.174629s | 0.153927s | 1.13x |
| `PhylogeneticRegression._format_result` zip-based JSON mappings | 250k taxon residual/fitted rows plus 4k coefficients, side-by-side previous index lookups | 0.077218s | 0.062417s | 1.24x |
| `PhylogeneticRegression._format_result` bulk residual/fitted conversion | 250k taxon residual/fitted rows plus 4k coefficients, side-by-side previous per-item NumPy scalar conversion | 0.137516s | 0.109884s | 1.25x |
| `PhyloPath._fit_gls_from_vcv` | 420 taxa SPD VCV x 3-predictor design matrix, post-lambda GLS step | 0.0078s | 0.0012s | 6.4x |
| `PhyloPath._fit_gls_from_vcv` combined Cholesky RHS solve | 120 repeated 420-taxon SPD VCV x 3-predictor design matrix fits, SciPy already warm | 0.053644s | 0.043705s | 1.23x |
| `PhyloPath._fit_gls_from_vcv` combined normal-equation RHS | 420 taxa SPD VCV x 6-predictor design matrix, side-by-side previous explicit normal-matrix inverse | 0.000462s | 0.000388s | 1.19x |
| `PhyloPath._dsep_test` | 260 taxa SPD VCV x 3 conditioning/predictor columns, including lambda estimation | 0.0783s | 0.0446s | 1.8x |
| `PhyloPath._design_matrix_from_z` preallocated PGLS design matrix | 10k repeated 260-taxon x 4-predictor design matrix builds, side-by-side previous `np.column_stack([ones] + predictors)` setup | 0.113720s | 0.021587s | 5.27x |
| `phylo_path` module import without `scipy.stats` | cold process import for path-analysis command module | 0.616865s | 0.486959s | 1.27x |
| `phylo_path` module import without eager SciPy linalg | cold process import for path-analysis command module | 0.326131s | 0.157336s | 2.1x |
| `phylo_path` module import without eager NumPy/PGLS helper | cold subprocess import after lazy NumPy proxy, postponed annotations, and localized PGLS imports | 0.157336s | 0.031335s | 5.02x |
| `phylo_path` module import without eager JSON/plot config helpers | median cold subprocess import after lazy JSON wrapper and localized `PlotConfig` import | 0.011623s | 0.004829s | 2.41x |
| `phylo_path` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006188s | 0.005581s | 1.11x |
| `phylo_path` cached lazy SciPy special helpers | 200k chi-square and 200k Student-t p-value helper calls, SciPy already warm | 0.419375s | 0.344857s | 1.22x |
| `phylo_path` even-df chi-square survival helper | 200k Fisher's C chi-square survival probabilities at df=6 | 0.149460s | 0.054022s | 2.77x |
| `PhyloPath` Fisher's C scalar p-value logs | 16 p-values in one basis set, side-by-side previous scalar `np.log` generator | 0.000004624s | 0.000001691s | 2.73x |
| `PhyloPath` model-weight relative likelihoods | 1000 model deltas, side-by-side previous scalar `np.exp` list loop | 0.000241713s | 0.000063842s | 3.79x |
| `PhyloPath._print_text` batched report output | 100k ranked model rows plus 100k path coefficient rows, captured stdout and identical text | 0.299066s | 0.269681s | 1.11x |
| `PhyloPath._print_text` row-template formatting | 100k ranked model rows plus 100k path coefficient rows, captured stdout and identical text, side-by-side previous f-string row formatter comparison | 0.286284s | 0.260257s | 1.10x |
| `PhyloPath._print_text` percent row formatting | 100k ranked model rows plus 100k path coefficient rows, captured stdout and identical text, side-by-side previous `.format()` row formatter comparison | 0.267401s | 0.174840s | 1.53x |
| `PhyloPath._plot_dag` node circle rendering | 4096 data-coordinate node circles, real Matplotlib Agg render | 3.024579s | 0.155128s | 19.50x |
| `PhyloPath._parse_trait_file` streaming valid-row parser | 300k-row multi-trait TSV, 3 numeric trait columns, 100k shared taxa | 0.654959s | 0.457753s | 1.43x |
| `PhyloPath._parse_trait_file` all-shared parser fast path | 300k-row multi-trait TSV, 3 numeric trait columns, all taxa shared | 0.449593s | 0.311079s | 1.45x |
| `PhyloPath._parse_models_file` streaming parser | 300k candidate path models, three edges per model, comments/blanks included | 1.162768s | 0.905995s | 1.28x |
| `PhyloPath._is_dag` queue cursor | 40k-variable wide acyclic DAG with 20k roots and 20k dependent nodes | 0.095966s | 0.011607s | 8.27x |
| `PhyloPath._topological_order` queue cursor | 40k-variable wide acyclic DAG with 20k roots and 20k dependent nodes | 0.099320s | 0.014319s | 6.94x |
| `PhyloPath._basis_set` precomputed parent sets | 350-variable sparse DAG, 60,782 d-separation basis statements | 1.211617s | 0.057178s | 21.19x |
| `PhyloPath._model_average` edge-indexed coefficient accumulation | 8k candidate models, 4k unique path coefficients, 12 coefficients per model | 1.562250s | 0.034478s | 45.31x |
| `PhyloPath._model_average` combined weight/coefficient sums | 12 coefficient entries for one edge, side-by-side previous three generator passes | 0.000007119s | 0.000002185s | 3.26x |
| `PhyloPath.run` model-variable trait-name index resolution | 40k parsed trait columns, 1200 model variables, first duplicate index preserved | 0.684638s | 0.002820s | 242.76x |
| `PhyloPath.run` multi-variable trait setup | 120k taxa x 12 parsed trait columns, five model variables extracted and copied | 0.062919s | 0.038572s | 1.63x |
| `PhyloPath.run` selected-variable standardization | 120k taxa x 64 parsed trait columns, 40 model variables, side-by-side previous raw-data dict and duplicate mean pass | 0.030258s | 0.025553s | 1.18x |
| `PhyloPath.run` selected-variable ndarray reductions | 120k taxa x 64 parsed trait columns, 40 model variables, side-by-side previous `np.mean`/`np.std` calls on copied columns | 0.108664s | 0.048551s | 2.24x |
| `PhyloPath.run` cached read-only setup | balanced 32768-tip cached tree, three traits for every tip, VCV/d-sep/coefficient/output mocked | 0.311719s | 0.074181s | 4.20x |
| `PhylogeneticSignal._blombergs_k` | 420 taxa SPD VCV, 1000 permutations | 0.1134s | 0.0338s | 3.4x |
| `PhylogeneticSignal._blombergs_k` Cholesky inverse construction | 450 taxa SPD VCV, 16 seeded permutations, side-by-side previous explicit inverse | 0.007445s | 0.003480s | 2.14x |
| `PhylogeneticSignal._kmult_permutations` | 320 taxa x 5 traits SPD VCV, 1000 permutations | 0.143512s | 0.049227s | 2.9x |
| `PhylogeneticSignal`/`NetworkSignal` permutation p-value counts | 1M permutation statistics, side-by-side previous `np.mean(permutations >= observed)` reduction | 0.000991s | 0.000390s | 2.54x |
| `PhylogeneticSignal._compute_r2_phylo` | 420 taxa SPD VCV, single continuous trait | 0.0047s | 0.0015s | 3.2x |
| `PhylogeneticSignal._compute_r2_phylo` combined RHS solve | 120 repeated 420-taxon SPD VCV R2 effect-size evaluations, SciPy already warm | 0.052471s | 0.039279s | 1.34x |
| `PhylogeneticSignal._compute_r2_phylo` white-noise variance reduction | 260 / 420 / 1000 trait values, side-by-side previous `np.var(x)` wrapper | 0.000006042s / 0.000006084s / 0.000006667s | 0.000005459s / 0.000005541s / 0.000006083s | 1.11x / 1.10x / 1.10x |
| `PhylogeneticSignal._log_likelihood` | 420 taxa SPD VCV, single continuous trait | 0.0074s | 0.0005s | 14.5x |
| `PhylogeneticSignal._log_likelihood` combined Cholesky RHS solve | 120 repeated 420-taxon SPD VCV likelihood evaluations, SciPy already warm | 0.059433s | 0.040428s | 1.47x |
| `PhylogeneticSignal._log_likelihood` cached SciPy Cholesky wrappers | 120 repeated 420-taxon SPD VCV likelihood evaluations, SciPy already warm, side-by-side previous import-on-call wrappers | 0.040999s | 0.039617s | 1.03x |
| `PhylogeneticSignal._pagels_lambda` | 260 taxa SPD VCV, single continuous trait, `max_lambda=1.0` | 0.8285s | 0.0721s | 11.5x |
| `PhylogeneticSignal._pagels_lambda`/`NetworkSignal._pagels_lambda` lambda-matrix diagonal restoration | 8 / 40 / 260 / 900 / 2000 taxa VCV transform with precomputed diagonal, side-by-side previous `np.fill_diagonal` restoration per lambda evaluation | 0.000008163s / 0.000009956s / 0.000119735s / 0.003229405s / 0.009533993s | 0.000004753s / 0.000001878s / 0.000031953s / 0.002818906s / 0.007760870s | 1.72x / 5.30x / 3.75x / 1.15x / 1.23x |
| `PhylogeneticSignal._parse_trait_file` streaming two-column parser | 500k two-column trait rows with comments/blanks, all taxa shared, randomized old/new measurement order | 0.457297s | 0.441401s | 1.04x |
| `PhylogeneticSignal._parse_trait_file` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.436331s | 0.239903s | 1.82x |
| `PhylogeneticSignal._parse_trait_file` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.235394s | 0.229988s | 1.02x |
| `PhylogeneticSignal.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, VCV/statistic calculation, and output mocked | 0.425021s | 0.000195s | 2181.46x |
| `phylogenetic_signal` module import without `scipy.stats` | cold process import for Pagel-lambda-capable command module | 0.608692s | 0.409510s | 1.5x |
| `phylogenetic_signal` module import without eager SciPy linalg/optimize | cold process import for phylogenetic-signal command module | 0.485930s | 0.214970s | 2.3x |
| `phylogenetic_signal` module import without eager NumPy/PGLS helper | cold subprocess import after lazy NumPy proxy, postponed annotations, and localized PGLS import | 0.079401s | 0.026265s | 3.02x |
| `phylogenetic_signal` module import without eager JSON/trait parsing helpers | median cold subprocess import after lazy helper wrappers | 0.007508s | 0.005705s | 1.32x |
| `phylogenetic_signal` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005558s | 0.004152s | 1.34x |
| `PhyloImpute._estimate_complete_case_stats` | 420 taxa SPD VCV x 8 traits | 0.0098s | 0.0005s | 18.0x |
| `PhyloImpute._estimate_complete_case_stats_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV x 8-trait complete-case stat estimates, SciPy already warm | 0.053690s | 0.039704s | 1.35x |
| `PhyloImpute._estimate_complete_case_stats_inverse` | 420 taxa SPD VCV x 700 traits, inverse fallback path | 0.062003s | 0.021273s | 2.9x |
| `PhyloImpute._impute_taxon` | 260 taxa SPD VCV x 8 traits, 3 missing traits for one taxon | 0.0148s | 0.0044s | 3.4x |
| `PhyloImpute._impute_taxon_cholesky` cached trait context | 260 taxa x 80 traits, 60 missing traits for one taxon | 0.037352s | 0.025250s | 1.5x |
| `PhyloImpute._impute_taxon_cholesky` combined phylogenetic RHS solve | 120 repeated 260-taxon SPD VCV x 8 traits, 3 missing traits for one taxon, SciPy already warm | 0.150362s | 0.144876s | 1.04x |
| `PhyloImpute._impute_taxon_inverse` cached trait context | 180 taxa x 80 traits, 60 missing traits for one taxon | 0.265825s | 0.187101s | 1.4x |
| `PhyloImpute._other_taxon_indices` vectorized index construction | 140 repeated all-but-one index vectors over 40k taxa, side-by-side previous list comprehension plus `np.asarray` setup | 2.176238s | 0.010129s | 214.86x |
| `PhyloImpute.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, VCV build, complete-case stats, and output mocked | 0.395725s | 0.000771s | 513.26x |
| `PhyloImpute.run` missing-trait row scan | 80k taxa x 80 traits, 2% rows with 3 missing traits, side-by-side previous per-row `np.isnan` scans and list conversions | 0.124522s | 0.010385s | 11.99x |
| `PhyloImpute.run` complete-case mask count | 80k-taxon complete-case boolean mask after missing-trait scan | 0.000053595s | 0.000005332s | 10.05x |
| `PhyloImpute._build_data_matrix` direct NumPy construction | 250k taxa x 12 traits with ordered row lookup and preserved `NaN` missing values | 0.479432s | 0.177250s | 2.70x |
| `PhyloImpute._parse_trait_file_with_na` single-pass parser | 200k-row multi-trait TSV, 8 numeric/NA trait columns, 100k shared taxa | 1.289345s | 1.186224s | 1.09x |
| `PhyloImpute._parse_trait_file_with_na` off-tree row retention skip | 200k-row multi-trait TSV, 8 numeric/NA trait columns, 100k shared taxa and 100k off-tree taxa | 0.780822s | 0.605037s | 1.29x |
| `PhyloImpute._parse_trait_file_with_na` numeric-row fast path | 200k-row multi-trait TSV, 8 numeric/NA trait columns, 100k shared taxa and 100k off-tree taxa | 0.605037s | 0.458054s | 1.32x |
| `PhyloImpute._parse_trait_file_with_na` all-shared parser fast path | 200k-row multi-trait TSV, 8 numeric/NA trait columns, all taxa shared | 0.510872s | 0.385042s | 1.33x |
| `PhyloImpute._parse_trait_file_with_na` off-tree numeric validation skip | 200k-row multi-trait TSV, 8 numeric/NA trait columns, 100k shared taxa and 100k off-tree taxa, side-by-side previous off-tree value-list allocation | 1.006604s | 0.912679s | 1.10x |
| `PhyloImpute._subset_to_shared_taxa` cached membership set | 20k ordered taxa, 10k shared taxa, NumPy matrix row subset | 3.554433s | 0.002021s | 1758.5x |
| `PhyloImpute._write_output_tsv` row-slice TSV formatting | 30k taxa x 80 imputed trait matrix, identical six-decimal TSV text | 0.762587s | 0.667688s | 1.14x |
| `PhyloImpute._write_output_tsv` chunked `savetxt` formatting | 30k taxa x 80 imputed trait matrix, identical six-decimal TSV text | 0.644563s | 0.396832s | 1.62x |
| `PhyloImpute._print_text` batched imputed rows | 100k imputed missing-value rows, captured stdout and identical text | 0.141269s | 0.128685s | 1.10x |
| `PhyloImpute._print_text` combined width scan | 100k imputed missing-value rows, captured stdout and identical text, side-by-side previous two-generator width scan | 0.135697s | 0.126773s | 1.07x |
| `phylo_impute` module import without eager SciPy linalg | cold process import for phylogenetic-imputation command module | 0.319911s | 0.148381s | 2.2x |
| `phylo_impute` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.148381s | 0.025825s | 5.75x |
| `phylo_impute` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.006538s | 0.004955s | 1.32x |
| `phylo_impute` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.037504s | 0.031000s | 1.21x |
| `TraitCorrelation._compute_correlation_matrices` | 420 taxa SPD VCV x 10 traits | 0.0099s | 0.0024s | 4.2x |
| `TraitCorrelation._compute_correlation_matrices_inverse` | 420 taxa SPD VCV x 700 traits, inverse fallback path | 0.738626s | 0.133975s | 5.5x |
| `TraitCorrelation._correlation_and_p_values` | 420 taxa x 900-trait synthetic covariance matrix | 11.083108s | 0.092105s | 120.3x |
| `TraitCorrelation._correlation_and_p_values` upper-triangle p-values | 420 taxa x 900-trait synthetic covariance matrix | 0.078657s | 0.045366s | 1.73x |
| `TraitCorrelation._t_two_tailed_p_values` cached SciPy special helper | 20k repeated p-value vector calls after SciPy warm import, 190 t-statistics | 0.376058s | 0.359423s | 1.05x |
| `TraitCorrelation._plot_heatmap` significance star rendering | 64 x 64 trait heatmap, 4032 off-diagonal significant cells, real Matplotlib Agg star render | 0.672074s | 0.023131s | 29.06x |
| `TraitCorrelation._draw_significance_stars` no-significance early exit | 2500 x 2500 all-nonsignificant p-value matrix, no scatter calls, side-by-side previous off-diagonal mask setup | 0.067631s | 0.002929s | 23.09x |
| `TraitCorrelation._draw_significance_stars` sparse significant coordinates | 2500 x 2500 p-value matrix, sparse symmetric significant pairs, fake scatter calls, side-by-side previous full-mask grouping | 0.218060s | 0.065593s | 3.32x |
| `TraitCorrelation._draw_significance_stars` flat-index sparse coordinates | 2500 x 2500 p-value matrix, sparse symmetric significant pairs, fake scatter calls, side-by-side previous `np.nonzero` coordinate path | 0.100309s | 0.012676s | 7.91x |
| `TraitCorrelation._draw_significance_stars` min-guard no-hit return | 2500 x 2500 all-nonsignificant p-value matrix, no scatter calls, side-by-side previous coordinate extraction path | 0.110337s | 0.001154s | 95.65x |
| `TraitCorrelation._print_text` joined matrix rows | 700-trait correlation matrix, captured stdout and identical text | 0.299372s | 0.264128s | 1.13x |
| `TraitCorrelation._print_text` inline star thresholds | 700-trait correlation matrix, captured stdout and identical text, optimized text baseline | 0.265591s | 0.231011s | 1.15x |
| `TraitCorrelation._print_text` local row-cell joins | 700-trait correlation matrix, captured stdout and identical text, side-by-side previous row append comparison | 0.232587s | 0.204576s | 1.14x |
| `TraitCorrelation._print_text` upper-triangle significance count | 2500-trait p-value matrix with sparse significant pairs, side-by-side previous `np.triu_indices` count path | 0.105923s | 0.005865s | 18.06x |
| `TraitCorrelation._print_json` row-slice payload construction | 700-trait correlation and p-value matrices, identical rounded payload | 0.314889s | 0.258569s | 1.22x |
| `TraitCorrelation._print_json` vectorized matrix rounding and upper-triangle pairs | 700-trait correlation and p-value matrices, print_json stubbed, identical rounded payload | 0.260765s | 0.014778s | 17.65x |
| `TraitCorrelation._build_significant_pairs` vectorized value gather | 2500-trait correlation and p-value matrices, 390385 significant upper-triangle pairs with identical JSON rows | 0.859738s | 0.579571s | 1.48x |
| `TraitCorrelation._build_significant_pairs` sparse row scan | 2500-trait correlation and p-value matrices, 3095 significant upper-triangle pairs with identical JSON rows | 0.024440s | 0.012631s | 1.93x |
| `trait_correlation` module import without `scipy.stats` | cold process import for trait-correlation command module | 0.621309s | 0.338620s | 1.8x |
| `trait_correlation` module import without eager SciPy linalg | cold process import for trait-correlation command module | 0.310193s | 0.162593s | 1.9x |
| `trait_correlation` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.162593s | 0.031299s | 5.19x |
| `trait_correlation` module import without eager JSON/plot/trait helpers | median cold subprocess import after localizing PlotConfig and lazy helper wrappers | 0.015338s | 0.007276s | 2.11x |
| `trait_correlation` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005491s | 0.003498s | 1.57x |
| `TraitCorrelation.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, VCV build, correlations, plotting, and output mocked | 0.231448s | 0.000038s | 6090.74x |
| `TraitCorrelation._compute_correlation_matrices_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV x 10-trait correlation computations, SciPy already warm | 0.065236s | 0.051288s | 1.27x |
| `NetworkSignal._validate_tree` direct traversal | balanced 65536-tip tree, min-tip, branch-length, and polytomy checks | 0.5475s | 0.0253s | 21.6x |
| `NetworkSignal._validate_tree` unordered validation scan | balanced 131072-tip tree, min-tip, branch-length, and polytomy checks, optimized helper baseline | 0.030270s | 0.020724s | 1.46x |
| `NetworkSignal._compute_network_vcv` vectorized covariance assembly | balanced 2048-tip tree DAG with one hybrid edge | 2.500761s | 0.045117s | 55.4x |
| `NetworkSignal._blombergs_k` | 420 taxa SPD VCV, 1000 permutations | 0.1238s | 0.0262s | 4.7x |
| `NetworkSignal._blombergs_k` Cholesky inverse construction | 450 taxa SPD VCV, 16 seeded permutations, side-by-side previous explicit inverse | 0.005925s | 0.002848s | 2.08x |
| `NetworkSignal._blombergs_k_permutations` algebraic batch ratios | 420 taxa SPD VCV, 1000 permutations | 0.021850s | 0.014678s | 1.49x |
| `NetworkSignal._blombergs_k_permutations` invariant trait sum | 260 / 420 / 1000 trait values, side-by-side previous `np.sum(x)` wrapper | 0.000001625s / 0.000001625s / 0.000004125s | 0.000000834s / 0.000000875s / 0.000002292s | 1.95x / 1.86x / 1.80x |
| `NetworkSignal._log_likelihood` | 420 taxa SPD VCV, single continuous trait | 0.0088s | 0.0005s | 17.6x |
| `NetworkSignal._log_likelihood_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV likelihood evaluations, SciPy already warm | 0.054446s | 0.039508s | 1.38x |
| `NetworkSignal._log_likelihood_cholesky` cached SciPy Cholesky wrappers | 120 repeated 420-taxon SPD VCV likelihood evaluations, SciPy already warm, side-by-side previous import-on-call wrappers | 0.090518s | 0.069456s | 1.30x |
| `NetworkSignal._pagels_lambda` | 260 taxa SPD VCV, single continuous trait, `max_lambda=1.0` | 0.7070s | 0.0638s | 11.1x |
| `NetworkSignal.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, DAG construction, VCV, signal statistic, and output mocked | 0.497249s | 0.003029s | 164.16x |
| `NetworkSignal._output` indexed hybrid donor lookup | 10k-tip network output with 2k hybrid recipients sharing a late donor parent, captured stdout | 1.392050s | 0.002664s | 522.62x |
| `NetworkSignal._infer_hybrid_edges` topology index helper | 200k hybrid quartet JSON records, three alternating dominant topologies | 0.647379s | 0.552106s | 1.17x |
| `NetworkSignal._infer_hybrid_edges` cached topology validation and gamma sums | 200k hybrid quartet JSON records, three alternating dominant topologies | 0.469039s | 0.190320s | 2.46x |
| `NetworkSignal._infer_hybrid_edges` direct topology-pair mapping | 200k hybrid quartet JSON records, three alternating dominant topologies, side-by-side previous set-based swap comparison | 0.374105s | 0.178060s | 2.10x |
| `NetworkSignal._infer_hybrid_edges` plain-dict pair counts | 500k quartet JSON records with repeated pair evidence, identical stable descending-frequency edge output | 0.837850s | 0.691837s | 1.21x |
| `NetworkSignal._infer_hybrid_edges` combined pair stats | 500k quartet JSON records with repeated pair evidence, side-by-side previous parallel count/sum dictionaries | 0.307262s | 0.269797s | 1.14x |
| `NetworkSignal._parse_trait_file` streaming two-column parser | 500k two-column trait rows with comments/blanks, all taxa shared, randomized old/new measurement order | 0.457297s | 0.441136s | 1.04x |
| `NetworkSignal._parse_trait_file` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.437650s | 0.239459s | 1.83x |
| `NetworkSignal._parse_trait_file` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.235204s | 0.223874s | 1.05x |
| `network_signal` module import without `scipy.stats` | cold process import for Pagel-lambda-capable command module | 0.585185s | 0.411641s | 1.4x |
| `network_signal` module import without eager SciPy linalg/optimize | cold process import for network-signal command module | 0.428415s | 0.170449s | 2.5x |
| `network_signal` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.070375s | 0.026032s | 2.70x |
| `network_signal` module import without eager JSON helpers | median cold subprocess import after lazy JSON output and load wrappers | 0.006586s | 0.005310s | 1.24x |
| `network_signal` module import without `typing` startup | median cold subprocess import after converting the remaining annotation-only typing alias to a built-in postponed annotation | 0.006474s | 0.004728s | 1.37x |
| `RateHeterogeneity` setup traversal pipeline | balanced tree with 768 tips, Fitch regimes and per-regime VCV setup | 0.1894s | 0.1775s | 1.1x |
| `RateHeterogeneity`/`OUwie` non-binary regime state-set merge | 200 repeated 2048-child multifurcation merges with empty and overlapping state sets, side-by-side previous child-list plus `child_sets[1:]` loop | 0.590826s | 0.423113s | 1.40x |
| `rate_heterogeneity` module import without eager `scipy.stats` | cold process import for LRT-capable command module | 0.607287s | 0.413648s | 1.5x |
| `rate_heterogeneity` module import without eager SciPy linalg/optimize | cold process import for rate-heterogeneity command module | 0.439255s | 0.158809s | 2.8x |
| `rate_heterogeneity` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.093604s | 0.032009s | 2.92x |
| `rate_heterogeneity` module import without eager JSON/pickle/plot helpers | median cold subprocess import after localizing PlotConfig, pickle, and plotting helpers | 0.020559s | 0.006999s | 2.94x |
| `rate_heterogeneity` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005381s | 0.003747s | 1.44x |
| `RateHeterogeneity.run` copied-tree prune setup | balanced 32768-tip tree, all tips shared with trait/regime data | 0.0675s | 0.0075s | 9.0x |
| `RateHeterogeneity.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, regime parsing, fitting, and output mocked; protective prune-copy retained | 1.725839s | 0.618584s | 2.79x |
| `RateHeterogeneity.run` all-shared read-only setup | balanced 32768-tip cached tree, trait/regime parsing, fitting/output mocked, ladderize off | 0.330438s | 0.084517s | 3.91x |
| `RateHeterogeneity._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.035051s | 0.024073s | 1.46x |
| `RateHeterogeneity._iter_postorder` reverse-preorder helper | balanced 131072-tip tree, postorder generator materialized as a list, side-by-side previous visited-tuple helper | 0.091678s | 0.045672s | 2.01x |
| `RateHeterogeneity._parse_trait_file` streaming parse | 500k two-column trait rows with comments/blanks | 0.729276s | 0.653401s | 1.12x |
| `RateHeterogeneity._parse_trait_file` valid-row parsing | 500k two-column trait rows with comments/blanks, all taxa shared | 0.474866s | 0.457581s | 1.04x |
| `RateHeterogeneity._parse_trait_file` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.438800s | 0.238584s | 1.84x |
| `RateHeterogeneity._parse_trait_file` strict split row parsing | 500k two-column trait rows with comments/blanks, all taxa shared, identical traits and error column counts | 0.404013s | 0.371445s | 1.09x |
| `RateHeterogeneity._parse_regime_file` all-shared parser fast path | 500k two-column regime rows with comments/blanks, all taxa shared, randomized old/new measurement order | 0.397329s | 0.232471s | 1.71x |
| `RateHeterogeneity._parse_regime_file` binary strict split row parsing | 500k two-column regime rows with comments/blanks, all taxa shared, identical regime assignments and error column counts | 0.487937s | 0.352139s | 1.39x |
| `RateHeterogeneity._print_text_output` batched regime rows | 100k regime sigma rows, captured stdout and identical text | 0.047003s | 0.036191s | 1.30x |
| `RateHeterogeneity._count_regime_tips` one-pass counts | 1M tip-regime assignments x 64 regimes, identical per-regime counts including absent regimes | 1.537367s | 0.051320s | 29.96x |
| `RateHeterogeneity._assign_branch_regimes` binary state-set merge | balanced 65536-tip tree, three alternating regimes, side-by-side previous generic child-set merge | 0.224014s | 0.182105s | 1.23x |
| `RateHeterogeneity` / `OUwie` regime ambiguity tiebreak | 1M small ambiguous regime sets, identical lexicographic selected regime without sorting | 0.188272s | 0.152076s | 1.24x |
| `RateHeterogeneity._plot_regime_tree` rectangular setup | balanced 32768-tip tree, precomputed branch regimes and parent map | 0.6891s | 0.1084s | 6.4x |
| `RateHeterogeneity._plot_regime_tree` circular setup | balanced 32768-tip tree, precomputed branch regimes and parent map | 0.8730s | 0.2226s | 3.9x |
| `RateHeterogeneity._plot_regime_tree` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.047597s | 0.037744s | 1.26x |
| `RateHeterogeneity._plot_regime_tree` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.057115s | 0.045826s | 1.25x |
| `RateHeterogeneity._plot_regime_tree` rectangular batched regime branches | balanced 2048-tip tree, 3 regimes, real Matplotlib Agg branch/label/legend render | 1.970551s | 0.531163s | 3.71x |
| `RateHeterogeneity._plot_regime_tree` circular batched regime branches/arcs | balanced 2048-tip tree, 3 regimes, real Matplotlib Agg branch/arc/legend render | 1.475253s | 0.261012s | 5.65x |
| `RateHeterogeneity._build_per_regime_vcv` branch accumulation | balanced 1024-tip synthetic root-to-tip paths x 3 regimes | 0.345525s | 0.019535s | 17.7x |
| `RateHeterogeneity._build_per_regime_vcv` single-tip diagonal updates | balanced 2048-tip prepared branch groups x 3 regimes | 0.067988s | 0.054495s | 1.25x |
| `RateHeterogeneity._chi2_sf` df=2 scalar p-values | 2k three-regime LRT chi-square survival probabilities | 0.049870s | 0.000256s | 194.8x |
| `RateHeterogeneity._fit_single_rate` | 420 taxa SPD VCV, single continuous trait | 0.0079s | 0.0005s | 15.8x |
| `RateHeterogeneity._fit_single_rate_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV single-rate BM likelihood evaluations, SciPy already warm | 0.053752s | 0.045064s | 1.19x |
| `RateHeterogeneity._fit_single_rate_cholesky` cached SciPy Cholesky wrappers | 120 repeated 420-taxon SPD VCV single-rate BM likelihood evaluations, SciPy already warm, side-by-side previous import-on-call wrappers | 0.070947s | 0.057375s | 1.24x |
| Tree Cholesky logdet diagonal extraction | 4 / 16 / 64 / 256 / 1024 / 2048 square factors, side-by-side previous `np.diag` copy before log-sum | 0.000001570s / 0.000002759s / 0.000002770s / 0.000004272s / 0.000005422s / 0.000009910s | 0.000001388s / 0.000001408s / 0.000001560s / 0.000002150s / 0.000005251s / 0.000009730s | 1.13x / 1.96x / 1.78x / 1.99x / 1.03x / 1.02x |
| `RateHeterogeneity._fit_multi_rate` | 120 taxa x 3 SPD regime covariance matrices | 0.5818s | 0.0843s | 6.9x |
| `RateHeterogeneity._multi_rate_log_likelihood_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV multi-rate likelihood evaluations, SciPy already warm | 0.053230s | 0.041447s | 1.28x |
| `OUShiftDetection` setup traversal pipeline | balanced tree with 4096 tips, parent map, lineages, edges, descendant counts | 0.2012s | 0.0344s | 5.9x |
| `TraitRateMap`/`RateHeterogeneity`/`OUShiftDetection` `_build_parent_map` no-preorder direct map | balanced 65536-tip tree, shared optional-preorder parent-map helper pattern, optimized helper baseline | 0.030767s | 0.017838s | 1.72x |
| `OUShiftDetection.run` copied-tree prune setup | balanced 32768-tip tree, all tips shared with trait data | 0.0675s | 0.0068s | 9.9x |
| `OUShiftDetection.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, lineage setup, fitting, and output mocked; protective prune-copy retained | 1.376328s | 0.549209s | 2.51x |
| `OUShiftDetection.run` all-shared read-only setup | balanced 32768-tip cached tree, trait parsing, lineage setup, fitting/output mocked | 0.310091s | 0.064237s | 4.83x |
| `OUShiftDetection._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.049777s | 0.036699s | 1.36x |
| `OUShiftDetection._parse_trait_file` streaming valid-row parser | 500k two-column trait rows with comments/blanks, all taxa shared | 0.472496s | 0.460270s | 1.03x |
| `OUShiftDetection._parse_trait_file` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.438335s | 0.237122s | 1.85x |
| `OUShiftDetection._parse_trait_file` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.237912s | 0.228169s | 1.04x |
| `OUShiftDetection._describe_edge` direct terminal-name traversal | balanced 32768-tip tree, labels for all non-root candidate edges | 0.962313s | 0.129219s | 7.45x |
| `OUShiftDetection._count_descendants` reverse-preorder postorder helper | balanced 8192-tip tree, descendant counts for all candidate edges | 0.010322s | 0.006827s | 1.51x |
| `OUShiftDetection._count_descendants` binary child count | balanced 32768-tip tree, descendant counts for all candidate edges, side-by-side previous generator-sum helper | 0.052036s | 0.043487s | 1.20x |
| `OUShiftDetection._print_text_output` batched detected shifts | 100k detected shifts, captured stdout and identical text | 0.058158s | 0.034934s | 1.66x |
| `OUShiftDetection._print_text_output` paired shift lines | 100k detected shifts, captured stdout and identical text, side-by-side previous two-append loop | 0.036603s | 0.032883s | 1.11x |
| `OUShiftDetection._fit_ou1_for_alpha` | balanced 256-tip shared-path covariance, single continuous trait | 1.1596s | 0.1786s | 6.5x |
| `OUShiftDetection._fit_and_score_config` | balanced 256-tip shared-path covariance, 3 synthetic shift columns | 0.0775s | 0.0117s | 6.6x |
| `OUShiftDetection._gls_profile_likelihood_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV x 6-column design matrix GLS likelihood evaluations, SciPy already warm | 0.062531s | 0.058617s | 1.07x |
| `OUShiftDetection` cached SciPy numerical wrappers | 1k Cholesky factor/solve calls, 1k triangular solves, and 100 bounded `minimize_scalar` calls, SciPy already warm, side-by-side previous import-on-call wrappers | 0.020695s | 0.016566s | 1.25x |
| `OUShiftDetection._build_indicator_design_matrix` lineage-row cache | 2048-tip balanced synthetic lineage, 80 eight-shift configs | 0.161885s | 0.001258s | 128.6x |
| `OUShiftDetection._build_shift_weight_matrix` baseline weight total | per-row shifted-regime weight vector with 2 / 3 / 4 / 8 / 16 / 32 / 128 columns, side-by-side previous `np.sum` wrapper | 0.000005610s / 0.000005376s / 0.000005409s / 0.000005934s / 0.000005992s / 0.000005708s / 0.000005718s | 0.000003521s / 0.000003595s / 0.000002428s / 0.000003187s / 0.000002629s / 0.000003032s / 0.000002396s | 1.59x / 1.50x / 2.23x / 1.86x / 2.28x / 1.88x / 2.39x |
| `OUShiftDetection._extract_lasso_configs` flat coefficient indices | 5000 shift coefficients x 1200 LASSO-path steps, sparse nonzero coefficients, side-by-side previous `np.where(...)[0]` extraction | 0.051099s | 0.034082s | 1.50x |
| `OUShiftDetection._extract_lasso_configs` column L2 normalization | 50x200 / 200x1000 / 1000x3000 / 5000x1000 residualized shift-design matrices, side-by-side previous `np.linalg.norm(..., axis=0)` | 0.000007166s / 0.000145125s / 0.001731750s / 0.010120416s | 0.000004417s / 0.000133042s / 0.000690542s / 0.001168125s | 1.62x / 1.09x / 2.51x / 8.66x |
| `OUShiftDetection._compute_pbic_from_vcv` Cholesky information matrix | 520 taxa SPD VCV x 7-column indicator design | 0.003160s | 0.000817s | 3.9x |
| `OUShiftDetection._compute_pbic_from_info` determinant-only correction | 120-parameter SPD information matrix, side-by-side previous scaled inverse determinant | 0.001747838s | 0.000093647s | 18.66x |
| `ou_shift_detection` module import without eager `sklearn` | cold process import for OU-shift command module | 1.241723s | 0.512947s | 2.4x |
| `ou_shift_detection` module import without eager SciPy linalg/optimize | cold process import for OU-shift command module | 0.449763s | 0.183488s | 2.5x |
| `ou_shift_detection` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.084523s | 0.025853s | 3.27x |
| `ou_shift_detection` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.008109s | 0.006566s | 1.24x |
| `ou_shift_detection` module import without eager pickle | median cold subprocess import after lazy pickle proxy | 0.006528s | 0.005131s | 1.27x |
| `ou_shift_detection` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006865s | 0.004855s | 1.41x |
| `PhyloAnova._cholesky_transform` | 420 taxa SPD VCV x 6 responses, 4-group design | 0.0062s | 0.0012s | 5.2x |
| `PhyloAnova._run_anova` | 500 taxa, 4 groups, 1000 RRPP permutations | 0.0383s | 0.0118s | 3.2x |
| `PhyloAnova._run_anova` projection residual SS | 500 taxa, 4 groups, 1000 RRPP permutations, side-by-side previous materialized full-model residuals | 0.016085s | 0.012779s | 1.26x |
| `PhyloAnova._run_manova` | 500 taxa x 6 responses, 4 groups, 1000 RRPP permutations | 0.0652s | 0.0386s | 1.7x |
| `PhyloAnova._run_manova` projection residual SSCP | 500 taxa x 6 responses, 4 groups, 1000 RRPP permutations, side-by-side previous materialized full-model residual matrices | 0.054512s | 0.048688s | 1.12x |
| `PhyloAnova` permutation summary reductions | 1M permutation statistics, side-by-side previous duplicate `std` reduction | 0.003426s | 0.002183s | 1.57x |
| `PhyloAnova` permutation p-value counts | 1M permutation statistics, side-by-side previous `np.mean(permutations >= observed)` reduction | 0.000991s | 0.000390s | 2.54x |
| `PhyloAnova._run_pairwise` residual mean permutations | 600 taxa x 5 responses, 6 groups, 1000 RRPP permutations | 0.2743s | 0.2302s | 1.2x |
| `PhyloAnova._run_pairwise` vectorized residual mean permutations | 600 taxa x 5 responses, 6 groups, 1000 RRPP permutations | 0.2299s | 0.1177s | 2.0x |
| `PhyloAnova._run_pairwise` contrast-weight residual permutations | 600 taxa x 5 responses, 6 groups, 1000 RRPP permutations | 0.119985s | 0.062728s | 1.9x |
| `PhyloAnova._run_pairwise` direct L2 distance reductions | permutation difference matrices shaped 1000x2 / 1000x6 / 10000x6 / 100000x6 / 10000x32 plus observed vectors with 2 / 6 / 32 responses, side-by-side previous `np.linalg.norm` calls | matrices: 0.000029080s / 0.000033346s / 0.000239950s / 0.003130689s / 0.000655159s; vectors: 0.000002851s / 0.000002613s / 0.000002659s | matrices: 0.000014866s / 0.000020804s / 0.000107190s / 0.001361021s / 0.000194944s; vectors: 0.000001434s / 0.000002033s / 0.000002128s | matrices: 1.96x / 1.60x / 2.24x / 2.30x / 3.36x; vectors: 1.99x / 1.29x / 1.25x |
| `PhyloAnova._run_pairwise` flat group-pair mask indices | 120 sparse boolean masks over 500k observations, side-by-side previous `np.where(mask)[0]` extraction used in pairwise setup | 0.208942s | 0.113959s | 1.83x |
| `PhyloAnova._build_design_matrix` cached group lookup | 500k taxa, 12 groups, treatment-coded design matrix | 0.152966s | 0.047562s | 3.22x |
| `PhyloAnova._run_manova` Pillai trace reduction | 2 / 3 / 4 / 8 / 16 trait eigenvalue vectors, side-by-side previous `np.sum` wrapper used in observed and permuted Pillai statistics | 0.000008343s / 0.000007664s / 0.000007563s / 0.000007925s / 0.000007421s | 0.000005759s / 0.000005310s / 0.000003918s / 0.000004917s / 0.000004665s | 1.45x / 1.44x / 1.93x / 1.61x / 1.59x |
| `PhyloAnova._permutation_p_value_and_z` ndarray mean/std reductions | 10 / 100 / 1000 / 10k / 100k permutation statistics, side-by-side previous `np.mean` and `np.std` wrappers | 0.000021215s / 0.000020485s / 0.000020951s / 0.000030437s / 0.000162023s | 0.000015491s / 0.000016347s / 0.000019552s / 0.000025833s / 0.000142907s | 1.37x / 1.25x / 1.07x / 1.18x / 1.13x |
| `PhyloAnova._prepare_phylomorphospace_overlay` direct traversal | balanced 65536-tip tree, parent map and ancestral coordinate setup | 0.958097s | 0.304198s | 3.15x |
| `PhyloAnova._prepare_phylomorphospace_overlay` single-pass parent map and child means | balanced 32768-tip tree, parent map and ancestral coordinate setup | 0.198770s | 0.067740s | 2.93x |
| `PhyloAnova._plot_boxplot` vectorized group masks | 500k taxa across 12 groups, univariate plot group-value preparation, side-by-side previous per-group Python list masks | 0.377527s | 0.053611s | 7.04x |
| `PhyloAnova._plot_phylomorphospace` branch rendering | 4096 phylomorphospace branch segments, real Matplotlib Agg line render | 0.699849s | 0.032767s | 21.36x |
| `PhyloAnova._plot_phylomorphospace` PCA variance total | 2 / 3 / 4 / 8 / 16 / 32 / 128 / 1024 singular values, side-by-side previous `np.sum(S ** 2)` | 0.000005887s / 0.000006054s / 0.000005545s / 0.000005994s / 0.000006021s / 0.000006877s / 0.000004525s / 0.000006631s | 0.000001281s / 0.000001132s / 0.000001423s / 0.000001156s / 0.000001312s / 0.000001259s / 0.000000879s / 0.000001432s | 4.60x / 5.35x / 3.90x / 5.19x / 4.59x / 5.46x / 5.15x / 4.63x |
| `PhyloAnova._plot_phylomorphospace` vectorized group masks | 500k taxa across 12 groups, scatter coordinate preparation, side-by-side previous per-group Python list masks | 0.396877s | 0.075103s | 5.28x |
| `PhyloAnova._print_results` batched pairwise text output | 100k pairwise comparison rows, captured stdout and identical text | 0.094132s | 0.077935s | 1.21x |
| `PhyloAnova.run` cached read-only setup | balanced 32768-tip cached tree, one response trait for every tip, VCV/cholesky/RRPP/output mocked | 0.463085s | 0.176876s | 2.62x |
| `PhyloAnova.run` group/response matrix setup | 120k taxa x 12 parsed trait columns, non-default group column, side-by-side previous per-row conditional response-column scan | 1.181767s | 0.711145s | 1.66x |
| `PhyloAnova`/`PhyloHeatmap`/`PhyloPath._needs_default_branch_lengths` unordered scan | balanced 131072-tip tree, complete branch lengths, optimized helper baseline | 0.023242s | 0.016890s | 1.38x |
| `PhyloAnova._parse_trait_file` single-pass parser | 200k-row mixed group/continuous TSV, comments/blanks before and within data, 100k shared taxa | 0.235966s | 0.220250s | 1.07x |
| `PhyloAnova._parse_trait_file` precomputed column metadata | 100k all-shared mixed group/continuous TSV rows, non-default group column, side-by-side previous row formatter and shared-taxa filter | 0.144863s | 0.084070s | 1.72x |
| `phylo_anova` module import without eager SciPy linalg | cold process import for phylogenetic-ANOVA command module | 0.318135s | 0.183414s | 1.7x |
| `phylo_anova` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.183414s | 0.031434s | 5.84x |
| `phylo_anova` module import without eager JSON/plot helpers | median cold subprocess import after localizing PlotConfig and lazy JSON wrapper | 0.012744s | 0.005728s | 2.22x |
| `phylo_anova` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.037702s | 0.033729s | 1.12x |
| `CharacterMap.run` taxon-set setup | balanced 32768-tip tree, all tips shared with character matrix | 0.0695s | 0.0093s | 7.5x |
| `CharacterMap.run` cached read-only tree setup | balanced 32768-tip cached tree, character parsing, parsimony reconstruction, plotting, and output mocked; protective working copy retained | 1.069391s | 0.600203s | 1.78x |
| `CharacterMap.run` all-clean read-only setup | balanced 32768-tip cached binary tree, all branch lengths present, character parsing, parsimony/plot/output mocked | 0.522463s | 0.206576s | 2.53x |
| `CharacterMap._assign_node_labels` direct preorder | balanced 65536-tip tree, node/tip label assignment with identical labels | 2.497245s | 0.370307s | 6.74x |
| `CharacterMap._summarize_character_states` | 1200 taxa x 5000 characters, wildcard-aware state counts | 0.724018s | 0.332562s | 2.2x |
| `CharacterMap._parse_character_matrix` streaming parser | 300k-row TSV character matrix, 16 characters, blank rows included | 0.564221s | 0.483253s | 1.17x |
| `CharacterMap._print_summary` batched summary output | 100k captured character-map text summaries, alternating numeric/N/A metrics, identical stdout text | 0.149204s | 0.078751s | 1.89x |
| `CharacterMap._print_verbose` grouped change rows | 1k characters, 20k classified changes, captured stdout and identical text | 0.957803s | 0.006996s | 136.91x |
| `CharacterMap._plot_character_map` phylogram layout setup | balanced 32768-tip tree, precomputed classified changes and parent map | 0.5485s | 0.0787s | 7.0x |
| `CharacterMap._plot_character_map` cladogram layout setup | balanced 32768-tip tree, precomputed classified changes and parent map | 0.7097s | 0.1002s | 7.1x |
| `CharacterMap._iter_preorder` binary-child fast path | balanced 131072-tip tree, preorder generator materialized as a list, side-by-side previous `reversed(children)` helper | 0.058023s | 0.041283s | 1.41x |
| `CharacterMap._plot_character_map` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, cladogram coordinate setup | 0.060839s | 0.045614s | 1.33x |
| `CharacterMap._plot_character_map` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.055037s | 0.045283s | 1.22x |
| `CharacterMap._plot_character_map` rectangular batched base branches | balanced 2048-tip tree, empty changes, real Matplotlib Agg branch/label render | 1.929078s | 0.528533s | 3.65x |
| `CharacterMap._plot_character_map` rectangular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 1.237836s | 0.041628s | 29.74x |
| `CharacterMap._plot_character_map` circular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 0.613073s | 0.029609s | 20.71x |
| `CharacterMap._plot_character_map` change markers | 4096 classified character-change markers with variable colors, real Matplotlib Agg scatter render | 9.393315s | 0.056376s | 166.62x |
| `character_map` module import without eager circular-layout NumPy | cold subprocess import after shared lazy circular-layout arc cache | 0.080137s | 0.032536s | 2.46x |
| `character_map` module import without eager pickle/JSON/parsimony/plot helpers | median cold subprocess import after localizing pickle, PlotConfig/layout helpers, circular/color helpers, and lazy JSON/parsimony wrappers | 0.015168s | 0.005650s | 2.68x |
| `character_map` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.005191s | 0.003692s | 1.41x |
| `StochasticCharacterMap._prune_tree_to_tip_states` | balanced 2048-tip tree, prune 1024 tips before SIMMAP setup | 0.6334s | 0.4164s | 1.5x |
| `StochasticCharacterMap._prune_tree_to_tip_states` batch standard-tree pruning | balanced 2048-tip tree, prune 1024 tips before SIMMAP setup | 0.4293s | 0.0240s | 17.9x |
| `StochasticCharacterMap._prune_tree_to_tip_states` target scan child push | balanced 262144-tip tree, collect 131072 missing-tip targets, side-by-side previous `reversed(children)` scan | 0.092651s | 0.081674s | 1.13x |
| `StochasticCharacterMap.run` all-tip-state read-only setup | balanced 32768-tip cached tree, trait data for every tip, fitting/simulation/output mocked | 0.248931s | 0.047217s | 5.27x |
| `SimmapSummary` branch label generation | balanced 2048-tip tree, labels for all clades in text/JSON/CSV output | 0.0395s | 0.0076s | 5.2x |
| `QuartetPie._plot_quartet_pie` rectangular setup | balanced 32768-tip tree, precomputed proportions and branch labels | 0.3998s | 0.1165s | 3.4x |
| `QuartetPie._plot_quartet_pie` circular setup | balanced 32768-tip tree, precomputed proportions and branch labels | 0.4671s | 0.1889s | 2.5x |
| `QuartetPie._plot_quartet_pie` node-position preorder reuse | balanced 32768-tip tree, parent map and preorder list already available, phylogram coordinate setup | 0.051277s | 0.037255s | 1.38x |
| `QuartetPie._plot_quartet_pie` circular coordinate clade-list reuse | balanced 32768-tip tree, node positions plus preorder/tip lists already available | 0.055361s | 0.045029s | 1.23x |
| `QuartetPie._plot_quartet_pie` rectangular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 1.231566s | 0.042586s | 28.92x |
| `QuartetPie._plot_quartet_pie` circular clade-color overlay rendering | balanced 2048-tip tree, all branches highlighted by color-file clade, real Matplotlib Agg overlay render | 0.614624s | 0.029184s | 21.06x |
| `QuartetPie` JSON/CSV branch label generation | balanced 2048-tip tree, labels for all internal nodes in JSON and CSV output | 0.1075s | 0.0413s | 2.6x |
| `QuartetPie` direct JSON/CSV output traversal | balanced 2048-tip tree, labels for all internal nodes in JSON and CSV output | 0.023108s | 0.010872s | 2.1x |
| `QuartetPie._collect_clade_tip_names` reverse-preorder postorder helper | balanced 4096-tip tree, descendant tip-name cache for JSON/CSV output | 0.011898s | 0.005105s | 2.33x |
| `QuartetPie._collect_clade_tip_names` binary child tuple cache | balanced 8192-tip tree, descendant tip-name cache for JSON/CSV output | 0.008510s | 0.004790s | 1.78x |
| `QuartetPie.run` native confidence map | balanced 2048-tip species tree, branch info for every internal node | 5.706448s | 0.001047s | 5450.7x |
| `QuartetPie._iter_preorder` order-preserving child push | balanced 131072-tip tree, direct preorder helper, side-by-side previous `reversed(children)` helper | 0.033571s | 0.020755s | 1.62x |
| `QuartetPie._collect_branch_confidence` preorder child push | balanced 131072-tip tree, confidence on every internal clade, side-by-side previous preorder helper | 0.047254s | 0.033501s | 1.41x |
| `quartet_pie` module import without unused NumPy | cold subprocess import after removing unused NumPy import | 0.083696s | 0.032468s | 2.58x |
| `quartet_pie` module import without eager helper modules | median cold subprocess import after lazy JSON wrapper plus localized plot, quartet, circular, and color helper imports | 0.017209s | 0.005564s | 3.09x |
| `quartet_pie` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.036311s | 0.030926s | 1.17x |
| `CharacterMap._summarize_character_states` full-column Counter | 1200 taxa x 5000 characters, wildcard-aware state counts | 0.361532s | 0.194653s | 1.9x |
| `CharacterMap._summarize_character_states` full summary byte-count reuse | 1200 taxa x 5000 single-character states, wildcard-aware state lists and counts | 0.206011s | 0.131361s | 1.57x |
| `CharacterMap._summarize_character_states` byte-matrix state-list transpose | 1200 taxa x 5000 single-character states, wildcard-aware state lists and counts | 0.135322s | 0.100202s | 1.35x |
| `CharacterMap._summarize_ascii_matrix` large-alphabet column counts | three repeated 5000-taxon x 2000-character ASCII matrices with 36 observed states plus wildcards, side-by-side previous per-symbol equality reductions | 2.852672s | 0.460152s | 6.20x |
| `retention_index` full-column Counter | 1200 taxa x 5000 characters, wildcard-aware RI state counts | 0.370389s | 0.149759s | 2.5x |
| `retention_index` ASCII single-character fast path | 1200 taxa x 5000 single-character states, wildcard-aware RI state counts | 0.125720s | 0.039702s | 3.17x |
| `CharacterMap.run` state summary plus RI counts | 1200 taxa x 5000 characters, wildcard-aware counts and synthetic observed steps | 0.394492s | 0.241499s | 1.63x |
| `CharacterMap.run` counts-only byte state summary | 1200 taxa x 5000 single-character states, wildcard-aware counts | 0.183195s | 0.051354s | 3.57x |
| `CharacterMap._summarize_character_ascii` counts-only streaming matrix | 1200 taxa x 5000 single-character states, identical wildcard-aware counts without row-list or giant text materialization | 0.238175s | 0.207201s | 1.15x |
| `ThresholdModel._run_mcmc` continuous/continuous | 180 taxa SPD VCV, 1500 generations, sample every 10 | 0.3091s | 0.0750s | 4.1x |
| `ThresholdModel._run_mcmc` discrete/continuous | 180 taxa SPD VCV, 700 generations, sample every 10 | 12.5360s | 1.3572s | 9.2x |
| `ThresholdModel._run_mcmc` one-sided truncated-normal fast path | 120 taxa SPD VCV, discrete/continuous traits, 300 generations, sample every 10 | 0.456215s | 0.253992s | 1.8x |
| `ThresholdModel._run_mcmc` incremental Gibbs residuals | 120 taxa SPD VCV, discrete/continuous traits, 300 generations, sample every 10 | 0.272798s | 0.218870s | 1.2x |
| `ThresholdModel._run_mcmc` cached Gibbs/stat products | 180 taxa SPD VCV, discrete/continuous traits, 700 generations, sample every 10 | 0.627373s | 0.546933s | 1.15x |
| `ThresholdModel._run_mcmc` cached truncated-normal float bounds | 180 taxa SPD VCV, discrete/continuous traits, 700 generations, sample every 10 | 0.558707s | 0.498435s | 1.12x |
| `ThresholdModel._run_mcmc` cached proposal scales/log sigma | 180 taxa SPD VCV, discrete/continuous traits, 700 generations, sample every 10 | 0.927273s | 0.860394s | 1.08x |
| `ThresholdModel._run_mcmc` scalar math and cached Gibbs context fast path | 180 taxa SPD VCV, discrete/continuous traits, 700 generations, sample every 10 | 0.860394s | 0.433986s | 1.98x |
| `ThresholdModel._run_mcmc` continuous-trait initial variance reuse | 4 / 10 / 100 / 1000 / 10k / 100k liability values, side-by-side previous duplicate `np.var` expression | 0.000022695s / 0.000014483s / 0.000016306s / 0.000016240s / 0.000034778s / 0.000212574s | 0.000010199s / 0.000010521s / 0.000009156s / 0.000010178s / 0.000019218s / 0.000071040s | 2.23x / 1.38x / 1.78x / 1.60x / 1.81x / 2.99x |
| `ThresholdModel._run_mcmc` continuous-trait ndarray reductions | 4 / 10 / 100 / 1000 / 10k liability values, side-by-side previous `np.var`/`np.mean` wrappers during initial setup | 0.000008083s / 0.000008042s / 0.000008166s / 0.000009084s / 0.000018125s | 0.000007000s / 0.000007000s / 0.000007084s / 0.000008041s / 0.000016791s | 1.15x / 1.15x / 1.15x / 1.13x / 1.08x |
| `ThresholdModel._run_mcmc` covariance diagonal mean setup | SPD VCV matrices with 120 / 180 / 500 / 1000 taxa, side-by-side previous `np.mean(np.diag(C))` copy | 0.000003209s / 0.000003167s / 0.000003458s / 0.000003666s | 0.000001583s / 0.000001584s / 0.000001708s / 0.000001875s | 2.03x / 2.00x / 2.02x / 1.96x |
| `ThresholdModel._sample_liabilities_gibbs` inline one-sided inverse-CDF draws | 700 dense 180-tip Gibbs sweeps, alternating binary states, cached Gibbs context | 0.252991s | 0.239303s | 1.06x |
| `ThresholdModel._sample_truncated_normal` cached float bounds | 126k one-sided inverse-CDF draws | 0.202700s | 0.138897s | 1.46x |
| `ThresholdModel._sample_truncated_normal` cached lazy SciPy special helpers | 126k one-sided inverse-CDF draws | 0.227829s | 0.172824s | 1.32x |
| `ThresholdModel._sample_truncated_normal` scalar math normal CDF | 126k alternating one-sided inverse-CDF draws, copied old sampler baseline, identical RNG streams | 0.103077s | 0.092541s | 1.11x |
| `ThresholdModel._vcv_inverse_and_logdet` Cholesky setup | 500-taxon SPD VCV, side-by-side previous explicit inverse plus `slogdet` | 0.023747s | 0.002966s | 8.01x |
| `ThresholdModel` full inverse-matrix scalar sum | 2 / 4 / 16 / 64 / 256 / 1024 / 2048 square inverse VCV matrices, side-by-side previous `np.sum(C_inv)` wrapper | 0.000001934s / 0.000002439s / 0.000001560s / 0.000002062s / 0.000055042s / 0.000314244s / 0.001257795s | 0.000000699s / 0.000000699s / 0.000000752s / 0.000001328s / 0.000014229s / 0.000237814s / 0.001173287s | 2.77x / 3.49x / 2.07x / 1.55x / 3.87x / 1.32x / 1.07x |
| `ThresholdModel._bivariate_quadratic_stats_from_products` vector sums | 40-taxon cached `C^-1x` product vectors, identical sufficient-stat tuple | 0.000012523s | 0.000005615s | 2.23x |
| `ThresholdModel._initialize_liabilities` vectorized discrete liabilities | 20k binary discrete taxa, scalar SciPy stream preserved | 2.035507s | 0.005099s | 399.2x |
| `ThresholdModel._summarize_posterior` single-sort median/HPD | five 500k-sample posterior traces, identical summary statistics | 0.317456s | 0.248927s | 1.28x |
| `ThresholdModel._summarize_posterior` small-trace mean reductions | five 1k-sample posterior traces, side-by-side previous `np.mean` wrapper after single-sort summary setup | 0.337282s | 0.166087s | 2.03x |
| `ThresholdModel._output_text` batched summary output | 100k captured threshold-model text summaries, identical stdout text | 0.414422s | 0.337860s | 1.23x |
| `ThresholdModel._output_text` single-report formatting | 100k captured threshold-model text summaries, identical stdout text, side-by-side previous newline-joined list comparison | 0.301403s | 0.285244s | 1.06x |
| `ThresholdModel._output_text` percent report template | 100k captured threshold-model text summaries, identical stdout text, side-by-side previous f-string report formatter comparison | 0.333694s | 0.258711s | 1.29x |
| `ThresholdModel._parse_multi_trait_file` streaming parser | 500k-row mixed discrete/continuous TSV with an extra column, all taxa shared | 0.750809s | 0.734721s | 1.02x |
| `ThresholdModel._parse_multi_trait_file` all-shared parser fast path | 500k-row mixed discrete/continuous TSV with an extra column, all taxa shared | 0.743827s | 0.714690s | 1.04x |
| `ThresholdModel._parse_multi_trait_file` bounded split and all-shared direct return | 500k-row mixed discrete/continuous TSV with an extra column, all taxa shared, side-by-side previous full-split/filter-copy comparison | 0.710943s | 0.493796s | 1.44x |
| `threshold_model` module import without eager `scipy.stats` | cold process import for threshold-model command module | 0.624578s | 0.362909s | 1.7x |
| `threshold_model` module import without eager SciPy special/shared linalg | cold process import for threshold-model command module | 0.331525s | 0.144341s | 2.3x |
| `threshold_model` module import after lazy shared VCV Bio.Phylo | cold subprocess import of threshold-model command module | 0.196941s | 0.121166s | 1.63x |
| `threshold_model` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and `sys.float_info` constants | 0.121166s | 0.032146s | 3.77x |
| `threshold_model` module import without eager JSON/plot/VCV helpers | median cold subprocess import after localizing PlotConfig, JSON output, and VCV helper import | 0.019158s | 0.005869s | 3.26x |
| `threshold_model` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006424s | 0.004589s | 1.40x |
| `ThresholdModel.run` all-shared read-only setup | balanced 32768-tip cached tree, two traits for every tip, VCV/MCMC/summary/output mocked | 0.326072s | 0.102394s | 3.18x |
| `ouwie` module import without eager SciPy linalg/optimize | cold process import for OUwie command module | 0.442738s | 0.169719s | 2.6x |
| `ouwie` module import without eager NumPy | cold subprocess import after lazy NumPy proxy and postponed annotations | 0.169719s | 0.025757s | 6.59x |
| `ouwie` module import without eager JSON helper | median cold subprocess import after lazy JSON wrapper | 0.007878s | 0.006670s | 1.18x |
| `ouwie` module import without eager pickle | median cold subprocess import after lazy pickle proxy | 0.006479s | 0.005000s | 1.30x |
| `ouwie` module import without `typing` startup | median cold subprocess import after converting annotation-only typing aliases to built-in postponed annotations | 0.006873s | 0.004510s | 1.52x |
| `OUwie.run` copied-tree prune setup | balanced 32768-tip tree, all tips shared with trait/regime data | 0.0672s | 0.0077s | 8.7x |
| `OUwie.run` cached read-only tree setup | balanced 32768-tip cached tree, trait parsing, regime parsing, model fitting, and output mocked; protective prune-copy retained | 1.174518s | 0.358352s | 3.28x |
| `OUwie.run` all-shared read-only setup | balanced 32768-tip cached tree, trait/regime parsing, model fitting/output mocked | 0.292509s | 0.073640s | 3.97x |
| `OUwie.run` shared trait/regime setup | 300k trait taxa x 300k regime taxa with 225k overlap, identical filtered maps and sorted names | 2.820324s | 2.602742s | 1.08x |
| `OUwie._parse_trait_file` streaming valid-row parser | 500k two-column trait rows with comments/blanks, all taxa shared | 0.470183s | 0.452638s | 1.04x |
| `OUwie._parse_trait_file` all-shared parser fast path | 500k two-column trait rows with comments/blanks, all taxa shared | 0.433626s | 0.241180s | 1.80x |
| `OUwie._parse_trait_file` two-column split fast path | 500k two-column trait rows with comments/blanks, all taxa shared, side-by-side previous partition parser comparison | 0.241349s | 0.224103s | 1.08x |
| `OUwie._parse_regime_file` streaming valid-row parser | 500k two-column regime rows with comments/blanks, all taxa shared | 0.408915s | 0.397772s | 1.03x |
| `OUwie._parse_regime_file` all-shared parser fast path | 500k two-column regime rows with comments/blanks, all taxa shared, side-by-side previous mismatch-set/filtering comparison | 0.979687s | 0.489047s | 2.00x |
| `OUwie._build_parent_map` direct traversal | balanced 65536-tip tree, object parent map for regime setup | 0.208800s | 0.022670s | 9.21x |
| `OUwie._build_parent_map` unordered child push | balanced 65536-tip tree, object parent map for regime setup, optimized helper baseline | 0.021098s | 0.017901s | 1.18x |
| `OUwie._assign_branch_regimes_and_root` combined setup | balanced 32768-tip tree, three alternating regimes, copied two-pass branch/root helpers baseline | 0.510217s | 0.084746s | 6.02x |
| `OUwie._assign_branch_regimes_and_root` binary state-set merge | balanced 32768-tip tree, three alternating regimes, side-by-side previous generic child-set merge | 0.153070s | 0.115853s | 1.32x |
| `OUwie._build_root_to_tip_paths` set-backed tip filter | balanced 8192-tip tree, all ordered tips included in lineage setup | 0.299838s | 0.033505s | 8.95x |
| `OUwie._build_lineage_info` set-backed tip filter | balanced 8192-tip tree, all ordered tips included in lineage setup, identical lineage output | 0.491299s | 0.072740s | 6.75x |
| `OUwie._print_text_output` batched model and parameter rows | 100k regimes, 300k best-model parameter rows, captured stdout and identical text | 0.093383s | 0.061903s | 1.51x |
| `OUwie._compute_model_comparison` scalar AICc weights | seven synthetic model-result dictionaries, identical normalized weights | 0.000018015s | 0.000006892s | 2.61x |
| `OUwie._compute_model_comparison` unweighted multi-regime sigma2 average | 25 model-result dictionaries x 100k regime sigma2 values, no regime-tip weights, identical R2 inputs | 0.023221s | 0.014002s | 1.66x |
| `OUwie._concentrated_ll_bm` | 220 taxa SPD VCV, single continuous trait | 0.0021s | 0.0002s | 10.5x |
| `OUwie._concentrated_ll_bm_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV concentrated BM1 likelihood evaluations, SciPy already warm | 0.053318s | 0.042459s | 1.26x |
| `OUwie._fit_bms` | 120 taxa x 3 SPD regime covariance matrices | 1.5649s | 0.1192s | 13.1x |
| `OUwie._bm_log_likelihood_from_vcv_cholesky` combined RHS solve | 120 repeated 420-taxon SPD VCV BMS likelihood evaluations, SciPy already warm | 0.054405s | 0.043169s | 1.26x |
| `OUwie._ou_profile_likelihood_cholesky` combined RHS solve | 128 taxa SPD VCV x 3-column weight matrix | 0.000100s | 0.000072s | 1.4x |
| `OUwie._ou_gls_log_likelihood` Cholesky solve | 120 taxa SPD VCV x 3-column weight matrix | 0.000447s | 0.000073s | 6.1x |
| `OUwie._ou_profile_likelihood_cholesky` residual solve reuse | 3000 repeated 256-taxon balanced-lineage likelihoods, 3-column weight matrix | 0.640801s | 0.554231s | 1.16x |
| `OUwie._ou_gls_log_likelihood_cholesky` residual solve reuse | 3000 repeated 256-taxon balanced-lineage GLS likelihoods, 3-column weight matrix | 0.715722s | 0.559679s | 1.28x |
| `OUwie._fit_ou1` | balanced 128-tip synthetic lineage, single continuous trait | 3.6601s | 2.0022s | 1.8x |
| `OUwie._fit_oum` | balanced 128-tip synthetic lineage x 3 regimes | 4.1251s | 2.1616s | 1.9x |
| `OUwie._build_ou_vcv_single_alpha` shared-path cache | balanced 256-tip synthetic lineage, alpha=0.7 | 0.024579s | 0.010688s | 2.3x |
| `OUwie._fit_ou1` single-alpha shared-path cache | balanced 256-tip synthetic lineage, single continuous trait | 7.616826s | 0.171313s | 44.5x |
| `OUwie._fit_oum` single-alpha shared-path cache | balanced 128-tip synthetic lineage x 3 regimes | 1.931321s | 0.188768s | 10.2x |
| `OUwie._build_weight_matrix_single_alpha` cache reuse | balanced 2048-tip synthetic lineage x 4 regimes, 50 alpha values | 0.655751s | 0.018331s | 35.8x |
| `OUwie._prepare_single_alpha_ou_cache` branch accumulation | balanced 2048-tip synthetic lineage | 0.652449s | 0.059365s | 11.0x |
| `OUwie._prepare_single_alpha_ou_cache` single-tip diagonal updates | balanced 2048-tip synthetic lineage | 0.081723s | 0.064831s | 1.26x |
| `OUwie._build_per_regime_vcv` branch accumulation | balanced 1024-tip synthetic lineage x 3 regimes | 0.307685s | 0.029079s | 10.6x |
| `OUwie._build_per_regime_vcv` single-tip diagonal updates | balanced 4096-tip synthetic lineage x 3 regimes | 0.267320s | 0.195444s | 1.37x |
| `OUwie.run` per-regime VCV lineage-cache reuse | balanced 2048-tip synthetic lineage x 3 regimes, lineage info plus per-regime VCV setup | 0.168087s | 0.099107s | 1.70x |
| `OUwie._build_ou_H_matrices` branch outer products | balanced 512-tip synthetic lineage x 3 regimes, alpha=0.7 | 0.116251s | 0.013122s | 8.9x |
| `OUwie._build_ou_H_matrices` single-tip diagonal updates | balanced 2048-tip synthetic lineage x 3 regimes, alpha=0.7 | 0.129839s | 0.092253s | 1.41x |
| `OUwie._build_ou_vcv_multi_alpha` branch outer products | balanced 512-tip synthetic lineage x 3 regimes, regime-specific alpha/sigma2 | 0.434786s | 0.009932s | 43.8x |
| `OUwie._build_ou_vcv_multi_alpha` single-tip diagonal updates | balanced 2048-tip synthetic lineage x 3 regimes, regime-specific alpha/sigma2 | 0.154475s | 0.098951s | 1.56x |

Profiling summary:

- `print_json` baseline time recursively copied every payload through
  `to_builtin_json_types`, even when it already contained only JSON-native
  Python values. The optimized helper serializes through `json.dumps` with a
  small `default` hook, preserving sorted-key output and existing NumPy
  scalar/array support while converting only unsupported NumPy values that the
  encoder encounters. The shared JSON helper now identifies NumPy values
  structurally instead of importing NumPy at module import time. A later startup
  pass also defers stdlib `json` until `print_json` is called, so command
  discovery that only imports the helper avoids the serializer import. Direct
  `to_builtin_json_types` conversion still scans nested payloads for non-builtin
  values but now reuses dict/list branches that do not change, avoiding a full
  duplicate of large JSON-native row blocks. JSON-native scalar leaves now
  return before the NumPy provenance checks, preserving NumPy scalar/array
  conversion while reducing recursive conversion overhead for mixed payloads
  dominated by builtin rows.
- `helpers.caching` baseline import time eagerly imported `json` and `pickle`
  and constructed global cache instances, which also initialized the default
  temporary cache directory. The optimized helper defers JSON and pickle until
  key generation or cache reads/writes, and creates global cache instances only
  on `get_result_cache()` or `get_alignment_cache()`, preserving singleton
  behavior and the module-level `pickle.loads` patch point. A later startup
  pass also defers `hashlib` until `_get_cache_key()` is called. Repeated
  cache-key generation now caches the resolved `hashlib.md5` function after the
  first key, avoiding repeated import dispatch while keeping import-only callers
  free of `hashlib` startup. Cache clearing now scans directory entries with
  `os.scandir()` and removes `entry.path` for `.pkl` files, avoiding a full
  `listdir()` filename list and repeated path joins during cleanup.
- `Phykit.__init__` previously built the full top-level help parser and
  dedented the long command listing before dispatching normal commands. The
  optimized constructor dispatches ordinary commands and aliases directly, while
  preserving the full top-level parser for `phykit`, `phykit -h`, and
  `phykit --help`. CLI startup also keeps `add_plot_arguments` behind a local
  forwarding wrapper, so non-plot command discovery avoids importing the shared
  plot configuration helper. A follow-up startup pass removed unused module-level
  logger/handler construction, avoiding an eager stdlib `logging` import during
  `phykit.phykit` import. Parser construction now also defers `argparse` until
  `_new_parser()` runs; the global `SUPPRESS` value is synchronized with
  `argparse.SUPPRESS` at that point, and `str2bool` imports `argparse` only for
  the invalid-value error path. Parser-description dedenting now goes through a
  lazy `_dedent()` helper, so import-only callers avoid `textwrap` until a help
  description or banner actually needs dedenting. The default version and
  invalid-alias banners now cache their dedented text after first use, preserving
  lazy `textwrap` startup and custom `help_header` version behavior while
  avoiding repeated dedent work. Geological timescale epoch
  selection now skips an overwritten preliminary interval list, preserving the
  same included epoch bands while scanning the epoch table once per call. Epoch
  selection now also stops at the first out-of-range epoch because the bundled
  epoch table is ordered by increasing end age, avoiding older epochs that
  cannot be included for shallower chronogram ranges.
  Period-range selection stops at the first out-of-range period because the
  bundled period table is ordered by increasing end age, preserving included
  bands while avoiding later table entries for shallower period-scale plots.
- `trait_parsing.parse_multi_trait_file` baseline time materialized the whole
  file with `readlines()` and then built a second stripped/comment-filtered
  `data_lines` list. The optimized parser finds the first non-comment header
  and parses data rows directly from the file iterator, preserving filtering,
  warnings, shared-taxa output, and logical data-row error numbering. The shared
  parser now converts valid numeric rows through a list comprehension and only
  falls back to the detailed per-column loop when conversion fails, preserving
  non-numeric trait messages while reducing valid-row overhead. A follow-up
  narrow-row pass converts common one-to-four-trait rows with direct float calls
  while still routing conversion failures through the detailed per-column error
  helper. The shared
  trait-matrix builder now passes parsed trait rows directly to NumPy instead of
  copying each row with an inner trait-column loop; TraitCorrelation,
  PhylogeneticOrdination, and multivariate PhylogeneticSignal use the helper.
  Regression-style design matrix setup now gathers the response and predictor
  columns in one selected-column pass with `itemgetter`, and builds the
  intercept matrix directly; one-to-four-predictor designs now fill each
  selected column with `np.fromiter()` to avoid a temporary tuple matrix, while
  wider designs retain the faster `itemgetter` matrix path. PhylogeneticRegression,
  PhyloLogistic, and PhylogeneticGLM use the helper. Their run paths now build
  a first-occurrence trait-name index map once before validation and column
  extraction, preserving `list.index()` behavior for duplicate names while
  avoiding repeated linear scans across wide trait headers. A later parser pass
  returns immediately for exact tree/trait taxon matches after all rows are
  validated and the configured minimum shared-taxa threshold is satisfied,
  avoiding shared/warning set construction while preserving too-few-shared-taxa
  errors.
  No-axis boolean guard reductions in tree/statistics helpers now use ndarray
  `any()`/`all()` methods, applying the shared mask-reduction timing win while
  leaving axis-based reductions on their existing `np.any` paths where
  benchmarks were mixed.
  Discordance-VCV trait subsetting now uses the sorted `shared_taxa` contract
  from `build_discordance_vcv()` and compares ordered taxon lists directly
  instead of allocating two sets before every filter; FitContinuous,
  TraitCorrelation, PhylogeneticOrdination, PhylogeneticSignal,
  PhylogeneticRegression, and PhylogeneticGLM use the helper.
- `IdentityMatrix._compute_identity_matrix` baseline time was dominated by
  pair-level valid-character list comprehensions and repeated `np.array`
  construction. The optimized path precomputes one sequence matrix and one
  validity mask matrix per alignment. A later taxon-heavy pass routes ASCII
  alignments through the same matrix-product strategy used by partition
  identities, computing compared-site and per-symbol match matrices once for
  all taxon pairs while retaining the pairwise path for non-ASCII sequences.
  The dense products now use float masks so NumPy routes the count matrices
  through optimized BLAS; counts remain exact for practical alignment lengths,
  and the final identity matrix/output format is unchanged. Longer clean ASCII
  alignments with moderate taxon counts now bypass validity-mask construction
  and use direct blockwise sequence comparisons, while shorter, taxon-heavy, or
  invalid-symbol alignments stay on the existing matrix-product path. Fully
  identical sequence dictionaries now fill the identity matrix directly from one
  sequence before byte-matrix construction, preserving zero off-diagonal
  identity when no valid sites are available. The identical-sequence dictionary
  guard now scans taxon names with an iterator instead of allocating
  `taxa_names[1:]`, preserving early mismatch exits without copying the taxon
  name tail. The non-ASCII pairwise fallback now counts valid and matching
  boolean masks with `np.count_nonzero` instead of summing booleans, preserving
  Unicode fallback behavior while avoiding a general reduction.
- `IdentityMatrix` clustered heatmap setup baseline time computed the same
  hierarchical linkage once for taxon ordering and again for dendrogram
  plotting. The optimized run path returns the linkage from ordering and passes
  it into the plot helper, while direct plot-helper calls still compute linkage
  when none is provided. A later startup pass replaced the eager SciPy
  clustering imports with same-name lazy wrappers, preserving existing patch
  points and deferring the clustering stack until ordering or dendrogram work
  actually needs it. The clustering wrappers now cache the resolved SciPy
  callables after first use, avoiding repeated import dispatch in clustered
  ordering and dendrogram rendering while preserving wrapper patch points. A
  subsequent startup pass removed the unused eager Bio.SeqIO import, moved
  Bio.Phylo import into the tree-sort branch, and moved the shared Bio.AlignIO
  import into the file-read helper so importing identity-matrix code does not
  load Biopython parsers before they are needed. The shared protein/nucleotide
  alignment detector now checks ASCII sequences with byte deletion instead of
  building an uppercased character set per record, retaining the original
  Unicode fallback behavior. Text summary output now keeps the single batched
  print call but uses a static `%s` template fed by
  the existing rounded values, preserving exact stdout text while reducing
  per-report formatting overhead.
  Another startup pass defers NumPy, JSON output, and plot configuration behind
  a module-level proxy/wrapper or argument-processing import; the existing
  SciPy wrapper patch points and tree-sort behavior are unchanged. A later
  startup pass localizes the shared tree-base helper import to the tree-sort
  branch and replaces annotation-only `typing` aliases with built-in postponed
  annotations, so command discovery no longer imports tree caching/hash helpers
  or `typing`.
- `IdentityMatrix._determine_order` tree sorting still parses Newick input
  through Bio.Phylo, but terminal order extraction now uses the shared direct
  traversal before falling back to `get_terminals()` for nonstandard trees.
- `IdentityMatrix._compute_partition_identities` baseline time scanned every
  upper-triangle taxon pair inside every partition. The optimized ASCII path
  builds one byte matrix, uses matrix products to compute all compared-site and
  per-symbol match counts for a partition at once, and retains the pairwise
  implementation as a fallback for non-ASCII sequences. A follow-up pass
  precomputes the float validity mask once and slices it per partition, avoiding
  repeated dense boolean-to-float conversions while preserving the existing
  per-symbol match products. Fully identical sequence dictionaries now compute
  each partition's constant pair identity from one sequence and tile it over
  upper-triangle pairs, preserving all-invalid partitions as zero identity.
- `IdentityMatrix._partition_identity_strip` baseline time rebuilt a
  pair-to-index dictionary and scanned all other taxa for every strip row and
  partition. The optimized helper accumulates pairwise partition identities onto
  both participating taxa with `np.add.at`, divides by `n_taxa - 1`, and then
  applies the existing taxon order. Text summary output now batches the
  ten-line report into one newline-joined print while preserving exact stdout
  text and the surrounding broken-pipe handling. Identity matrix summaries now
  reduce the condensed `squareform` vector with ndarray methods and reuse
  argmin/argmax locations for min/max values, preserving reported taxon pairs
  while avoiding generic NumPy reduction dispatch and extra scans.
- `IdentityMatrix._parse_partitions` baseline time split valid partition rows
  into temporary lists for comma, equals, and dash separators. The optimized
  parser uses `str.partition()` in the same delimiter order, preserving comment
  skips, comma/no-comma forms, whitespace handling, and 1-indexed-to-slice
  coordinate conversion while reducing per-row allocation.
- `PairwiseIdentity._process_pair_batch` baseline time was dominated by
  repeated `np.isin` calls for the same taxa when `--exclude-gaps` was enabled.
  The optimized path precomputes each taxon's gap mask once and reuses it
  across pair comparisons. ASCII sequence arrays are now stored as byte arrays,
  with a Unicode fallback for non-ASCII sequence content. A later pass added an
  equal-length ASCII matrix path for `--exclude-gaps`, computing all pairwise
  valid-symbol matches by matrix products while preserving the existing
  full-length denominator semantics. Clean ASCII alignments with no gap or
  ambiguous symbols now route `--exclude-gaps` through the faster default matrix
  path for both full pair rows and summary-only stats, preserving the existing
  mask path for gapped data. Verbose `run` now batches pair rows into one
  newline-joined print while preserving stdout text, JSON payloads, summary
  output, plot reporting, and broken-pipe handling. A later startup pass deferred SciPy
  hierarchical-clustering imports until heatmap plotting, so normal non-plot
  text and JSON runs no longer import the clustering stack. A subsequent startup
  pass defers the annotation-only Bio.Align import and optional tqdm progress
  import until type checking or the multiprocessing progress path needs them.
  A later startup pass builds the class-level byte gap lookup tables lazily while
  preserving the module-level `np.isin` patch point used by fallback-path tests.
  Another startup pass defers multiprocessing, summary-statistics, JSON, and
  plot-config imports behind module-level wrappers while preserving the
  `mp.Pool`, `print_json`, and `print_summary_statistics` patch points. A
  follow-up startup pass removes the runtime `TYPE_CHECKING` dependency and
  converts annotation-only typing aliases to postponed built-in annotations. The
  default equal-length ASCII path now uses blocked row comparisons, counting
  gap and ambiguous-symbol matches exactly as the pairwise fallback does while
  avoiding per-pair NumPy call overhead. Non-verbose, non-plot runs now use a
  stats-only matrix path that summarizes upper-triangle identities directly,
  avoiding pair-label and dictionary construction while preserving the full
  result path for verbose, JSON-row, and heatmap outputs. Fully identical
  normalized alignments now compute the constant identity from one sequence
  before matrix construction, preserving the full-length denominator and
  double-gap exclusion behavior for `--exclude-gaps`; summary-only runs return
  constant statistics directly while full-result runs still emit every pair row.
  The identical-sequence guard now scans the sequence list by index instead of
  materializing `sequences[1:]`, preserving early mismatch exits without a
  temporary list copy. Alignments with fewer than two records now return the
  existing no-pair result and no-values summary before matrix probing or
  fallback sequence-array construction.
  Mixed-length or non-ASCII sequential fallback runs now compute the pair count
  arithmetically and stream `itertools.combinations` into the batch worker,
  reserving full pair-list materialization for the multiprocessing chunking path.
  Summary-only runs no longer materialize taxa IDs before scoring, because the
  taxa label list is only needed for heatmap plotting.
- `VariableSites`, `ParsimonyInformative`, and `AlignmentEntropy` baseline time
  was dominated by per-site filtering/counting. The optimized paths uppercase
  each sequence once, identify valid symbols once per alignment, and count each
  valid symbol across all sites with vectorized NumPy operations. Variable
  sites, parsimony-informative sites, and alignment entropy now also use
  byte-backed matrix setup for normal ASCII alignments, retaining a Unicode
  fallback for non-ASCII sequence content. Later variable-sites,
  parsimony-informative, and entropy passes replaced their ASCII `np.isin` gap
  masks with byte lookup tables, retaining the same uppercase normalization and
  Unicode fallback behavior. A later pass builds those masks by OR-ing per-gap
  byte comparisons for variable-site, parsimony-informative, and entropy ASCII
  paths, avoiding the lookup-gather mask allocation while preserving
  valid-column semantics. A later variable-sites startup pass builds the
  byte lookup tables lazily and defers the annotation-only Bio.Align import
  while preserving the module-level `np.isin` patch point used by fallback
  tests. Variable-sites now tests ASCII columns for variability with valid-symbol
  min/max reductions instead of scanning once per observed symbol, which reduces
  repeated full-matrix passes for protein alphabets while retaining the Unicode
  fallback. Clean protein variable-site alignments now skip the masked min/max
  temporaries after the existing gap-code scan finds no invalid symbols,
  preserving the gappy path. Clean DNA variable-site alignments use the same
  direct min/max path once the gap-code scan finds no invalid symbols. Fully
  identical normalized alignments now return zero variable sites before
  byte-matrix construction, covering conserved multi-symbol sequences while
  leaving non-identical alignments on the existing min/max path. Raw-identical
  alignments now test equality before uppercasing every row, avoiding the
  normalization pass for already identical inputs while preserving mixed-case
  equivalence. A follow-up pass reuses the shared no-slice equality helper,
  avoiding `sequences[1:]` materialization while preserving early exit for
  heterogeneous alignments. The single-record path now returns zero variable
  sites after reading the alignment length, avoiding sequence materialization
  because one taxon cannot vary by column. The Unicode fallback now counts final
  variable columns with `np.count_nonzero` instead of summing a boolean vector.
  Parsimony-informative clean ASCII alignments pass blocks directly into the
  256-bin counter after the same gap-code scan proves there is no validity mask
  to apply, preserving the masked path for gappy alignments. Fully identical
  normalized alignments now return zero parsimony-informative sites before
  byte-matrix construction, covering conserved multi-symbol sequences while
  leaving non-identical alignments on the existing block-count path. A follow-up
  identical-row pass scans the existing sequence list directly instead of
  materializing `sequences[1:]`, reducing temporary allocation for high-taxon
  conserved alignments. A later setup pass combines sequence-list construction
  with identical-row detection, removing an extra pass over conserved inputs
  while preserving the existing block-count and Unicode fallback paths. The
  single-record path now returns zero parsimony-informative sites after reading
  the alignment length because one taxon cannot provide two recurrent states.
  Parsimony-informative ASCII columns now use chunked 256-bin column
  histograms to count symbols appearing at least twice, avoiding one full
  matrix scan per observed protein symbol while keeping memory bounded by the
  column chunk size. The public parsimony-informative per-column occurrence
  helper now counts directly from records, avoiding BioPython column-slice
  construction while preserving `Counter` output, and the predicate exits as
  soon as two recurrent states are found instead of scanning every state count.
  The Unicode fallback now counts final recurrent-state columns with
  `np.count_nonzero` instead of summing a boolean vector.
  A follow-up parsimony-informative startup pass applies the same lazy lookup
  construction and Bio.Align annotation deferral while preserving its
  module-level `np.isin` patch point. A later alignment-entropy
  startup pass applies lazy lookup construction while preserving its module-level `np.isin`
  patch point. Alignment-entropy byte histograms now encode block counts in
  column-major order, avoiding a tiled site-offset array for each block while
  preserving the same entropy values. Identical-row entropy shortcuts now scan
  the existing sequence list directly instead of materializing `sequences[1:]`,
  reducing temporary allocation for high-taxon conserved alignments. A
  follow-up startup pass keeps JSON output for the variable-sites and
  parsimony-informative commands behind module-level
  forwarding wrappers, preserving their patch points while avoiding JSON helper
  startup during command discovery. A later alignment-entropy startup pass keeps JSON output
  behind a module-level forwarding wrapper and localizes `PlotConfig` to
  argument processing, preserving patch points while avoiding helper startup
  during command discovery. A later variable-sites startup pass removes the
  runtime `TYPE_CHECKING` dependency and converts annotation-only typing aliases
  to postponed built-in annotations. A follow-up parsimony-informative startup
  pass applies the same runtime `TYPE_CHECKING` removal and built-in annotation
  conversion. A later shared-helper startup pass localizes
  `helpers.files` cache-key hashing so alignment command discovery does not
  import `hashlib`, and a later cache-key pass returns the raw path/size/mtime_ns
  key directly instead of hashing it. It converts `AlignmentEntropy`
  annotation-only typing aliases to built-in annotations so it no longer imports
  `typing`. Shared single-column list parsing first moved to a one-shot
  splitlines/strip path; a later real-file benchmark pass now streams lines
  directly, preserving whitespace stripping and blank-line behavior while
  avoiding whole-file line-list allocation for large taxa/tree/alignment lists.
  Shared alignment format detection now splits a PHYLIP-like first line at most
  once, preserving the two-column numeric-header check while avoiding duplicate
  `split()` calls on common alignment startup.
  `AlignmentEntropy.run` now batches verbose per-site text rows into one
  newline-joined print while preserving summary,
  JSON, plot, and stdout text behavior. A later nonverbose summary pass defers
  per-site row dictionaries until verbose or plot output needs them and
  computes the mean from the existing Python list with `sum() / len()` instead
  of materializing a NumPy array. A follow-up verbose text pass formats rows
  directly from the entropy list when JSON/plot output do not need structured
  row dictionaries. Clean protein ASCII entropy blocks now bypass the per-block
  valid-mask boolean index, preserving the existing gappy path while reducing
  clean 1200 x 12000 protein helper time from 0.100458s to 0.090770s. A later
  clean-protein entropy pass also skips constructing the full valid mask before
  calling that helper when no protein gap bytes are present, leaving the DNA and
  Unicode paths unchanged.
  Entropy term calculation now uses masked `np.log2` into a scratch array
  instead of `np.where`, so sparse protein count matrices avoid computing logs
  for zero-probability symbols while preserving exact entropy values.
  Mask-presence checks in shared alignment helpers now use ndarray
  `any()`/`all()` reductions instead of top-level `np.any`/`np.all` dispatch,
  preserving boolean results while reducing dispatch overhead across empty,
  small, and large masks. Alignment plotting and PhyloGWAS contingency guards
  now use the same no-axis mask method calls where no public NumPy patch point
  is required.
  Entropy matrices now reduce probability/log-probability terms with an
  `einsum` column dot for both small DNA-sized and larger protein-sized
  alphabets, avoiding the previous in-place product plus `np.sum` path.
  DNA-sized ASCII entropy counts now use `np.count_nonzero` instead of summing
  boolean equality masks while preserving the faster `np.sum` path for Unicode matrices.
  Conserved alignments with one valid observed symbol now return zero site
  entropies directly after the existing `np.unique` check, skipping count,
  probability, and log calculations while leaving multi-symbol inputs on the
  existing path. Fully identical normalized alignments now return zero site
  entropies before byte-matrix construction, covering conserved multi-symbol
  sequences while leaving non-identical alignments on the existing count path.
  The ASCII entropy path now defers the Unicode fallback's gap-character set
  construction until a `UnicodeEncodeError` occurs, so common byte-matrix inputs
  use the cached byte gap codes without also building an unused Python set.
  `_entropy_from_ascii_codes` now converts the final NumPy array with
  `ndarray.tolist()`, preserving the Python float list output while avoiding one
  Python-level conversion loop over every site.
- `SumOfPairsScore.determine_number_of_matches_and_total_pairs` baseline time
  compared every taxon pair independently. The optimized complete equal-length
  path counts taxa whose query and reference residues match at each site, then
  converts those per-site counts to matching pair counts with `n choose 2`;
  current code builds byte-backed alignment matrices for ASCII FASTA content
  and uses `np.count_nonzero` for the per-site matching-taxon counts.
  `run` now tries that complete-record path before materializing the full
  `n choose 2` reference pair list, preserving the mixed-length and incomplete
  pair-list fallback while avoiding quadratic setup for ordinary FASTA inputs.
  Exact equal-length query/reference records now return the full-match pair
  count after the existing length validation, avoiding the matrix stack entirely
  for unchanged alignments. Complete pair-set validation now streams over the
  common ordered `itertools.combinations` layout before falling back to the
  previous unordered set comparison, preserving shuffled pair-list behavior
  without allocating the expected pair set on the hot path. The complete-record
  path now totals the per-site matching-pair vector with its ndarray `sum()`
  method, avoiding an extra lazy NumPy reduction dispatch. A later ordered
  validation pass compares directly against the combinations iterator instead
  of walking nested slices of the record ID list, preserving the same unordered
  fallback while avoiding repeated slice allocation. Incomplete pair sets whose
  involved taxa have unchanged raw query/reference sequences now return
  full-match totals from the requested pair lengths before building arrays or
  entering the multiprocessing fallback.
- `SumOfPairsScore._process_pair_batch` fallback baseline time converted every
  sequence into Unicode character arrays and rebuilt trimmed arrays for
  mixed-length comparisons. The optimized path uses byte-backed sequence arrays
  for ASCII FASTA content and slices those cached arrays for trimmed
  comparisons, retaining a Unicode fallback for non-ASCII strings. The small
  sequential mixed-length path now uses the same cached array slicing instead
  of a Python character loop.
- `SumOfPairsScore._read_fasta` baseline time materialized `SeqRecord`
  objects through `SeqIO.to_dict(SeqIO.parse(...))`. The optimized loader uses
  `SimpleFastaParser` to keep plain sequence strings while preserving
  first-token IDs, mixed case, and duplicate-key `ValueError` behavior. A later
  pass reuses the shared unique first-token parser, preserving multiline
  sequence joining, internal whitespace removal, mixed case, and duplicate-key
  errors while avoiding parser tuple construction. A startup pass defers NumPy
  behind a module-level proxy so import-only callers do not pay for vectorized
  comparison setup until the scoring path needs it. Another startup pass imports
  the FASTA parser inside the read helper, so importing the sum-of-pairs command
  module does not load Bio.SeqIO.FastaIO. A follow-up startup pass keeps JSON
  output behind a module-level forwarding wrapper, preserving the patch point
  while avoiding JSON helper startup during command discovery. A later startup
  pass postpones annotations and converts annotation-only `typing` aliases to
  built-in annotations, so command discovery
  no longer loads `typing`.
- `CompositionPerTaxon` and `GCContent` baseline time was dominated by NumPy
  array construction for simple per-symbol counts. The optimized paths use
  uppercase sequence strings and C-level `str.count` calls over the small
  sequence alphabet. `CompositionPerTaxon` now uses a byte-backed alignment
  matrix for the current larger-alignment path, counting each observed symbol
  across all taxa at once while retaining a Unicode fallback. A later pass stores
  the ASCII matrix as `uint8` values and uses per-row `np.bincount` so each taxon
  row is counted once instead of comparing the full matrix once per observed
  symbol. A subsequent pass replaces ASCII invalid-character `np.isin` calls
  with a cached 256-entry lookup table while preserving the Unicode fallback.
  The current ASCII path uses symbol-wise matrix counts for small observed
  alphabets, which is faster for DNA-style inputs, while retaining per-row
  `bincount` for protein-sized alphabets where it remains cheaper. A later pass
  also routes large, short protein alignments through symbol-wise matrix counts,
  avoiding thousands of tiny per-row `bincount` calls while keeping long protein
  alignments on the row-count path. A later ASCII pass detects alignments with
  no invalid symbols and uses the observed byte codes plus full sequence lengths
  directly, avoiding full valid-mask construction for no-gap DNA and protein
  inputs while keeping the filtered-mask path for gap-bearing alignments.
  Identical alignments with multiple valid symbols now compute one composition
  vector from the normalized first sequence and copy it for each taxon, avoiding
  full alignment-matrix counting while preserving per-row arrays and output
  ordering. The identical-frequency normalization now uses the count vector's
  ndarray `sum()` method, avoiding a lazy NumPy reduction dispatch on that
  shortcut. A follow-up identical-row pass scans the existing sequence list
  directly instead of materializing `sequences[1:]`, reducing temporary
  allocation for high-taxon conserved alignments. A later identical-row pass
  scans the normalized record tuples directly before building the separate
  sequence list, preserving the matrix path for heterogeneous alignments while
  reducing conserved-input setup work.
  Raw-identical conserved rows now avoid uppercasing every record before the
  identical shortcut, while mixed-case identical rows still compare through the
  same normalized sequence check and heterogeneous rows still uppercase before
  the matrix path.
  Conserved alignments with one valid observed symbol now fill the one-column
  frequency output directly from valid sequence lengths, skipping the count
  matrix while preserving `0.0` frequencies for taxa with no valid sites.
  A later startup pass defers NumPy behind a module-level proxy and postpones
  annotations while preserving the `np` spy and patch points used by tests.
  `GCContent` now also uses byte lookup tables for ASCII sequence
  validity and GC counts, and a later pass made those tables case-insensitive
  to avoid allocating uppercase copies for ASCII rows while retaining the
  string-count path for non-ASCII sequence content. A later per-sequence ASCII
  pass scans the joined byte buffer for invalid symbols and uses full sequence
  lengths directly for no-gap rows, avoiding the separate validity lookup count
  while preserving GC byte counting. The total-GC ASCII flat path now uses the
  same no-gap valid-length shortcut after sampling both ends of the byte buffer
  for ambiguous symbols, falling back to the lookup count when invalid symbols
  are present. A variable-length ASCII
  fallback pass counts valid and GC bytes with `np.count_nonzero` instead of
  reducing boolean arrays with `np.sum`; the mixed Unicode/ASCII total-GC
  fallback now uses the same count primitive for ASCII rows. A later mixed
  fallback pass batches all ASCII fallback rows into one byte buffer and counts
  valid/GC bytes once, while retaining the existing uppercase string counting
  for non-ASCII rows. Verbose GC text output now batches per-sequence rows into
  one newline-joined print while preserving the same stdout text, summary
  output, JSON payloads, and broken-pipe handling. Identical sequences now
  count valid and GC characters from the normalized first sequence and reuse
  that result for both verbose per-sequence rows and aggregate GC, avoiding
  full byte-matrix construction while preserving invalid-symbol handling. A
  later helper pass keeps that normalized identical-row count in bytes for
  ASCII sequences and avoids re-uppercase work, while retaining the existing
  non-ASCII string-count fallback semantics. A
  follow-up identical-row pass checks raw string equality before normalizing
  each row, avoiding repeated uppercase allocations for already-identical rows
  while preserving case-insensitive matching. A
  later startup pass builds lookup tables lazily and
  defers the annotation-only Bio.Align import so import-only callers do not
  load NumPy. A follow-up startup pass keeps JSON output behind a module-level
  forwarding wrapper, preserving the patch point while avoiding JSON helper
  startup during command discovery. A later GC-content startup pass removes the
  runtime `TYPE_CHECKING` dependency and converts annotation-only typing aliases
  to postponed built-in annotations. `GCContent` now also preserves its
  module-level alignment-reader patch point with a lazy forwarding wrapper, so
  command discovery avoids `phykit.helpers.files` once the shared alignment base
  import is also deferred. `CompositionPerTaxon` now uses the same
  forwarding-wrapper pattern for JSON output, preserving its patch point while
  avoiding JSON helper startup during command discovery. A later composition
  startup pass converts its annotation-only `typing` aliases to built-in
  postponed annotations, so command discovery no longer loads `typing`.
- `AlignmentLength.run` baseline time materialized the full Biopython
  alignment only to read its width. FASTA inputs now validate row lengths with
  `SimpleFastaParser` and return the shared alignment length directly, falling
  back to the existing parser for non-FASTA, empty, or inconsistent files. A
  later startup pass imports the FASTA parser inside the fast-path helper,
  avoiding Bio.SeqIO.FastaIO startup for import-only callers. A follow-up
  startup pass keeps JSON output behind a module-level forwarding wrapper,
  preserving the patch point while avoiding JSON helper startup during command
  discovery. A later startup pass removes the annotation-only `typing` import by
  using postponed built-in annotations. A later direct-scanner pass removes
  `SimpleFastaParser` from the common FASTA length path entirely, preserving its
  wrapped-sequence and internal-space length behavior while avoiding tuple
  yields and full sequence string construction. A follow-up scanner pass reads
  the same fast path as bytes, using C-level ASCII checks and falling back to
  the full parser for non-ASCII sequence lines to preserve compatibility.
- `AlignmentLengthNoGaps.get_sites_no_gaps_count` baseline time constructed a
  Unicode character matrix and called `np.isin` once per alignment column. The
  optimized path uppercases sequences once, builds a byte-backed alignment
  matrix, computes one alignment-wide invalid-character mask, and counts
  columns without invalid symbols from that mask. A later pass replaced the
  ASCII `np.isin` mask with a case-insensitive byte lookup table so gap
  detection avoids per-row uppercase copies and per-symbol membership checks.
  A later pass detects gap-bearing columns by OR-ing per-gap byte comparisons,
  avoiding construction of a full taxa-by-sites boolean lookup mask while
  preserving DNA/protein ambiguity semantics. Clean ASCII alignments now first
  scan the byte stream for gap/ambiguous byte codes and return the full
  alignment length immediately when none are present, while gap-bearing inputs
  fall through to the existing column-reduction path. A follow-up no-gap pass
  performs that precheck with module-level byte constants before resolving lazy
  NumPy gap-code arrays, and now checks those byte constants with direct integer
  byte membership instead of allocating one-byte probes. Fully identical
  alignments now count gap-free sites from the first sequence alone, preserving
  DNA/protein gap and ambiguity semantics while avoiding alignment-wide byte
  joins and matrix construction for unchanged inputs. The identical-sequence
  guard now uses the same index-based scan as `PairwiseIdentity`, avoiding
  `sequences[1:]` allocation while retaining early mismatch exits. The
  identical-sequence Unicode fallback now counts uppercase and lowercase gap
  codes directly instead of allocating an uppercased copy before counting.
  A subsequent startup pass builds those lookup tables lazily and defers the
  annotation-only Bio.Align import. The shared `alignment.base` RCV lookup
  tables now use the same lazy NumPy construction, preserving the module-level
  `np.isin` patch point used by fallback-path tests. A follow-up startup pass
  keeps JSON output behind a module-level forwarding wrapper, preserving the
  patch point while avoiding JSON helper startup during command discovery. A
  later startup pass removes the runtime `TYPE_CHECKING` dependency and converts
  annotation-only typing aliases to postponed built-in annotations. The shared
  alignment reader now uses a lazy forwarding wrapper, preserving the base-class
  helper path while avoiding `phykit.helpers.files` during import-only command
  discovery.
- `MaskAlignment.calculate_keep_mask` baseline time constructed a Unicode
  character matrix and, when entropy filtering was enabled, called `np.unique`
  once per column. The optimized path uses a byte-backed alignment matrix and
  computes entropy from one symbol-by-site count matrix while preserving
  all-invalid sites as zero-entropy columns. A later pass replaced the ASCII
  `np.isin` validity mask with a byte lookup table while preserving the Unicode
  fallback path. For protein-sized ASCII alphabets, a later pass computes
  entropy from block-wise byte `bincount` tables while retaining the simpler
  comparison path for DNA-sized alphabets. Protein entropy block counts now use
  column-major encoded byte indices, avoiding a tiled site-offset array for each
  block while preserving the same entropy thresholds. When entropy is disabled
  and the thresholds accept every possible column, `calculate_keep_mask` now
  returns an all-true mask directly instead of building occupancy arrays. A
  follow-up moves that all-pass return before sequence materialization because
  the result only depends on the alignment length. Clean ASCII alignments with
  no entropy threshold now also return that all-true mask
  after a byte scan proves every column has occupancy 1.0 and gap fraction 0.0,
  preserving the slower occupancy path for gappy or entropy-filtered inputs.
  Clean ASCII alignments with an entropy threshold skip occupancy averaging and
  feed unmasked blocks into the entropy counter, preserving the validity-mask
  path for gap-bearing entropy filters. Entropy block calculations now use
  masked `np.log2` with `out`/`where` instead of boolean-indexing the
  probability array, avoiding temporary positive-probability slices.
  Protein-sized entropy matrices now reduce probability/log-probability terms
  with an `einsum` column dot while DNA-sized matrices keep the previous
  in-place product path. DNA-sized
  ASCII entropy counts now use `np.count_nonzero` instead of summing boolean
  equality masks while preserving the faster `np.sum` fallback for Unicode
  matrices. When only one valid symbol is observed, entropy-filtered runs now
  keep the zero-entropy columns directly instead of constructing count and
  probability arrays.
  Fully identical normalized alignments now derive occupancy and zero entropy
  from one sequence before matrix construction, preserving gap, occupancy, and
  entropy thresholds for conserved multi-symbol inputs. A follow-up pass uses
  the shared no-slice equality helper for that shortcut, avoiding
  `sequences[1:]` materialization while preserving early exit for heterogeneous
  alignments. The identical-sequence valid-site mask now uses the shared ASCII
  invalid-byte lookup, with the prior Python generator retained for Unicode
  fallback.
  `MaskAlignment.apply_mask` also filters a single alignment matrix now, instead
  of converting and masking each taxon's sequence independently. When every site
  is kept, `apply_mask` now returns uppercased sequence strings directly instead
  of materializing and slicing that matrix. When no sites are kept and the mask
  length matches all sequence lengths, it returns empty sequence strings
  directly while preserving the matrix path for mismatched masks. Text-mode `run` now emits the masked FASTA records as one
  newline-joined print while preserving record order, JSON output, and stdout
  text. A later startup pass builds the lookup tables lazily and postpones
  annotations so import-only callers do not load NumPy while fallback tests can
  still patch module-level `np.isin`. A follow-up startup pass keeps JSON output
  behind a module-level forwarding wrapper, preserving the patch point while
  avoiding JSON helper startup during command discovery. A later startup pass
  removes the annotation-only `typing` import under postponed annotations.
- `EvolutionaryRatePerSite.calculate_evolutionary_rate_per_site` baseline time
  filtered each column and called `np.unique` per site. The optimized path
  builds one byte-backed alignment matrix, counts each valid symbol across all
  sites once, and computes PIC values from vectorized per-site frequency sums.
  A later pass replaced the ASCII `np.isin` gap mask with a byte lookup table
  while preserving the Unicode fallback path. A subsequent startup pass builds
  those lookup tables lazily and defers the annotation-only Bio.Align import
  while preserving the module-level `np.isin` patch point used by fallback
  tests. For protein-sized ASCII alphabets, a later pass switches from one full
  matrix comparison per valid symbol to block-wise byte `bincount` totals and
  sum-of-squares, while retaining the simpler comparison path for DNA-sized
  alphabets. The protein block counter now encodes columns in site-major order,
  avoiding a tiled site-offset vector for each block while preserving identical
  per-site totals and sum-of-squares. Text-mode `run` now batches per-site rows
  into one newline-joined print, preserving the same stdout text and leaving
  JSON/plot payloads unchanged. A later terminal-output pass skips the intermediate row
  dictionaries unless JSON or plot output needs them. A follow-up startup pass
  defers the JSON helper behind a local wrapper and imports `PlotConfig` only
  while processing command arguments, so import-only callers avoid those helper
  modules. A later startup pass removes the runtime `TYPE_CHECKING` dependency
  and converts annotation-only typing aliases to postponed built-in annotations.
  Conserved alignments with one valid observed symbol now return zero PIC values
  directly after the existing valid-symbol discovery step, skipping count and
  frequency-sum work while preserving all multi-symbol paths. Fully identical
  normalized alignments now return zero PIC values before byte-matrix
  construction, covering conserved multi-symbol sequences while leaving
  non-identical alignments on the existing count path. Raw-identical alignments
  now test equality before uppercasing every row, avoiding the normalization pass
  for already identical inputs while preserving mixed-case equivalence. A
  follow-up pass reuses the shared no-slice equality helper so this shortcut no
  longer materializes `sequences[1:]` and still exits early for heterogeneous
  alignments.
  The public gap-removal helper now reuses cached `str.translate` deletion
  tables for single-character gap lists, preserving the previous fallback for
  unusual multi-character gap tokens. The public per-column occurrence helper
  now counts directly from records for single-character gap lists, avoiding
  BioPython column-slice construction while preserving the existing
  multi-character gap-token fallback. A follow-up extends the direct record-wise
  occurrence counter to multi-character gap-token configurations too, avoiding
  column-slice construction while preserving the same token semantics. A later
  ASCII validity pass derives DNA
  valid symbols from the observed byte-code set and lets no-gap protein block
  counts skip the full valid mask, while protein alignments containing gap
  symbols keep the filtered-mask path. Column count sum-of-squares reductions
  now use an `einsum` column dot, avoiding temporary squared count matrices in
  ASCII block counters and Unicode fallback counters.
- `CompositionalBiasPerSite.calculate_compositional_bias_per_site` baseline
  time filtered each column, called `np.unique`, and ran `chisquare` per site.
  The optimized path uses the same symbol-by-site count matrix pattern to
  compute chi-square statistics and p-values in bulk while preserving `"nan"`
  correction slots for all-invalid and single-symbol columns. A later pass
  replaced the ASCII `np.isin` gap mask with a byte lookup table while
  preserving the Unicode fallback path. A subsequent pass removed the eager
  `scipy.stats` import by using local integer-degree chi-square survival
  calculations and a local one-dimensional Benjamini-Hochberg correction. A
  later chi-square p-value pass vectorizes the square-root step before applying
  the scalar complementary-error function loop, preserving the no-SciPy
  calculation while reducing df=1 p-value overhead. The corrected p-value
  reconstruction now preallocates the full result instead of repeatedly
  inserting `"nan"` slots, preserving interleaved site positions in linear time.
  When every site has a valid corrected p-value, the restore helper now returns
  the corrected list directly instead of allocating and filling an identical
  no-`"nan"` result.
  Column count sum-of-squares reductions now use an `einsum` column dot,
  avoiding temporary squared count matrices in ASCII block counters and Unicode
  fallback counters.
  Text-mode `run` now batches per-site rows into one
  newline-joined print while preserving stdout text, JSON payloads, and plot
  reporting. A later terminal-output pass skips the intermediate row
  dictionaries unless JSON or plot output needs them. A later startup pass
  builds the byte lookup tables lazily and defers the annotation-only Bio.Align
  import while preserving the module-level `np.isin` patch point used by
  fallback-path tests. A follow-up startup pass keeps JSON output behind a
  module-level forwarding wrapper and localizes `PlotConfig` to argument
  processing, preserving patch points while avoiding helper startup during
  command discovery. Small p-value correction lists now use the same
  Benjamini-Hochberg ordering and reverse cumulative-minimum semantics in
  Python, avoiding NumPy startup for direct helper use while retaining the
  vectorized NumPy path for large site-wide correction arrays. A later startup
  pass converts remaining annotation-only `typing` names to built-in postponed
  annotations, so command discovery no longer imports `typing` for this module.
  Protein ASCII count blocks now encode valid observations in site-major order,
  avoiding the tiled site-offset vector while preserving the same category
  counts, per-site totals, and sum-of-squares. The public per-column occurrence
  helper now counts directly from records, avoiding BioPython column-slice
  construction while preserving first-seen count order. A later ASCII validity
  pass derives DNA valid symbols from the observed byte-code set and lets
  protein block counts skip the full valid mask when a byte-level scan finds no
  protein gap symbols, while gap-containing protein alignments keep the filtered
  mask path. Conserved alignments with one valid observed symbol now return zero
  statistics and `"nan"` corrected p-values directly after the existing
  valid-symbol discovery step, skipping count, statistic, p-value, and FDR work.
  Fully identical normalized alignments now take the same zero-statistic and
  `"nan"` p-value result before matrix construction, covering conserved
  multi-symbol sequences while preserving the existing non-identical paths. A
  follow-up identical-row pass scans the existing sequence list directly instead
  of materializing `sequences[1:]`, reducing temporary allocation for high-taxon
  conserved alignments. All-valid p-value correction now converts the p-value
  array directly to a list instead of first applying an all-true boolean mask,
  preserving mixed-NaN filtering for gappy/all-invalid sites. Per-site
  `Power_divergenceResult` rows are now built
  with a localized list comprehension, preserving float coercion and namedtuple
  output while avoiding the explicit Python append loop.
- `RelativeCompositionVariabilityTaxon.calculate_rows` baseline time built a
  Unicode character matrix and populated per-taxon composition counts through
  nested taxon/symbol loops. The optimized path uses a byte-backed alignment
  matrix and computes each valid symbol's counts across all taxa in one vector
  operation before applying the unchanged RCVT formula. A later pass replaced
  the ASCII `np.isin` validity mask with a byte lookup table while preserving
  the Unicode fallback path. A protein-alphabet follow-up uses per-row byte
  `bincount` histograms when the ASCII valid-symbol set is large, avoiding a
  full alignment-matrix comparison for every amino-acid symbol while leaving
  the small-alphabet comparison path in place. The observed-symbol validity
  pass derives valid ASCII symbols from byte codes present in the alignment and
  uses full sequence lengths when no invalid symbols are observed, avoiding the
  full validity mask for no-gap inputs and reducing setup work for gap-bearing
  alignments. Gap-bearing paths now count row valid lengths with
  `np.count_nonzero` instead of summing boolean masks before float conversion.
  Many-short no-gap protein alignments now build the per-taxon symbol count
  table with one row-offset `bincount`, avoiding tens of thousands of tiny
  row-level `bincount` calls while keeping the per-row path for medium or long
  alignments where it remains faster. Count-matrix column totals now use the
  ndarray reduction method, avoiding generic `np.sum` dispatch while leaving
  the mixed-performance row sums unchanged. Identical alignments now return
  zero-valued per-taxon rows before matrix construction because each taxon's
  composition equals the average composition after case normalization.
  Raw-identical alignments now test equality before uppercasing every row,
  avoiding the normalization pass for already identical inputs while preserving
  mixed-case equivalence. A follow-up pass reuses the shared no-slice sequence
  equality helper to avoid materializing `sequences[1:]` while preserving
  early exit for heterogeneous alignments. Single-record alignments now return
  the known zero-valued row
  before sequence materialization. Text-mode `run` now batches terminal
  rows into one print, and plot setup now sorts an index array over extracted
  RCVT values instead of sorting full row dictionaries while preserving stable
  descending order for tied values.
  per-taxon output rows into one newline-joined print while preserving the same
  stdout text; JSON and plot reporting are unchanged. A subsequent startup pass
  builds the lookup tables lazily while preserving the module-level `np.isin`
  patch point used by fallback-path tests. A follow-up startup pass keeps JSON
  output behind a module-level forwarding wrapper and localizes `PlotConfig` to
  argument processing, preserving patch points while avoiding helper startup
  during command discovery.
- `AlignmentOutlierTaxa.calculate_outliers` baseline time combined a Unicode
  character matrix with nested composition-count loops, per-taxon/per-taxon
  distance loops, and per-column entropy counting. The optimized path uses a
  byte-backed alignment matrix, vectorized symbol counts for composition and
  entropy burden, and per-taxon vectorized pairwise distance scans while
  preserving no-overlap distances as `None` in public rows. A follow-up stores
  ASCII alignments as `uint8` and builds the validity mask with a cached
  256-entry invalid-character lookup, retaining the Unicode fallback. A later
  pass computes all pairwise overlap and match counts with BLAS-backed matrix
  products for moderate taxon counts, keeping the original row scan as a memory
  guard for oversized pairwise matrices. The matrix path now counts comparable
  pair rows with `np.count_nonzero` instead of summing boolean rows; the full
  public method moved from 0.734118s to 0.008124s on the 180 x 1200 benchmark
  while preserving no-overlap distances as `None`. Row assembly now zips the feature arrays and
  reuses each scalar for threshold checks, rounded public rows, and nested
  reason payloads, avoiding repeated NumPy indexing while preserving reason
  order and `None` no-overlap distances. Text-mode `run` now batches the header,
  thresholds, and outlier rows into one newline-joined print while preserving
  the same stdout text and JSON output. A later startup pass defers NumPy and
  JSON helper imports behind module-level proxy/wrapper objects while
  preserving existing test patch points. A follow-up startup pass converts
  annotation-only `typing` aliases to built-in postponed annotations, so command
  discovery no longer loads `typing`. Conserved alignments with one valid
  observed symbol and equal valid lengths now return the constant feature result
  directly after valid-symbol discovery, skipping composition, pairwise distance,
  RCVT, entropy, threshold, and reason-building work while leaving uneven-length
  single-symbol alignments on the normal path. Alignments with fully identical
  normalized sequences now compute the same constant feature result from the
  first sequence before alignment matrix construction, covering conserved
  multi-symbol alignments while leaving all-invalid identical alignments on the
  existing normal path to preserve `None` no-overlap distance rows. That
  identical-sequence guard now uses the shared index-based scanner instead of
  allocating `sequences[1:]`, preserving early mismatch short-circuiting without
  list-copy overhead. The identical-sequence Unicode valid-length fallback now
  counts invalid characters with `str.count` instead of looping over every
  character in Python. Single-record alignments now use the same constant-result
  formatter after computing the one sequence's valid length, skipping matrix
  construction and threshold work while preserving `None` long-branch proxy
  rows.
  Variable-composition ASCII alignments with no invalid symbols now use a
  byte-scan guard to skip invalid lookup construction, validity-mask allocation,
  and validity-weighted entropy averaging while preserving the existing mask
  path for ambiguous or gapped data. Gapped paths now count row valid lengths
  with `np.count_nonzero` instead of summing boolean masks before float
  conversion. A later protein-oriented pass counts ASCII symbols by row/site
  with `bincount` when the symbol-count setup is large enough to beat repeated
  equality reductions. The row-count helper now also keeps that encoded
  histogram path for large taxon counts with short alignments, avoiding one
  Python `bincount` loop per taxon while leaving longer matrices on the existing
  row loop. The site-count helper now keeps encoded histograms through 32M
  cells, avoiding repeated equality scans on moderate-large protein matrices.
  The all-valid ASCII long-branch proxy now reuses per-site
  symbol counts to compute each taxon's exact mean distance to all other taxa,
  avoiding the previous taxon-by-taxon pairwise match matrix while preserving
  the gapped and Unicode fallback paths. Composition-distance setup now computes
  row L2 norms with an `einsum` reduction instead of `np.linalg.norm`, avoiding
  linalg dispatch while preserving the same rounded row distances. Entropy
  burden now computes per-site probability/log-probability products with an
  `einsum` column dot, avoiding a temporary product matrix while preserving
  entropy values.
- `PlotAlignmentQC` composition-distance scatter panel baseline time called
  Matplotlib `scatter` once per taxon, creating thousands of artists for large
  alignments. The optimized path batches normal and flagged taxa into two
  scatter calls and still annotates flagged taxa individually. Plot setup now
  extracts taxa, feature arrays, long-branch NaNs, and flagged masks once and
  reuses them across occupancy/gap bars, the composition-distance panel, and
  heatmap z-score setup. The flagged mask is now filled during that row pass,
  removing the previous second `np.fromiter` pass over every taxon. Occupancy
  and gap bar colors are now selected with vectorized mask mapping, and the two
  bar panels reuse one x-position array. A later
  startup pass defers NumPy, JSON output, plot config, and the alignment-outlier
  service behind lazy proxy/wrapper functions while preserving the existing
  `AlignmentOutlierTaxa` monkeypatch point. A follow-up startup pass converts
  annotation-only `typing` aliases to postponed built-in annotations, so command
  discovery no longer loads `typing`.
- `ColumnScore.get_columns_from_alignments` baseline time constructed Unicode
  character matrices for query and reference alignments, then joined each
  matrix column into a string. The optimized path uppercases sequence strings
  once and transposes them with `zip`, preserving the existing list-of-column
  strings API with far less allocation overhead. A later pass lets
  `ColumnScore.run` compare ASCII alignment columns as byte-matrix records
  directly, preserving the set-intersection column-score semantics and falling
  back to the public column-list helper for non-ASCII or ragged inputs. A later
  same-object pass handles self-comparisons by uniquing one byte matrix, and
  `run` reuses a single parsed alignment when query and reference paths match.
  Separate query/reference alignments with identical normalized sequences now
  use the same one-matrix unique-column path, preserving duplicate-column
  semantics while avoiding the second unique set and intersection. Repeated-row
  ASCII alignments now count the unique symbols in the repeated sequence, or
  the symbol-set intersection for separate repeated-row alignments, preserving
  the same unique-column semantics without building transposed byte matrices.
  A follow-up helper pass scans repeated rows without materializing
  `sequences[1:]`, which removes a large temporary list for high-taxon repeated
  alignments. Separate alignments with different nonzero taxon counts now return
  the existing zero-match result from the query alignment length before sequence
  materialization.
  A startup pass defers NumPy, Bio.AlignIO, annotation-only Bio.Align, and JSON
  helper imports behind module-level proxies/wrappers while preserving existing
  test patch points. A follow-up startup pass removes the runtime
  `TYPE_CHECKING` dependency and converts annotation-only typing aliases to
  postponed built-in annotations.
- `DNAThreader.normalize_n_seq` baseline time sliced a Biopython `Seq` object
  for every codon, creating many temporary `Seq` instances. The optimized path
  converts the nucleotide sequence to a plain string once and slices that string
  for codon lookup, preserving the existing gap and missing-codon behavior.
- `DNAThreader._thread_sequence` baseline time rebuilt Unicode NumPy arrays and
  boolean masks for every sequence. The optimized path applies the same masking
  and stop-codon restoration rules directly over strings and the existing keep
  mask, avoiding per-sequence array allocation.
- `DNAThreader._thread_sequence` still materialized the full normalized
  nucleotide sequence even though masking can only inspect the prefix bounded by
  `keep_mask`. The optimized codon-aligned path builds only the reachable
  normalized prefix, keeps the non-triplet DNA fallback, and preserves the same
  stop-codon restoration behavior.
- `DNAThreader._thread_sequence` default no-ClipKIT runs use an all-true
  `keep_mask`. The optimized path now returns the bounded normalized prefix
  directly for that case after validating the inspected mask, avoiding a second
  per-position output assembly pass. A later pass has `thread()` compute the
  all-true mask status once and pass it into each `_thread_sequence` call,
  avoiding repeated mask scans across many records.
- `DNAThreader._thread_sequence` plain all-keep protein prefixes without gap,
  unknown, or stop symbols map directly onto the nucleotide prefix. The
  optimized path returns that nucleotide prefix directly when the DNA is long
  enough, while falling back to normalized-prefix assembly for gapped,
  stop-containing, or short-DNA cases. A later masked-prefix pass uses
  `itertools.compress` for the same no-gap/no-stop direct nucleotide prefix
  branch, avoiding Python generator filtering over each nucleotide. A follow-up
  pass replaces repeated Python-level gap-symbol generator scans with
  C-backed string membership checks over the inspected protein prefix.
- A subsequent all-keep `_thread_sequence` pass uses a dedicated gappy fallback
  that emits bounded nucleotide spans or gap spans by 9-site groups, avoiding
  the older codon-chunk loop when every site is kept. `thread()` now also skips
  precomputing the shared mask plan when the mask is all true; in the benchmark
  shape this avoids an additional 4.518145s per 5k plan builds that no-ClipKIT
  runs do not need. The all-keep helper now groups consecutive gap and non-gap
  amino-acid runs, emitting dash spans and nucleotide spans in chunks while
  preserving the same short-DNA padding behavior.
- `DNAThreader._create_thread_mask_plan` now detects full-keep and full-trim
  9-site groups directly from their bit tuples, avoiding per-group keep-index
  list construction and `itemgetter` setup. Full-keep emitters append the
  nucleotide segment directly, preserving the existing shared-plan output for
  mixed ClipKIT masks. Full-trim non-gap groups now advance the nucleotide
  offset without slicing the 9-character segment that will be discarded,
  preserving codon consumption before later kept groups.
- A follow-up repeated-mask pass attaches uniform-selector metadata to shared
  group-emitter plans when every 9-site amino-acid group uses the same ClipKIT
  mask. `_thread_sequence` then groups consecutive gap and non-gap protein runs,
  filters whole nucleotide spans with the repeated selector, and emits gap spans
  in chunks while preserving the existing varied-mask emitter path and short-DNA
  padding fallback. A subsequent prefix-uniform pass applies the same run
  grouping to the fully uniform prefix when only the final emitter group is
  partial, preserving the existing tail-group padding and stop-handling behavior.
- `DNAThreader._thread_sequence` full-length masks previously copied the active
  mask window while checking whether all sites were kept, then copied it again
  for masked output assembly. The optimized path uses the existing full-length
  mask directly and indexes the mask in the output loop, preserving the shorter
  mask fallback. A later repeated-mask pass has `thread()` precompute grouped
  per-codon keep-bit plans once for the shared ClipKIT mask and pass them into
  `_thread_sequence`, avoiding repeated mask indexing and repeated amino-acid
  lookup in gapped protein paths while preserving direct helper behavior
  without a plan.
- A later `_thread_sequence` pass reuses a class-level gap-character set and
  skips kept-index tracking unless terminal stop-codon restoration can run,
  avoiding one list append per kept nucleotide in the common no-terminal-stop
  masked-output path.
- A subsequent `_thread_sequence` pass filters the nucleotide prefix directly
  for no-gap, no-stop mixed-mask inputs and builds no-stop gapped output from
  codon groups without first materializing the normalized nucleotide prefix.
- A follow-up no-stop `_thread_sequence` pass replaces the amino-acid loop plus
  partial-repeat tail with direct reachable codon-chunk iteration, preserving
  the same gapped masked output while trimming repeated branch work.
- A terminal-stop `_thread_sequence` pass now emits the reachable codon-aligned
  chunks directly, avoiding normalized protein/nucleotide prefix strings and
  redundant kept-index tracking in the stop-restoration branch.
- `DNAThreader.create_mask` baseline time read and split the ClipKIT log twice:
  once for the empty-log guard and again while expanding keep/remove rows. The
  optimized path stores the parsed rows locally and expands the same data once.
- A later `DNAThreader.create_mask` pass removes the intermediate parsed-row
  list entirely. It streams the ClipKIT log, checks the second status token's
  boundary directly, and expands each status into three nucleotide-mask booleans.
- A follow-up `create_mask` pass reads the ClipKIT log as bytes so the streaming
  parser avoids per-line text decoding while preserving the same second-token
  `keep` boundary checks, including final rows without trailing newlines.
- A later `create_mask` pass expands each row with shared true/false codon-mask
  triplets instead of allocating a fresh `(keep, keep, keep)` tuple for every
  ClipKIT row, preserving the returned plain boolean list. `thread()` now uses a
  private mask builder that returns the all-sites-kept flag while expanding the
  ClipKIT log, avoiding a second full `all(mask)` scan for all-keep logs while
  preserving the public `create_mask()` return type and monkeypatch fallback.
  A later parser pass compares the bounded four-byte status field directly,
  avoiding per-row `startswith` dispatch while preserving the exact existing
  rule that `keep` must start immediately after the first space.
  The compatibility `clipkit_log_data` property now streams rows from the file
  handle instead of materializing `readlines()` before splitting, preserving the
  returned parsed row lists with lower peak memory.
- `DNAThreader.run` text output now joins threaded FASTA records into one
  printed block, preserving the same header/sequence ordering, JSON payload, and
  pipe behavior while avoiding two print calls per record. A later startup pass
	  defers Bio.SeqIO and annotation-only Bio.Seq imports behind a small module
	  proxy, preserving `SeqIO.parse` and `SeqIO.to_dict` behavior for runtime calls
	  and tests that patch those functions. A follow-up startup pass keeps JSON
	  output behind a module-level forwarding wrapper, preserving the patch point
	  while avoiding JSON helper startup during command discovery. A later startup
	  pass removes the runtime `TYPE_CHECKING` dependency and converts
	  annotation-only typing aliases to postponed built-in annotations.
- `CreateConcatenationMatrix._process_alignment_file` baseline time scanned
  the full `SeqRecord` list once for every requested taxon. The optimized path
  uses `SimpleFastaParser`, builds a first-seen record-ID sequence map once,
  and fills requested taxa from that map, preserving first-token IDs and the
  previous duplicate-ID behavior. A later all-present pass returns the parsed
  first-seen sequence map directly when every requested taxon is present,
  avoiding a redundant requested-taxa fill pass.
- `CreateConcatenationMatrix.get_list_of_taxa_and_records` baseline time
  materialized `SeqRecord` objects in the sequential two-alignment path. The
  optimized helper uses `SimpleFastaParser` and lightweight record-like objects,
  preserving first-token IDs, sequence strings, and duplicate record order. A
  later pass stores those parsed records in a slots-backed lightweight object,
  reducing per-record allocation and attribute lookup overhead while keeping the
  `record.id` and `record.seq` helper contract. A later startup pass imports the
  FASTA parser inside the FASTA-reading helpers, avoiding Bio.SeqIO.FastaIO
  startup for import-only callers. A later FASTA parser pass reuses the shared
  first-token parser helper and constructs the slots-backed records directly,
  preserving duplicate record order while avoiding `SimpleFastaParser` overhead.
- `MemoryEfficientAlignmentProcessor.calculate_column_stats_streaming` baseline
  reparsed the FASTA stream once for dimensions and then once per alignment
  column. The optimized helper accumulates per-column unique characters and gap
  counts during a single stream pass, preserving the existing variable-site,
  gap-count, and conservation outputs.
- `StreamingFastaReader.get_sequence_count` baseline time iterated every FASTA
  line through Python. The optimized helper scans the memory map for header-line
  markers, preserving the same count semantics for non-empty FASTA files while
  avoiding per-line iteration.
- `NumpyParallel.parallel_pairwise_operation` with explicit `num_workers=1`
  now fills the result matrix directly instead of materializing all pair payloads
  and result tuples before writing the matrix. The default parallel path still
  dispatches through `ParallelProcessor.parallel_map`. Worker selection now
  caches the process CPU count after the first lookup, preserving lazy
  multiprocessing startup while avoiding repeated `cpu_count()` calls across
  repeated batch setup. `ParallelProcessor.parallel_reduce` without an initial
  value now consumes the mapped results through an iterator after taking the
  first value, avoiding a near-full-list `results[1:]` copy before reduction.
  `NumpyParallel.parallel_apply_along_axis` now passes `range` index streams to
  `parallel_map` instead of eagerly materializing `list(range(...))`, preserving
  ordered row/column application while avoiding large temporary index lists.
- `CreateConcatenationMatrix._compute_effective_occupancy` baseline time
  checked every concatenated character in Python when threshold filtering was
  enabled. The optimized path counts each fixed invalid symbol with C-level
  `str.count` calls and derives informative positions from the total length.
- `CreateConcatenationMatrix._compute_effective_occupancy` still counted invalid
  symbols independently for every gene block. The optimized path joins each
  taxon's parts once and performs one invalid-symbol count pass per taxon,
  reducing repeated block-loop overhead for many-gene concatenations. A later
  pass uses ASCII byte deletion to count informative positions in C while
  retaining the character-wise fallback for non-ASCII sequence content.
- `CreateConcatenationMatrix._plot_concatenation_occupancy` baseline time filled
  the plot state matrix one character at a time. The optimized path builds state
  values for each sequence block with a byte lookup table, retains a Unicode
  fallback, and leaves matplotlib rendering behavior unchanged.
- A later occupancy-matrix pass batches each gene's equal-length ASCII sequence
  blocks into one byte lookup and scatters the resulting state rows into the
  matrix, while preserving the per-sequence fallback for Unicode and unexpected
  block lengths. Gene-boundary lines in the occupancy plot now render through
  one `LineCollection` with the same x-data/y-axes transform as `axvline`,
  preserving boundary placement while avoiding one Matplotlib line artist per
  gene boundary. The represented-occupancy sort now counts matching state cells
  with `np.count_nonzero(..., axis=1)` instead of summing the boolean mask,
  preserving row order while reducing the NumPy reduction cost.
- `CreateConcatenationMatrix.add_to_occupancy_info` rebuilt and sorted the full
  taxa set for every occupancy row. The optimized main concatenation workflow
  reuses the already sorted taxa list and total taxon count while keeping direct
  helper calls compatible with the previous sorted unique missing-taxa output.
- `CreateConcatenationMatrix.process_taxa_sequences` rebuilt the global taxa set
  for every sequentially processed alignment. The optimized main concatenation
  workflow reuses one cached taxa set while preserving direct helper-call
  compatibility. A later helper pass builds the present-taxa set while appending
  present sequences, avoiding a separate record scan while preserving duplicate
  record handling. For the common unique sorted-taxa list, a later pass scans
  that existing order to append missing blocks instead of materializing a set
  difference; duplicate-taxa direct calls retain the unique set-difference
  fallback. A follow-up pass reuses the parser's already-computed present-taxa
  set in the sequential workflow and returns immediately when all taxa are
  present, avoiding duplicate set insertions and missing-taxa scans. A later
  pass reuses the same slots-backed parsed record representation in the
  sequential workflow, reducing per-record append-loop attribute overhead without
  changing duplicate record handling or direct `SeqRecord` helper compatibility.
  The append loop now uses that parser-created record type to skip per-record
  sequence coercion on the internal path while retaining string conversion for
  external or mixed record-like helper inputs. A later workflow pass computes
  each alignment's ordered missing-taxa list once and reuses it for both missing
  sequence insertion and occupancy-row formatting. A follow-up multi-alignment
  path caches the per-taxon sequence-list handles once, avoiding repeated
  `defaultdict` lookups while appending each ordered parsed alignment result.
  Threshold filtering now batches the exclusion summary and per-taxon effective
  occupancy rows into one newline-joined print while preserving exact stdout
  text. FASTA output now accumulates bounded row chunks and writes each chunk as
  one string, preserving record order and text while reducing per-record file
  write overhead without materializing the full output. Command-module import now
  defers `textwrap.dedent` until the verbose start-message path is used.
- `CreateConcatenationMatrix._get_taxa_from_alignment` baseline time
  materialized `SeqRecord` objects while collecting only FASTA IDs for the
  global taxon list. The optimized path uses `SimpleFastaParser` and preserves
  the first-token ID behavior for headers with descriptions. A later pass uses
  the shared header-only parser, avoiding sequence assembly entirely when only
  taxon IDs are needed. A startup pass makes the occupancy-state NumPy lookup
  table lazy and defers concurrency, JSON, and plot config helper imports until
  their command paths need them, while preserving module-level patch points. A
  follow-up startup pass converts annotation-only `typing` aliases to built-in
  postponed annotations, so command discovery no longer loads `typing`. The
  alignment-list reader now uses a local forwarding wrapper, preserving the
  module-level patch point while avoiding the last eager service-level
  `phykit.helpers.files` import during command discovery.
- `TaxonGroups._extract_taxa` FASTA mode baseline time materialized full
  `SeqRecord` objects before reading only their IDs. The optimized path uses
  Biopython's lightweight `SimpleFastaParser` and preserves SeqIO's first-token
  record ID behavior for headers with descriptions. A later pass uses the shared
  header-only parser, preserving the same first-token IDs while skipping
  sequence assembly.
- `TaxonGroups._extract_taxa` tree mode still parses Newick input through
  Bio.Phylo, but terminal-name collection now reuses the shared direct
  traversal instead of materializing terminal clade lists through
  `get_terminals()`. Text-mode `run` now builds the group report as one
  newline-joined string while preserving the same stdout text and JSON output. A
  follow-up pass groups each file directly after taxon extraction, avoiding the
  intermediate path-to-taxa dictionary and second grouping pass. A later startup
  pass keeps JSON output behind a module-level forwarding wrapper, preserving
  the patch point while avoiding JSON helper startup during command discovery.
  List-file parsing now streams over the list file instead of materializing a
  full `splitlines()` list first, while preserving comment/blank filtering and
  relative path resolution. The current list parser also resolves relative
  paths with a precomputed parent string, `os.path.isabs`, and a small
  Path-compatible string normalizer with a fast path for already-normalized
  rows, avoiding a `Path` object for every file-list row while preserving
  absolute paths and relative path output. A follow-up common-case pass checks
  for path separators once before looking for separator-dependent
  normalization patterns, avoiding extra substring scans for simple file names
  with extensions. Group-size ordering now sorts by `len(files)` with
  `reverse=True` instead of negating each key, preserving stable order for
  equal-sized groups while reducing key-call overhead.
  A later startup pass localizes tree-only helpers to tree extraction, so
  FASTA-mode command discovery avoids tree caching/hash helpers. `_extract_taxa`
  now checks input existence with `os.path.exists()` instead of allocating a
  `Path` object per file, preserving the same missing-file error while reducing
  overhead when grouping many small files.
- `OccupancyFilter._extract_taxa` FASTA mode baseline time also materialized
  full `SeqRecord` objects before counting only taxon IDs. The optimized path
  uses `SimpleFastaParser` for extraction. A later pass uses the shared
  header-only parser, preserving first-token IDs while skipping sequence
  assembly. A startup pass imports the FASTA parser inside FASTA
  extraction/filtering helpers, avoiding Bio.SeqIO.FastaIO startup for
  import-only callers. A follow-up startup pass keeps JSON output behind a
  module-level forwarding wrapper, preserving the patch point while avoiding
  JSON helper startup during command discovery. A later startup pass localizes
  tree-only helpers to tree extraction/filtering paths and removes
  annotation-only `typing` startup, so FASTA-mode command discovery avoids tree
  caching/hash helpers and `typing`. List-file parsing now streams directly
  over the list file instead of materializing a `splitlines()` list, preserving
  comment/blank filtering and relative path resolution. The current parser
  also uses string-based absolute-path checks, a precomputed parent prefix, and
  a small Path-compatible string normalizer with a fast path for
  already-normalized rows, avoiding per-line `Path` allocation while preserving
  resolved output. A follow-up common-case pass checks for path separators once
  before looking for separator-dependent normalization patterns, avoiding extra
  substring scans for simple file names with extensions. `_extract_taxa` now
  mirrors `TaxonGroups` by using `os.path.exists()` for the per-file existence
  guard, preserving missing-file behavior while avoiding `Path` allocation in
  large small-file occupancy runs.
- `OccupancyFilter._extract_taxa` tree mode mirrors the TaxonGroups tree-mode
  optimization by reusing the shared direct terminal-name traversal after
  Newick parsing and falling back to `get_terminals()` for nonstandard trees.
- `OccupancyFilter._filter_fasta` baseline time streamed kept `SeqRecord`
  objects through `SeqIO.write`. The optimized path uses `SimpleFastaParser`
  and writes retained FASTA records directly with 60-character wrapping while
  preserving full header descriptions and retained/removed counts. The wrapping
  helper now emits a whole wrapped sequence with one joined write while
  preserving empty-sequence behavior and line widths. A later pass builds the
  wrapped chunks as a list before joining, avoiding generator overhead in large
  retained FASTA outputs. A later streaming parser pass only accumulates
  sequence lines for retained records, preserving wrapped output while skipping
  sequence joining and cleanup for records that will be dropped.
- `OccupancyFilter._print_text` now batches the occupancy summary, removed-taxa
  table, per-taxon status rows, and output-file rows into one newline-joined
  print while preserving the same stdout text and broken-pipe handling.
- `OccupancyFilter.run` no longer sorts all occupied taxa before threshold
  classification; that list was unused, and downstream reports still perform
  their own deterministic sorting for displayed taxa. A follow-up pass batches
  per-file taxon sets through `Counter.update()`, preserving identical
  occupancy counts while avoiding one Python-level increment per taxon.
- `OccupancyFilter._filter_tree` baseline time deep-copied the parsed tree and
  pruned by tip-name strings. The optimized path prunes terminal clade objects
  directly on the isolated parsed tree before writing the filtered Newick output.
- `OccupancyFilter._filter_tree` batch standard-tree pruning keeps the same
  parsed-tree filtering semantics but collects terminal objects with a direct
  traversal and removes all filtered tips with the shared standard-tree batch
  prune helper, falling back to per-tip pruning for nonstandard tree objects.
  The terminal scan now pushes binary children right-then-left and indexes
  multifurcations backward, preserving terminal order while avoiding
  `reversed(children)` iterator setup and localizing prune-target bookkeeping.
- `RenameFastaEntries.replace_ids_in_file_and_write` baseline time streamed
  renamed `SeqRecord` objects through one `SeqIO.write` call. The optimized
  command path uses `SimpleFastaParser` and writes FASTA records directly with
  60-character wrapping, preserving renamed headers, unrenamed descriptions,
  and total/renamed counts. A later pass combines each streamed record's header
  and wrapped sequence into one write call in the main file-renaming loop while
  preserving falsey mapped identifiers through explicit id-map membership
  checks. The wrapping helper now uses one joined write per sequence instead of
  one write per wrapped line. A later pass uses list-backed chunk construction
  for the same wrapping semantics. The main streaming writer now uses a
  sentinel-safe `idmap.get()` and writes sequences that already fit within the
  60-character FASTA line width directly, avoiding chunk-list construction while
  preserving falsey mapped identifiers. ID-map loading now uses an explicit
  split loop instead of feeding split rows through `dict()`, preserving
  whitespace tokenization and duplicate-key overwrite behavior. A subsequent
  startup pass imports Bio.SeqIO and the FASTA parser inside the writer/read
  helpers, avoiding Biopython parser startup for import-only callers. A
  follow-up startup pass keeps JSON
  output behind a module-level forwarding wrapper, preserving the patch point
  while avoiding JSON helper startup during command discovery. A later startup
  pass removes the runtime `TYPE_CHECKING` dependency and converts
  annotation-only typing aliases to postponed built-in annotations.
- `Faidx.run` baseline text output called `print` once per requested FASTA
  entry. The optimized path formats all requested blocks in order and emits one
  joined print call, preserving the same FASTA text stream while reducing write
  overhead for large multi-entry requests. JSON output now also reuses each
  indexed record lookup while building rows. A later pass replaced
  `SeqIO.index` with a `SimpleFastaParser` scan that stores only requested
  records while preserving request order, first-token FASTA names, missing-entry
  `KeyError`s, and duplicate-ID `ValueError`s. A later pass uses the shared
  selected-entry parser to skip sequence assembly for unrequested records while
  still scanning the full file for later duplicate IDs. A startup pass imports
  the FASTA parser inside `_fetch_entries`, avoiding Bio.SeqIO.FastaIO startup
  for import-only callers. A follow-up startup pass keeps JSON output behind a
  module-level forwarding wrapper, preserving the patch point while avoiding
  JSON helper startup during command discovery. A later startup pass removes the
  annotation-only `typing` import under postponed built-in annotations. Shared
  FASTA parsing now skips `"".join(...)` for single-line sequence records while
  preserving whitespace cleanup and multiline record behavior, benefiting
  `Faidx`, `SumOfPairsScore`, taxon grouping, subsampling, and other direct
  first-token FASTA readers.
- `AlignmentSubsample._read_alignment` baseline time materialized
  `SeqRecord` objects before extracting IDs and sequence strings. The optimized
  path uses `SimpleFastaParser`, preserving first-token IDs and last duplicate
  record wins behavior. A later pass reuses the shared case-preserving
  first-token parser, preserving multiline sequence joining and internal
  whitespace removal without uppercasing sampled sequences. A startup pass
  imports the FASTA parser inside the read helper, avoiding Bio.SeqIO.FastaIO
  startup for import-only callers. A follow-up startup pass keeps JSON output
  behind a module-level forwarding wrapper, preserving the patch point while
  avoiding JSON helper startup during command discovery. A later startup pass
  localizes the `random` import to command execution and the partition/site
  helper imports to their specific modes, so import-only command discovery does
  not load the sampler. A follow-up parser pass handles common RAxML-style
  partition rows with split operations while preserving comment skips,
  invalid-row skips, whitespace tolerance, and trailing text tolerance. A later
  startup pass converts annotation-only `typing` aliases to postponed built-in
  annotations.
- `AlignmentSubsample._run_sites` baseline time built every sampled sequence
  with a Python generator indexing selected columns one at a time. The
  optimized path builds one `itemgetter` for the selected index list and applies
  it to each taxon sequence, preserving single-site and bootstrap duplicate
  index behavior. A full-site non-bootstrap follow-up now reuses the parsed
  sequence mapping directly when the selected count equals the alignment length,
  avoiding random sampling, `itemgetter` construction, and per-taxon sequence
  rebuilds while leaving bootstrap and partial-site sampling on the existing
  path.
- `AlignmentSubsample._run_partitions` baseline time walked every selected
  partition and then every taxon, repeatedly looking up taxon sequences while
  collecting partition slices. The optimized path builds selected slice ranges
  once, then assembles each taxon's subsampled sequence with the source
  sequence held locally, preserving partition ordering, duplicate partition
  names, and coordinates.
- `AlignmentSubsample._print_summary` baseline time emitted the fixed summary
  header and each output path through separate `print()` calls. The optimized
  path builds the same lines once and emits one joined report, preserving
  captured stdout exactly.
- `AlignmentSubsample` FASTA and partition output baselines wrote each record
  or partition row separately. The optimized paths batch bounded groups of
  rows into joined writes, preserving insertion order and exact output text
  while avoiding one file write per taxon or partition on large subsamples.
- `AlignmentRecoding.recode_alignment` baseline time performed a Python
  dictionary lookup and gap check for every character. The optimized path builds
  a `str.translate` table once per recoding operation, including lowercase
  equivalents for recoded symbols while preserving gap and ambiguous symbols.
  `run` now uses a string-backed internal recoding path for text and JSON
  output, avoiding per-sequence list materialization while keeping
  `recode_alignment` list-compatible for callers.
  Text-mode `run` now batches recoded FASTA records into one newline-joined
  print while preserving record order, stdout text, and JSON output. A later
  startup pass postpones the annotation-only Bio.Align import and wraps
  `print_json` lazily so command discovery no longer initializes Biopython
  alignment internals or JSON helpers. A follow-up startup pass removes the
  runtime `TYPE_CHECKING` dependency and converts annotation-only typing aliases
  to postponed built-in annotations. `read_recoding_table` now bounds
  whitespace splitting to the first three fields, preserving first-two-token
  semantics and ignored trailing columns while reducing token allocation for
  large custom recoding tables. `_build_translation_table` now inserts
  uppercase and lowercase translation entries directly instead of allocating a
  two-item set for every recoding symbol, while preserving skipped uppercase and
  lowercase gap symbols.
- `PhyloGwas._read_fasta` baseline time materialized `SeqRecord` objects before
  extracting IDs and uppercase sequence strings. The optimized path uses
  `SimpleFastaParser`, preserving first-token IDs, uppercase conversion, and
  last duplicate record wins behavior. A later pass reuses the shared streaming
  first-token parser, preserving multiline sequence joining and internal
  whitespace removal without constructing parser tuples. A startup pass defers
  both the FASTA parser import and optional Bio.Phylo tree parser import until
  FASTA loading or phylogenetic classification actually needs them. A subsequent
  startup pass defers NumPy, JSON helper, and plot config imports behind a
  module-level proxy/wrapper or argument-processing import while preserving
  existing test patch points. A later shared-startup pass removes
  annotation-only `typing` imports from the error, file, and alignment-base
  helpers loaded by `phylo_gwas` command discovery. Phenotype TSV parsing now
  bounds row splitting to the first two tab-separated columns, preserving
  comment/blank skipping, malformed-row skipping, ignored trailing columns, and
  last duplicate taxon wins behavior while avoiding allocation of unused extra
  fields. A follow-up parser pass reads rows in binary mode and decodes only
  the retained first two fields, preserving the same parsed phenotype mapping
  while avoiding text decoding and string allocation for ignored trailing fields.
- `PhyloGwas` allele extraction baseline time performed one dictionary lookup
  per taxon per tested alignment column. The optimized path builds one ASCII
  byte matrix for shared taxa, extracts columns through `tobytes().decode()`,
  and falls back to the original string path for non-ASCII alignments. A later
  site-scan pass precomputes non-ambiguous ASCII columns once with the same byte
  lookup and only iterates valid columns; non-ASCII alignments keep the original
  per-column ambiguity check. A subsequent pass samples unambiguous ASCII
  columns and, when invariant or multiallelic columns are present, builds a
  stricter biallelic mask so those columns never enter the per-site association
  tests; all-biallelic scans keep the lighter non-ambiguous mask.
- `PhyloGwas._benjamini_hochberg` baseline time walked sorted p-values in a
  Python reverse loop to enforce monotonic adjusted values. The optimized path
  computes ranked p-values in NumPy and applies a vectorized reverse cumulative
  minimum, preserving scalar BH results including ties.
- `PhyloGwas` partition annotation baseline time scanned every partition for
  every tested site. The optimized path builds a binary-search lookup for
  sorted, non-overlapping partition files and keeps the original linear lookup
  as a fallback for overlapping, unsorted, or invalid ranges.
- `PhyloGwas._parse_partition_file` now lazily compiles and reuses the
  RAxML-style row matcher within the parser, preserving whitespace tolerance,
  trailing text tolerance, skipped malformed rows, and 0-based range conversion
  while avoiding a regex lookup/compile-cache path for every partition row.
- `PhyloGwas._test_site_categorical` baseline time built a contingency table
  and then rescanned each site once per phenotype group to compute derived
  allele frequencies. The optimized path computes frequencies directly from
  a vectorized `bincount` table, computes the multi-group chi-square p-value
  directly, and lets `run` reuse per-taxon group codes, group row counts, and
  the group-index map. ASCII biallelic sites now build one byte vector for both
  allele counting and derived-state encoding, retaining the Counter fallback for
  non-ASCII alleles. A later categorical pass lets `run` hand ASCII byte-matrix
  columns directly to the categorical tester, avoiding column decode, list
  construction, and re-encoding before the same byte-counting logic. A
  subsequent byte path uses a fixed 256-bin histogram for allele counts instead
  of sorting unique byte values. A minor follow-up reuses precomputed
  `group_codes * 256` offsets for multi-group byte columns and passes the
  chi-square survival function from `run`, avoiding one mask allocation and
  repeated import dispatch per categorical site. A later byte path identifies
  major/minor alleles with min/max counts instead of a 256-bin histogram while
  still rejecting invariant and multiallelic byte columns. The two-group path
  now computes the same two-sided Fisher exact p-value locally from the 2x2
  table, using a log-space hypergeometric recurrence to avoid repeated
  `scipy.stats.fisher_exact` calls and repeated log-combination work at GWAS
  site scale. The categorical run path now also caches those local Fisher
  p-values by 2x2 count table for the duration of a GWAS run, preserving direct
  helper behavior while avoiding repeated hypergeometric summations for common
  count patterns. A later pass reuses the same Fisher cache for repeated
  log-combination terms across distinct 2x2 tables, preserving the exact
  two-sided summation while reducing repeated `lgamma` calls. A later prepared
  two-group byte path counts the minor allele in group 1 directly and derives
  group 0 from the known minor count, avoiding the per-site 512-bin grouped
  histogram and generic frequency dictionary used by multi-group tests. Byte
  major/minor helpers now use ndarray `min()`/`max()` directly, avoiding lazy
  NumPy proxy dispatch in both categorical and continuous ASCII paths. The
  Benjamini-Hochberg correction now scales and cumulative-mins the sorted
  p-value buffer in place, then scatters into an uninitialized result array,
  preserving adjusted p-values while reducing temporary allocation pressure.
  Obvious non-ASCII multiallelic sites now skip after seeing a third unique
  allele instead of building a full `Counter` and sorting every unique allele;
  true binary Unicode sites retain the previous sorted-label tie behavior and
  major/minor frequency ordering.
- `PhyloGwas._test_site_continuous` baseline time called SciPy's
  `pointbiserialr` and rebuilt phenotype arrays for every biallelic site. The
  optimized path computes the equivalent Pearson point-biserial statistic and
  Student-t p-value directly, while `run` reuses centered phenotype values and
  builds each ASCII allele indicator vector from a byte buffer with a non-ASCII
  fallback. A later pass routes ASCII alignment columns directly into the
  continuous test as `uint8` arrays, avoiding per-site list decoding while
  preserving allele labels and the non-ASCII fallback path. A subsequent byte
  path computes the same binary correlation from allele counts and the summed
  centered phenotype values for derived alleles, avoiding per-site float binary
  vector allocation. The same min/max byte allele detector removes the fixed
  histogram setup from ordinary biallelic byte columns and now reduces through
  the byte array's own `min()`/`max()` methods. A later scalar cleanup
  uses `math` for per-site square-root and NaN checks and delays importing
  `stdtr` until a non-skipped site actually needs a Student-t p-value. A
  follow-up byte path computes the centered phenotype dot product directly
  against the derived-allele boolean mask, avoiding a masked phenotype-value copy
  for every continuous byte-column site.
- `PhyloGwas._prune_tree_to_taxa` baseline time rescanned all current
  terminals once for every removed tip before phylogenetic pattern
  classification. The optimized path materializes terminals once and prunes
  non-shared terminal clades directly. A later no-op shared-taxa setup pass
  collects terminals through direct standard-tree traversal instead of
  `get_terminals()`, retaining the generic fallback for nonstandard tree
  objects. The direct terminal helper now preserves tip order while pushing
  binary children right-then-left and indexing multifurcations backward,
  avoiding `reversed(children)` iterator setup.
- `PhyloGwas._classify_phylo_pattern` baseline time recomputed an MRCA and
  descendant tip list for every significant-site classification. The optimized
  run path builds all monophyletic descendant taxon sets once in postorder and
  classifies each derived-taxon set by hash lookup, retaining the original MRCA
  fallback for direct helper calls without a precomputed index. A follow-up
  pass builds that monophyletic-set index with a direct stack-based postorder
  traversal for standard trees, retaining the generic `find_clades` fallback
  for nonstandard tree objects. A later traversal pass builds that postorder
  by collecting preorder once and iterating it in reverse, preserving the
  monophyletic-set index while avoiding visited-flag stack tuples. Balanced
  binary nodes now combine their two child tipsets directly, while unary nodes
  and multifurcations retain their existing paths. The run path now also builds
  each significant site's minor-allele taxon list in one direct shared-taxon
  scan, avoiding the previous temporary per-site allele list before
  phylogenetic-pattern classification.
- `PhyloGwas._create_manhattan_plot` partition shading baseline time added one
  `axvspan` artist per shaded partition. The optimized path batches alternating
  gray spans into one `PatchCollection` with the same x-data/y-axes transform,
  preserving partition labels while avoiding one patch artist per shaded
  partition. Point series preparation now fills positions, log-transformed
  p-values, point colors, and the FDR threshold in one pass over result rows,
  avoiding separate list comprehensions for each series while preserving the
  same scatter inputs.
- `PhyloGwas._print_text_output` baseline time fully sorted all significant
  sites before displaying only the top ten and emitted report lines one at a
  time. The optimized text path selects the same top ten with partial selection
  for shuffled site-order results and emits the joined report in one `print()`,
  preserving captured stdout exactly. A follow-up pass reuses a cached
  `itemgetter("p_value")` key for the partial top-hit selection, avoiding one
  Python lambda call per significant site while preserving the same selected
  rows and report text. JSON output now passes the existing result-row list into
  the payload instead of shallow-copying every result dictionary before
  serialization, preserving identical JSON text while avoiding per-row dict
  allocation. Significant-site summary accounting now builds the significant-row
  list and counts mono/polyphyletic rows in one pass, preserving the same row
  objects and summary counts. CSV output now uses direct required-field access
  for result dictionaries produced by the GWAS run path while retaining the
  previous blank-default `get()` fallback for incomplete external result rows.
- `Alignment.calculate_rcv` baseline time built a Unicode character matrix and
  filled the per-sequence count matrix through nested sequence/symbol loops.
  The optimized path uses a byte-backed alignment matrix and computes each
  valid symbol's counts across all taxa in one vector operation before applying
  the unchanged aggregate RCV formula. A later pass replaced the ASCII
  `np.isin` validity mask with a byte lookup table while preserving the Unicode
  fallback path. Protein-scale ASCII alphabets now use per-row byte bincounts
  for the count matrix, avoiding one full alignment scan per amino-acid symbol
  while leaving the faster small-alphabet path in place. The observed-symbol
  validity pass derives the valid ASCII symbol set from byte codes present in
  the alignment and uses full sequence lengths when no invalid symbols are
  observed, avoiding full validity-mask construction for no-gap inputs and
  reducing filtered-symbol setup for gap-bearing alignments. Gap-bearing paths
  now count row valid lengths with `np.count_nonzero` instead of summing boolean
  masks before float conversion. Identical alignments now return zero RCV
  before NumPy matrix construction, preserving case-insensitive behavior by
  comparing uppercased sequence strings. A
  follow-up pass scans those sequence strings by index instead of evaluating the
  shortcut over `sequences[1:]`, avoiding a large temporary list while
  preserving early exit for heterogeneous alignments. The final per-taxon RCV
  total now uses the ndarray reduction method directly, avoiding the generic
  `np.sum` dispatch on the realistic taxon-vector sizes covered by the RCV
  benchmarks. Clean large-short ASCII matrices now use a single encoded
  `bincount` count matrix pass, matching the existing row-wise count output
  while avoiding one Python loop iteration per taxon. Single-record alignments
  now return zero RCV before sequence materialization, preserving the existing
  zero result while avoiding unnecessary uppercase copies. The `rcv` command module
  now keeps JSON output behind a
  module-level forwarding wrapper,
  preserving the patch point while avoiding JSON helper startup during command
  discovery.
- `OccupancyPerTaxon.calculate_occupancy_per_taxon` baseline time allocated a
  NumPy character array for each sequence before checking invalid symbols. The
  optimized path uppercases each sequence once and uses C-level `str.count`
  calls over the small invalid-symbol set. The current ASCII path uses a byte
  validity lookup table for larger alignments; a later pass made that table
  case-insensitive so ASCII rows no longer allocate uppercase copies, with a
  string-count fallback for non-ASCII sequences. A variable-length ASCII
  fallback pass reuses the sequence strings gathered during the matrix attempt
  and counts valid bytes with `np.count_nonzero` instead of reducing boolean
  arrays with `np.sum`; a later pass replaced that per-record NumPy setup with
  `bytes.translate(..., invalid_bytes)` for variable-length ASCII rows while
  preserving the matrix path and Unicode fallback. A later equal-length ASCII
  pass scans the joined byte buffer for invalid symbols first and returns
  occupancy `1.0` for all no-gap taxa, avoiding full matrix lookup construction
  for all-valid DNA and protein alignments. Identical equal-length ASCII records
  now count valid symbols once and repeat the occupancy for each taxon, avoiding
  byte-matrix construction for conserved gappy alignments while preserving DNA
  and protein invalid-symbol rules. A follow-up identical-row pass scans the
  existing sequence list directly instead of materializing `sequences[1:]`,
  reducing temporary allocation for high-taxon conserved alignments. A later
  helper pass combines equal-length validation and identical-row
  detection before building a separate sequence list, speeding conserved
  alignments and late variable-length rejection while preserving matrix
  behavior for non-identical equal-length inputs. Single-record alignments now
  count occupancy directly with the same ASCII/Unicode helper used by the
  fallback path, avoiding record-data and matrix-helper setup. Text-mode
  `run` now batches per-taxon rows into one newline-joined print, preserving
  the same stdout text and leaving JSON output unchanged. A subsequent startup
  pass defers NumPy, validity lookup
  construction, and the JSON output helper until occupancy calculation or JSON
  output actually needs them. A follow-up startup pass removes the annotation-only
  `typing` import under postponed annotations.
- `Dstatistic._run_alignment_mode` baseline time scanned every site in Python
  once for total ABBA/BABA counts and again for block jackknife counts. The
  optimized path builds byte-backed arrays for the four focal sequences,
  computes ABBA/BABA boolean masks once, and derives total and per-block counts
  from those masks, with a scalar fallback for non-ASCII sequence data. A later
  clean-ASCII pass samples and scans for skip codes, bypassing the validity mask
  entirely for all-valid quartets while preserving the mask path for ambiguous
  or gapped data. Quartet inputs where P1/P2 are identical or P3/outgroup are
  identical now return zero totals and zero-valued block arrays before NumPy
  byte-array setup because ABBA/BABA patterns are impossible in those cases.
  A later partial-pattern pass detects aligned quartets where only ABBA or only
  BABA can occur and skips the impossible boolean mask and block aggregation,
  using guarded sequence identity checks to keep mixed-pattern inputs on the
  two-mask path.
- `Dstatistic._jackknife_d_values` baseline time computed every leave-one-block
  D value in a Python loop. The optimized helper computes all leave-one-out
  ABBA/BABA totals and D values with vectorized NumPy arithmetic, preserving
  zero-denominator blocks as zero. The helper now totals each block vector with
  the ndarray `sum()` method, avoiding two lazy NumPy reduction dispatches.
  Alignment-mode jackknife setup now computes the jackknife mean through the
  ndarray `mean()` method, avoiding the generic `np.mean` wrapper.
  Alignment-mode jackknife standard-error setup now computes the centered
  sum-of-squares with a dot product, avoiding the temporary squared deviation
  reduction while preserving the same standard-error formula.
- `Dstatistic._read_fasta` and `Dfoil._read_fasta` baseline time materialized
  `SeqRecord` objects while alignment mode only needed IDs and uppercase
  sequences. The optimized path uses `SimpleFastaParser`, preserving first-token
  IDs, uppercase conversion, and last duplicate record wins behavior. A later
  pass replaces the two command-local `SimpleFastaParser` readers with a shared
  streaming first-token parser that preserves multiline sequence joining,
  internal space/carriage-return removal, uppercase conversion, and duplicate
  replacement without constructing parser tuples. A startup pass imports
  SimpleFastaParser inside the read helpers and imports Bio.Phylo inside
  D-statistic gene-tree parsing, avoiding Biopython parser startup for
  import-only callers. Follow-up startup passes also defer NumPy and JSON output
  imports behind module-level wrappers while preserving the existing vectorized
  site-pattern paths and JSON output behavior. Gene-tree mode now builds the
  direct preorder list once and reuses its reverse for the descendant-taxon pass,
  preserving nonterminal preorder while avoiding visited-flag stack tuples. A
  later DFOIL startup pass removes annotation-only `typing` aliases under
  postponed built-in annotations. A follow-up D-statistic startup pass applies
  the same conversion, so command discovery no longer loads
  `typing`. The scalar fallback for non-ASCII alignments now uses direct skip
  checks and computes block membership only for counted sites, preserving
  ABBA/BABA totals and jackknife block arrays while avoiding per-site list and
  generator allocation.
- `Dfoil._count_site_patterns` baseline time scanned every five-taxon site in
  Python, built temporary allele sets, and assembled string pattern keys. The
  optimized path builds byte-backed arrays, filters skipped and non-biallelic
  sites with boolean masks, and counts pattern codes with `np.bincount`, while
  keeping a scalar fallback for non-ASCII sequence data. A later clean-ASCII
  pass samples and scans for skip codes, bypassing the validity mask entirely
  for all-valid sequence sets while preserving the mask path for ambiguous or
  gapped data. The scalar fallback now mirrors the pattern-code idea directly:
  it uses direct skip checks, verifies all derived symbols match the first
  derived state, and indexes the 16-pattern table from bit flags instead of
  allocating temporary allele sets and pattern strings. Fully invariant
  five-taxon inputs now return the same all-zero pattern dictionary before
  NumPy array setup, matching the existing skipped-invariant scalar semantics.
  Short ASCII alignments with skip codes now build the validity mask through a
  cached byte lookup table, while longer alignments retain the previous
  six-code vector loop that remains faster at larger sizes.
- `Dstatistic` alignment-mode and gene-tree-mode text report baselines emitted
  each report line through separate `print()` calls. The optimized helpers build
  the same lines and emit each report with one `print()`, preserving captured
  stdout exactly. A follow-up pass formats both fixed reports directly as one
  string, avoiding the temporary line lists while preserving the same text for
  significant and non-informative branches.
- `Dfoil._print_text_output` baseline time emitted each report section through
  separate `print()` calls. The optimized path builds the same report lines and
  emits them with one `print()`, preserving captured stdout exactly. A
  follow-up pass caches the fixed informative-pattern groups and formats the
  same report directly, avoiding per-report pattern-list and line-list
  construction.
- `Dstatistic._get_quartet_topology` baseline time called
  `clade.get_terminals()` for every internal clade while classifying each gene
  tree. The optimized path computes descendant taxon sets once in postorder and
  then evaluates bipartitions in the same `get_nonterminals()` order as before,
  preserving support-threshold behavior. A follow-up pass avoids building the
  full `all_taxa - tips` complement for every internal bipartition; after the
  four requested quartet taxa are known to be present, topology classification
  only needs the two quartet taxa found on the current side of the split.
- `TipToTipDistance` and `PatristicDistances` baseline time was dominated by
  repeated Biopython `tree.distance` calls, each of which searches paths and
  common ancestors. The optimized path caches root-to-tip depths and ancestor
  paths once per tree, then computes pair distances from those cached values.
  All-pairs tip-name setup now also uses the shared direct terminal-name
  traversal for parsed trees before entering the cached pairwise distance core.
  TipToTipDistance all-pairs text output now batches pairwise rows into one
  newline-joined print while preserving plot status output and stdout text.
  Large all-pairs heatmap matrix construction now detects rows already in sorted
  upper-triangle order and fills the symmetric matrix from a dense distance
  vector, retaining the taxon-index dictionary fill for smaller matrices or
  arbitrary row orders.
  PatristicDistances verbose text output now batches pairwise rows into one
  newline-joined print while preserving empty-output behavior and stdout text.
  A later startup pass defers PatristicDistances' annotation-only Bio.Phylo
  import and optional tqdm progress import until the progress path actually
  needs it. A follow-up startup pass keeps module-level `mp` and `pickle`
  patch points as lazy proxies and imports `functools.partial` only in the
  large-tree multiprocessing branch. The latest startup pass preserves the
  module-level summary-statistics and JSON helper patch points as lazy
  forwarding functions, so import-only callers no longer initialize those
  helper modules. A later startup pass converts annotation-only `typing` aliases
  to built-in postponed annotations so command discovery no longer loads
  `typing`. Non-verbose PatristicDistances text and JSON runs now compute
  distance summaries with a stats-only path that skips pair-label tuple
  allocation while preserving verbose pair output. Deep-tree stats-only
  PatristicDistances now switches from materialized root paths for every tip to
  an indexed binary-lifting LCA path, reducing pectinate-tree memory pressure
  while keeping the shallow balanced-tree path unchanged. The stats-only
  PatristicDistances setup now mirrors the shared all-pairs helper's
  order-preserving binary-child push instead of extending a `reversed(children)`
  iterator for every internal node. The shared all-pairs distance helper now
  uses the same deep-tree LCA-index fallback for verbose PatristicDistances,
  `tip_to_tip_distance --all-pairs`, and other callers that need combo labels
  as well as distances. Cached read-only PatristicDistances run setup now uses
  the explicit unmodified tree read helper, avoiding a defensive copy of the
  cached parsed tree before dispatching to the pairwise distance/statistics
  routines.
- `SpuriousSequence.get_branch_lengths_and_their_names` baseline time
  materialized every terminal clade through Bio.Phylo before collecting terminal
  branch lengths. The optimized path uses a direct stack traversal for standard
  parsed trees and retains the original terminal-list fallback for nonstandard
  tree objects. A later startup pass defers its annotation-only Bio.Phylo
  import until tree I/O is actually requested by the shared tree base. The
  direct terminal helper now pushes binary children right-then-left and indexes
  multifurcations backward, preserving terminal order while reducing balanced
  131072-tip terminal traversal median time from 0.045395s to 0.040772s. The
  branch-length collection path now records terminal branch lengths and names
  during that standard-tree traversal, avoiding a separate terminal-list
  materialization pass while preserving the fallback route for nonstandard tree
  objects.
- `SpuriousSequence.run` text-output baseline rounded threshold/median for
  every flagged row, built JSON payload rows even for text output, and called
  `print` once per reported row. The optimized text path rounds shared values
  once, builds only the text rows it needs, and emits the same stdout text in a
  single output call while preserving JSON payload shape. Cached read-only
  `SpuriousSequence.run` now also uses the explicit unmodified tree read helper
  to avoid copying the cached parsed tree before collecting terminal branches.
  A follow-up startup pass keeps JSON output behind a module-level forwarding
  wrapper, preserving the patch point while avoiding JSON helper startup during
  command discovery. A later startup pass replaces the module-level
  `statistics.median` import with a local numeric median helper, preserving
  odd/even median behavior while avoiding `statistics` startup during command
  discovery. A follow-up startup pass removes the annotation-only `typing`
  import by using postponed built-in collection annotations. Large branch-length
  medians now use lazy NumPy partition after collection, avoiding a full Python
  sort while preserving the lightweight sorted path for small inputs and eager
  import behavior for command discovery.
- `DensityMap._terminal_clades` applies the same direct terminal traversal to
  plot y-order setup and the final tip count, avoiding Bio.Phylo terminal-list
  materialization while preserving terminal order. The direct preorder and
  terminal helpers now push binary children right-then-left and index
  multifurcations backward, avoiding `reversed(children)` iterator setup while
  preserving Biopython traversal order.
- `DensityMap._plot_density_map` layout setup baseline time rebuilt
  rectangular coordinates, circular parent maps, and branch drawing traversal
  with separate preorder/postorder `find_clades()` scans. The optimized path
  uses the shared direct `compute_node_positions` helper, reuses the existing
  id-to-parent map for circular coordinates, and iterates one direct preorder
  list for rectangular and circular branch drawing. A later setup pass passes
  that preorder list into `compute_node_positions` and passes the same preorder
  list plus the prepared tip list into `compute_circular_coords`, avoiding
  additional direct tree walks while preserving circular tip order. A later branch-rendering
  pass batches rectangular posterior segments and vertical connectors into
  `LineCollection`s, and batches circular radial posterior segments plus gray
  internal arcs into collections while vectorizing state-color blending per
  branch. Branch posterior aggregation now preselects the branch histories once
  and advances through each history across sorted segment midpoints, preserving
  the previous missing-history denominator while avoiding a full history scan
  per segment. A later startup pass
  defers NumPy behind a lazy proxy, leaving array startup until density-map
  simulation setup or posterior color blending first needs numeric arrays.
  A subsequent startup pass localizes JSON output, `PlotConfig`, shared layout
  helpers, circular-layout drawing helpers, and color annotations to output or
  plot paths, avoiding those helper imports for import-only callers. A later
  `run` setup pass reads the cached tree without copying, counts missing tip
  states directly, and only creates a working copy when pruning or ladderizing
  would mutate the tree. Text summary output now batches the five-line report
  into one newline-joined print while preserving exact stdout text. A follow-up
  startup pass converts annotation-only `typing` aliases to built-in postponed
  annotations, so command discovery no longer imports `typing`.
- `RobinsonFouldsDistance.run` now finds the first terminal for outgroup rooting
  with direct preorder traversal instead of materializing every terminal clade.
  Nonstandard tree objects retain the previous terminal-list fallback. A later
  helper pass descends the leftmost child chain directly instead of maintaining
  a traversal stack, preserving the same rooting tip while reducing depth-17
  balanced-tree lookup time from 2.650us to 1.039us. A later
  startup pass postpones annotations and removes the annotation-only
  `Bio.Phylo.Newick` import, leaving tree parser startup to the shared tree
  reader instead of module import. A follow-up startup pass keeps JSON output
  behind a module-level forwarding wrapper, preserving the patch point while
  avoiding JSON helper startup during command discovery.
- `KuhnerFelsensteinDistance.run` uses the shared tree helper for the same
  rooting-tip setup, avoiding full terminal-list materialization while retaining
  the original fallback for nonstandard tree objects. A later startup pass
  defers its annotation-only Bio.Phylo import until tree I/O is actually
  requested by the shared tree base. A follow-up startup pass keeps JSON output
  behind a module-level forwarding wrapper, preserving the patch point while
  avoiding JSON helper startup during command discovery.
- `TipToTipDistance.calculate_tip_to_tip_distance` applies the same idea to the
  single-pair command path: it finds both tips while recording parent links and
  root depths in one preorder traversal, then computes the distance from the
  two root depths and LCA depth. Nonstandard tree objects retain the previous
  `check_leaves()`/`TreeMixin.distance()` fallback. A later scan pass uses each
  clade's child list directly instead of calling `is_terminal()` for every
  visited clade, preserving the same standard-tree traversal and fallback
  behavior while reducing large-tree single-pair setup overhead. A later startup pass
  deferred SciPy hierarchical-clustering imports until heatmap plotting, so
  normal non-plot single-pair and all-pairs runs no longer import the clustering
  stack. Cached read-only `TipToTipDistance.run` now uses the explicit
  unmodified tree read helper to avoid copying the cached parsed tree before
  dispatching to single-pair or all-pairs distance calculation. A subsequent
  startup pass defers NumPy and Bio.Phylo TreeMixin imports behind module-level
  proxies, preserving `TreeMixin.find_any` and `TreeMixin.distance` patch
  points while avoiding those imports for import-only callers. A follow-up
  startup pass keeps JSON output behind a module-level forwarding wrapper and
  localizes `PlotConfig` to argument processing, preserving patch points while
  avoiding JSON and plot-config helper startup during command discovery. A later
  startup pass converts annotation-only `typing` aliases to built-in postponed
  annotations so command discovery no longer loads `typing`. Same-tip requests
  now use a validation-only terminal scan that returns zero once the shared
  terminal is found, preserving missing-tip validation while skipping branch
  depth maps and ancestor-set construction.
- `TipToTipNodeDistance.calculate_tip_to_tip_node_distance` baseline time
  searched for both tips and then used Biopython `trace()`, which performs
  additional tree searches. The optimized path finds both tips and records
  parent links in one preorder traversal, reconstructs the two root paths, and
  falls back to the previous `check_leaves()`/`trace()` behavior for nonstandard
  tree objects. A follow-up traversal pass checks `clades` directly while
  searching for terminals, avoiding one `is_terminal()` method call per visited
  clade on standard parsed trees. A later startup pass postpones annotations
  and replaces the eager Bio.Phylo TreeMixin import with a module-level lazy
  proxy, preserving the existing `TreeMixin.find_any` and `TreeMixin.trace`
  patch points while keeping Bio.Phylo out of import-only callers. Cached read-only
  `TipToTipNodeDistance.run` now uses the explicit unmodified tree read helper
  to avoid copying the cached parsed tree before node-distance calculation. A
  follow-up startup pass keeps JSON output behind a module-level forwarding
  wrapper, preserving the patch point while avoiding JSON helper startup during
  command discovery. A later startup pass removes the annotation-only `typing`
  import so command discovery no longer loads `typing`. A subsequent path
  calculation pass avoids building and reversing two root-path lists by recording
  distances from one tip to each ancestor and climbing from the other tip until
  the LCA is found; the same pass preserves preorder traversal without allocating
  `reversed(children)` for binary parsed trees. Same-tip requests now use a
  validation-only terminal scan and return zero once the shared terminal is
  found, preserving missing-tip validation while skipping parent-map and
  ancestor-distance work.
- `RobinsonFouldsDistance.calculate_robinson_foulds_distance` baseline time
  compared each non-root internal clade by calling `common_ancestor()` on the
  other tree and then collecting descendant tips again. The optimized path
  extracts rooted descendant tip sets for both trees in one postorder traversal
  per tree and compares the resulting sets, preserving the historical rooted RF
  semantics and normalization. A later pass avoids the remaining
  `count_terminals()` traversal during normalization by counting terminal clades
  directly for standard Bio.Phylo trees, with the original count method retained
  as fallback. Split extraction now also uses an iterative direct postorder walk
  for standard Bio.Phylo trees, retaining the previous `find_clades()` path for
  nonstandard tree-like objects.
- `KuhnerFelsensteinDistance.calculate_kf_distance` baseline time collected
  descendant terminal names with `clade.get_terminals()` for every non-root
  branch while building branch-score split maps. The optimized path reuses
  child descendant tip sets in one postorder traversal per tree while preserving
  terminal branches, internal branches, and zero treatment for missing branch
  lengths. A follow-up traversal pass replaces the remaining generic Bio.Phylo
  postorder iterator with a direct clade list for standard trees, retaining the
  generic fallback for nonstandard objects. A later split-map pass special-cases
  binary child unions with direct frozenset `|` operations and keeps explicit
  single-child and multifurcating paths, preserving split-key overwrite behavior
  while reducing balanced 32768-tip split-map time from 0.049567s to 0.035310s.
  A later startup pass converts annotation-only `typing` aliases to built-in
  postponed annotations, so command discovery no longer loads `typing`.
  Same-object method calls now return `(0.0, 0.0)` before split extraction;
  equal-but-separate tree objects still take the normal branch-score path. A
  later accumulation pass avoids allocating the union of split keys by consuming
  a copied second split map while accumulating branch-score deltas and total
  branch length together, preserving common and one-sided split behavior.
- `TreeSpace._build_distance_matrix` baseline time deep-copied every tree even
  when all taxa were already shared, then filled RF/KF pairwise distances in
  nested Python loops. The optimized path only copies trees when pruning is
  needed and constructs RF/KF distance matrices from split-presence or
  branch-length matrices with NumPy operations. TreeSpace auto-K detection now
  builds the normalized Laplacian by row/column scaling the affinity matrix
  directly instead of materializing a dense diagonal scaling matrix. Its
  degree-vector setup now uses the matrix's own row-sum reduction, preserving
  the same normalized Laplacian while avoiding generic NumPy reduction dispatch.
  `einsum` was not used for the row sums because the largest 5000-row affinity
  benchmark regressed. A startup pass defers direct NumPy imports and Bio.Phylo
  loading behind lazy proxies, so command discovery avoids array and
  tree-format parser startup until tree parsing, distance-matrix construction,
  clustering, or plotting runs. A
  follow-up startup pass keeps JSON output behind a forwarding wrapper and
  localizes `PlotConfig` to argument processing, avoiding those helpers during
  import-only command discovery. A later startup pass defers pickle until a tree
  copy is required for pruning. Text output now batches cluster summary rows
  into one newline-joined print while preserving exact stdout text. Tree-source
  cleanup now streams over the input file and strips/filters comments in one
  pass instead of materializing `splitlines()`, preserving multi-Newick and
  path-list behavior. Path-list resolution now uses string parent-prefix,
  absolute-path, and existence checks instead of allocating a `Path` object for
  every listed tree, preserving normal relative/absolute path behavior while
  reducing overhead before tree parsing. Shared-taxa setup now walks gene-tree
  indices directly instead of materializing `gene_trees[1:]`, preserving
  intersection behavior while avoiding one large list copy for large tree sets.
- `SpectralDiscordance._build_bipartition_matrix` baseline time deep-copied
  every tree even when all taxa were already shared, then checked every known
  bipartition against every gene tree while filling the matrix. The optimized
  path only copies trees when pruning is needed and populates the matrix from
  each tree's extracted splits through a bipartition-to-column lookup.
  Shared-taxa setup now uses the shared terminal-name traversal for every gene
  tree and the optional species tree, and walks gene-tree indices directly
  instead of materializing `gene_trees[1:]` during the intersection loop. Split extraction now also uses a direct
  postorder clade list for standard trees in both NRF and WRF paths, retaining
  the Bio.Phylo traversal fallback for nonstandard objects. Equal-size
  bipartition canonicalization now compares the smallest taxon on each disjoint
  half instead of sorting both halves, preserving the documented sorted
  lexicographic tiebreak. Spectral clustering uses the same direct
  normalized-Laplacian scaling as TreeSpace, including ndarray row sums for the
  degree vector. Eigenvector row normalization now computes row L2 norms with
  an `einsum` reduction before the square root, avoiding `np.linalg.norm`
  dispatch while preserving normalized rows. K-means++
  initialization now computes each new center's squared distances with one
  `einsum` over a single difference buffer, preserving the fixed-seed center
  choices while avoiding an extra squared temporary. PCA variance explained now
  computes singular-value total variance with a dot product, avoiding the
  temporary squared array reduction after SVD. Gene-tree file-list rows now use
  a bound string existence check before `Phylo.read`, preserving the parser's
  current path interpretation while avoiding a `Path` object per listed file.
  Top-loading reporting now uses guarded partial selection for the common
  no-tie case and falls back to the previous full sort whenever tied loadings
  could affect exact reported membership or ordering.
  The command run path now reuses
  the cached parsed species tree for read-only setup; species-tree pruning in
  bipartition matrix construction remains copy-isolated by `_copy_prune_if_needed`.
  normalized-Laplacian row/column scaling as TreeSpace, preserving eigengaps
  without allocating `diag(d^-1/2)`. A later startup pass defers NumPy behind a
  lazy proxy, postpones annotations, and imports Bio.Phylo only inside
  gene-tree parsing, so command discovery avoids numerical and tree parser
  startup. A follow-up startup pass keeps JSON output behind a forwarding
  wrapper and localizes `PlotConfig` to argument processing, avoiding those
  helper modules during import-only command discovery. Another startup pass
  removes the annotation-only `typing` import by using postponed built-in
  collection annotations. JSON score payload construction now precomputes PC
  labels and iterates score rows directly, preserving nested `gene_tree_N`
  payloads and cluster assignments while avoiding repeated two-dimensional
  indexing.
- `QuartetNetwork._compute_quartet_cfs` baseline time represented quartet and
  split membership as Python sets and repeatedly intersected them while
  scanning every tree's bipartitions for every quartet. The optimized path
  keeps public bipartition extraction behavior but uses integer bitmasks for
  the internal all-quartet count loop, and bipartition extraction now reuses
  postorder descendant-tip sets instead of calling `get_terminals()` per
  internal clade. A later pass precomputes each quartet's combined mask and
  three topology pair masks once, so the per-tree split scan no longer rebuilds
  those masks for every topology check. The public frozenset bipartition helper
  now uses a direct standard-tree preorder plus reversed postorder setup while
  preserving the previous `get_nonterminals()` output order and generic
  fallback behavior. NANUQ distance and quartet report setup now use small
  stable comparison helpers for dominant and top-two topology indices, avoiding
  per-quartet `max()+index()` and three-item `sorted()` calls while preserving
  first-index tie behavior.
- `QuartetNetwork._extract_bipartition_masks` baseline time made a postorder
  pass to build clade masks, then a second Bio.Phylo nonterminal traversal to
  filter split masks. The optimized path emits eligible masks during the same
  postorder pass. A later pass uses a direct stack-based postorder traversal for
  standard Bio.Phylo trees, retaining `find_clades(order="postorder")` as the
  fallback for nonstandard tree objects. A follow-up mask-union pass handles
  binary and unary clades directly before the generic multifurcating loop,
  preserving polytomy skipping while reducing balanced 32768-tip mask extraction
  time from 0.053298s to 0.049024s.
- `QuartetNetwork._compute_circular_split_weights` equal-size canonical split
  tiebreaks now compare the minimum taxon on each complementary side instead of
  sorting both sides, preserving the sorted-lexicographic choice because the
  two sides are disjoint.
- `QuartetNetwork._build_splits_graph` baseline time enumerated every
  `2 ** n_splits` sign vector and then compared every valid-node pair to find
  graph edges. The optimized path reuses the existing forbidden-pair rules while
  assigning signs incrementally, pruning invalid partial assignments early, and
  finding edges by flipping one split at a time. A later startup pass replaces
  the eager Bio.Phylo import with a lazy `Phylo.read` proxy, preserving the
  tree-parser patch point while avoiding Biopython startup for import-only
  callers. A follow-up startup pass keeps JSON output behind a forwarding
  wrapper and localizes `PlotConfig` to argument processing, avoiding JSON and
  plotting-helper startup during command discovery. A later startup pass converts
  annotation-only `typing` aliases to built-in annotations so command discovery
  no longer loads `typing`. A traversal follow-up replaces the visited-flag
  postorder stack with a preorder collection plus one reverse step, preserving
  direct postorder while avoiding per-node tuple allocation during bitmask
  extraction.
- `QuartetNetwork._draw_quartet_network` now batches split-graph edges and
  pendant taxon edges into two `LineCollection`s, preserving the black
  linewidth styling while avoiding one Matplotlib artist per edge. Pendant-edge
  scaling now computes the split-graph x/y extent in one pass without temporary
  coordinate lists.
- `QuartetNetwork.run` text output now batches the summary header and
  per-quartet rows into one newline-joined print while preserving exact stdout
  text.
- `quartet_utils.compute_gcf_per_node` baseline time called `get_terminals()`
  repeatedly while extracting gene-tree splits and decomposing each species-tree
  branch into quartet groups. The optimized path computes descendant tip sets
  once per tree in postorder, then reuses those cached sets for split matching
  and species-node C1/C2/S/D decomposition. A later pass uses direct
  standard-tree preorder/postorder traversals instead of Bio.Phylo terminal and
  nonterminal materialization while retaining the generic traversal fallback.
  Another pass counts gene-tree support for each canonical split once, then
  looks up the concordant and discordant split frequencies per species branch
  instead of scanning every gene-tree split set three times per branch. A
  helper-level traversal pass now builds descendant tip sets by collecting
  preorder once and iterating it in reverse, preserving filtered tip-set output
  while avoiding visited-flag stack tuples. A follow-up helper pass combines
  binary child tipsets directly, reuses unary child tipsets, and retains the
  generic union path for multifurcations. Canonical split construction now uses
  a tuple container instead of a temporary list for the two frozenset members.
- `quartet_utils.parse_astral_annotations` and `parse_astral_branch_info`
  baseline time used Bio.Phylo's generic preorder iterator while scanning
  ASTRAL/wASTRAL labels. Both parsers now reuse the shared direct preorder
  traversal for standard parsed trees, keeping the fallback path inside
  `_preorder_clades` for nonstandard tree-like objects. The direct preorder
  helper now pushes binary children explicitly and indexes larger child lists
  backward, preserving Biopython preorder while avoiding `reversed(children)`.
- `ConcordanceAsr._compute_gcf_per_node` used the same repeated
  descendant-tip scans in its weighted-ASR-specific gCF implementation. The
  optimized path caches descendant tip sets once per species tree and gene
  tree, while keeping `_get_four_groups` compatible with callers that do not
  pass a cache. A later iterator pass uses the service's direct preorder helper
  for gene-tree internal split extraction and species-tree node iteration,
  retaining the generic Bio.Phylo iterators as fallbacks. Equal-size canonical
  split tiebreaks now compare the minimum taxon on each complementary side
  instead of sorting both sides, preserving the
  sorted-lexicographic choice because the two sides are disjoint.
- `ConcordanceAsr._run_distribution` result assembly baseline time recomputed
  descendant taxa with `get_terminals()` while matching species-tree nodes to
  per-gene-tree ASR estimates. The optimized weighted, distribution, uncertainty
  plot, and ASR re-keying paths compute descendant tip sets once per tree and
  reuse them for frozenset lookups and sorted output labels.
  `_run_asr_on_tree` setup now uses the shared terminal-name traversal and skips
  the second terminal-name scan when no pruning is needed. A later descendant
  cache pass replaces the remaining `find_clades(order="postorder")` helper
  with a direct reverse-preorder traversal and binary-child merge fast path,
  and the gene-tree list parser now streams cleaned rows directly into inline
  Newick/path parsing instead of first materializing a cleaned list.
  preserving the same clade-id-to-frozenset mapping with a generic fallback.
  A later no-prune
  pass also skips the protective tree copy entirely when all tree tips have
  trait values, while still copying before pruning missing trait taxa.
  `_normalize_taxa` and `run` species-copy pruning setup also use the shared
  terminal-name traversal for species and gene tree taxon extraction. A later
  run-level pass skips the second species-tree copy when every retained species
  tip has a trait value and ladderizing is disabled, while preserving copy
  isolation before trait pruning or ladderizing. Gene-tree source cleanup now
  streams over the input file and strips/filters comments in one pass instead
  of materializing `splitlines()`, preserving inline-Newick and path-list
  behavior. The path-list branch now resolves relative rows with a precomputed
  parent string prefix and a bound absolute-path check instead of constructing
  `Path` objects per listed tree. Weighted-ASR variance decomposition now uses
  scalar accumulation for the usual small topology-weight lists instead of
  allocating NumPy arrays. Distribution-ASR result assembly now computes
  per-node estimate means, population variances, and positive-sigma2 averages
  with scalar Python loops on the non-CI path, avoiding repeated NumPy dispatch
  while preserving CI percentile calculations when requested.
- `ConcordanceAsr._collect_uncertainty_node_data` baseline time then scanned
  every ancestral-estimate entry for every species-tree internal node while
  preparing uncertainty plots. The optimized path indexes entries once by their
  descendant tuple and preserves first-match behavior for duplicate descendant
  keys.
- `ConcordanceAsr._plot_concordance_contmap` setup baseline time delegated node
  labels and parent maps to helper methods using generic traversal, rebuilt
  rectangular coordinates with terminal/preorder/postorder scans, and scanned
  preorder again for branch and gCF-dot drawing. The optimized path derives
  labels, parent links, coordinates, tips, and plot loops from one direct
  preorder list plus the shared direct coordinate helper. A later setup pass
  passes that preorder list into `compute_node_positions` and passes the same
  preorder list plus the prepared tip list into `compute_circular_coords`,
  avoiding additional direct tree walks while preserving circular tip order. A
  later traversal-helper pass keeps preorder output identical while using a
  one-/two-child fast path instead of `reversed(children)` for direct preorder
  materialization. A later rectangular
  rendering pass batches horizontal and vertical gray tree branches into two
  `LineCollection`s while leaving gCF markers, labels, and colorbars unchanged.
  A later marker-rendering pass batches variable-size, variable-color gCF dots
  into one `scatter` call per plot while preserving black marker edges, labels,
  and colorbars. Uncertainty plot mean markers now compute marker positions
  with Python arithmetic over each per-node estimate list instead of dispatching
  to NumPy once per plotted row.
  A later startup pass
  defers direct NumPy imports and top-level Bio.Phylo loading behind lazy
  proxies, so command discovery avoids scientific-array and tree-format parser
  startup until gene-tree parsing or ASR calculations run. A follow-up startup
  pass also localizes the `AncestralReconstruction` helper import to `run()`,
  defers `pickle` behind a module-level proxy, and imports plot configuration,
	  circular-layout, and color-annotation helpers only during argument processing
	  or plot rendering. A later startup pass keeps JSON output behind a lazy
	  forwarding wrapper so command discovery avoids the shared JSON helper.
	  Another startup pass converts annotation-only `typing` aliases to built-in
	  postponed annotations so command discovery no longer loads `typing`. Text
	  output now batches the summary header and ancestral-estimate rows into one
	  newline-joined print while preserving exact stdout text for mixed CI/non-CI
	  reports.
- `EvoTempoMap._classify_gene_trees` baseline time repeatedly called
  `clade.get_terminals()` while extracting gene-tree branch-length splits,
  decomposing species-tree branches into C1/C2/S/D groups, and building branch
  labels. The optimized path computes clade-to-descendant-taxa maps once per tree
  and reuses them throughout classification. Shared-taxa/no-prune setup now uses
  the direct terminal-name traversal for species and gene-tree tip sets and
  reuses those gene-tree tip sets before deciding whether terminal-object
  pruning is needed. Gene-tree branch-length validation now also uses the
  existing direct preorder traversal for standard parsed trees while retaining
  the generic Bio.Phylo fallback. A later clade-taxa pass builds the direct
  postorder helper by reversing root-right-left preorder, preserving Biopython
  postorder while avoiding visited-flag stack tuples; split extraction remains
  dominated by canonical split construction rather than traversal overhead.
  The parent-map helper now builds the child-to-parent map directly instead of
  first materializing an ordered preorder list. A later direct-preorder helper
  pass avoids `reversed(children)` by pushing binary children right-then-left
  and indexing multifurcations backward. Gene-tree source cleanup now streams
  over the input file and strips/filters comments in one pass instead of
  materializing `splitlines()`, preserving inline-Newick and path-list behavior.
  The path-list branch now uses the same precomputed parent-prefix resolver,
  avoiding `Path` joins for every listed gene-tree file. A later startup pass
  defers `pathlib.Path` behind the path-list wrapper, so command discovery avoids
  loading `pathlib` until gene-tree path parsing is needed.
  Equal-size canonical split tiebreaks now compare the minimum taxon on each
  complementary side instead of sorting both sides, preserving the
  sorted-lexicographic choice because the two sides are disjoint.
- `EvoTempoMap._compute_global_treeness` still used repeated descendant-tip
  scans while building species-tree and gene-tree split sets for the global
  concordant/discordant treeness comparison. The optimized path reuses the same
  postorder clade-to-taxa maps already used by branch classification. A later
  pass moved those clade-to-taxa maps and cached split-set construction onto
  direct standard-tree traversals, retaining generic Bio.Phylo traversal as a
  fallback for nonstandard trees. A follow-up summary pass computes concordant
  and discordant treeness means and medians directly from the collected Python
  lists, avoiding NumPy startup when one group is too small for Mann-Whitney
  testing.
- `EvoTempoMap._compute_treeness` baseline time duplicated the previous
  two-traversal treeness calculation for every gene tree in global treeness
  comparisons. The optimized path reuses `Tree`'s one-pass internal/total branch
  length helper while preserving this helper's `0.0` result for zero-length trees.
- `EvoTempoMap._test_branch` eagerly converted concordant and discordant branch
  length lists into NumPy arrays before determining whether both groups had
  enough observations for Mann-Whitney and permutation tests. The optimized path
  computes early-return summary statistics directly from the lists and allocates
  arrays only after both groups are testable, preserving summary values while
  avoiding NumPy and SciPy startup for insufficient-data branches.
- `EvoTempoMap._fdr` baseline time sorted p-values into Python tuples and walked
  them in a reverse Python loop. The optimized path keeps the same
  Benjamini-Hochberg adjustment semantics, including ties, while using NumPy
  ordering and a reverse cumulative minimum for large arrays. A later small-list
  branch uses the same ordering and reverse cumulative-minimum semantics in
  Python, avoiding NumPy startup for direct helper use while retaining the
  vectorized path for large correction sets. A later pass deferred the
  `scipy.stats.mannwhitneyu` import until a statistical test has enough
  observations to run, preserving Mann-Whitney U behavior while avoiding
  `scipy.stats` on normal module import. A follow-up startup pass defers direct
  NumPy imports and Bio.Phylo loading behind lazy proxies, so command discovery
  avoids array and tree-format parser startup until statistical calculations,
  plotting, or gene-tree parsing run. The command run path now reuses the
  cached parsed species tree for read-only branch classification and global
  treeness setup. A later startup pass keeps JSON output behind a lazy forwarding
  wrapper and localizes `PlotConfig` construction to argument processing, avoiding
  those helper imports during command discovery. Another startup pass replaces
  annotation-only `typing` aliases with built-in postponed annotations, so command
  discovery no longer imports `typing`. Plot strip points now batch all
  concordant values into one scatter collection and all discordant values into
  another while preserving the deterministic jitter seed, point styling, boxplots,
  and significance labels. Significant branch labels now render through one text-star
  scatter collection, preserving their placement above each significant box group
  while avoiding one text artist per branch. Text output now batches the branch
  table, global summary, and optional verbose concordant/discordant length
  details into one newline-joined print while preserving exact stdout text.
- `Hybridization._count_topologies` and
  `DiscordanceAsymmetry._count_topologies` already used one-pass gene-tree split
  extraction, but still decomposed species-tree branches with repeated
  descendant-tip scans. The optimized path caches species-tree descendant taxa
  once and reuses that cache for C1/C2/S/D groups and branch labels. Species
  tree taxon-set setup now uses the shared terminal-name traversal. Equal-size
  canonical split tiebreaks now compare the minimum taxon on each complementary
  side instead of sorting both sides, preserving the sorted-lexicographic choice
  because the two sides are disjoint. Topology support counting now checks
  concordant, alt1, and alt2 splits in one pass over the gene-tree split list
  for each species branch instead of scanning that list three times. Plot setup
  now reuses the same cached
  descendant-taxa strategy for branch-result lookup instead of calling
  `clade.get_terminals()` for every internal species-tree node, while also using
  direct terminal and parent-map traversals for standard Bio.Phylo trees. A
  later mirrored terminal-helper pass preserves tip order while pushing binary
  children explicitly and indexing multifurcations backward, avoiding the
  `reversed(children)` iterator on common bifurcating
  trees. Gene-tree source cleanup now streams over the input file and
  strips/filters comments in one pass instead of materializing `splitlines()`,
  preserving inline-Newick and path-list behavior. Their path-list branches now
  resolve relative gene-tree rows through a precomputed parent string prefix and
  bound absolute-path check, preserving absolute-row behavior while avoiding
  per-row `Path` joins. `EvoTempoMap._output_text`
  now prepares verbose length details during the table-row pass and uses a
  cached table-row formatter, preserving exact stdout text while avoiding a
  second verbose traversal. Their FDR helpers now use the same small-list Python
  Benjamini-Hochberg path and large-array NumPy reverse cumulative-minimum path
  as `EvoTempoMap._fdr`, preserving scalar results including tied p-values while
  avoiding NumPy startup for direct small correction sets. A later pass removed
  the eager `scipy.stats` import by evaluating the symmetric two-sided binomial
  test with lazy `scipy.special.bdtr`. Follow-up startup
  passes for `hybridization` and `discordance_asymmetry` defer NumPy behind a
  proxy and import Bio.Phylo only inside gene-tree parsing, leaving command
  discovery free of array and tree-format parser startup. Later startup work
  keeps JSON output behind forwarding wrappers and localizes `PlotConfig`,
  cladogram layout, circular layout, and color annotation helper imports to
  argument parsing and plotting paths. A follow-up startup pass converts
  annotation-only `typing` aliases to postponed built-in annotations, so command
  discovery no longer loads `typing`. `Hybridization.run`
  now reuses the cached parsed species tree for read-only topology setup; when
  plotting with ladderization enabled it still reads a copied tree before
  mutating branch order. `DiscordanceAsymmetry.run` uses the same cached
  read-only species-tree path with the same copied-tree exception for
  ladderized plots. Their mirrored parent-map helpers now localize stack
  operations and push child lists directly because the result is a map rather
  than an ordered traversal. Their mirrored clade-taxon helpers now collect
  preorder once and build descendant tipsets from reversed preorder, combining
  binary child sets directly while retaining unary, multifurcating, and generic
  fallback behavior. `Hybridization._output_text` now batches the summary,
  significant-branch table, and non-significant count into one joined report
  while preserving stdout text and broken-pipe handling. It now renders
  significant rows in the same pass that counts non-significant rows and uses a
  cached row formatter, preserving exact table text while avoiding the
  intermediate significant-branch list. Later rendering passes batch
  rectangular horizontal
  branches, rectangular gray connectors, circular radial branches, and circular
  internal arcs into `LineCollection`s for both plotters while preserving
  per-branch colors and line widths. `Hybridization._plot` now builds direct
  preorder and postorder clade lists once and reuses them across coordinate
  setup, result matching, branch rendering, and star placement, retaining the
  `find_clades()` fallback for nonstandard tree-like objects.
  `DiscordanceAsymmetry._plot` now reuses the same direct traversal-list pattern
  across its rectangular and circular setup. A later startup pass converts
  annotation-only `typing` aliases to built-in postponed annotations, so
  command discovery no longer imports `typing`. `Hybridization` now batches
  significant branch star markers into one `scatter` call per plot, preserving
  red star styling while avoiding one marker collection per significant branch.
  A later `Hybridization._plot` pass caches colormap results by hybrid score
  inside each rectangular/circular plot call, preserving colors while avoiding
  repeated normalization and colormap calls when many branches share the same
  score. `DiscordanceAsymmetry` uses the same marker batching for
  FDR-significant branches with a favored alternative. A matching
  `DiscordanceAsymmetry._plot` pass caches colormap results by asymmetry ratio
  inside each rectangular/circular plot call, preserving colors while avoiding
  repeated normalization and colormap calls when branches share the same ratio.
  Both plotters now compute tiny child-y coordinate means with Python arithmetic
  during postorder setup, preserving node positions while avoiding one NumPy
  dispatch per internal branch.
  `DiscordanceAsymmetry._output_text` now batches the branch table, summary, and
  optional verbose per-branch details into one newline-joined print while
  preserving exact stdout text. It also prepares verbose details during the
  table-row pass and uses a cached table-row formatter, avoiding a second
  traversal in verbose output.
- `PolytomyTest._evaluate_tree_triplets_fast` baseline time called
  Biopython's `common_ancestor()` for every group triplet, repeatedly searching
  root-to-tip paths. The optimized path builds root-to-tip path and descendant
  clade-tip caches once per tree, then finds each triplet LCA from cached paths.
  The legacy sister-counter path now checks whether a triplet is a polytomy
  once before iterating tip pairs instead of repeating that traversal for each
  pair, and skips all pairwise distances for polytomy triplets. A later startup
  pass replaces the eager Bio.Phylo import with a lazy `Phylo.read` proxy,
  preserving the existing patch point while avoiding Biopython tree parser
  startup for import-only callers. A subsequent startup pass avoids importing
  `unittest.mock` solely to detect test doubles, replacing it with a lightweight
  module-name predicate. A follow-up startup pass keeps `mp.Pool`,
  `mp.cpu_count`, and `ThreadPoolExecutor` patch points as lazy wrappers,
  localizes `functools.partial` to the multiprocessing branch, and defers shared
  tree-copy `pickle`/`copy` imports until a tree object is actually copied. A
  path-cache LCA follow-up compares the three triplet path nodes directly
  instead of slicing `nodes[1:]` for every depth row, preserving first-divergence
  behavior while avoiding per-depth generator/list overhead. A
  later legacy-triplet pass replaces repeated
  `len(list(tree.get_terminals())) == 3` checks with a direct standard-tree
  terminal counter that exits after a fourth tip and retains the generic fallback
  for nonstandard tree objects. The triplet polytomy detector now likewise uses
  a direct standard-tree internal-node count that exits after the second
  internal clade while retaining the generic `get_nonterminals()` fallback for
  mocks and nonstandard tree objects. A follow-up fallback pass short-circuits
  those nonstandard `get_terminals()` and `get_nonterminals()` scans as soon as
  the exact-three-tip or single-internal-node result is known, avoiding full list
  materialization and full generator scans in false cases. Legacy triplet-tree
  preparation now builds
  the prune list by scanning `tips` against a three-item triplet set, avoiding a
  full `set(tips)` allocation for every triplet and making prune order
  deterministic. The fast triplet evaluator now reuses one triplet set for both
  tip-presence checks and sister-pair assignment, and precomputes group
  membership frozensets once per group block instead of rebuilding them for
  every represented triplet. Branch-length reset now uses a direct
  standard-tree stack traversal, with the generic `find_clades()` fallback kept
  for mocks and nonstandard tree objects. A later startup pass keeps JSON
  output behind a lazy forwarding wrapper so command discovery avoids the shared
  JSON helper. A subsequent startup pass converts annotation-only `typing`
  aliases to built-in postponed annotations so command discovery no longer loads
  `typing`. Gene-support frequency text output now batches the summary into one
  newline-joined print while preserving exact stdout text and broken-pipe
  handling. Triplet setup now gets the first group identifier through dictionary
  iteration instead of materializing all group keys, preserving insertion-order
  behavior while avoiding an O(n) allocation for large group maps. The tree-list
  reader now uses a local forwarding wrapper, preserving the module-level patch
  point while avoiding `phykit.helpers.files` during import-only command
  discovery.
- `TransferAnnotations.run` taxa validation setup baseline time collected
  target and source taxon sets through `get_terminals()`. The optimized path
  uses the shared direct terminal-name traversal while retaining fallback
  behavior through `get_tip_names_from_tree`.
- `TransferAnnotations._extract_annotations` and `_transfer` baseline time
  computed each internal clade's descendant tip names independently while
  building and matching bipartition keys. The optimized path computes
  descendant taxon sets once per source or target tree and reuses them for all
  bipartition lookups, while keeping `_get_bipartition` compatible with
  uncached callers. Clade-taxa collection now handles unary and binary child
  nodes directly with localized dictionary lookups, preserving multifurcation
  behavior while avoiding the generic child-combine helper on common internal
  shapes. A later startup pass replaces the eager Bio.Phylo import with a lazy
  `Phylo.read`/`Phylo.write` proxy, preserving those module patch points while
  avoiding Biopython startup for import-only callers. The run path now reads the
  source tree through the shared cache without copying because it
  is read-only; the target tree still uses the normal copied read before
  annotations are mutated onto it. A follow-up startup pass keeps JSON output
  behind a lazy forwarding wrapper so plain imports avoid the shared JSON helper.
  A later traversal pass uses a direct standard-tree postorder list, combines
  binary child taxon frozensets directly, and delays complement construction
  until the clade side is larger than half of the tree. The annotated-tree
  writer now sends Bio.Phylo output to an in-memory text handle, normalizes
  Biopython comment syntax there, and writes the final Newick file once instead
  of writing, reading, and rewriting the same path. Text summary output now
  batches the three-line report into one newline-joined print while preserving
  exact stdout text and broken-pipe handling. A later startup pass postpones
  annotations and converts annotation-only `typing` aliases to built-in
  annotations, so command discovery no longer loads `typing`.
- `color_annotations` now uses direct terminal traversal for MRCA validation,
  highlighted clade tip ids, and rectangular/circular range drawing. MRCA
  validation stops once all requested taxa have been found, avoiding full-tree
  terminal-name materialization for common small annotation lookups while
  retaining the original fallback for nonstandard clade objects. A later pass
  also collects highlighted clade branch ids with a direct descendant-node
  traversal, keeping `find_clades(order="preorder")` only as the fallback. A
  follow-up branch-id pass localizes stack operations and pushes child lists
  directly because the helper returns a set rather than an ordered traversal.
  A later ordered-traversal pass avoids `reversed(children)` iterator allocation
  in terminal-list extraction and MRCA taxa validation while preserving
  left-to-right tip order.
  Circular range wedge drawing now derives the minimum and maximum angular gaps
  in one pass over sorted tip angles instead of building two gap lists.
  Rectangular range drawing now computes independent x/y bounds in one tip scan
  without temporary coordinate lists. Label
  recoloring now scans Matplotlib text artists once and applies matching colors
  from the parsed label map. Color-file parsing now streams TSV rows directly
  from the file handle instead of first materializing `readlines()`, preserving
  label/range/clade parsing, comment filtering, and permissive malformed-row
  skipping. A later parser pass consolidates exact-case and mixed-case type
  dispatch for range and clade rows, preserving unknown-type skipping while
  avoiding duplicate row parsing branches. `Chronogram` first adopted the
  helper for rectangular and circular color-file label annotations; a follow-up
  pass routes the same helper through the remaining tree plotters with copied
  label-recolor loops.
- `circular_layout.compute_circular_coords` now computes tip index ranges with
  one direct postorder traversal and emits coordinates with one preorder
  traversal, replacing repeated Bio.Phylo terminal and clade scans while
  preserving the existing tip angle and internal midpoint semantics. Circular
  direct terminal and preorder helper traversals now push binary children
  right-then-left and index multifurcations backward, preserving Biopython order
  without allocating `reversed(children)` iterators. Circular
  branch drawing now reuses one direct preorder traversal for radial segments
  and internal arcs, caches arc interpolation fractions, scales child
  coordinates for radial segments instead of recomputing trigonometry, and
  scans child angles without temporary lists. A later pass computes each arc's
  61 polyline points with NumPy vectorized trigonometry while preserving the
  same interpolation fractions and Matplotlib plot call. Radial branch point
  interpolation now hoists sine and cosine outside the per-point loop. A later
  simple-axis fallback pass pushes binary children right-then-left and
  multifurcations backward, avoiding `reversed(children)` iterator allocation
  while preserving preorder plot-call order. The real-Matplotlib drawing path
  still batches radial segments and arc polylines into two `LineCollection`s,
  preserving the existing per-plot fallback for simple axes and nonstandard tree
  objects while avoiding thousands of Matplotlib `Line2D` artists during actual
  rendering. Gradient and discrete multi-segment
  circular branch helpers now use collection fast paths on real Matplotlib axes,
  with the original per-segment plot loops retained for simple axes. A follow-up
  whole-tree gradient helper batches many radial gradient branches into one
  `LineCollection`; `ContMap` and `AncestralReconstruction` use it for circular
  contMap rendering to avoid one collection per branch. Colored internal arcs now
  also have a whole-tree `LineCollection` helper, preserving the existing arc
  interpolation and sweep normalization while avoiding one `Line2D` artist per
  arc. Tip-label placement uses the same direct terminal traversal before placing
  text and localizes repeated math, coordinate, id, and axis-text lookups inside
  the placement loop. The legacy fallbacks are retained for nonstandard tree
  objects. A startup pass replaces local annotation-only `typing` aliases with
  built-in generics, so circular plotting helper import no longer loads `typing`.
- `plot_config.build_parent_map` now builds child-to-parent maps with direct
  traversal for standard Bio.Phylo trees, reducing shared setup time before
  rectangular and circular tree plotting helpers. A startup pass replaces
  annotation-only `typing` aliases with built-in generics, so plot argument
  processing no longer loads `typing` through this shared helper. A later pass
  localizes stack operations and pushes children directly because the helper
  returns a parent map rather than an ordered traversal.
- `plot_config.compute_node_positions` applies the same direct traversal pattern
  to shared rectangular phylogram and cladogram layout setup, preserving terminal
  y-order, branch-length x-coordinates, cladogram depth scaling, and internal
  midpoint y-coordinates. A follow-up pass computes internal-node midpoint
  y-coordinates with inline sum/count reductions instead of allocating a
  temporary child-y list for every internal node. Rectangular branch and
  tip-label drawing now also
  avoid generic tree traversal on standard trees. A later real-Matplotlib
  drawing pass batches horizontal and vertical rectangular branch segments into
  two `LineCollection`s while keeping the existing per-plot fallback for simple
  axes. The branch drawing loops now hoist parent, coordinate, color-callable,
  segment-append, and plot lookups out of the per-branch hot path for both
  simple axes and the `LineCollection` builder. A follow-up pass collects
  preorder once and computes internal y
  positions in reverse preorder, avoiding visited-flag stack entries while
  preserving terminal y-order and cladogram depth scaling. The direct terminal
  and preorder helpers now preserve Biopython traversal order while pushing
  binary children right-then-left and indexing multifurcations backward, avoiding
  `reversed(children)` iterator setup on standard clades.
- `BipartitionSupportStats.get_bipartition_support_vals` baseline time called
  `get_terminals()` for every supported internal node when preparing verbose
  bipartition output. The optimized path caches descendant terminal names in
  postorder and still iterates `get_nonterminals()` for output-order
  compatibility. A later pass preserves that preorder internal-node output with
  a direct stack traversal, removing the remaining generic `get_nonterminals()`
  scan on standard Bio.Phylo trees while retaining fallback behavior. Verbose
  text output now batches support rows into one newline-joined print while
  preserving empty-output behavior and stdout text. A later startup pass defers
  annotation-only Bio.Phylo imports and avoids importing NumPy for JSON scalar
  cleanup unless a summary calculation actually needs the shared statistics
  helper. Non-verbose text and JSON runs now collect support values directly
  into a NumPy array and skip descendant terminal-name construction entirely,
  while verbose output still uses the name-preserving path. Cached read-only
  `BipartitionSupportStats.run` now also uses the explicit unmodified tree read
  helper to avoid copying the cached parsed tree before extracting support
  values. A later startup pass defers stdlib JSON until JSON output is actually
  serialized. The latest startup pass also keeps summary-statistics helpers
  behind forwarding wrappers, avoiding `phykit.helpers.stats_summary` during
  command module import. Threshold summary rows now print as one newline-joined
  block while preserving exact stdout text and broken-pipe handling. A later
  startup pass removes the remaining annotation-only `typing` import by using
  built-in postponed annotations. The direct support-value array helper now
  pushes binary children explicitly and indexes larger child lists backward,
  preserving support order while avoiding the per-node `reversed(children)`
  iterator. The verbose direct helper now gets postorder data by reversing the
  preorder list, avoiding visited-flag stack tuples while preserving support
  and terminal-name output order. Verbose JSON builtin conversion now still
  scans nested payloads for NumPy values but reuses dict/list branches that are
  already builtins, avoiding a second allocation of large row payloads. Builtin
  scalar leaves now return before the NumPy scalar check, preserving NumPy
  conversion while reducing recursive conversion overhead for large JSON
  payloads that are already plain Python values.
- `BipartitionSupportStats.calculate_threshold_stats` previously counted values
  below each threshold twice. The optimized path counts once for a single
  threshold, and for multiple thresholds sorts support values once before using
  binary search for each cutoff. When support values already arrive as a NumPy
  array from the non-verbose extraction path, threshold counting now uses
  NumPy sort and `searchsorted`, preserving the list/bisect fallback. The
  single-threshold NumPy path now counts the boolean comparison with
  `np.count_nonzero`, avoiding a Python scalar loop for the common
  `--thresholds 90` style case. Large Python support-value lists now also use
  NumPy counting for one cutoff and the same NumPy sort/searchsorted path for
  many thresholds, preserving the Python fallbacks for small inputs to avoid
  unnecessary NumPy startup and conversion.
- `InternalBranchStats.get_internal_branch_lengths` used the same repeated
  terminal-name extraction for verbose internal-branch rows. The optimized path
  caches descendant terminal names once in postorder and keeps
  `get_nonterminals()` iteration for output-order compatibility. A later
  startup pass defers the annotation-only Bio.Phylo import; NumPy is still
  loaded through the shared statistics helper used for summaries.
- `NeighborNet._compute_distances_from_alignment` baseline time scanned every
  taxon pair and every site in Python to count comparable sites and residue
  matches. The optimized equal-length FASTA path builds one byte-backed
  alignment matrix, computes comparable-site counts with matrix products, and
  accumulates per-symbol match matrices before applying p-distance, identity,
  or Jukes-Cantor transforms. Uneven or non-ASCII sequences retain the legacy
  pairwise fallback. A later pass switched those dense count products from
  integer masks to float masks so NumPy can use the faster BLAS path while
  preserving the same count values. Clean, longer equal-length ASCII alignments
  now use direct blockwise sequence comparisons after the existing validity mask
  proves every cell is comparable; alignments with gaps or ambiguous symbols
  stay on the per-symbol match-product path.
- `NeighborNet._read_alignment` baseline time materialized `SeqRecord` objects
  before immediately converting them into a taxon list and uppercase sequence
  dictionary. The optimized loader uses `SimpleFastaParser`, preserving taxa
  order, first-token IDs, uppercase conversion, and last duplicate wins. A later
  pass reuses the shared ordered first-token parser, preserving duplicate IDs in
  the taxa list, multiline sequence joining, and internal whitespace removal
  while avoiding parser tuple construction. A startup pass imports the FASTA
  parser inside the read helper, so importing the neighbor-net command module
  does not load Bio.SeqIO.FastaIO. A subsequent startup pass postpones
  annotations and defers NumPy behind a module-level proxy while preserving the
  existing `np` patch point used by tests. A
  follow-up startup pass keeps JSON output behind a forwarding wrapper and
  localizes `PlotConfig` to argument processing, avoiding those helpers during
  import-only command discovery. Precomputed CSV distance matrices now parse
  each data row to floats and assign the row slice into the preallocated matrix
  in one operation, preserving labeled/unlabeled headers and short-row zero-fill
  behavior while avoiding per-cell NumPy assignments. Simple unquoted CSV
  matrices now parse numeric row text with `np.fromstring`, falling back to the
  existing `csv.reader` path for quoted fields or malformed simple rows. NJ
  circular ordering now extracts terminal names from standard Biopython trees
  with a direct order-preserving stack traversal after NJ construction,
  retaining the `get_terminals()` fallback for nonstandard tree-like objects.
- `NeighborNet._enumerate_circular_splits` now computes the full taxa set once
  per ordering and reuses it while canonicalizing each circular split, avoiding
  full-ordering materialization in the nested start/size loop. A later pass
  maintains each contiguous side as a rolling set for a fixed start position,
  avoiding per-size index-list construction while preserving split order and
  canonicalization. A follow-up pass avoids full complement-set construction for
  non-tie split sizes by deriving the shorter complement as its own circular
  arc, preserving the same canonical split order. A later complement-reuse pass
  builds the full taxa set once with indexed access and derives the complement
  from the rolling side set, preserving the no-full-order-iteration guard and
  the previous canonical split order. Equal-size split canonicalization now
  compares the minimum taxon on each disjoint half instead of sorting both
  halves; this preserves the sorted lexicographic tiebreak because equal-size
  disjoint sorted lists differ at their smallest element. Circular split
  validation now scans adjacent membership states directly and checks the
  wrap-around edge once, preserving the two-gap criterion while avoiding a
  modulo operation for every boundary. Non-circular splits now return as soon
  as the third boundary is seen, preserving the two-gap criterion while avoiding
  the rest of the ordering scan for rejected splits. Split-direction gap
  detection uses the same direct boundary scan and appends the wrap-around gap
  explicitly, preserving direction vectors while avoiding one modulo per taxon
  per split. Network drawing now carries accepted split gap positions from the
  filter step into direction setup, avoiding a second boundary scan for each
  accepted split while preserving filtered split order and direction vectors.
  Split-direction setup now also caches each taxon's circular cos/sin
  coordinates once per ordering instead of recomputing them for every accepted
  split membership centroid.
  The same circular split validation, direction-scan, and accepted-gap reuse
  pattern is mirrored in `ConsensusNetwork` and `QuartetNetwork`.
- `NeighborNet._estimate_split_weights` baseline time built the NNLS design
  matrix with Python loops over taxon pairs and circular splits. The optimized
  path builds a taxon-by-split membership matrix once and derives pairwise split
  separation with vectorized boolean XOR before calling the unchanged dense
  NNLS solver.
- `NeighborNet._build_splits_graph` baseline time enumerated all `2 ** n_splits`
  sign vectors before rejecting invalid combinations and then checked every
  valid-node pair for edges. The optimized path assigns signs incrementally and
  prunes forbidden split-pair states as soon as they are introduced, then builds
  edges by flipping one split sign and checking set membership. Network plotting
  now batches internal split-graph edges and pendant taxon edges into
  `LineCollection`s while preserving per-pendant linewidth and alpha. Pendant
  edge scaling now computes the split-graph x/y extent in one pass without
  temporary coordinate lists. A later
  fallback-rendering pass batches unlabeled no-split taxon points into one
  `scatter` collection instead of one marker artist per taxon. Text output now
  batches the analysis header, positive split rows, and output path into one
  newline-joined print while preserving exact stdout text.
- `RelativeRateTest._run_single` baseline time called the scalar Tajima test
  for every ingroup taxon pair, scanning every alignment site each time. The
  optimized equal-length ASCII path builds valid-site, differs-from-outgroup,
  and matches-outgroup boolean matrices, then obtains all pairwise `m1`/`m2`
  counts with matrix products before applying the same chi-square and multiple
  testing corrections. Uneven or non-ASCII sequences retain the legacy pairwise
  fallback. Its multiple-testing helpers keep list outputs while using NumPy for
  large correction arrays: Bonferroni correction applies the cap in one
  vectorized operation, and FDR computes the same Benjamini-Hochberg adjustment
  with sorting and a reverse cumulative minimum instead of tuple sorting plus a
  Python reverse loop. A later small-list path handles typical small correction
  sets with Python arithmetic, avoiding NumPy startup while preserving the large
  vectorized path. The standalone Tajima helper now uses a thresholded ASCII
  byte-array path for long sequences, preserving the scalar path for short or
  Unicode inputs. Clean ASCII triplets now sample and scan for skip codes before
  bypassing the validity mask, while ambiguous or gapped triplets keep the
  validity-mask path. The large vectorized pairwise path now computes df=1
  chi-square survival probabilities with the direct
  `erfc(sqrt(chi2 / 2))` formula via `scipy.special`, avoiding the heavier
  `scipy.stats.chi2.sf` distribution machinery while matching values to
  floating-point precision. The same large-pairwise path now scans upper
  triangle rows directly while building result dictionaries, avoiding the
  large `np.triu_indices` arrays and preserving `itertools.combinations` pair
  order. A later setup pass builds the ingroup byte matrix from one joined
  uppercase buffer and reshape, avoiding one `np.frombuffer` array per ingroup
  taxon plus the `np.vstack` join while preserving the Unicode fallback through
  the existing `UnicodeEncodeError` path. Small pair counts retain the scalar
  `_tajima_result` path to avoid importing SciPy for typical small alignments.
- `RelativeRateTest._identify_outgroup` now scans each root child only until it
  can prove whether that child is a singleton, avoiding full terminal-list
  materialization for a large ingroup clade when the rooted tree has one
  outgroup taxon. Nonstandard tree-like objects retain the generic fallback,
  and that fallback now applies the same singleton-only scan instead of
  materializing and sorting every root-child terminal list.
  Cached read-only `RelativeRateTest.run` now avoids copying the cached parsed
  tree because outgroup identification, alignment analysis, output, and optional
  plotting do not mutate it. Batch alignment-list cleanup now streams the list
  file and strips/filters comments in one pass instead of materializing
  `splitlines()`, preserving relative-path resolution and empty-list handling.
  Batch path resolution now computes the list-file parent prefix once and uses
  a bound absolute-path check for each row instead of constructing `Path`
  objects per alignment path. Module import now defers `pathlib.Path` until the
  batch alignment-list path is used.
  Text output for single-alignment and batch summaries now builds the full report
  and prints it once, preserving exact row formatting while avoiding one `print`
  call per result row. Batch output also summarizes each taxon-pair's gene
  results in one loop, collecting chi-square values and rejection counts together
  before applying the unchanged median and percentage formatting.
- `RelativeRateTest._plot_heatmap` baseline setup filled the heatmap matrix,
  then scanned every pairwise result again for each populated cell to decide
  whether to draw a significance marker. The optimized setup fills a symmetric
  FDR matrix alongside the heatmap values and finds significant cells directly
  from that matrix. Significance marker rendering now uses one scatter collection
  with a text-star marker instead of one Matplotlib text artist per significant
  cell. No-significance heatmaps now check the finite minimum FDR before
  extracting marker coordinates, avoiding a full empty `np.where` result when
  there is nothing to draw. A later startup pass defers NumPy behind a
  module-level proxy, so command
  import avoids numerical startup until vectorized tests or heatmap setup
  actually run. A follow-up startup pass keeps JSON output behind a forwarding
  wrapper and localizes `PlotConfig` to argument processing,
  avoiding those helper modules during import-only command discovery. A later
  startup pass converts annotation-only `typing` aliases to built-in annotations
  so command discovery no longer loads `typing`. The alignment reader now uses
  a same-name forwarding wrapper, preserving the module-level patch point while
  avoiding `phykit.helpers.files` during import-only command discovery.
- `CovaryingEvolutionaryRates.correct_branch_lengths` baseline time pickled
  trees, launched process batches, and repeatedly called `common_ancestor()` to
  match reference-tree branches to gene-tree branches. The optimized path first
  tries an exact-topology split map: descendant tip sets are computed once per
  tree and used to look up matching branch lengths directly. If any reference
  branch is not an exact split in either gene tree, the existing MRCA-based
  fallback is retained. A later pass removed the eager `scipy.stats` import by
  replacing z-score and Pearson coefficient calculations locally while retaining
  the beta-distribution p-value through lazy `scipy.special.betainc`. A later
  Pearson helper pass computes the correlation directly from centered dot
  products instead of normalizing two temporary vectors with `np.linalg.norm`.
  The exact
  split path now preserves the existing terminal-first/nonterminal row order
  with a direct preorder traversal instead of calling generic
  `get_terminals()`/`get_nonterminals()` after the postorder split-map pass.
  The branch-length map used for each exact-matching tree now also uses direct
  standard-tree postorder traversal before falling back to Bio.Phylo's generic
  `find_clades(order="postorder")`. Balanced binary nodes now combine the two
  child tip frozensets directly, while multifurcations retain the generic
  union path. The reference-tree split and tip-name
  metadata pass now uses the same direct postorder traversal, so the entire
  exact-split path avoids generic traversal for standard parsed trees. A later
  traversal pass builds both postorder metadata helpers by reversing
  root-right-left preorder, removing visited-flag stack entries while
  preserving split maps and terminal/nonterminal reference order. A later
  reference metadata pass avoids `reversed(children)` iterator allocation while
  preserving terminal/nonterminal reference order and combines binary child
  tipsets directly. A later
  startup pass defers NumPy behind a module-level proxy while keeping the same
  numerical call sites. A subsequent startup pass preserves the module-level
  `ProcessPoolExecutor` and `as_completed` patch points through lazy proxies,
  so import-only callers avoid concurrent-futures startup. A follow-up startup
  pass keeps module-level `pickle.dumps`/`pickle.loads` patch points behind a
  lazy proxy and imports `PlotConfig` only while processing command arguments.
  A later startup pass keeps JSON output behind a forwarding wrapper, removing
  the last eager helper import from command discovery. Verbose text output now
  batches branch rows into one newline-joined print while preserving exact row
  text and the separate plot-status message. Run setup now builds the shared-tip
  set once and scans the existing unique tree-tip lists to prepare deterministic
  prune lists without constructing three full temporary tip sets. Outlier
  filtering now filters NumPy branch-length arrays directly and uses normalized
  set membership for ordinary Python row lists, avoiding a temporary NumPy mask
  followed by a Python enumerate loop for every input type.
- `LastCommonAncestorSubtree.run` baseline time performed an extra
  pickle/unpickle copy after `read_tree_file()` had already returned a copied
  tree from the cache. The optimized path calls `common_ancestor()` directly on
  that already-isolated tree before writing the subtree.
- The second-stage `LastCommonAncestorSubtree.run` LCA lookup benchmark measured
  the remaining MRCA selection cost after copy removal. The optimized path maps
  terminal names once, builds root-to-node paths in preorder, and finds the
  deepest shared prefix for the requested taxa. A later pass combines terminal
  lookup and root-path construction into one direct stack traversal, retaining
  `common_ancestor()` as the fallback for nonstandard trees or missing taxa.
  A follow-up helper pass replaces per-node root-path list copies with parent
  pointers and depth-aligned ancestor climbs, preserving the direct traversal
  fallback behavior while reducing allocation during repeated MRCA lookups.
  A later helper pass only stores terminal clades whose names were requested,
  preserving the full-tree traversal and duplicate-name overwrite behavior for
  selected taxa while avoiding dictionary entries for unrelated terminal names.
  A subsequent localized-query pass stops the direct traversal once every
  requested terminal has been found, preserving the `common_ancestor()` fallback
  for missing taxa while avoiding unrelated subtrees. A later stack-push pass
  keeps left-to-right traversal order but handles binary children without
  allocating a `reversed()` iterator, retaining the indexed multifurcation path
  for non-binary nodes. The parent-depth target merge now walks target indices
  directly instead of slicing `targets[1:]`, preserving the ancestor-climb
  behavior while avoiding one per-query target-list copy. Single-unique-taxon
  requests now return the matching
  terminal clade from a direct traversal without building parent/depth maps,
  while missing taxa and nonstandard trees still fall back through
  `common_ancestor()`.
  JSON output now also counts subtree terminals with the shared direct terminal
  count helper before falling back to `count_terminals()`.
  The cached read-only setup follow-up switches `run()` to
  `read_tree_file_unmodified()`, avoiding an unnecessary copy of the cached
  parsed tree because the command only locates and writes an MRCA clade. A
  larger-tree side-by-side check confirms that avoiding the old pickle-copy
  core remains the dominant win even for opposite terminal MRCA lookups. A
  startup pass keeps JSON output behind a forwarding wrapper, avoiding JSON
  helper startup during command discovery. A later startup pass removes the
  annotation-only `typing` import by postponing annotations and using a built-in
  annotation. The taxa-list reader now uses a local forwarding wrapper,
  preserving the module-level patch point while avoiding `phykit.helpers.files`
  during import-only command discovery.
- `MonophylyCheck.run` exact-clade baseline time rooted the tree with every
  non-target tip as the outgroup before computing the MRCA. The optimized path
  first computes descendant taxon sets in postorder and, when the requested
  taxa already exactly match a clade, uses that clade directly. Non-exact cases
  retain the previous root-with-outgroup and MRCA behavior. Bootstrap support
  summaries now collect nonterminal confidence values with a direct stack
  traversal for standard clades, retaining `get_nonterminals()` fallback for
  nonstandard tree-like objects. A later startup pass postpones annotations and
  removes the annotation-only `Bio.Phylo.Newick` import, leaving tree parser
  startup to the shared tree reader instead of module import. The run path now
  reuses the cached parsed tree for exact-clade and insufficient-taxa checks;
  non-exact fallback rerooting first copies the tree so cached state remains
  isolated. A follow-up startup pass keeps JSON output behind a forwarding
  wrapper, avoiding JSON helper startup during command discovery. A later
  exact-clade lookup pass replaces per-node descendant `frozenset` construction
  with selected/terminal descendant counts for unique-tip standard trees, while
  falling back to the set-based path if duplicate terminal names make count
  semantics ambiguous. A later startup pass defers the summary-statistics helper
  used for bootstrap support summaries until those summaries are computed. Text
  output now batches mixed support/status rows into one newline-joined print
  while preserving exact stdout text and the existing offending-taxa sorting. A
  later startup pass removes annotation-only `typing` aliases by using built-in
  postponed annotations. The bootstrap support-value scan now pushes child lists
  directly because the values are immediately reduced to summary statistics, so
  traversal order is not observable. The exact-clade count lookup now uses a
  direct reverse-preorder traversal for standard parsed trees, retaining the
  generic postorder fallback for nonstandard tree-like objects and the
  duplicate-tip set fallback. A follow-up pass carries exact-clade descendant
  counts upward in explicit stack frames, avoiding per-node count dictionaries
  while keeping the same duplicate-tip fallback and nonstandard-object failure
  handling. The taxa-list reader now uses a local forwarding wrapper, preserving
  the module-level patch point while avoiding `phykit.helpers.files` during
  import-only command discovery.
- `HiddenParalogyCheck.run` sequential exact-clade baseline time reread and
  rerooted the tree for each requested clade. The optimized sequential path
  builds exact descendant taxon sets once from the master tree and classifies
  exact clades immediately; non-exact clades still use the previous fresh-tree,
  root-with-outgroup, and MRCA path. Cached read-only master-tree setup now
  avoids copying the cached parsed tree used only for tip extraction and
  exact-clade indexing. A follow-up exact-index pass builds the postorder clade
  list with a direct standard-tree stack traversal before falling back to
  `find_clades(order="postorder")`. Binary exact-index nodes now combine the two
  child descendant-taxon sets with direct set union, keeping the variadic union
  path for polytomies. Requested clade filtering now intersects the existing
  master-tip set with each requested clade iterable directly, avoiding a
  temporary `set(clade)` allocation while preserving duplicate and off-tree taxon
  semantics in both sequential and batch paths.
- `HiddenParalogyCheck._process_clade_batch` exact-clade baseline time reread
  and rerooted the tree for each clade in a multiprocessing batch. The
  optimized batch path reads the tree once per batch to build the exact-clade
  index, skips rereads for exact clades, and still falls back to fresh-tree
  rerooting for non-exact clades.
- `HiddenParalogyCheck._process_clade_batch` shared exact-clade index baseline
  time still read and indexed the tree once per batch. The multiprocessing run
  path now passes the master exact-clade index into each batch worker, avoiding
  per-batch tree reads for exact-clade-only batches while retaining the direct
  helper fallback that builds an index when none is provided.
- `HiddenParalogyCheck._process_clade_batch` non-exact fallback now collects
  MRCA terminal names with a direct standard-clade stack traversal before
  falling back to Bio.Phylo `get_terminals()` for nonstandard clade-like
  objects. A follow-up pass avoids order-preserving child reversal in that helper
  because it returns a set, preserving the same terminal-name result while
  reducing stack traversal overhead. Text output now batches clade status rows
  into one newline-joined print. Clade-file parsing now bulk-reads rows before
  applying the same whitespace split, preserving empty clade rows while reducing
  per-line iterator overhead on large clade lists. A later parser pass iterates
  the file handle directly with the same whitespace split, preserving blank-line
  compatibility while avoiding the intermediate `splitlines()` list.
  into one newline-joined print while preserving empty-output behavior and stdout
  text. Clade-list parsing now iterates directly over the file object instead
  of materializing `readlines()` first, while preserving blank-line entries as
  empty clades. A later startup
  pass wraps `Phylo.read` in a lazy module-level proxy, preserving existing
  patch points while avoiding Bio.Phylo parser startup until a tree read is
  actually needed. A subsequent startup pass keeps `mp.cpu_count` and `mp.Pool` behind
  a lazy proxy and imports `functools.partial` only when the large-clade
	  multiprocessing branch is used.
	  A follow-up startup pass keeps JSON output behind a
	  forwarding wrapper, avoiding JSON helper startup during command discovery.
	  A later startup pass converts annotation-only `typing` aliases to postponed
	  built-in annotations, so command discovery no longer loads `typing`.
- `BranchLengthMultiplier.run` baseline time also made a second full-tree
  pickle/unpickle copy before mutating branch lengths. The optimized path uses
  the isolated tree returned by `read_tree_file()` directly, then writes that
  scaled tree. A later traversal pass avoids building separate terminal and
  nonterminal lists or using Bio.Phylo's generic traversal before multiplying
  branch lengths, collecting standard clades with a direct stack traversal
  instead. The command module now also defers annotation-only Bio.Phylo imports
  and benefits from the shared JSON helper's lazy NumPy handling at startup.
  A follow-up startup pass keeps the module-level `print_json` patch point as a
  forwarding wrapper, avoiding JSON helper startup during command discovery.
  A later factor-one no-op pass reads the cached tree without copying or
  multiplying branches; JSON output still counts branch-bearing nodes
  read-only so the `scaled_branches` field remains compatible. A
  later traversal pass multiplies and counts branches while popping the direct
  traversal stack, avoiding a temporary all-node list while preserving fallback
  behavior for nonstandard tree objects. A later startup pass removes the
  remaining annotation-only `typing` import by using a built-in postponed
  annotation. The factor-one JSON count helper now localizes stack operations
  inside the read-only direct traversal, preserving branch-count semantics while
  reducing balanced 131072-tip median count time from 0.015458s to 0.014903s.
  The mutating branch-scaling helper now uses the same localized stack
  operations during direct traversal, reducing balanced 131072-tip median scale
  time from 0.068762s to 0.055513s.
- `DVMC.determine_dvmc` baseline time traversed terminals once to count species
  and again to build the terminal list used for depths. The optimized path uses
  the length of the terminal list already required for distance extraction,
  preserving the same depth-based fast path and distance fallback.
- `DVMC.determine_dvmc` one-pass terminal-distance pass replaces the remaining
  terminal-list plus `tree.depths()` traversals for parsed Biopython trees with
  a direct stack traversal that carries root-relative distances to terminal
  nodes. The existing depth/distance logic remains as a fallback for nonstandard
  tree-like objects. The command module now imports NumPy only when calculating
  DVMC, avoiding NumPy and annotation-only Bio.Phylo startup for command import.
  Cached read-only `DVMC.run` now also uses the explicit unmodified tree read
  helper to avoid copying the cached parsed tree before running the statistic.
  A follow-up startup pass keeps JSON output behind the same module-level
  forwarding wrapper used by other lightweight tree-output commands. A later
  scalar-statistics pass accumulates count, sum, and sum of squared
  root-to-tip distances during the standard-tree traversal, so ordinary parsed
  trees no longer import NumPy for the DVMC calculation. A subsequent startup
  pass removes the remaining annotation-only `typing` import with a built-in
  postponed annotation, so command discovery no longer loads `typing`.
- `RobinsonFouldsDistance.calculate_robinson_foulds_distance` now compares
  rooted descendant split sets using compact per-tip integer ids for standard
  trees, preserving the public name-based `get_all_bipartitions()` API and the
  legacy rooted split semantics. The multiple-tree batch helper routes through
  the same optimized calculation. Same-object comparisons now skip split
  extraction and return zero RF through the existing terminal-count
  normalization denominator, preserving the three-tip denominator edge case.
  A later startup pass preserves the module
  `ProcessPoolExecutor` patch point through a lazy proxy so import-only callers
  avoid concurrent-futures startup. A follow-up startup pass defers pickle until
  the multiple-tree batch path serializes work. A later split-helper pass
  special-cases binary child unions with direct frozenset `|` operations in both
  named and compact-id direct helpers, while keeping explicit unary and
  multifurcating paths; this reduces balanced 32768-tip compact split-id
  extraction from 0.048783s to 0.039282s. A follow-up direct-helper pass builds
  a preorder clade list and processes it in reverse, avoiding per-edge visited
  stack tuples while preserving the same fallback signal for malformed standard
  tree objects; compact split-id extraction drops further from 0.037405s to
  0.033784s. The direct terminal-count helper used for RF normalization now
  localizes stack operations, shaving balanced 262144-tip count time from
  0.027284s to 0.026463s. Another startup pass removes the annotation-only
  `typing` import by using postponed built-in collection annotations.
- `TreenessOverRCV` now defers importing the alignment RCV helper until `run()`
  computes the alignment component, reducing cold command-module import cost
  while preserving text and JSON output. A follow-up startup pass keeps JSON
  output behind a module-level forwarding wrapper, preserving the patch point
  while avoiding JSON helper startup during command discovery. A later startup
  pass removes the annotation-only `typing` import by using postponed built-in
  annotations. A later run-path pass reads the tree through the unmodified
  cached-tree path and passes it into `calculate_treeness`, avoiding the
  previous pickle copy for the read-only treeness calculation while leaving the
  alignment RCV calculation unchanged.
- `Tree.get_tip_names_from_tree` baseline time materialized terminal clade
  objects through Bio.Phylo before extracting names. The optimized helper walks
  parsed tree or clade objects directly while preserving terminal order and
  falling back to `get_terminals()` for test doubles and nonstandard objects.
  First-tip lookup now descends the leftmost child chain directly instead of
  maintaining a traversal stack, preserving the same first terminal while
  reducing depth-17 balanced-tree lookup time from 2.652us to 1.049us.
  A later pass keeps the same order-preserving traversal but appends children
  directly instead of extending from a `reversed()` iterator in terminal-name and
  terminal-clade collection. Follow-up terminal-name and terminal-clade passes
  special-case binary clades with explicit right-then-left stack pushes,
  preserving output order while avoiding a per-node `reversed()` iterator on the
  common bifurcating-tree path. A later terminal-count pass binds stack
  operations in the shared count-only and total-length-plus-count helpers,
  shaving loop overhead for commands that only need tip counts or count/length
  summaries.
  `TipLabels.run` now uses the same helper before text or JSON output instead
  of materializing terminal clade objects itself. The shared tree base module
  now also defers importing Bio.Phylo and NumPy until tree I/O or terminal root
  distance arrays are actually requested. A later startup pass removes its
  annotation-only `typing` import using postponed built-in annotations. A
  follow-up startup pass localizes cache-key hashing to `_get_file_hash`, so
  import-only callers avoid `hashlib` until tree-file reads need cache
  invalidation. A later cache-key pass removes hashing from tree reads too,
  using the raw path/size/mtime_ns key for the same cache invalidation inputs. For
  cached read-only runs,
  `TipLabels.run` now also uses the explicit unmodified tree read helper to
  avoid copying the cached parsed tree while preserving the defensive-copy
  default for tree-mutating commands. A follow-up startup pass keeps JSON
  output behind a module-level forwarding wrapper, preserving the patch point
  while avoiding JSON helper startup during command discovery. A later startup
  pass removes the annotation-only `typing` import, so command discovery no
  longer loads `typing`. JSON output row construction now uses literal row
  dictionaries instead of `dict(taxon=...)` calls while preserving the same
  `rows` and `tips` payload.
  longer loads `typing`.
- `tree_paths.build_root_path_map` now builds root-to-node path lists with a
  direct stack traversal for standard trees rooted at `tree.root`, avoiding
  Bio.Phylo's generic preorder iterator in ancestral-reconstruction path setup.
  Explicit alternate-root calls and nonstandard tree objects retain the legacy
  traversal behavior. The shared object-parent-map helper now localizes stack
  operations and pushes child lists directly because the result is a mapping,
  not an ordered traversal.
- `LTT` baseline time computed terminal lists and root-depth maps separately
  for the gamma statistic and lineage-through-time table. The optimized run path
  builds that setup context once and passes it into both helpers, while the
  helper defaults retain the existing standalone API behavior. Terminal clade
  setup now uses direct standard-tree traversal before falling back to
  `get_terminals()`. A later terminal helper pass preserves tip order while
  pushing binary children explicitly and indexing multifurcations backward,
  avoiding the `reversed(children)` iterator on common bifurcating trees.
  Root-depth setup now uses a direct stack traversal for
  standard parsed trees before falling back to generic `tree.depths()`. With
  the shared depth map available, `_compute_gamma` and `_compute_ltt` now
  collect internal node depths with direct stack traversals instead of separate
  `find_clades(order="level")` scans, while retaining those scans as a fallback
  for nonstandard tree-like objects. Cached read-only `LTT.run` now also uses
  the explicit unmodified tree read helper to avoid
  copying the cached parsed tree before computing terminal, depth, gamma, and
  LTT data. A later startup pass keeps JSON output behind a local forwarding
  wrapper and localizes `PlotConfig` to argument processing, avoiding those
  helper imports for import-only callers. A subsequent startup pass removes the
  annotation-only `typing` import under postponed annotations, so command
  discovery no longer imports `typing`. Verbose text output now batches the
  branching-time and lineage-through-time rows into one newline-joined print
  while preserving exact stdout text. A later run-path pass computes gamma and
  LTT output from one shared internal-depth list, avoiding the second internal
  traversal and sort while preserving the standalone helper APIs. A later
  internal-depth helper pass avoids the `reversed(children)` iterator by pushing
  binary children right-then-left and indexing multifurcations backward. The
  standalone gamma helper now accumulates the cumulative-statistic sum directly,
  avoiding temporary partial and cumsum lists while preserving the returned
  gamma, p-value, branching times, and internode intervals. A follow-up gamma
  pass combines the `ST` total and cumulative-statistic accumulation in one
  shared helper, avoiding a second pass over the internode intervals in both
  standalone and combined gamma/LTT paths. The combined gamma/LTT helper now
  builds LTT rows with an iterator helper instead of slicing
  `internal_depths[1:]`, preserving the same row sequence without copying the
  sorted internal-depth tail.
- `CollapseBranches.run` baseline time counted internal nodes before and after
  collapsing even for normal text output, where the collapsed branch count is
  not reported. The optimized path performs those count traversals only for
  JSON output and keeps the collapse predicate and tree-writing behavior
  unchanged. JSON count setup now counts internal nodes on standard Bio.Phylo
  trees with a direct stack traversal and falls back to `get_nonterminals()`
  for nonstandard tree objects. The no-op collapse path now pre-scans standard
  trees for weak support values and skips Bio.Phylo `collapse_all()` when no
  branch can be collapsed, falling back to the original collapse call whenever
  the scan finds a candidate or cannot inspect the tree directly. A later
  read-only no-op pass scans the cached tree before copying it and only copies
  when a branch will actually be collapsed, preserving cache safety for
  mutating paths while avoiding the pickle copy for unchanged trees.
  JSON no-op setup now combines the standard-tree internal-node count and
  weak-support scan in one traversal while preserving the early-exit scan used
  by non-JSON runs. A later non-JSON pass routes standard trees directly to
  the boolean weak-support scan instead of the JSON count-and-scan helper, so
  runs with an early weak branch can proceed to copy/collapse without counting
  the entire tree first. A follow-up startup pass keeps JSON output behind a
  module-level forwarding wrapper, preserving the patch point while avoiding
  JSON helper startup during command discovery. Another startup pass removes
  the annotation-only `typing` import with a built-in postponed annotation, so
  command discovery no longer loads `typing`.
- `RenameTreeTips.run` baseline time made a second full-tree pickle/unpickle
  copy before renaming terminal labels. The optimized path renames the isolated
  tree returned by `read_tree_file()` directly before writing it. A later
  startup pass defers its annotation-only Bio.Phylo import until tree I/O is
  actually requested by the shared tree base. Its
  terminal-renaming helper now collects standard Bio.Phylo terminal clades with
  a direct stack traversal and falls back to `get_terminals()` for nonstandard
  tree objects. A follow-up startup pass keeps JSON output behind a module-level
  forwarding wrapper, preserving the patch point while avoiding JSON helper
  startup during command discovery. A no-match pass scans the cached tree
  read-only before copying; when no tip names match the id map, it writes the
  unchanged tree and reports zero renamed tips without mutating cached state.
  A later default-output pass short-circuits after finding the first matching
  tip before copying, preserving the exact renamed-tip count for JSON runs while
  avoiding a full pre-count traversal for normal text output. The terminal
  renaming helper now applies mapped tip names during direct traversal instead
  of materializing terminal clades for a second rename pass. A later read-only
  matching pass lets the JSON count and non-JSON boolean preflight scans push
  children without preserving left-to-right order, since those helpers only
  report a count or existence check and do not expose traversal order. The same
  unordered child-push pattern is now used by the direct terminal-renaming pass,
  because the exposed result is only the renamed tree and count. Empty id maps
  now return immediately from the has/count/replace helper paths, avoiding a
  full-tree scan when a mapping file has no usable entries.
- `RootTree.run` baseline time made a second full-tree pickle/unpickle copy
  before rerooting. The optimized path reroots the isolated tree returned by
  `read_tree_file()` directly before writing it. A later startup pass preserves
  `Phylo.BaseTree.Tree.root_with_outgroup` as a patchable module path through a
  lazy proxy, so command import no longer initializes Bio.Phylo. A follow-up
  startup pass keeps JSON output behind a module-level forwarding wrapper,
  preserving the patch point while avoiding JSON helper startup during command
  discovery. The outgroup-list reader now uses the same forwarding-wrapper
  pattern, avoiding `phykit.helpers.files` import during root-tree command
  discovery while preserving the module-level patch point used by tests.
- `PruneTree.run` baseline time made a second full-tree pickle/unpickle copy
  before pruning. The optimized path prunes the isolated tree returned by
  `read_tree_file()` directly before writing it. A later no-op pass computes
  the prune list against the cached tree first and only copies/prunes when that
  list is non-empty, so keep-all runs write the unchanged tree read-only.
- `PruneTree.run` keep-mode complement baseline time checked every tree tip
  against the original taxa list, making complement construction scale with
  `tips * input_taxa`. The optimized path builds a taxa set for membership
  checks while preserving tree-tip order in the generated prune list.
- `PruneTree.run` keep-mode and `--ignore-branch-labels` selection baseline
  time still materialized terminal clade objects through Bio.Phylo before
  extracting names. The optimized path reuses the shared direct terminal-name
  traversal while preserving tree-tip order and branch-label stripping behavior.
  JSON output now uses the shared direct terminal count helper for
  `remaining_tips`, retaining `count_terminals()` as the fallback for
  nonstandard tree-like objects. A follow-up startup pass keeps JSON output
  behind a module-level forwarding wrapper, preserving the patch point while
  avoiding JSON helper startup during command discovery. A later startup pass
  removes the annotation-only `typing` import by using postponed built-in
  annotations. Branch-label stripping now checks for `{` before running the
  compiled label regex, avoiding regex work for ordinary tip names while
  preserving labeled-tip cleanup. The taxa-list reader now uses a local
  forwarding wrapper, preserving the module-level patch point while avoiding
  `phykit.helpers.files` during import-only command discovery.
- `InternodeLabeler.run` baseline time made a second full-tree
  pickle/unpickle copy before assigning internal-node labels. The optimized
  path labels the isolated tree returned by `read_tree_file()` directly before
  writing it.
- `InternodeLabeler.add_labels_to_tree` baseline time materialized all
  internal nodes through Bio.Phylo `get_nonterminals()`. The optimized path
  labels standard Bio.Phylo trees with a direct preorder traversal and keeps the
  original traversal as a fallback for nonstandard tree objects. The command
  module now defers annotation-only Bio.Phylo imports and inherits the lazy
  shared JSON helper startup behavior. A follow-up startup pass keeps JSON
  output behind a module-level forwarding wrapper, preserving existing patch
  points while avoiding JSON helper startup during command discovery. A later
  startup pass removes the annotation-only `typing` import, so command
  discovery no longer loads `typing`. A later direct-helper pass avoids
  `reversed(children)` iterator allocation while preserving preorder label
  assignment by pushing binary children right-then-left and using indexed reverse
  order for multifurcations.
- `CollapseBranches.run` baseline time made a second full-tree pickle/unpickle
  copy before counting and collapsing supported branches. The optimized path
  collapses the isolated tree returned by `read_tree_file()` directly before
  writing it.
- `Spr.run` single-taxon setup baseline time materialized terminal clade lists
  for both taxon validation and subtree lookup. The optimized path builds the
  terminal-name map once with direct standard-tree traversal, keeping the
  generic terminal fallback for nonstandard tree objects. A later helper pass
  preserves preorder while pushing binary children explicitly and indexing
  multifurcations backward, avoiding the `reversed(children)` iterator in the
  clade-list plus parent-map setup.
- `Spr._print_output_summary` now batches the fixed output-file summary into one
  newline-joined print while preserving exact stdout text. Stdout mode now batches
  generated SPR Newick strings into one newline-joined print while preserving the
  same stream text and broken-pipe handling.
- `Spr._generate_spr_trees` baseline time matched the copied subtree and each
  copied regraft target by recomputing descendant taxon sets after every full
  tree copy. The optimized generator indexes original preorder clades once,
  uses the same preorder positions in each copied tree before mutation, and
  precomputes regraft descriptions from cached descendant taxa. A later pass
  replaces remaining Bio.Phylo `find_clades()` traversals in the generator with
  direct preorder/postorder helpers and builds copied clade lists plus parent
  maps in one traversal. A subsequent startup pass defers Bio.Phylo writing and
  Newick `Clade` construction behind lazy module-level proxies. A later startup
  pass removes an unused `json` import, keeps JSON output behind a lazy wrapper,
  and defers pickle until tree copies are generated. A follow-up startup pass
  removes annotation-only `typing` aliases by using postponed built-in
  annotations. A later regraft-description pass builds the descendant-taxon
  cache with a reverse-preorder postorder helper and a binary-node union fast
  path, avoiding visited stack tuples and repeated empty-set unions while
  preserving Biopython postorder. A follow-up subtree preorder pass avoids
  `reversed(children)` when collecting subtree node IDs, using the same
  right-then-left binary push and indexed multifurcation traversal as the
  clade-list setup.
- `NearestNeighborInterchange._build_parent_map` baseline time called
  `tree.get_path()` for every clade while building parent links for NNI
  generation. The optimized helper assigns each child to its parent during one
  tree traversal. A later pass moves standard Bio.Phylo trees onto a direct
  stack traversal and keeps `find_clades()` only as the fallback for
  nonstandard tree-like objects. A startup pass also postpones annotations and
  defers `Phylo.write` behind a lazy proxy, preserving the write patch point. A
  follow-up startup pass keeps JSON output behind a module-level forwarding
  wrapper, preserving the patch point while avoiding JSON helper startup during
  command discovery. A later startup pass also keeps pickle behind a lazy proxy
  until NNI tree copies are generated. A follow-up startup pass converts the
  runtime branch-spec type alias and annotation-only `typing` aliases to built-in
  generics, keeping command discovery free of `typing`. A later helper pass
  localizes stack operations and pushes children directly because the result is
  a parent map rather than an ordered traversal. NNI neighbor generation now
  also uses a direct standard-tree level-order nonterminal scan, retaining
  `get_nonterminals(order="level")` as the fallback for nonstandard trees.
- `NearestNeighborInterchange._generate_targeted_nnis` terminal-name setup
  baseline time materialized terminal clade objects before validating targeted
  branch taxa. The optimized path uses the shared direct terminal-name
  traversal and retains `get_terminals()` as a fallback for nonstandard trees.
- `IndependentContrasts.run` baseline time made a second full-tree
  pickle/unpickle copy before validation, optional pruning, polytomy
  resolution, and PIC traversal. The optimized path performs those mutations on
  the isolated tree returned by `read_tree_file()` directly. Tip-name setup now
  also uses the shared direct terminal-name traversal before matching trait
  taxa. JSON output now builds the contrast-node list with a comprehension and
  serializes the all-Python payload directly with `json.dumps(sort_keys=True)`,
  avoiding the generic recursive NumPy-to-builtin conversion copy while
  preserving sorted-key JSON output. Polytomy resolution now uses a direct stack
  traversal for standard Bio.Phylo trees and lazily imports the Newick clade
  helper only when an actual multifurcation must be resolved, avoiding generic
  traversal overhead for already-binary trees. A follow-up no-op binary-tree
  pass resolves polytomies during the direct traversal instead of first
  materializing the full clade list, and skips child-iterator setup for terminal
  clades. The no-op resolution scan now pushes child lists directly because
  resolution is local to each node; nested multifurcation output remains
  identical to the previous left-to-right scan. Text output now batches the
  header, contrast rows, and summary statistics into one newline-joined print
  while preserving stdout text. A later direct-postorder pass builds the PIC
  traversal list with reverse preorder instead of stack entries carrying
  visited flags, preserving left-to-right postorder while reducing helper
  overhead. A later startup pass defers NumPy behind a module-level proxy, so
  import-only callers avoid numerical startup. A follow-up startup pass defers
  stdlib JSON until JSON output is actually serialized. A later trait setup
  pass streams the trait file, selects shared traits in tree-tip order without
  a set intersection, and skips the prune-membership set when all tree tips
  have trait values. A text-output pass uses compact descendant-tip previews
  and counts for deep trees, avoiding full descendant-label materialization
  when the text renderer only displays the first three labels and total count;
  JSON output and balanced text output keep the full-label path. A subsequent
  run-setup pass reads the cached tree directly for binary all-shared inputs and
  copies only before default branch-length assignment, pruning, or polytomy
  resolution would mutate the tree. A later parser pass returns the parsed
  trait dictionary directly when its taxon order exactly matches tree-tip order,
  preserving dict/JSON order while avoiding the tree-tip-order filtering copy
  for all-shared files. A later startup pass postpones annotations and removes
  the annotation-only `typing` import, so command discovery no longer loads
  `typing`. A boolean-scan pass pushes child lists directly in the default
  branch-length and polytomy checks, since those helpers only return whether a
  condition exists and do not expose traversal order.
- `ParsimonyScore.run` baseline time made a second full-tree pickle/unpickle
  copy before validation, polytomy resolution, optional pruning, and Fitch
  traversal. The optimized path performs those mutations on the isolated tree
  returned by `read_tree_file()` directly. Taxon-set setup now also uses the
  shared direct terminal-name traversal before matching alignment taxa. Verbose
  JSON output now serializes the all-Python payload directly with
  `json.dumps(sort_keys=True)`, avoiding the generic recursive JSON conversion
  over large per-site score lists while preserving sorted-key output. Verbose
  text output now batches the total score, blank separator, and per-site rows
  into one newline-joined print while preserving nonverbose integer output,
  JSON output, and stdout text. A follow-up startup pass defers stdlib JSON
  until JSON output is actually serialized. A subsequent run-setup pass reads
  the cached tree directly for binary all-shared inputs and copies only before
  polytomy resolution or pruning would mutate the tree. The Fitch scorer now
  returns zero total and per-site scores immediately for identical shared
  sequences after confirming every terminal has sequence data, preserving the
  previous missing-terminal `KeyError` path while avoiding the internal-node
  pass on conserved alignments.
- `ParsimonyScore._resolve_polytomies` baseline time scanned every clade with
  Bio.Phylo's postorder traversal, even for already-binary trees. The optimized
  setup path uses direct standard-tree postorder traversal and retains the
  original traversal fallback for nonstandard tree objects. A later startup
  pass defers NumPy behind a module-level proxy and imports the Newick clade
  helper only when a multifurcation is actually resolved. A follow-up direct
  postorder pass builds that traversal list with reverse preorder instead of
  visited-flag stack entries, preserving left-to-right postorder while reducing
  no-op dichotomy scan overhead. A later no-op binary-tree pass resolves during
  direct traversal instead of first materializing the full clade list and skips
  child-iterator setup for terminal clades. The preliminary polytomy existence
  scan now also pushes child lists directly because it returns only a boolean.
  The no-op resolution scan likewise pushes child lists directly because each
  node's resolution is local; nested multifurcation output remains identical to
  the previous left-to-right scan.
- `parsimony_utils.build_parent_map` baseline time used Bio.Phylo's generic
  preorder iterator to build parent links for downstream ACCTRAN, DELTRAN, and
  change-classification helpers. The optimized path walks standard parsed
  trees directly with an explicit stack and keeps the generic preorder fallback
  for nonstandard tree objects. A later pass localizes stack operations and
  pushes child lists directly because the helper returns a parent map rather
  than an ordered traversal.
- `parsimony_utils.fitch_downpass`, `fitch_uppass_acctran`,
  `fitch_uppass_deltran`, and `detect_changes` now use direct standard-tree
  preorder/postorder helpers before falling back to Bio.Phylo traversal. The
  measured rows report the downpass and ACCTRAN phases; DELTRAN and change
  detection share the same direct preorder traversal path. A later downpass pass
  adds a binary-node branch that avoids building per-character child-set lists on
  resolved trees while keeping the generic path for multifurcations. A follow-up
  direct-postorder pass builds postorder lists with reverse preorder instead of
  stack entries carrying visited flags, reducing traversal overhead for
  `fitch_downpass`, `resolve_polytomies`, and other direct-postorder callers.
  The direct postorder helper now also localizes stack and output-list methods
  inside that reverse-preorder loop, preserving identical clade order while
  reducing per-node attribute lookups for all callers.
  `classify_changes` now counts `(character, new_state)` transitions with a
  plain dictionary instead of `Counter`, avoiding the subclass overhead while
  preserving the same lookups for transitions collected from the first pass.
  `retention_index` now counts large ASCII state alphabets with one per-character
  `bincount` histogram instead of repeated full-matrix equality scans, while
  retaining the equality-count path for smaller DNA-sized alphabets where it is
  competitive.
  `resolve_polytomies` now resolves during a direct standard-tree scan and
  imports `Newick` lazily only when it finds a multifurcation; nested
  multifurcation output remains identical to the previous postorder-list path.
  A later bitmask pass caches terminal mask vectors for repeated character
  state patterns, avoiding redundant per-character tip remapping while keeping
  the public node-state return values unchanged.
- `PhyloLogistic._build_logistic_vcv` baseline time copied the full tree for
  each nonzero-alpha VCV evaluation, mutated branch lengths on the copy, and
  then built the VCV. The optimized path computes transformed root-to-node
  depths directly from the original tree for the standard VCV builder, avoiding
  repeated tree serialization while preserving the copy/mutate fallback for
  custom builders.
- `discrete_models.felsenstein_pruning` caches transition matrices by branch
  length during each likelihood evaluation. This avoids repeated matrix
  exponentials on trees with many equal branch lengths while keeping the cache
  scoped to the current Q matrix. `fit_q_matrix` now prepares the invariant
  postorder traversal, child indices, branch lengths, and tip-state indices once
  per model fit, then reuses that context across optimizer likelihood
  evaluations. The prepared-context setup now uses a direct postorder traversal
  for standard Bio.Phylo trees before falling back to the generic traversal for
  nonstandard tree objects. A later likelihood pass precomputes the fixed tip
  likelihood rows and the internal-node work list in that context, so optimizer
  evaluations copy the tip matrix and iterate only internal nodes while retaining
  the same transition matrices and pruning arithmetic. Generic-state root
  likelihood totals now use `np.dot` for the prior-weighted conditional
  likelihood vector, preserving the scalar log-likelihood while avoiding the
  temporary product array reduction. `FitDiscrete.run` now prepares the same
  pruning context once and passes it into each selected model fit instead of
  rebuilding identical tree/state metadata for ER, SYM, and ARD.
  A later run-setup pass reads the cached tree without copying when branch
  lengths are complete, copies before default branch-length assignment, and only
  copies for pruning when the trait file omits tree tips. The complete-branch
  preflight now pushes child lists directly because it only tests whether any
  branch length is missing. Discrete trait parsing now streams both two-column
  and named multi-column TSV files, avoiding
  `readlines()` materialization while preserving comment/blank filtering,
  logical error line numbers, and shared-taxon validation.
  Two-column discrete trait streams now parse valid rows with a tab partition
  and reserve full tab counting for the invalid-column error path, avoiding a
  temporary split list per valid row while preserving logical row numbering and
  column-count validation.
  Multi-column discrete trait parsing now counts tabs to preserve exact row-width
  validation and uses a bounded split to extract the taxon plus selected trait,
  avoiding allocation of unused trailing trait fields when an early column is
  selected.
  Two-state ER transition matrices now use the closed-form CTMC solution inside
  `matrix_exp`, avoiding SciPy matrix exponentials for the common binary ER
  model while leaving larger-state matrices on the existing SciPy path. A later
  binary-transition pass extends that closed-form path to unequal-rate two-state
  CTMCs, avoiding SciPy `expm` for binary ARD/SYM transition matrices while
  preserving the generic fallback for larger-state matrices. A later two-state
  ER fit pass evaluates the optimizer objective
  directly from the scalar ER rate, avoiding per-evaluation `Q` matrix
  allocation and equal-rate matrix-shape detection while preserving the final
  fitted `Q` and log-likelihood. A later pass routes that one-parameter
  two-state ER fit through bounded scalar optimization instead of multi-start
  vector optimizers, preserving the fitted log-likelihood and leaving
  multi-parameter SYM/ARD fits on the existing optimizer path. A follow-up
  traversal pass builds the direct postorder list with reverse preorder instead
  of visited-flag stack entries, preserving postorder while reducing setup
  overhead. Two-column and multi-column discrete trait parsing now stream
  non-comment rows directly from the file instead of first materializing
  `readlines()` and a stripped `data_lines` list, preserving logical data-line
  numbering and shared-taxa validation. A later validation pass returns
  immediately for exact tree/trait taxon matches after at least three taxa are
  shared, avoiding shared/warning set construction while preserving
  too-few-shared-taxa errors.
  overhead for both plain pruning and
  prepared pruning contexts.
  Import-time NumPy startup is deferred behind a lazy proxy, so command modules
  that import the discrete helper do not pay NumPy cost until Q-matrix,
  pruning, or fitting code actually runs.
  A later startup pass replaced eager
  `scipy.linalg.expm` and `scipy.optimize.minimize` imports with same-name lazy
  wrappers, so import-only callers avoid SciPy linalg/optimize startup while
  matrix exponential and fitting paths still call the same SciPy implementations.
  The lazy `expm` wrapper now caches the resolved SciPy callable after first use,
  avoiding repeated import dispatch during pruning and stochastic mapping paths
  with many transition-matrix exponentials, and the lazy `minimize` wrapper does
  the same for repeated optimizer calls.
  The prepared-pruning likelihood now has a scalar two-state path that reuses
  the same closed-form transition probabilities without constructing tiny
  transition matrices or running NumPy matrix-vector products for every internal
  branch; larger-state models continue through the existing generic pruning
  loop. The two-state scalar transition and log-likelihood paths now use
  `math.exp` and `math.log` instead of scalar NumPy ufunc calls, reducing
  dispatch overhead in repeated binary pruning evaluations.
  A later `fit_discrete` startup pass
  keeps the public helper patch points as thin lazy wrappers and stores the
  model-name constant locally, so importing the command module no longer imports
  NumPy or the shared discrete-model helper until traits are parsed or a model is
  fitted. A follow-up startup pass keeps JSON output behind the same lazy-wrapper
  pattern, avoiding the JSON helper during import-only command discovery.
  A later startup pass postpones annotations and removes the annotation-only
  `typing` import, so command discovery no longer loads `typing`.
  Text output now batches the model-comparison header and rows into one
  newline-joined print while preserving exact stdout text.
- `Tree.validate_tree` baseline time materialized terminals and then traversed
  clades again when validating or filling branch lengths. The optimized path
  counts tips and records missing non-root branch lengths in one direct
  standard-tree traversal, preserving fallback behavior for nonstandard tree
  objects and keeping root branch-length handling unchanged. A later pass
  localizes stack operations and pushes child lists directly because validation
  does not expose traversal order, reducing the remaining optimized helper cost
  for both required-length checks and default-length assignment. The generic
  fallback now counts terminals only until `min_tips` is reached, avoiding full
  terminal-list materialization when branch-length fallback traversal is not
  needed. `shared_tips`
  now computes the set intersection once before checking emptiness and returning
  the list, preserving unique-tip semantics while avoiding a duplicate
  intersection for overlapping tree/taxon lists.
- `Tree.prune_tree_using_taxa_list` baseline time pruned named taxa directly,
  making Bio.Phylo resolve each string target during every prune. The optimized
  base helper resolves named terminal clades once and prunes objects directly
  for normal Bio.Phylo trees, retaining the original string-target fallback for
  nonstandard tree doubles and missing taxa. A later setup pass builds that
  name-to-terminal map with direct standard-tree traversal instead of
  `get_terminals()`. The batch pruning path now removes unique resolved
  terminal objects from standard Bio.Phylo trees in one postorder traversal,
  preserving repeated `prune()` output and branch-length collapse semantics;
  duplicate prune requests, all-tip prune requests, missing taxa, and
  nonstandard tree objects keep the original repeated-prune fallback. A later
  pass builds that postorder mutation order by collecting root-right-left
  preorder and reversing it, avoiding visited-flag stack entries while
  preserving branch-length collapse output. A target-setup pass now stores only
  requested terminal objects while still counting all terminals for the all-tip
  guard, avoiding dictionary entries for unrelated tips. A later setup pass
  localizes stack operations and pushes children directly because the helpers
  only expose name-indexed maps, selected targets in caller order, and terminal
  counts; traversal order is not externally visible.
- `ThresholdModel._prune_tree_to_taxa` baseline time rebuilt the full terminal
  list once per taxon being pruned. The optimized path materializes terminals
  once and prunes terminal clade objects directly before VCV and MCMC setup.
  A later pass reuses the shared standard-tree batch prune helper for resolved
  terminal objects, pruning the keep complement in one postorder traversal while
  preserving Bio.Phylo branch-length collapse semantics and the original
  fallback for nonstandard tree objects.
- `TreeSpace._build_distance_matrix` baseline time copied each tree that
  needed pruning, then pruned copied tips by name. The optimized path scans the
  copied tree's terminals once and prunes terminal clade objects directly
  before RF/KF split extraction.
- `TreeSpace._build_distance_matrix` copied-tree batch pruning setup preserves
  the same copied-tree pruning semantics but removes all missing copied tips
  with the shared standard-tree batch prune helper, falling back to per-tip
  pruning for nonstandard tree objects. A later RF-only pass skips copied-tree
  pruning entirely and extracts shared-taxa-filtered RF splits directly from
  each original tree; KF distance keeps the copied-prune path because pruning
  can collapse branch lengths that affect branch-score distances.
- `TreeSpace._extract_splits` and `_extract_splits_with_lengths` now use a
  direct postorder traversal for standard Bio.Phylo trees, falling back to
  `find_clades(order="postorder")` for nonstandard tree objects. This preserves
  the existing canonical split and branch-length overwrite semantics while
  reducing per-tree RF/KF split extraction overhead. A later traversal pass
  builds that direct postorder list by reversing preorder, removing visited-flag
  stack entries while preserving left-to-right postorder. A follow-up child-union
  pass special-cases binary and unary clades before falling back to the generic
  multifurcating union, preserving canonical split and branch-length overwrite
  behavior while shaving RF/KF split extraction time on shared-taxa-filtered
  balanced trees. Equal-size canonical split tiebreaks now compare the minimum
  taxon on each complementary side instead of sorting both sides, preserving the
  sorted-lexicographic choice because the two sides are disjoint. A startup pass
  converts annotation-only `typing` aliases to built-in postponed annotations,
  so command discovery no longer loads `typing`.
- `TreeSpace._write_distance_matrix` now iterates each NumPy row directly when
  formatting CSV output, preserving labels and six-decimal values while avoiding
  repeated two-dimensional indexing in the inner formatting loop.
- `TreeSpace._parse_trees_from_source` now streams path-list rows directly
  instead of building a complete cleaned row list before parsing. Leading
  inline-Newick rows are buffered until the input type is known, preserving the
  previous mixed-input behavior while still avoiding pre-cleaning for path-list
  files.
- `ConsensusNetwork._prune_to_taxa` baseline time collected terminal names,
  then pruned by name for missing-taxa shared mode. The optimized path prunes
  terminal clade objects from the initial terminal pass directly. It now uses
  the shared standard-tree batch pruning path to avoid a repeated path lookup
  for every removed tip while preserving the repeated-prune fallback.
- `QuartetNetwork._prune_to_taxa` uses the same direct terminal-object pruning
  and batch-pruning approach for missing-taxa shared mode, avoiding per-tip
  path lookups before quartet concordance calculations.
- `build_discordance_vcv` baseline time copied each gene tree, then pruned copied
  tips by name. The optimized path scans copied terminals once and prunes terminal
  clade objects directly before VCV averaging. A later pass routes copied
  standard trees through the shared batch-pruning helper, collapsing missing-tip
  removal into one postorder traversal while keeping the object-prune fallback.
  When a gene tree already contains exactly the shared taxa, a follow-up pass
  builds its VCV directly and skips the copy/prune helper entirely.
  The helper itself now also scans the original tree first and returns it
  unchanged when every terminal is shared, while preserving copy-before-prune
  behavior whenever a terminal must be removed.
  Extra-tip standard gene trees now also avoid copy/prune when possible by
  building the retained-subset VCV directly and skipping branches ancestral to
  every retained taxon, matching the root-collapse semantics of pruning. A later
  pass writes retained terminal branches directly to the diagonal instead of
  constructing one-element index arrays and `np.ix_` tuples.
  Follow-up VCV accumulation work replaces multi-tip `np.ix_` block construction
  with broadcast advanced indexing and caches NumPy dtype/conversion lookups in
  the branch loop, preserving identical matrices while reducing per-branch
  indexing overhead in both standard and pruned-subset VCV builders.
  Branch-length validation also walks standard gene trees directly before
  falling back to `find_clades()`. Shared-taxa setup now uses the direct
  terminal-name traversal for species and gene trees. A later shared-setup pass
  collects each gene tree's terminal names and missing-branch-length flag in
  one traversal, preserving the existing shared-taxa error before any cached
  branch-length error is raised. The PSD correction step now computes only the
  smallest eigenvalue first and returns unchanged
  already-PSD matrices without a full eigendecomposition; the full
  eigendecomposition is still used when correction is actually needed. Corrected
  matrix reconstruction now scales eigenvector columns by clipped eigenvalues
  instead of materializing a dense diagonal eigenvalue matrix and reads the
  smallest eigenvalue through the eigenvalue array's own `min()` method. A later
  startup pass replaced the eager `scipy.linalg.eigvalsh` import with a same-name lazy
  wrapper, so import-only callers avoid linalg startup while PSD checks still
  call SciPy's Hermitian eigensolver. A repeated-check pass caches the resolved
  SciPy `eigvalsh` function after first use, keeping the import deferral while
  avoiding repeated import machinery on later PSD checks. Later startup passes
  defer Bio.Phylo and NumPy behind module-level proxies and postpone annotations,
  preserving the existing `Phylo.read` and `np` patch points used by tests.
- `SpectralDiscordance._copy_prune_if_needed` uses the same direct terminal-object
  and batch-pruning approach for shared-taxa bipartition extraction, avoiding repeated
  name-based target lookup after copying. A later pass avoids generic
  `get_terminals()` traversal for standard trees, using direct terminal scans for
  both no-prune detection and copied-tree removal target collection. A follow-up
  startup pass defers pickle until copied-tree pruning is required. The all-shared
  preflight scan now localizes stack operations and pushes child lists directly,
  while copied-tree prune target collection pushes binary children right-then-left
  and indexes multifurcations backward to preserve terminal order without
  `reversed(children)` iterator setup.
- `PrintTree.run` with `--remove` baseline time made the same redundant
  pickle/unpickle copy before clearing branch lengths. The optimized path
  clears branch lengths on the isolated tree returned by `read_tree_file()`
  before drawing or formatting it. Its branch-clearing loop now also walks all
  clades with a direct stack traversal instead of Bio.Phylo's generic
  traversal. A later helper pass clears branch lengths during that traversal
  instead of storing all nodes and looping over them again. A later startup pass
  defers Bio.Phylo behind a lazy
  `Phylo.draw_ascii` proxy while preserving the existing module patch point.
  No-remove `PrintTree.run` now uses the explicit unmodified tree read helper
  to avoid copying the cached parsed tree; `--remove` continues to use the
  defensive-copy reader because it mutates branch lengths before output. A
  follow-up startup pass keeps JSON output behind a module-level forwarding
  wrapper, preserving the patch point while avoiding JSON helper startup during
  command discovery. A later startup pass postpones annotations and removes the
  annotation-only `typing` import, so command discovery no longer loads
  `typing`. The branch-clearing helper now also localizes stack operations in
  the direct traversal, reducing balanced 131072-tip median clear time from
  0.031891s to 0.026640s.
- `ConsensusTree._prune_to_taxa` baseline time collected terminal names once,
  then asked Biopython to find each named tip again while pruning to shared
  taxa. The optimized path prunes terminal clade objects from the initial
  terminal pass directly. A later pass moved consensus tree, consensus network,
  and quartet network pruning onto a shared direct terminal-object traversal for
  standard parsed trees, retaining the `get_terminals()` fallback for
  nonstandard tree-like objects.
- `ConsensusTree._tips` and `ConsensusNetwork._tips` baseline time materialized
  terminal clade objects through Bio.Phylo before creating taxon sets for
  missing-taxa normalization. The optimized helpers reuse the shared direct
  terminal-name traversal and fall back to `get_terminals()` for nonstandard
  tree-like objects. A later startup pass replaces the eager Bio.Phylo and
  Bio.Phylo.Consensus imports with lazy module-level proxies, preserving
  `Phylo.read`, `Consensus.strict_consensus`, and
  `Consensus.majority_consensus` patch points while keeping consensus imports
  light for command discovery and import-only callers. A follow-up startup pass
  keeps JSON output behind a module-level forwarding wrapper, preserving the
  patch point while avoiding JSON helper startup during command discovery. A
  later startup pass converts annotation-only `typing` aliases to built-in
  postponed annotations, so consensus tree command discovery no longer imports
  `typing`. Consensus tree, consensus network, and quartet network source
  cleanup now streams over the input file and strips/filters comments in one
  pass instead of materializing `splitlines()`, preserving multi-Newick and
  path-list behavior. The path-list branch now resolves relative paths with a
  precomputed parent string and bound `os.path` helpers instead of per-row
  `Path` construction, preserving normal missing-file errors while reducing
  large tree-list setup overhead. Identical-taxon normalization in consensus
  tree, consensus network, and quartet network now compares the extracted taxon
  sets before building the shared intersection; the shared set is only needed
  for pruning modes after a mismatch is found. Consensus tree, consensus
  network, and quartet network now perform that identical-set check with an
  iterator loop instead of allocating `tip_sets[1:]`, preserving early mismatch
  exits without copying the remaining taxon-set list.
- `ConsensusNetwork._count_splits` baseline time extracted each tree's splits
  by calling `clade.get_terminals()` for every internal clade. The optimized
  split extractor computes descendant tip sets once in postorder and applies the
  existing polytomy, trivial-split, and canonicalization rules from cached
  clade tip sets. Allow-mode taxon setup now also uses the shared direct
  terminal-name traversal before extracting splits against each tree's own taxon
  set. A later shared-taxa pass counts integer split masks internally and
  converts each unique mask back to the public `frozenset` split key only once,
  preserving rooted split counts and canonical equal-size split behavior. A
  later traversal pass builds the direct postorder mask extractor by reversing a
  root-right-left preorder list and specializes binary child mask unions,
  preserving the same public split keys while avoiding visited-tuple stack
  entries. A follow-up split-set counting pass batches extracted split sets
  through `Counter.update()` and builds the uniform possible counter from
  `dict.fromkeys`, preserving the same returned `Counter` objects while
  avoiding per-split Python increment loops in allow-mode and fallback paths.
  Equal-size frozenset canonicalization now compares the minimum taxon
  on each complementary side instead of sorting both sides, preserving the
  sorted-lexicographic choice because the two sides are disjoint. Equal-size
  mask canonicalization now compares the lowest set bit on each disjoint half
  instead of materializing both sorted name lists, preserving the same
  lexicographic side choice because mask bit order follows sorted taxon names.
- `ConsensusNetwork._build_splits_graph` baseline time enumerated every
  `2 ** n_splits` sign vector and then compared every valid-node pair to find
  graph edges. The optimized path reuses the existing forbidden-pair rules while
  assigning signs incrementally, pruning invalid partial assignments early, and
  finding edges by flipping one split at a time. A later edge-heavy pass keeps
  the recursive pruning behavior but tracks valid sign vectors as integer masks
  for edge discovery, avoiding tuple slicing for every node/split flip while
  returning the same public tuple sign vectors. Circular-order setup now also
  reads the majority-consensus terminal order through the shared direct
  terminal-name traversal instead of materializing terminal clade objects.
  Split-direction setup now caches each taxon's circular cos/sin coordinates
  once per ordering instead of recomputing them for every accepted split
  membership centroid.
  Network plotting now batches internal split-graph edges and pendant taxon
  edges into `LineCollection`s while preserving per-pendant linewidth and alpha.
  Pendant-edge scaling now computes the split-graph x/y extent in one pass
  without temporary coordinate lists.
  A later fallback-rendering pass batches unlabeled no-split taxon points into
  one `scatter` collection instead of one marker artist per taxon. Text output
  now batches the summary header and filtered split rows into one
  newline-joined print while preserving exact stdout text.
  Startup passes keep Bio.Phylo consensus helpers behind lazy proxies, keep
  JSON output behind a forwarding wrapper, and localize `PlotConfig` to argument
  processing so command discovery avoids parser, JSON, and plotting startup. A
  later startup pass converts annotation-only `typing` aliases to built-in
  postponed annotations so command discovery no longer loads `typing`.
- `LBScore.calculate_lb_score` builds the cached tree-distance table once and
  reuses it for both global average tip distance and per-taxon average distance.
  The optimized path now accumulates per-tip distance sums directly from the
  cached pair list, avoiding a frozenset-keyed lookup table and nested per-tip
  pair scans while preserving the existing `set(tip)` compatibility behavior
  in per-taxon distance denominators. A later standard-tree path computes both
  the total pairwise distance and each tip's sum of distances to all other tips
  with one postorder traversal and one reroot traversal, falling back to the
  pairwise cache for nonstandard tree objects or duplicate tip names. Verbose
  text output now batches taxon score rows into one newline-joined print while
  preserving empty-output behavior and stdout text. The optional tqdm progress
  import is now deferred until the large multiprocessing distance path needs it;
  per-taxon score calculation now uses equivalent Python list arithmetic,
  avoiding NumPy startup on normal LB-score runs. A later denominator pass keeps
  the historical `set(tip)` compatibility rule but computes the denominator by
  counting matching tip-name characters directly, avoiding a full `tip_set`
  copy for every tip in the linear helper and pairwise-cache path. The fallback
  per-taxon distance path now also builds the tip set once before constructing
  legacy `set(tip)` other-taxon lists. A later
  startup pass postpones annotations and
  removes the annotation-only `Bio.Phylo.Newick` import, leaving tree parser
  startup to the shared tree reader instead of module import. A subsequent
  startup pass preserves the module-level `ProcessPoolExecutor`, `as_completed`,
  and `mp.cpu_count` patch points through lazy proxies, so importing the command
  no longer initializes multiprocessing or concurrent futures. The latest
  startup pass keeps the module-level `pickle`, summary-statistics, and JSON
  helper patch points as lazy forwarding objects/functions, so import-only
  callers no longer initialize those helpers. A later startup pass converts
  annotation-only `typing` aliases to built-in postponed annotations, so command
  discovery no longer imports `typing`. Cached read-only `LBScore.run` now uses
  the explicit unmodified tree read helper to avoid copying the cached parsed
  tree before score calculation and output dispatch. The linear component
  helper's postorder setup now pushes binary children right-then-left and
  indexes multifurcations backward, preserving traversal order while reducing
  balanced 32768-tip median helper time from 0.159546s to 0.132146s.
- `Saturation.loop_through_combos_and_calculate_pds_and_pis` reuses the cached
  tree-distance helper for patristic distances while preserving combo order and
  the existing uncorrected-distance calculations. When cached pairwise
  distances are returned in the requested combo order, the optimized loop now
  consumes that distance list directly instead of building and indexing a
  `frozenset` dictionary; arbitrary combo orders keep the dictionary fallback.
  No-gap uncorrected-distance matrix calculation now uses direct blockwise
  sequence comparisons instead of constructing an all-true validity matrix and
  accumulating per-symbol match products; gap-excluding runs keep the
  validity-mask path. A later guard lets `exclude_gaps=True` clean alignments
  with no observed gap positions reuse that no-gap matrix path after the usual
  sequence validation, preserving the gappy validity-mask path as soon as any
  selected taxon has a gap. The gappy matrix path now stacks raw gap masks and
  inverts the single stacked matrix, avoiding one temporary inverted boolean
  array per taxon. Standard upper-triangle pair requests now extract
  matrix distances row-wise, avoiding per-pair taxon-index dictionary lookups
  while preserving the fallback for custom pair orders. When cached tree distances
  are available and every
  requested taxon has the same normalized sequence, uncorrected distances now
  return as a constant vector before sequence-array and matrix construction;
  gap-excluding all-ambiguous identical sequences still return `NaN` distances.
  That identical-sequence guard now uses a local index-based scanner instead of
  allocating `sequences[1:]`, retaining early mismatch exits without list-copy
  overhead. The identical-sequence Unicode valid-site fallback now counts gap
  characters with `str.count` instead of looping over every character in Python.
  Raw-identical requested tips now compare raw sequence strings before falling
  back to uppercase comparisons, preserving case-insensitive matching while
  avoiding repeated uppercase allocations for already-normalized conserved rows.
  Cached read-only `Saturation.run` now uses the explicit unmodified tree read
  helper to avoid copying the cached parsed tree before tip-pair setup and
  pairwise distance calculation. Verbose text output now batches pairwise rows
  into one newline-joined print while preserving the same stdout text.
- The `saturation` command now keeps the module-level `mp.cpu_count` and
  `mp.Pool` access points through a lazy proxy and imports `functools.partial`
  only when the large-workload multiprocessing branch is used, so command
  startup no longer initializes `multiprocessing`. A later startup pass
  localizes `PlotConfig` to CLI argument processing, keeping plotting helper
  setup out of import-only command discovery. A follow-up startup pass keeps JSON
  output behind a lazy forwarding wrapper, so plain imports avoid the shared JSON
  helper. Another startup pass removes the runtime `TYPE_CHECKING` dependency and
  converts annotation-only typing aliases to postponed built-in annotations. The
  alignment reader now uses a same-name forwarding wrapper, preserving the
  module-level patch point while avoiding `phykit.helpers.files` during
  import-only command discovery.
- `Saturation.loop_through_combos_and_calculate_pds_and_pis` byte sequence array
  setup stores normal ASCII alignments as byte arrays instead of Unicode scalar
  arrays, preserving uppercasing and gap masking while retaining a Unicode
  fallback for non-ASCII sequence content. A later pass stores those byte
  arrays as `uint8` values and counts valid matches with boolean masks instead
  of slicing both sequences for every taxon pair. Mock and nonstandard tree
  objects still fall back when cached pairwise tree-distance setup is
  unavailable. A subsequent path handles the common ordered-combo case with
  one set of all-pairs matrix products over valid sequence masks, guarded by
  a maximum taxa count and falling back to the pair loop for non-ASCII,
  uneven-length, or larger inputs. Gap masks for byte-array alignments now use
  a cached 256-entry lookup table instead of rebuilding `np.isin` membership
  arrays for every sequence; non-byte arrays keep the existing fallback. A
  later startup pass defers NumPy and annotation-only Biopython imports while
  preserving the module-level `np` patch point used by tests.
- `vcv_utils.build_vcv_matrix` baseline time was dominated by repeated
  root-to-tip and pairwise `tree.distance` calls. The optimized path caches
  root depths and ancestor paths once, subtracts any explicit root branch
  offset to match Biopython `tree.distance(tree.root, tip)` semantics, and
  fills shared-path entries from cached MRCA depths.
- The second-stage `Tree.calculate_pairwise_tip_distances_fast`,
  `vcv_utils.build_vcv_matrix`, and `vcv_utils.build_transformed_vcv_matrix`
  path setup benchmarks measured the remaining per-tip `tree.get_path()`
  overhead on 1024-tip trees. The optimized paths build one object-keyed
  parent map and derive root-to-tip paths by walking parent links, retaining
  fallback behavior when parent links cannot be resolved. The shared parent-map
  helper now builds standard Bio.Phylo parent maps with direct traversal before
  falling back to `find_clades()`. A later startup pass defers pickle until tree
  transformation or pruning helpers need a protective copy. A follow-up startup
  pass converts annotation-only `typing` names in `vcv_utils` and
  `tree_paths` to built-in postponed annotations, and localizes `Path` and
  `StringIO` to gene-tree parsing. Shared gene-tree source cleanup now streams
  over the input file and strips/filters comments in one pass instead of
  materializing `splitlines()`, preserving inline-Newick and path-list behavior.
  The shared path-list branch now resolves relative rows with a precomputed
  parent string prefix and bound absolute-path check, preserving absolute-row
  behavior while avoiding per-row `Path` joins.
  Gene-tree branch-length validation now pushes child lists directly because it
  only returns whether a missing branch length exists. The combined terminal-name
  and branch-length scan preserves tree tip order while pushing binary children
  explicitly and indexing multifurcations backward. The all-shared prune
  preflight uses the same unordered direct scan before returning the original
  gene tree; pruning paths still keep their existing copied-tree collection order.
- A later `Tree.calculate_pairwise_tip_distances_fast` pass replaces
  Bio.Phylo `depths()`/`get_terminals()` setup and per-pair ancestor-set scans
  with one direct depth/parent traversal plus root-path prefix comparison for
  the selected tips. A follow-up setup pass localizes stack operations and
  pushes binary children explicitly while indexing multifurcations backward,
  preserving caller-ordered pair output and distances while avoiding the
  per-node `reversed(children)` iterator on common bifurcating trees.
- The branch-accumulation VCV pass removes the remaining pairwise MRCA scan for
  standard and transformed VCV matrices. For every ordered tip, it walks the
  root-to-tip path once and records descendant-tip indices per branch; each
  branch then adds its length to the descendant-tip covariance submatrix in one
  NumPy update. This preserves explicit-root-branch semantics by skipping the
  root itself, matching `tree.distance(tree.root, tip)`. A later startup pass
  replaces the eager Bio.Phylo import with a lazy `Phylo.read` proxy, preserving
  parse-gene-tree behavior while avoiding Biopython tree parser startup for
  import-only callers.
- `LTT._compute_gamma` and `LTT._compute_ltt` baseline time was dominated by
  repeated root-to-tip `tree.distance` calls for terminal and internal clades.
  The optimized paths compute tree depths once per calculation and reuse
  root-relative depths while preserving fallback behavior for nonstandard tree
  objects. Tip-height normalization now streams the maximum terminal height from
  the shared terminal/depth context without rebuilding a generator expression at
  each gamma/LTT call site.
- `Chronogram._compute_root_to_tip` baseline time called `tree.get_path()` for
  every node before summing branch lengths. The optimized path walks the tree
  once, carrying each parent's root-relative distance to its children.
- `Chronogram._print_json` baseline time called `clade.get_terminals()` for
  every internal node and `tree.count_terminals()` for the tip count. The
  optimized path builds sorted descendant tip-name tuples once in postorder and
  reuses them for all node-age payload entries. A later cache-construction pass
  keeps exact Biopython postorder order while building it as root-right-left
  preorder followed by an in-place reverse, avoiding per-node visited stack
  tuples in the direct standard-tree helper.
- `Chronogram` plot setup baseline time materialized terminal clades and used
  separate preorder/postorder scans for y coordinates, branch drawing, HPD bars,
  and node-age labels. The optimized rectangular/circular plot paths reuse one
  direct preorder list, the precomputed root-to-node distances, and the shared
  direct coordinate helper for y positions. A later setup pass passes that
  preorder list into `compute_node_positions` for rectangular y coordinates and
  passes the same preorder list plus prepared tips into `compute_circular_coords`
  for circular coordinates, avoiding additional direct tree walks while
  preserving tip order. A later traversal-helper pass keeps preorder output
  identical while using a one-/two-child fast path instead of
  `reversed(children)` for the generator used by rectangular and circular setup.
  The list-returning direct preorder helper now uses the same child-push pattern
  for plot overlay setup. A later rectangular plotting pass
  batches horizontal and vertical base branches into two `LineCollection`s,
  preserving branch styling while avoiding one Matplotlib `Line2D` artist per
  segment. Rectangular HPD bars now render through one translucent
  `PatchCollection`, preserving bar dimensions and styling while avoiding one
  patch artist per HPD node. Rectangular clade-color overlays now batch highlighted
  horizontal and vertical branch segments into two `LineCollection`s per clade,
  preserving color-file styling while avoiding one line artist per overlay segment.
  A later circular plotting pass batches radial branch segments and internal arcs
  into two `LineCollection`s while leaving timescale guide lines on their existing
  drawing path. A later HPD rendering pass batches circular radial confidence
  intervals into one translucent `LineCollection`, preserving the round-capped blue
  styling while avoiding one line artist per HPD node. A later startup pass defers
  NumPy behind a lazy proxy, leaving array startup until
  circular plotting first computes angle and radius guides. The run path now
  reuses the cached parsed tree for read-only chronogram setup; optional ladderization copies before
  mutating branch order. A follow-up startup pass keeps JSON, plot, timescale,
  circular-layout, and color-annotation helpers behind forwarding wrappers or
  localized imports so command discovery avoids those helper modules. A later
  startup pass also removes annotation-only `typing` imports from `Chronogram`
  and the shared tree base loaded by it. A later
  HPD parsing pass reuses cached compiled annotation regexes and walks standard
  trees with direct preorder traversal, retaining the generic `find_clades`
  fallback for nonstandard tree objects and ignoring clades without comments.
- `TotalTreeLength.run` and `EvolutionaryRate.run` baseline time used
  Bio.Phylo's generic total-branch-length traversal; evolutionary rate then
  traversed again to count terminals. The optimized shared helpers use direct
  stack traversals over parsed Biopython trees, with fallback to the original
  tree-like APIs for test doubles and nonstandard objects. The
  total-tree-length command module now also defers annotation-only Bio.Phylo
  imports and inherits the shared JSON helper's lazy NumPy handling at startup.
  A follow-up startup pass keeps JSON output behind a module-level forwarding
  wrapper, preserving the existing output path while avoiding JSON helper
  startup during command discovery. The evolutionary-rate command module now
  uses the same module-level forwarding wrapper for JSON output, preserving the
  patch point while avoiding JSON helper startup during command discovery.
  For cached read-only total-tree-length and evolutionary-rate runs, these
  commands now use an explicit unmodified tree read helper to avoid copying the
  cached parsed tree while leaving the defensive-copy default in place for
  mutation-prone services. A later helper pass localizes the stack pop/extend
  methods in the total-branch-length traversal, preserving fallback behavior
  while shaving repeated attribute lookups from shared total-length callers. A
  later total-tree-length startup pass removes annotation-only `typing` imports
  so command discovery no longer loads `typing`. A later evolutionary-rate
  startup pass applies the same annotation-only `typing` cleanup to that
  command module. The shared terminal-root-distance helper now pushes binary
  children explicitly and indexes larger child lists backward, preserving
  terminal distance order while avoiding the per-node `reversed(children)`
  iterator on common bifurcating trees.
- `Tree.calculate_treeness` baseline time walked nonterminals, then traversed
  the full tree again through Bio.Phylo's total-branch-length API. The
  optimized helper computes internal and total branch lengths in one direct
  stack traversal, preserving fallback behavior for tree-like test doubles.
  Cached read-only `Treeness.run` now also uses the explicit unmodified tree
  read helper to avoid copying the cached parsed tree before calculating the
  scalar. A follow-up startup pass keeps JSON output behind a module-level
  forwarding wrapper, preserving the patch point while avoiding JSON helper
  startup during command discovery. A later startup pass removes the
  annotation-only `typing` import, so command discovery no longer loads
  `typing`.
- `build_vcv_matrix` baseline time built a parent map and walked from every
  ordered terminal back to the root to collect clade membership. The optimized
  standard-tree path computes descendant ordered-tip indices in one postorder
  pass, adds each branch to the corresponding VCV block, and falls back to the
  parent-map route for duplicate or nonstandard trees. A later pass updates
  one-tip terminal branches with direct diagonal writes instead of constructing
  one-element index arrays and `np.ix_` tuples for every terminal branch.
- `DVMC.determine_dvmc` baseline time was dominated by repeated terminal
  `tree.distance` calls. The optimized path computes depths once and reuses
  root-relative terminal depths for the same variance calculation.
- `PhyloHeatmap.run` initial taxon setup and post-prune tip-order setup
  baseline time materialized terminal clade objects even though only names were
  needed. The optimized path uses the shared direct terminal-name traversal and
  retains fallback behavior through `get_tip_names_from_tree`. A later
  cached-tree setup pass reads the parsed tree without copying, copies before
  default branch-length assignment, pruning, or ladderizing, and reuses the
  initial tip order when no topology mutation occurs. Matrix parsing now streams
  the first non-comment header and data rows instead of materializing both
  `readlines()` and a filtered `data_lines` list, preserving logical data-row
  numbering in validation messages while lowering peak memory. A later parser
  pass returns immediately for exact tree/data taxon matches after all rows are
  validated and at least three taxa are shared, avoiding shared/warning set
  construction while preserving too-few-shared-taxa errors. The branch-length
  preflight now pushes child lists directly because it only tests whether any
  branch length is missing.
- `PhyloHeatmap._plot_phylo_heatmap_circular` baseline setup materialized
  terminal clades to map names to node ids and used extra preorder scans for
  clade-color overlays. The optimized path reuses a direct preorder list for
  terminal lookup and overlay branch iteration, while keeping circular layout
  and heatmap ring geometry unchanged. A later setup pass feeds the same
  preorder list into `compute_node_positions` for rectangular and circular tree
  panels and passes the preorder plus prepared tip list into
  `compute_circular_coords`, avoiding additional direct tree walks while
  preserving tip order. A later traversal-helper pass keeps preorder output
  identical while using a one-/two-child fast path instead of
  `reversed(children)` for direct preorder materialization. A later
  color-overlay rendering pass
  batches rectangular clade horizontal/vertical branches and circular clade
  radial branches into `LineCollection`s, preserving color-file clade styling
  while avoiding one `Line2D` artist per highlighted branch segment. A later
  circular heatmap rendering pass batches ring-cell `Wedge` patches into one
  `PatchCollection`, preserving per-cell colormap colors and no-edge styling
  while avoiding one Matplotlib patch artist per tip-trait cell. A startup
  pass defers NumPy behind a lazy proxy, so import-only command discovery avoids
  array startup until matrix standardization, clustering, or heatmap plotting
  runs. A follow-up startup pass keeps JSON, plot-config, plotting, and
  color-annotation helpers behind forwarding wrappers/local imports so plain
  imports avoid those helper modules. A later startup pass converts
  annotation-only `typing` aliases to built-in postponed annotations, so command
  discovery no longer imports `typing`. Heatmap matrix construction now builds
  the ordered row list through typed `np.asarray`, preserving float matrix
  output while avoiding the slower inline `np.array` conversion on large
  parser-shaped inputs.
- `Cophylo._rotate_tree` baseline time recomputed descendant terminal lists
  for every internal child while rotating nodes. The optimized path computes
  mapped descendant target positions once in postorder and reuses each clade's
  mean target position for the same swap decisions.
- `Cophylo._validate_tree` baseline time materialized terminal clade objects,
  then traversed all clades again to fill missing branch lengths. The optimized
  path counts tips and fills missing non-root branch lengths in one direct
  standard-tree traversal with the generic traversal retained as fallback. A
  later helper pass pushes child lists directly because validation only returns a
  tip count and fills branch lengths, without exposing traversal order.
- `Cophylo._rotate_tree` aggregate position means avoid storing full descendant
  target-position lists for every clade after the postorder optimization. The
  rotation decision now carries only the sum and count needed for each child
  mean, preserving the same terminal-scan reference behavior with less memory
  churn. A later pass builds one direct postorder clade list for standard
  Bio.Phylo trees and reuses it for both rotation phases, retaining
  `find_clades(order="postorder")` as the fallback for nonstandard tree-like
  objects. A follow-up helper pass collects root-right-left preorder and
  reverses it to preserve Bio.Phylo left-to-right postorder without visited-flag
  stack entries. The postorder helper now also localizes stack pop, extend, and
  clade-append operations inside the traversal loop. A later pass rotates binary
  children while computing each clade's aggregate target-position sum/count,
  removing the second scan and the mean-position cache while preserving
  unmapped-tip handling. Because postorder guarantees child aggregate entries
  are already populated, the same loop now uses direct dict lookups for child
  sums/counts instead of defaulted `.get()` calls.
- `Cophylo._draw_phylogram` rectangular plotting setup now reuses direct
  standard-tree terminal, preorder, and parent-map traversals when assigning
  node coordinates and drawing branches/labels, with generic Bio.Phylo
  traversal retained as fallback for nonstandard tree-like objects. A later
  parent-map pass localizes stack operations and pushes child lists directly
  because the helper returns a map rather than an ordered traversal. A later
  Matplotlib pass batches vertical and horizontal base branch segments into two
  `LineCollection`s on real axes, while preserving the per-plot fallback for
  lightweight axes. A later rectangular plotting pass batches middle-panel
  association connectors into one `LineCollection`, preserving their gray alpha
  styling while avoiding one Matplotlib artist per mapped taxon. Text summary
  output now batches the five summary lines into one newline-joined print while
  preserving exact stdout text. A later
  color-overlay rendering pass batches rectangular
  clade horizontal/vertical branches and circular clade radial branches into
  `LineCollection`s, preserving color-file styling while avoiding one artist per
  highlighted branch segment on each tree. A later circular setup pass reuses
  the preorder and terminal lists built for each tree when computing circular
  coordinates, clade overlays, and matched-tip id lookups, avoiding additional
  direct tree walks while preserving tip order. The shared preorder helper now
  pushes child lists directly instead of allocating `reversed(children)` iterator
  objects, preserving the same Bio.Phylo preorder for binary and multifurcating
  trees. The same order-preserving push pattern now applies to the terminal-clade
  helper used by rectangular drawing and circular matched-tip setup. A startup pass removes the only
  direct NumPy use by computing internal y-position means with Python arithmetic,
  so command import avoids NumPy entirely. The run path now routes tree2 through
  the shared cached tree reader while still receiving a mutable copy for
  validation, ladderization, and rotation. Later startup work keeps JSON output
  behind a forwarding wrapper and localizes `PlotConfig` plus
  cladogram layout helper imports to argument processing and plotting paths.
  A follow-up import pass removes the remaining annotation-only `typing` import
  with postponed built-in annotations. Mapping-file parsing now streams rows
  directly from the file handle instead of materializing `readlines()`, preserving
  column-count validation while reducing peak memory for large tanglegram maps.
- `SimmapSummary._summarize_per_branch` single-segment branch histories are now
  accumulated directly as full-branch dwelling time, while multi-segment
  histories compute dwelling times and transition counts in one pass. This
  avoids the generic segment loop and separate transition scan for branches
  without within-branch transitions. Standard parsed trees now collect branch
  clades with the class's direct preorder stack traversal instead of
  `find_clades(order="preorder")`, retaining the generic traversal fallback for
  nonstandard tree objects. A follow-up setup pass stores each branch's length,
  dwelling accumulator, and transition count in one branch-data table instead of
  building three separate dictionaries over the same clade list.
- `SimmapSummary._summarize_node_posteriors` now collects internal node ids with
  a direct standard-tree stack traversal before accumulating per-mapping state
  counts, retaining generic traversal fallback for nonstandard tree objects. A
  later helper pass preserves preorder while avoiding `reversed(children)`
  iterator allocation and reuses the posterior item table during per-mapping
  accumulation. A
  startup pass defers direct NumPy imports behind a lazy proxy, so command
  discovery avoids array startup until SIMMAP fitting, summarization, or plotting
  runs. Text and JSON output now build branch and node payload rows from one
  direct preorder list instead of running separate generic preorder traversals for
  branch summaries and node posteriors. CSV output now uses the same direct
  preorder reuse, collecting node posterior rows while streaming branch rows so it
  avoids a second generic tree traversal while preserving row order and text. A
  later startup pass converts annotation-only `typing` aliases in SIMMAP summary
  and its stochastic-map parent to built-in postponed annotations so command
  discovery no longer loads `typing`. A text-output pass precomputes percent
  row format strings for Q-matrix, transition, branch, and node-posterior rows,
  preserving byte-identical summary text while avoiding per-cell f-string joins
  in large branch/node tables.
- `StochasticCharacterMap._build_parent_map` baseline time used Bio.Phylo's
  generic `find_clades()` traversal to build id-to-parent mappings for SIMMAP
  setup. The optimized helper walks standard trees directly with a stack and
  retains the generic traversal as a fallback for nonstandard tree objects. A
  later helper pass localizes stack operations and pushes children directly
  because the result is a parent map rather than an ordered traversal.
- `StochasticCharacterMap._summarize_simulations` direct setup traversal keeps
  the prior direct accumulation strategy but builds branch-length and internal
  node metadata with a standard-tree stack walk before falling back to
  Bio.Phylo traversal for nonstandard tree objects. A later pass adds a
  single-segment branch-history fast path, accumulating full branch dwelling
  time immediately when no within-branch transition occurred and using one
  combined dwelling/transition pass for multi-segment histories.
  Text summary output now batches the fitted-Q table, dwelling rows, and
  transition rows into one newline-joined print while preserving exact stdout
  text. A later text-output pass precomputes the fitted-Q percent row format,
  preserving byte-identical dense matrix output while avoiding per-cell f-string
  joins for large state tables.
- `StochasticCharacterMap._build_simulation_metadata` baseline time used
  Bio.Phylo's generic preorder traversal and per-node terminal checks. The
  optimized helper reuses the same direct stack traversal as the parent-map
  setup, preserving preorder metadata and retaining the generic traversal as a
  fallback for nonstandard tree objects.
- `StochasticCharacterMap._run_single_simulation` now samples tiny categorical
  distributions with one uniform draw plus a cumulative-probability search
  instead of `Generator.choice(..., p=...)`. This matches NumPy's fixed-seed
  choice stream for the tested probability vectors while avoiding per-draw
  `choice` overhead in ancestral-state and branch-history simulation.
  Ancestral-state sampling now precomputes each internal node's conditional
  probability row for every possible parent state from the fixed Q matrix,
  branch length, and child likelihood vector, preserving random draw order while
  avoiding repeated transition lookup, likelihood multiplication, and
  normalization across simulations. A later pass stores cumulative probability
  rows for those precomputed ancestral-state distributions and for branch-state
  transition distributions, preserving the same one-uniform draw stream and
  `searchsorted(..., side="right")` thresholding while avoiding per-draw
  cumulative-sum construction. The CDF draw helper now uses direct comparisons
  for 2-, 3-, and 4-state distributions, preserving `side="right"` boundary
  behavior while avoiding `np.searchsorted` overhead in the common small-state
  simulation paths. A later pass extends the same direct-comparison path through
  8-state CDFs for larger discrete-state models. A follow-up pass orders the
  5- through 8-state comparisons as small binary searches, reducing the average
  number of scalar threshold checks while preserving the same boundary behavior.
  A subsequent pass prepares the
  branch-history rates and transition CDFs once for the fixed fitted `Q` matrix
  and reuses that context across simulations, avoiding repeated context
  construction inside `_run_single_simulation`. A later two-state branch-history
  pass stores deterministic next-state targets in that shared context; the
  simulator still consumes the same uniform draw as the CDF sampler but skips
  generic CDF thresholding when only one outgoing transition is possible.
  A startup pass defers direct NumPy imports behind a lazy proxy, so command import
  avoids array startup until Q-matrix fitting, ancestral sampling, or mapping
  simulation runs. A later startup pass converts annotation-only `typing` aliases
  to built-in postponed annotations, so command discovery no longer loads
  `typing`.
- `StochasticCharacterMap._plot_stochastic_map` layout setup baseline time
  rebuilt terminal order, rectangular coordinates, and branch drawing traversal
  with `get_terminals()` plus separate preorder/postorder `find_clades()` scans.
  The optimized path uses the shared direct `compute_node_positions` helper and
  one direct preorder list. A later pass also passes that already-built
  preorder list into `compute_node_positions`, avoiding another direct tree
  walk during stochastic-map coordinate setup. A later traversal-helper pass
  keeps preorder output identical while using a one-/two-child fast path instead
  of `reversed(children)` for direct preorder materialization. The circular plot path now
  passes its existing preorder and terminal lists into `compute_circular_coords`
  as well, preserving the previous left-to-right tip angle order while avoiding
  another circular-layout tree traversal. A later branch-rendering pass batches
  rectangular vertical connectors and mapped history segments into
  `LineCollection`s, and batches circular radial history segments plus internal
  state-colored arcs into collections while preserving branch-history colors
  and segment boundaries.
	  A later import-time pass benefits from the shared `discrete_models` lazy SciPy
	  wrappers, avoiding linalg/optimize startup until Q-matrix fitting or
	  transition matrix evaluation actually runs. A follow-up startup pass converts
	  annotation-only `typing` aliases to built-in postponed annotations, keeping
	  import-only callers free of `typing`. A subsequent startup pass localizes JSON
	  output, `PlotConfig`, discrete-model wrappers, circular-layout drawing helpers,
	  and color annotations to the methods that need them, leaving import-only
	  callers free of those helper modules.
- `Phenogram._fast_anc` and `ContMap._fast_anc` baseline time used repeated
  Bio.Phylo `find_clades()` traversals for node labeling, parent-map
  construction, ML passes, sigma2 contrasts, and result construction. The
  optimized path uses direct stack traversals and reuses preorder/postorder
  clade lists inside `_fast_anc`. A later pass keeps the same direct traversal
  but caches clade ids within the hot loops, switches repeated membership
  checks to single `dict.get()` lookups, and builds internal-node output from
  child lists instead of calling `is_terminal()` for every clade. Later
  fallback passes build `Phenogram` and `ContMap` `_iter_postorder` helpers by
  reversing root-right-left preorder, preserving Biopython postorder while
  reducing helper-only sigma2 contrast setup when callers do not supply a
  clade list. Their standalone parent-map helper now walks child lists directly
  because it returns only a map and does not expose traversal order.
- `Phenogram.run`, `ContMap.run`, and `TraitRateMap.run` copied-tree pruning
  setup baseline time materialized terminal clade objects where only tip names
  were needed. The optimized path uses the shared direct terminal-name
  traversal for copied trees before pruning. Cached read-only `Phenogram.run`
  now uses the unmodified tree reader before creating its own protective copy
  for pruning, removing one redundant cached parsed-tree copy while preserving
  prune isolation. Cached read-only `ContMap.run` applies the same reader
  switch while keeping its copied-tree prune and optional ladderize operations
  isolated from the cached tree. A later Phenogram run-setup pass skips that
  protective copy entirely when all tree tips have trait values, copying only
  before pruning missing trait taxa. ContMap now uses the same read-only
  all-shared path when ladderize is disabled, and copies only before pruning or
  ladderizing would mutate the cached tree. Their single-trait parsers now
  stream input rows, validate the two-column shape with `partition("\t")`, and
  avoid `readlines()` plus temporary split lists on valid rows. A later parser
  pass returns immediately for exact tree/trait taxon matches after all rows
  are validated and at least three taxa are shared, avoiding shared/warning set
  construction while preserving the minimum-shared-taxa error.
- `Phenogram._plot_phenogram` layout and branch setup baseline time rebuilt
  all-node estimates, terminal order, root-distance coordinates, and drawable
  branch iteration with `get_terminals()` plus three preorder `find_clades()`
  scans. The optimized path computes estimates, tips, and x-coordinates in one
  direct preorder pass and reuses that list for branch drawing. Plot color
  normalization now computes the estimate value range in one pass without
  materializing a temporary values list. A later Matplotlib pass batches the 50
  gradient subsegments for every branch into one
  `LineCollection`, preserving the sampled trait-gradient colors while avoiding
  one `Line2D` artist per subsegment. A later startup pass postpones annotations
  and defers NumPy behind a module-level proxy. A subsequent startup pass
	  localizes the protective-copy `pickle` import, keeps JSON output behind a
	  local forwarding wrapper, and localizes
	  `PlotConfig` to argument processing. A later startup pass converts
	  annotation-only `typing` aliases to built-in postponed annotations so command
	  discovery no longer loads `typing`. Text summary output now batches the
	  four-line report into one newline-joined print while preserving exact stdout
	  text.
- `ContMap._plot_contmap` layout setup baseline time rebuilt tip order,
  x-coordinates, internal y-coordinates, and branch drawing traversal with
  `get_terminals()` and separate `find_clades()` scans. The optimized path uses
  the shared direct `compute_node_positions` helper and a single direct preorder
  list for estimate collection plus rectangular and circular drawing loops. A
  later setup pass passes the `_prepare_contmap_plot_data` preorder list into
  `compute_node_positions` and passes that same preorder list plus the prepared
  tip list into `compute_circular_coords`, avoiding two additional direct tree
  walks while preserving the existing tip-angle order. A later traversal-helper
  pass keeps preorder output identical while using a one-/two-child fast path
  instead of `reversed(children)` for direct preorder materialization. Plot
  color normalization now computes the estimate value range in one pass without
  materializing a temporary values list. A later rectangular plotting pass
  batches the 50 per-branch gradient segments
  into one `LineCollection` and the vertical connectors into a second
  `LineCollection`, preserving the sampled gradient colors while avoiding one
  Matplotlib `Line2D` artist per tiny segment. Circular contMap rendering now
  uses the shared whole-tree radial gradient helper, avoiding one
  `LineCollection` per branch, and the shared colored-arc helper, avoiding one
  `Line2D` artist per internal arc. Text summary output now batches the
  four-line report into one newline-joined print while preserving exact stdout
  text.
- `cont_map` module import baseline paid NumPy startup directly and through the
  shared circular-layout helper's module-level arc array. The optimized path
  defers NumPy behind lazy proxies and materializes the circular arc fractions
  array only when circular drawing actually runs.
- Later `cont_map` startup work localizes the protective-copy `pickle` import,
  `PlotConfig`, layout helpers, circular/color helpers, and keeps JSON output
  behind the command module's forwarding wrapper. A follow-up startup pass
  removes the remaining annotation-only `typing` import by using built-in
  postponed annotations.
- `TraitRateMap` baseline branch-rate computation traversed the same tree
  repeatedly via Bio.Phylo for node labels, tip values, ancestral values,
  parent maps, and branch rates. The optimized path uses direct stack traversal
  helpers and shares preorder/postorder clade lists across the run-path
  computation steps. A later traversal-helper pass keeps preorder output
  identical while using a one-/two-child fast path instead of
  `reversed(children)` for direct preorder materialization. A later branch-rate
  pass keeps the same direct preorder
  list but caches clade ids, switches repeated membership checks to
  `dict.get()`, and labels terminal children from direct child lists. The
  parent-map helper now builds maps directly when no preorder list is supplied,
  keeping supplied-preorder behavior unchanged while avoiding generator and
  ordered-child overhead for map-only setup. A later ancestral-reconstruction
  pass builds postorder by reversing root-right-left preorder, preserving
  Bio.Phylo postorder while removing visited-flag stack entries. Cached read-only
  `TraitRateMap.run` now avoids the extra cached tree
  copy before creating its own protective copy for pruning and optional
  ladderizing. A later run-setup pass skips that protective copy entirely when
  all tree tips have trait values and ladderize is disabled, copying only before
  pruning or ladderizing would mutate the cached tree. A later startup pass
  postpones annotations and converts annotation-only `typing` aliases to
  built-in annotations, so command discovery no longer imports `typing`. Its
  single-trait parser now uses the same streaming, `partition("\t")`-based
  validation path and avoids `readlines()` plus temporary split lists on valid
  rows. The multi-trait parser now bounds row splitting to the selected column,
  preserving selected-column parsing and ignored trailing fields while avoiding
  tokenizing unused tail columns in wide trait tables.
- `TraitRateMap._plot_rate_map` layout setup baseline time rebuilt tip order,
  x-coordinates, and internal y-coordinates with `get_terminals()` and separate
  preorder/postorder `find_clades()` scans. The optimized path uses the shared
  direct `compute_node_positions` helper and reuses a single direct preorder
  list for rectangular and circular branch drawing loops. A later pass lets
  `compute_node_positions` consume that existing preorder list directly,
  avoiding a second coordinate-layout traversal in the TraitRateMap plot setup
  while producing identical node coordinates. A later circular setup pass passes
  the same preorder list plus the prepared tip list into
  `compute_circular_coords`, avoiding another direct tree walk while preserving
  the existing circular tip order. A later rectangular plotting pass
  batches gray vertical connectors and rate-colored horizontal branches into
  two `LineCollection`s, preserving per-branch colors and linewidths while
  avoiding one Matplotlib `Line2D` artist per segment. A later circular plotting
  pass batches rate-colored radial branches and gray internal arcs into two
  `LineCollection`s, preserving per-branch colors and widths while leaving the
  colorbar and annotations unchanged. A later rate-color pass caches colormap
  results by exact branch rate inside each plot call, preserving colors while
  avoiding repeated normalization and colormap calls when branches share a rate.
  Multi-trait parsing now streams the selected-column table directly from the
  file handle instead of materializing both `readlines()` and a stripped
  `clean_lines` list, preserving comment/blank filtering, extra-column
  tolerance, and logical data-row error numbering.
  A startup pass removes
  the only direct NumPy use by computing the branch-rate summary mean with
  Python arithmetic, so command import avoids NumPy entirely. Later startup work
  localizes the protective-copy `pickle` import, `PlotConfig`, layout helpers,
  circular/color helpers, and keeps JSON output behind the command module's
  forwarding wrapper. Text summary output now batches the conditional summary
  lines into one newline-joined print while preserving exact stdout text.
- `PhyloLogistic._build_logistic_vcv` baseline time still used repeated
  original-tree root-to-tip `tree.distance` calls for OU diagonal correction
  after building the transformed VCV. The optimized path reuses cached
  root-relative terminal depths for that correction and the starting mean
  root-to-tip distance helper. `_root_tip_distances` now uses a direct
  name-to-distance traversal for standard parsed trees before falling back to
  Bio.Phylo `depths()`/`get_terminals()` or per-tip `distance()`. A startup pass
  defers direct NumPy imports behind a lazy proxy, so command discovery avoids
  array startup until regression fitting begins. A later startup pass converts
  annotation-only `typing` names to built-in postponed annotations so command
  discovery no longer imports `typing`. A later likelihood-loop pass computes
  root-to-tip distances once per fit and reuses them across optimizer VCV
  builds, and uses scalar `math.exp` for per-branch OU transforms while retaining
  NumPy vector exponentials for diagonal correction. The saturated-starting-value
  fallback now counts validated binary response classes with one
  `np.count_nonzero` call and subtraction instead of two equality-mask
  reductions. OU diagonal correction now updates the VCV diagonal through a flat
  matrix view, avoiding `np.diag` allocation and `np.fill_diagonal` dispatch
  after the transformed VCV has already been built.
- `PhyloLogistic._compute_info_matrix` baseline time formed an explicit inverse
  of the OU-transformed VCV and allocated a dense diagonal weight matrix for
  every Firth-correction likelihood evaluation. The optimized path scales the
  design matrix directly and solves against one Cholesky factorization, while
  retaining the inverse-based implementation as a fallback. Logistic IRLS
  starting values now also row-scale the design matrix instead of materializing
  a dense diagonal weight matrix before cross-products, and use a direct ndarray
  max reduction for the small coefficient convergence vector. The inverse
  fallback uses the same row-scaled weighted design matrix to avoid allocating a
  dense diagonal weight matrix when Cholesky factorization cannot be used. A
  later pass removed the eager `scipy.stats` import by computing two-tailed
  z-test p-values with the standard-library complementary error function. A
  later p-value pass keeps the `scipy.stats` import out of this path but uses cached
  `scipy.special.erfc` as a vectorized ufunc for coefficient z statistics,
  preserving the same two-tailed normal probabilities with far less Python loop
  overhead. A subsequent startup pass replaced eager `scipy.linalg` and
  `scipy.optimize` imports with same-name lazy wrappers, so import-only callers
  avoid linalg/optimize startup while fitting still calls the same Cholesky and
  optimizer implementations. A repeated-call pass caches the resolved SciPy
  linalg/optimizer callables after first use, keeping import deferral while
  avoiding repeated import machinery during fitting iterations.
  Follow-up startup work defers pickle to the fallback VCV-copy path and keeps
  JSON output and trait parsing behind local forwarding wrappers, preserving
  test patch points while avoiding those helper imports for import-only callers.
  A cached-tree setup pass switches `run()` to `read_tree_file_unmodified()`,
  avoiding an unnecessary copy of the parsed species tree because the command
  validates, reads tips, prepares response/design arrays, and fits against
  non-mutating VCV builders. Text report output now batches the header,
  coefficient rows, and footer into one newline-joined print while preserving
  exact stdout text. JSON coefficient mapping now zips coefficient names and
  statistic arrays directly, preserving the same payload while avoiding repeated
  integer indexing. Final coefficient standard errors now extract the inverse
  diagonal from a Cholesky solve for positive-definite information matrices
  while preserving the explicit inverse fallback for non-Cholesky cases.
- `FitContinuous._vcv_ou` baseline time filled the OU-transformed VCV with a
  nested upper-triangle Python loop. The optimized path applies the same
  Martins-Hansen elementwise formula to the full matrix using broadcasted tip
  depths and shared-path values.
- `FitContinuous._build_transformed_vcv` baseline time transformed every
  root-to-tip path and then compared path prefixes for every tip pair. The
  optimized path accumulates each transformed branch length directly into the
  descendant-tip submatrix, preserving the same shared-prefix semantics for EB,
  Delta, and Kappa model VCVs. A later pass writes one-tip terminal branches
  directly to the diagonal instead of constructing one-element index arrays and
  `np.ix_` tuples for terminal-only branch updates.
- `FitContinuous._build_parent_map` baseline setup used Bio.Phylo's generic
  preorder `find_clades()` traversal. The optimized path walks standard parsed
  trees with a direct stack traversal and retains the original traversal as a
  fallback for nonstandard tree-like objects. Cached read-only
  `FitContinuous.run` now uses the explicit unmodified tree read helper to
  avoid copying the cached parsed tree before VCV/path setup and model fitting.
  A later helper pass localizes stack operations and pushes children directly
  because the result is an id-to-parent map rather than an ordered traversal.
  `_print_text_output` now builds the model-comparison report as one
  newline-joined string and tracks the best BIC model while formatting rows,
  preserving exact stdout text while avoiding per-row `print` calls and a
  separate BIC pass. A follow-up text-output pass reuses one row-format
  template and tracks the best BIC scalar separately, reducing per-row
  formatting and dictionary lookup overhead while preserving identical stdout
  text. A later model-comparison pass computes the tiny AIC weight vector with
  scalar `math.exp` calls, avoiding a NumPy array allocation for the fixed
  seven-model comparison. A later formatter pass switches the fixed-width model
  rows from `.format()` to percent formatting, preserving the same padded text
  with less per-row formatter overhead. Trait parsing now streams directly over the file handle,
  validates two-column rows with a single tab partition, and builds the trait
  taxa set directly from the parsed dictionary, avoiding full-file
  materialization and temporary split lists on valid rows. A later parser pass
  returns immediately for exact tree/trait taxon matches after all rows are
  validated and at least three taxa are shared, avoiding shared/warning set
  construction while preserving the minimum-shared-taxa error. A startup pass
  converts annotation-only `typing` aliases to built-in postponed annotations,
  so command discovery no longer loads `typing`.
- `PhylogeneticGLM._make_ultrametric` baseline time was dominated by repeated
  root-to-tip `tree.distance` calls. The optimized path computes tree depths
  once and reuses root-relative terminal depths to derive `D`, `Tmax`, and mean
  height; the root-height maximum is now reduced directly on the ndarray to
  avoid generic NumPy dispatch. `_root_tip_distances` now uses a direct
  name-to-distance traversal for standard parsed trees before falling back to Bio.Phylo
  `depths()`/`get_terminals()` or per-tip `distance()`. Poisson/logistic IRLS
  and Poisson GEE information assembly now row-scale the design matrix instead
  of materializing dense diagonal weight matrices before cross-products. The
  Poisson GEE score multiply now applies `R^-1` to the score vector first,
  avoiding a larger left-associated `X' R^-1` intermediate. A later pass
  avoids generic NumPy reduction dispatch in the Poisson GEE convergence check
  by reducing the small coefficient-delta vector directly. A later pass removed
  the eager `scipy.stats` import by computing two-tailed z-test p-values with
  the standard-library complementary error function. A later p-value pass keeps
  the `scipy.stats` import out of this path but uses cached
  `scipy.special.erfc` as a vectorized ufunc for coefficient z statistics,
  preserving the same two-tailed normal probabilities with far less Python loop
  overhead. A subsequent startup pass replaced the eager `scipy.optimize` import
  with a same-name lazy wrapper, so import-only callers avoid optimizer startup
  while fitting still calls the same SciPy implementation. A follow-up startup
  pass defers NumPy behind a module-level proxy and postpones annotations while
  preserving the normal module-level patch point. Cached read-only
  `PhylogeneticGLM.run` now uses the unmodified tree reader because validation,
  trait parsing, ultrametric corrections, and fitting read from the tree without
  changing the parsed species tree. Gene-tree discordance VCV setup now builds
  the shared-taxa matrix once and reuses it as `_precomputed_vcv` after
  metadata-driven trait subsetting, preserving the ordering returned by
  `build_discordance_vcv()`. Final logistic MPLE and Poisson GEE standard errors
  now extract the inverse diagonal from a Cholesky solve for
  positive-definite information matrices, preserving the explicit inverse
  fallback for non-Cholesky cases. The saturated logistic starting-value
  fallback now counts validated binary response classes with one
  `np.count_nonzero` call and subtraction instead of two equality-mask
  reductions. Logistic and Poisson starting-value loops now use direct ndarray
  max reductions for their small coefficient convergence vectors. The SciPy
  linalg and optimizer wrappers now cache their imported callables after first
  use, matching the logistic command's wrapper pattern and avoiding repeated
  import-on-call overhead during iterative fits while preserving module-level
  patch points. A later startup pass keeps
  JSON output behind a lazy forwarding wrapper so plain
  imports avoid the JSON helper module. Text report output now batches the
  header, coefficient rows, and footer into one newline-joined print while
  preserving exact stdout text. A follow-up text pass inlines the existing
  significance-code thresholds in the coefficient-row loop, preserving
  byte-identical report text while avoiding one method call per coefficient. A
  subsequent startup pass wraps
  multi-trait parsing in the same lazy-forwarding style and converts
  annotation-only `typing` aliases to built-in postponed annotations, so command
  discovery no longer imports the shared trait parser or `typing`.
- `Dtt._compute_dtt` baseline time combined repeated root/node `tree.distance`
  calls with repeated clade traversal and terminal extraction inside nested
  loops. The optimized path caches root-relative depths, the preorder clade
  list, and each clade's matched tip names. A later pass replaces pairwise
  nested-loop disparity scans with closed-form average squared-distance and
  Manhattan-distance formulas. The observed avg-Manhattan formula now computes
  the sorted weighted sum with `einsum`, preserving the closed-form result while
  avoiding a broadcasted product matrix. DTT trapezoidal integration now uses
  `np.trapezoid` when present and falls back to `np.trapz` for older NumPy
  environments. A later context-setup pass builds preorder clades, terminal
  clades, parent links, and root-relative depths in one direct standard-tree
  traversal, retaining the generic `depths()`/`find_clades()` path for
  nonstandard tree-like objects. Terminal tree height now streams the maximum
  root-relative depth from the prepared terminal list without materializing a
  temporary height list. `_simulate_null` now prepares the tree-only DTT context once
  and reuses the same time grid, lineage memberships, and clade index arrays
  across every simulated trait matrix. A later `avg_sq` simulation pass
  preserves the same fixed-seed random stream while generating Brownian-motion
  trait matrices in one batch and applying the average-squared disparity
  formula across all simulations and clades, avoiding one full `_compute_dtt`
  pass per simulation; non-default disparity indices retain the scalar loop.
  Copied-tree pruning setup now uses the shared direct terminal-name traversal
  instead of materializing terminal clade objects when only names are needed. A
  later pass replaces the remaining branching-time by clade nested scan with
  start/end lineage events swept over sorted node times, preserving the same
  active-lineage rule at each time point. A later null-simulation pass skips
  row-wise interpolation only when observed and simulated grids match exactly
  and are strictly increasing, and computes simulated MDI reductions with one
  vectorized trapezoid call. A follow-up startup pass defers NumPy behind a
  module-level proxy and postpones annotations while preserving runtime
  `np.trapezoid`/`np.trapz` selection. Later startup work localizes `PlotConfig`
  to argument processing and keeps JSON output and trait parsing behind local
  forwarding wrappers, avoiding those helper imports for import-only callers. A
  follow-up startup pass converts annotation-only `typing` aliases to built-in
  postponed annotations, so command discovery no longer loads `typing`.
  Cached read-only `Dtt.run` now uses the
  unmodified tree reader before making its own protective copy for pruning,
  removing one cached tree copy while preserving isolation for the mutating prune
  step. A later run-setup pass skips that protective copy entirely when all tree
  tips have trait data, copying only before pruning missing trait taxa. Trait
  matrix setup now uses the shared selected-column vector helper for
  `--trait` mode and the shared row-matrix helper for multivariate mode,
  preserving ordered taxa while avoiding nested Python row-copy loops. The
  batched simulation path now
  preallocates the optional terminal time/value extension instead of using
  `np.append` and `np.column_stack`, and then assigns the terminal column in
  place before masking positive values to zero instead of constructing an
  intermediate `np.where` result. A later pass reuses the prepared DTT
  context's postorder tree metadata to compute simulated clade disparities from
  subtree sums and counts, avoiding repeated fancy-indexed subtree slices in
  the avg-squared batched simulation path. A later run-level pass passes the
  already-computed observed DTT curve into null simulation, so MDI calculation
  avoids recomputing the observed curve from the same prepared context while
  preserving direct `_simulate_null()` callers that only provide observed times.
  MDI p-values now count extreme simulated MDI values with `np.count_nonzero`
  before dividing by the simulation count, avoiding boolean mean reductions.
  The observed DTT event sweep now averages each already-materialized lineage
  disparity list with Python `sum`/`len`, avoiding one NumPy array conversion
  per lineage event while preserving identical relative-disparity values.
  Observed average-squared disparity now uses flattened dot-product
  sum-of-squares for the trait matrix and column-sum vector, preserving the
  closed-form pairwise-distance formula while avoiding temporary squared arrays.
  Batched average-squared simulation reductions now compute row sum-of-squares
  through flattened `einsum` row dots, preserving the closed-form disparity
  values while avoiding temporary squared trait cubes for total and clade
  disparity calculations.
  Text output now batches the header
  and per-time-point rows into one newline-joined print while preserving exact
  stdout text.
- `Phylomorphospace._reconstruct_ancestral_scores` and
  `PhylogeneticOrdination._reconstruct_ancestral_scores` baseline time was
  dominated by repeated root-to-node `tree.distance` calls after pruning/copying
  the plotted tree. The optimized paths compute depths once and fill
  root-relative node distances from that cache before running the unchanged
  ancestral score reconstruction. Copied-tree pruning setup now also converts
  `ordered_names` to a set once and uses the shared terminal-name traversal,
  avoiding quadratic list membership checks and terminal clade materialization.
  A later phylomorphospace pass skips its protective reconstruction-tree copy
  entirely when all tree tips are present in the score matrix, while still
  copying before pruning missing taxa.
  A matching ordination pass applies the same copy-on-prune rule to its
  ancestral-score reconstruction helper.
  Cached read-only `Phylomorphospace.run` now uses the unmodified tree reader
  before reconstruction makes its own protective copy, removing a redundant
  cached parsed-tree copy without changing pruning isolation.
  Run setup now also builds the full ordered trait matrix once with the shared
  row-matrix helper and derives the selected x/y plot axes from that matrix,
  avoiding a second nested Python row-copy loop before reconstruction.
  Ordination tree-color file resolution now collects numeric tip values and
  covered taxa in one file pass before ancestral edge-color reconstruction,
  preserving missing-taxon and nonnumeric-value fallbacks while avoiding a
  second scan of large external color tables.
  Phylomorphospace branch rendering setup now also uses a direct preorder clade
  traversal when building plot branch segments and branch color values, with the
  generic traversal retained for nonstandard tree objects. That direct preorder
  helper now pushes binary children right-then-left and indexes multifurcations
  backward, preserving clade order while reducing balanced 131072-tip median
  helper time from 0.051598s to 0.047193s. Phylomorphospace tree plotting now
  computes the maximum root distance directly from the values view without
  materializing a temporary list. Phylogenetic
  ordination tree overlays now use the same direct preorder branch traversal for
  segment and edge-color setup, and the tree plot now computes the maximum root
  distance directly from the values view without materializing a temporary list.
  PCA proportion setup now totals eigenvalues with the ndarray reduction method,
  avoiding generic `np.sum` dispatch for the bounded component vector while
  preserving the same variance proportions.
  A later `phylomorphospace` startup pass defers NumPy behind a module-level
  proxy, leaving array startup until trait matrices, ancestral estimates, or
  plotting color arrays are built. A follow-up startup pass localizes pickle to
  ancestral-score reconstruction, imports `PlotConfig` only during argument
  processing, and keeps JSON output and trait parsing behind local forwarding
  wrappers so import-only callers avoid those helper modules. A later import pass
  converts annotation-only `typing` aliases to built-in postponed annotations, so
  command discovery no longer loads `typing`.
- `PhylogeneticOrdination._multi_trait_log_likelihood` baseline time formed an
  explicit inverse of every candidate lambda VCV matrix. The optimized path
  uses Cholesky factorization and triangular solves for positive-definite
  matrices, while retaining the original inverse-based calculation as a
  fallback. A later startup pass replaced eager `scipy.linalg` and
  `scipy.optimize` imports with same-name lazy wrappers, preserving Cholesky and
  bounded-optimizer behavior while avoiding linalg/optimize startup on module
  import. A follow-up startup pass defers NumPy behind a module-level proxy,
  postpones annotations, and imports the shared PGLS lambda-bound helper only
  inside the lambda-correction path that uses it. Lambda matrix transforms
  inside `_multi_trait_lambda` now copy the cached source diagonal through
  ndarray access and restore it via a flat diagonal stride, avoiding
  `np.fill_diagonal` dispatch for each bounded-search likelihood evaluation.
- `PhylogeneticOrdination._multi_trait_log_likelihood_inverse` keeps the
  fallback explicit inverse but multiplies it by `[ones, Y]` once, then reuses
  `C_inv_Y` and `C_inv_ones` to center traits and form `Z.T @ C_inv_Z`. This
  removes per-trait `ones @ C_inv @ Y[:, j]` products and the later dense
  `C_inv @ Z` multiply. A repeated-likelihood Cholesky pass now applies the
  same reuse to the positive-definite path, solving `[1, Y]` once and deriving
  `C_inv_Z` from those solved columns.
- `PhylogeneticOrdination` PCA setup also formed a full VCV inverse to center
  traits and compute `Z.T @ C_inv @ Z`. The optimized run path uses the same
  Cholesky factor to solve for centered traits and the weighted centered matrix,
  preserving the inverse fallback for nonstandard covariance inputs. PCA
  correlation mode now scales the covariance matrix and centered scores with
  one-dimensional diagonal factors instead of materializing a dense diagonal
  matrix for `D^-1 R D^-1` and `Z D^-1`. Cached read-only
  `PhylogeneticOrdination.run` now uses the unmodified tree reader because
  validation, tip extraction, VCV construction, centering, and ordination output
  do not mutate the parsed species tree. Text output for PCA and dimensional
  reduction now builds the full report and prints it once, preserving row
  formatting while avoiding one `print` call per trait or taxon row. PCA JSON
  payload construction now iterates eigenvector and score matrix rows directly,
  preserving nested labels and float conversion while avoiding repeated
  two-dimensional indexing. Dimensionality-reduction JSON output now converts
  embedding rows once and zips them with dimension labels, preserving rounded
  nested payload values while avoiding repeated two-dimensional indexing.
  A follow-up startup pass localizes pickle to
  ancestral-score reconstruction, imports `PlotConfig` only while processing
  arguments, and keeps JSON output and trait parsing behind local forwarding
  wrappers so import-only callers avoid those
  helper modules. A later startup pass converts annotation-only `typing`
  aliases to built-in postponed annotations, so import-only callers no longer
  import `typing`.
- `PhylogeneticOrdination._center_traits_by_vcv_inverse` applies the same
  combined RHS multiply in the PCA centering fallback, deriving `C_inv_Z` as
  `C_inv_Y - outer(C_inv_ones, a_hat)` instead of multiplying the inverse by
  centered traits a second time. A later pass uses NumPy broadcasting for both
  centered trait subtraction and weighted-centered fallback products, avoiding
  dense `outer(...)` temporaries with the same centered matrices. The Cholesky
  centering path now also derives the weighted centered traits from solved
  `[1, Y]` columns instead of issuing a second triangular solve for `Z`.
- `AncestralReconstruction._anc_ml` baseline time was dominated by computing
  cross-covariances with repeated per-tip `tree.get_terminals()`,
  `tree.common_ancestor()`, and root-distance calls for every labeled internal
  node. The optimized path caches terminal paths and root-relative depths once,
  then finds each node-tip MRCA from cached ancestor paths while preserving the
  same conditional mean and CI formulas.
- The second-stage `AncestralReconstruction._anc_ml` root-path setup benchmark
  measured the remaining `tree.get_path()` calls inside the optimized
  cross-covariance route. The new path builds root-to-node paths for all clades
  in one preorder traversal and reuses them for terminal and internal node
  ancestor lookups. Continuous and discrete `run` copied-tree pruning setup now
  uses the shared terminal-name traversal instead of materializing terminal
  clade objects. Cached read-only `AncestralReconstruction.run` now avoids the
  extra cached tree copy before dispatching to continuous or discrete handlers,
  which retain their own protective copies for pruning and optional ladderizing.
  A later run-setup pass skips those protective copies entirely for all-shared,
  non-ladderized continuous and discrete reconstructions while still copying
  before pruning missing taxa or ladderizing the tree.
- `AncestralReconstruction._fast_anc` keeps the same direct traversal path for
  standard trees but now builds the parent map from the already materialized
  preorder list. A later traversal-helper pass keeps preorder output identical
  while using a one-/two-child fast path instead of `reversed(children)` for
  direct preorder materialization, caches clade ids inside the ML passes, and replaces repeated
  membership checks plus small list allocations with streaming sums and
  `dict.get()` lookups. Its sigma2 contrast helper now streams valid children
  and accumulated contrast sums directly, preserving the existing binary and
  polytomy contrast order while avoiding per-node child and contrast lists.
  The BM log-likelihood helper now solves `[1, x]` from one Cholesky
  factorization for positive-definite VCVs, preserving the old inverse-based
  fallback while avoiding full inverse materialization in normal continuous ASR
  runs. A later parent-map helper pass localizes stack operations and pushes
  children directly because the result is an id-to-parent map rather than an
  ordered traversal.
- `AncestralReconstruction._print_text_output` now batches the continuous ASR
  report header and internal-node estimate rows into one newline-joined print
  while preserving exact stdout text. Discrete text output now precomputes
  percent row formats for the Q matrix and posterior table rows, preserving
  exact table text while avoiding per-cell f-string concatenation for large
  internal-node tables.
- `AncestralReconstruction._anc_ml` now precomputes inverse-weighted trait,
  intercept, and residual vectors once, then batches all internal-node
  cross-covariance inverse products for conditional means and CIs. This keeps
  the same scalar ML formulas while avoiding one dense vector-matrix multiply
  per labeled internal node. A later residual-product pass derives
  `C_inv @ residuals` from the already available `C_inv @ x` and `C_inv @ 1`
  vectors, avoiding another dense matrix-vector multiply before conditional
  means are formed.
- `AncestralReconstruction` continuous and discrete plot setup baseline time
  rebuilt parent maps, tip order, x positions, and internal y positions with
  separate generic tree traversals, then scanned preorder again for branch
  drawing, node pies/arcs, and CI lookup. The optimized plot paths reuse one
  direct preorder list plus the shared direct coordinate helper. A later setup
  pass also feeds that preorder list into `compute_node_positions` and feeds the
  same preorder list plus the prepared tip list into `compute_circular_coords`
  for both contMap and discrete-ASR plots, avoiding additional direct tree
  walks while preserving rectangular and circular tip order. A later contMap
  normalization pass computes the estimate value range in one pass without
  materializing a temporary values list. A later contMap rendering pass batches
  the 50 sampled gradient subsegments for every
  rectangular branch into one `LineCollection` and vertical connectors into a
  second `LineCollection`, preserving sampled trait-gradient colors while
  avoiding one artist per subsegment. Circular contMap rendering now uses the
  shared whole-tree radial gradient helper, avoiding one `LineCollection` per
  branch, and the shared colored-arc helper, avoiding one `Line2D` artist per
  internal arc. A later contMap CI overlay pass batches rectangular vertical
  bars, rectangular caps, and circular tangential bars into `LineCollection`s and
  draws point estimates with one `scatter` call per plot, preserving black CI
  styling while avoiding several Matplotlib artists per CI node. A later
  discrete-ASR rendering pass batches
  rectangular gray base branches into
  horizontal/vertical `LineCollection`s and batches rectangular clade
  horizontal/vertical branches plus circular clade radial branches into
  `LineCollection`s, preserving color-file styling while avoiding one artist per
  highlighted branch segment. A later
  discrete-ASR pie rendering pass batches rectangular and circular posterior
  state `Wedge` patches into one `PatchCollection` per plot, preserving state
  colors, black edges, and linewidth while avoiding one Matplotlib patch artist
  per nonzero posterior slice. A later
  import-time pass removes unused direct SciPy imports and benefits from the
  shared `discrete_models` lazy SciPy wrappers, so continuous/default imports
  avoid discrete-model linalg/optimize startup. A later startup pass also
  defers direct NumPy imports in both the command module and shared discrete
  helper. A follow-up startup pass keeps pickle, VCV, discrete-model, plot,
  circular-layout, and color-annotation helpers behind lazy module-level
  call-throughs, so import-only callers avoid those runtime-only helper modules.
  A later startup pass keeps JSON output behind a lazy forwarding wrapper so
  command discovery avoids the shared JSON helper. A subsequent startup pass
  removes the annotation-only `typing` import under postponed annotations, so
  command discovery no longer imports `typing`. A later startup pass keeps
  `heapq.merge` behind a lazy proxy so command discovery avoids `heapq` while
  preserving the descendant-name merge patch point. Text and JSON descendant-tip
  labels now use the shared direct terminal-name traversal for parsed clades,
  preserving sorted descendant names while avoiding Bio.Phylo terminal object
  materialization for each internal output row. Discrete text output now batches
  the Q matrix and posterior rows into one newline-joined print, and uses each
  posterior array's own `argmax()` to avoid repeated lazy NumPy proxy lookups in
  the per-node MAP-state loop. Discrete marginal-posterior normalization now
  uses each state vector's own `sum()` method, avoiding repeated lazy NumPy
  proxy lookups in the per-node upward and posterior reductions. Multi-trait parsing now finds the first
  non-comment header and parses subsequent data rows in one file pass instead
  of building both a full `readlines()` list and a stripped `data_lines` list,
  preserving comment skipping, missing-column validation order, and data-row
  error numbering. The discrete multi-trait parser now uses the same streaming
  shape, preserving categorical state values and validation order while avoiding
  both full-file and stripped-line buffers for large state tables. Two-column
  continuous and discrete trait parsers now stream rows and use `partition("\t")`
  for the valid-row path, avoiding full-file buffering while preserving
  column-count diagnostics, non-numeric continuous-trait errors, warnings, and
  shared-taxon filtering. A later parser pass returns immediately for exact
  tree/trait taxon matches after all rows are validated and at least three taxa
  are shared, avoiding shared/warning set construction while preserving
  too-few-shared-taxa errors.
- `FaithsPD.calculate_faiths_pd` baseline time used Biopython
  `common_ancestor()` and searched each root-to-tip path to exclude branches at
  or above the community MRCA. The optimized path caches each community
  root-to-tip path once, derives the MRCA from shared path prefixes, and skips
  excluded prefix clades by ID while preserving include/exclude-root semantics.
  A later pass removes the remaining per-community-tip `get_path()` calls by
  counting selected descendants in one tree traversal, finding the MRCA from
  those counts. When the deduplicated community includes every tree tip, a
  follow-up pass now returns the accumulated total branch length immediately,
  skipping selected-descendant counting and MRCA setup. A default-root pass sums
  unique root-to-tip branches for selected taxa directly after the tree scan,
  avoiding full selected-descendant counting for partial communities when no
  MRCA branch exclusion is needed. A later pass accumulates each unique selected
  branch length as it is first marked, removing the follow-up full-tree branch
  summation while preserving shared-branch de-duplication. A pectinate-tree
  pass stops each default-root upward walk as soon as it reaches a branch already
  seen from an earlier selected tip, because all ancestors above that branch
  have already been marked. The `include_root=False` selected-descendant pass
  now handles binary nodes with direct child indexing instead of a generator
  `sum`, preserving the generic loop for multifurcations. Cached read-only
  `FaithsPD.run` now also uses the explicit unmodified tree read helper to
  avoid copying the cached parsed tree before computing PD. A later startup pass
  keeps JSON output behind a lazy
  forwarding wrapper so plain command imports avoid the JSON helper module.
  Another startup pass converts annotation-only `typing` aliases to built-in
  postponed annotations, so command discovery no longer loads `typing`. The
  tree scan now stores clade objects only for selected terminal names, while
  retaining the unique terminal-name set for all-tip compatibility, and the
  selected-descendant pass checks `clades` directly instead of calling
  `is_terminal()` at each node. Single-tip `include_root=False` communities now
  return PhyKIT's documented `0.0` result immediately after taxon validation,
  skipping the selected-descendant count pass. Community taxa deduplication now uses
  order-preserving dictionary keys in both the file loader and calculation
  entry point, preserving blank filtering and first-occurrence order while
  reducing large community-list setup time. The initial full-tree scan now
  pushes binary children right-then-left and indexes multifurcations backward,
  preserving clade order while reducing balanced 32768-tip sparse-community
  median time from 0.043449s to 0.032935s. The taxa-list reader now uses a
  local forwarding wrapper, preserving the module-level patch point while
  avoiding `phykit.helpers.files` during import-only command discovery.
- `ParsimonyScore._fitch_parsimony` baseline time performed a full postorder
  tree traversal for every site and used Python set intersections at each
  internal node. The optimized path encodes state sets as integer bitmasks and
  performs the Fitch downpass across all sites at once per clade, with the
  original set-based implementation retained as a fallback for very large
  state alphabets. A later pass maps terminal sequences to bitmasks through a
  reusable byte lookup table and uses the direct postorder traversal helper,
  while retaining the per-character path for wider Unicode states.
  A subsequent pass routes alignments with 16 or fewer sites through a scalar
  bitmask downpass, avoiding one tiny NumPy array allocation per clade while
  preserving the vectorized path for wider alignments. A follow-up pass caches
  the read-only per-site mask list for repeated short sequence strings,
  avoiding duplicate terminal mask conversion in large low-site alignments. A
  subsequent postorder-helper pass reduces traversal overhead in both the scalar
  short-alignment path and the vectorized path. The large-state set fallback now
  also materializes postorder clades once and reuses that order for every site,
  preserving the fallback semantics while avoiding one Bio.Phylo traversal per
  site. Non-verbose scoring now requests a score-only path, avoiding unused
  per-site score array/list materialization while preserving verbose output and
  the direct helper default.
- `ParsimonyScore._parse_alignment` baseline time materialized `SeqRecord`
  objects before immediately converting them into uppercase strings. The
  optimized loader streams titles and sequences with `SimpleFastaParser`,
  preserving first-token IDs, uppercase conversion, and last duplicate wins. A
  later pass reuses the shared streaming first-token parser, preserving
  multiline sequence joining and internal whitespace removal while avoiding
  parser tuple construction. A startup pass imports the FASTA parser inside the
  parser helper, avoiding Bio.SeqIO.FastaIO startup for import-only callers. A
  subsequent startup pass also defers NumPy behind a module-level proxy. Later
  startup passes keep stdlib JSON and annotation-only `typing` imports out of
  module import while preserving JSON output and public annotations.
- `IndependentContrasts._compute_pic` baseline time called
  `clade.get_terminals()` for every internal node to build output labels,
  repeatedly walking subtrees that had already been processed. The optimized
  path propagates descendant tip labels through the same postorder pass used
  for PIC values. A later traversal pass handles standard Bio.Phylo trees with
  one direct postorder traversal, avoiding separate `get_terminals()` and
  `find_clades()` generator overhead. A follow-up pass uses scalar `math.sqrt`
  and merges already-sorted child descendant labels directly instead of
  re-sorting concatenated labels at every internal node.
  `find_clades()` setup while retaining the original fallback.
- `FitContinuous._concentrated_ll` baseline time formed an explicit inverse and
  determinant of each candidate VCV matrix. The optimized path uses Cholesky
  factorization and triangular solves for positive-definite matrices, while
  retaining the original inverse-based calculation as a fallback. BM, OU, EB,
  Lambda, Delta, Kappa, and White fits all call this likelihood kernel.
  Root-to-tip path setup now builds the relevant terminal map with a direct
  tree traversal instead of materializing all terminals through Bio.Phylo before
  walking parent links for EB, Delta, and Kappa transforms. A later startup pass
  replaced eager `scipy.linalg` and `scipy.optimize` imports with same-name lazy
  wrappers, so command import avoids loading SciPy linalg/optimize until model
  fitting actually needs Cholesky solves or bounded optimization. A follow-up
  startup pass defers direct NumPy use behind a lazy proxy and localizes the
  PGLS `max_lambda` import to Lambda-model runs, avoiding PGLS/NumPy startup for
  import-only command discovery. A later startup pass keeps JSON output behind a
  lazy forwarding wrapper so plain command imports avoid the JSON helper module.
  A repeated-likelihood pass now solves `[1, x]` as one combined Cholesky right-hand
  side and derives the residual solve from those columns, avoiding two extra
  triangular solves per positive-definite concentrated likelihood evaluation.
  SciPy linalg and optimizer wrappers now cache their imported callables after the
  first use, avoiding import-on-call overhead during repeated model likelihood
  evaluation while preserving the module-level wrapper patch points used by tests.
  Lambda model matrix transforms now copy the cached source diagonal through
  ndarray access and restore it via a flat diagonal stride, avoiding
  `np.fill_diagonal` dispatch for each bounded-optimizer likelihood evaluation.
- `pgls_utils.max_lambda` baseline time used repeated root-to-tip and
  root-to-node `tree.distance` calls to detect ultrametricity and compute the
  lambda upper bound. The optimized path computes `tree.depths()` once, derives
  root-relative heights from that cache, and retains the original distance-based
  implementation as a fallback for nonstandard tree objects. A later pass
  computes tip heights, ultrametricity, and max parent height in one direct
  traversal for standard Bio.Phylo trees before falling back to the depths path.
  A follow-up direct-path pass streams minimum and maximum tip heights during
  that traversal instead of allocating a full tip-height list and scanning it
  again.
- `pgls_utils.pgls_log_likelihood` baseline time formed an explicit inverse of
  each candidate VCV matrix before estimating GLS coefficients. The optimized
  path uses Cholesky factorization and triangular solves for positive-definite
  VCV matrices, while retaining the original inverse-based implementation as a
  fallback. `estimate_lambda` benefits directly because each bounded optimizer
  evaluation calls this likelihood kernel. A later startup pass replaced eager
  `scipy.linalg` and `scipy.optimize` imports with same-name lazy wrappers, so
  import-only callers avoid SciPy linalg/optimize startup while likelihood and
  lambda-estimation paths still call the same SciPy implementations. Follow-up
  startup work defers NumPy behind a module-level proxy and replaces
  annotation-only `typing` aliases with built-in postponed annotations. A
  repeated-call pass caches the resolved SciPy Cholesky/optimizer functions
  after first use, preserving import deferral while avoiding repeated import
  machinery inside optimizer likelihood loops. Cholesky log-determinant paths
  now reduce `log(diag(L))` through the ndarray's `sum()` method, avoiding
  generic `np.sum` dispatch in shared likelihood kernels. Lambda matrix
  transforms inside `estimate_lambda` now copy the cached source diagonal through
  ndarray access and restore it via a flat diagonal stride, avoiding
  `np.fill_diagonal` dispatch for every bounded-search likelihood evaluation.
- `pgls_utils.fit_gls` baseline time applied `C^-1` separately while building
  `X' C^-1 X`, `X' C^-1 y`, and the residual sum of squares. The optimized helper
  applies `C^-1` once to the combined `[X, y]` right-hand side and reuses those
  products. Coefficient estimation now also multiplies `X' C^-1 y` before the
  small coefficient solve, avoiding a larger left-associated intermediate.
  A follow-up pass obtains the returned coefficient covariance by solving the
  small identity system instead of calling `inv()` on `X' C^-1 X`, preserving
  singular-design error handling and the public covariance output. A later RHS
  construction pass preallocates the combined matrix and fills `X` and `y`
  directly, avoiding `np.column_stack` overhead in repeated small and mid-sized
  GLS fits while remaining neutral for large matrix-multiply-dominated fits.
- `PhylogeneticRegression._fit_model` baseline time formed an explicit inverse
  of the VCV matrix and reused it across GLS fitting, model statistics, and
  likelihood calculation. The optimized path computes the same quantities from
  one Cholesky factorization and triangular solves for positive-definite VCV
  matrices, while preserving the inverse-based path as a fallback. A later pass
  removed the eager `scipy.stats` import by evaluating coefficient Student-t
  p-values and model F p-values through lazy `scipy.special` calls. A subsequent
  startup pass replaced eager `scipy.linalg` imports with same-name lazy wrappers,
  so import-only callers avoid linalg startup while model fitting still calls the
  same Cholesky implementations. The lazy SciPy special helpers now cache
  `stdtr` and `fdtrc` after first use, avoiding repeated import dispatch when
  many regression coefficient and model p-values are evaluated. A later
  p-value pass routes single-predictor model F tests through the equivalent
  two-tailed Student-t survival probability (`F(1, df) = t(df)^2`), avoiding the
  slower `fdtrc` helper while retaining it for multi-predictor F tests. A follow-up
  cached-tree setup pass switches `run()` to `read_tree_file_unmodified()`
  because the command validates, reads tips, builds VCV matrices, and formats
  regression output without mutating the parsed species tree. Another startup
  pass defers NumPy behind a module-level proxy, postpones annotations, and
  imports shared PGLS helpers only inside the lambda and inverse-GLS paths that
  need them. A subsequent startup pass defers JSON output and trait-file parsing
  behind local forwarding wrappers, preserving existing patch points while
  avoiding those helper imports for import-only callers. Another startup pass
  converts annotation-only `typing` aliases to built-in postponed annotations, so
  command discovery no longer loads `typing`. Text output now batches coefficient
  rows and summary statistics into one newline-joined print while preserving exact
  stdout text. A repeated-fit pass now solves `[X, y, 1]` as one combined
  Cholesky right-hand side and derives residual and centered-response solves from
  those columns, preserving the fallback inverse implementation while avoiding
  four extra triangular solves per successful Cholesky fit. A later fit-kernel
  pass solves the coefficient target and identity matrix as one normal-equation
  right-hand side, returning coefficients and coefficient covariance without
  explicitly inverting `X' C^-1 X`. JSON result formatting now converts the
  residual and fitted NumPy vectors to Python lists once before dictionary
  construction, avoiding per-item NumPy scalar conversion while preserving the
  same JSON-compatible payload values.
- `PhyloPath._fit_gls_from_vcv` replaces the post-lambda explicit VCV inverse
  used by d-separation tests and path coefficient fits with Cholesky
  factorization and triangular solves. End-to-end `_dsep_test` improves more
  modestly because lambda estimation remains part of every test. Design-matrix
  setup for d-separation tests and path coefficient fits now preallocates the
  intercept-plus-predictor matrix directly from standardized trait vectors,
  avoiding repeated `np.column_stack` list and ones-array setup. A later pass
  removed the eager `scipy.stats` import by using lazy `scipy.special` calls for
  chi-square survival probabilities and two-tailed Student-t p-values. A
  subsequent startup pass replaced eager `scipy.linalg` imports with same-name
  lazy wrappers, so import-only callers avoid linalg startup while path fitting
  still calls the same Cholesky implementations. The lazy SciPy special helpers
  now cache the resolved `chdtrc` and `stdtr` callables after first use, avoiding
  repeated import dispatch across many d-separation and model-comparison
  p-values. A follow-up startup pass defers NumPy behind a module-level proxy,
  postpones annotations, and imports shared PGLS helpers only inside the GLS
  paths that need them, so import-only callers do not load NumPy through the
  helper module. A later startup pass keeps JSON output behind a lazy forwarding
  wrapper and localizes `PlotConfig` construction to argument processing, avoiding
  those helper imports during command discovery. A follow-up p-value pass uses
  the closed-form chi-square survival sum for positive even degrees of freedom,
  matching Fisher's C model-comparison calls (`df = 2k`) while retaining the
  cached SciPy-special fallback for non-even helper calls. Fisher's C now uses
  scalar `math.log` for each already-clamped p-value, avoiding scalar NumPy
  dispatch in the model-comparison loop. Model relative likelihoods now likewise
  use scalar `math.exp` for each delta value, avoiding scalar NumPy dispatch in
  the model-weight loop. DAG node circles are
  now batched into one `PatchCollection`, preserving data-coordinate circle
  geometry and labels while avoiding one Matplotlib patch artist per variable
  node. Text output now batches ranked model rows and path coefficient rows into
  one newline-joined print while preserving exact report formatting. A follow-up
  output pass caches the two fixed row formatters, preserving report text while
  avoiding repeated f-string formatter setup in large model/coefficient tables. A
  later formatter pass switches those fixed-width rows from `.format()` to percent
  formatting, preserving byte-identical output with less per-row overhead. A later
  run-setup pass reads the cached tree without copying when branch lengths are
  complete and copies before validation when default branch lengths would be
  assigned. The branch-length preflight now pushes child lists directly because
  it only returns whether any branch length is missing. Multi-variable trait
  setup now builds the ordered full trait matrix once with the shared row-matrix
  helper and copies requested variable columns from it, avoiding one row scan per
  model variable. Model-variable validation
  and column extraction now reuse one first-occurrence trait-name index map,
  preserving `list.index()` behavior for duplicate names while avoiding repeated
  linear scans across wide trait headers. Selected-variable standardization now
  standardizes each copied column in place, reusing the computed mean and
  avoiding the intermediate `raw_data` dictionary while preserving
  constant-column centering behavior. The same standardization loop now uses
  each copied column's ndarray `mean()` and `std()` methods, avoiding generic
  NumPy reduction dispatch. Trait parsing now streams directly over
  the file handle, converts valid numeric rows through a list comprehension, and
  falls back to the detailed nonnumeric-value loop only when conversion fails,
  preserving existing filtered taxa and error messages. Exact all-shared trait
  files now return the parsed trait mapping after the same row validation,
  avoiding the shared-taxa intersection and filtered dictionary copy.
  Candidate model-file parsing now streams rows directly from the file handle
  instead of materializing `readlines()`, preserving malformed-line and edge
  validation while reducing large model-set parse time and peak memory.
  A later startup pass converts annotation-only `typing` aliases to built-in
  postponed annotations, so import-only callers no longer import `typing`. A
  repeated-fit pass now solves `[X, y]` as one combined Cholesky right-hand side
  and derives residual solves from those columns, avoiding two extra triangular
  solves per successful post-lambda GLS fit while retaining the inverse fallback.
  A follow-up fit-kernel pass solves the coefficient target and identity matrix
  together, avoiding an explicit inverse of `X' C^-1 X` while preserving the
  returned coefficient covariance matrix.
  DAG validation and topological ordering now keep their existing queue ordering
  but use an index cursor instead of deleting from the front of the list,
  avoiding quadratic queue shifts on wide candidate model graphs. Basis-set
  construction now builds parent sets once while collecting adjacent edges,
  preserving statement order and conditioning sets while avoiding two full
  parent scans for every non-adjacent variable pair. Conditional model averaging
  now accumulates coefficients by edge in one pass over candidate models,
  preserving sorted coefficient output and the same weighted standard-error
  formula while avoiding one model scan per unique edge. A follow-up averaging
  pass combines weight and weighted-coefficient accumulation for each edge,
  avoiding one generator scan before the standard-error pass.
- `PhylogeneticSignal._blombergs_k` baseline time recomputed each permutation's
  GLS-centered residual ratio one shuffle at a time. The optimized path still
  uses the same seeded permutation stream, but evaluates batches of permutations
  with matrix operations. The K and K_mult setup paths now build the required
  covariance inverse from a Cholesky identity solve for positive-definite VCV
  matrices while preserving the explicit inverse fallback for non-Cholesky
  cases. `_kmult_permutations` applies the same batching to multivariate K,
  preserving trait-row correlations within each shuffled taxon order while
  evaluating the GLS-centered residual ratios in bulk. K and K_mult p-values
  now count extreme permutation statistics with `np.count_nonzero` before
  dividing by the permutation count, avoiding boolean mean reductions.
- `PhylogeneticSignal._compute_r2_phylo` baseline time formed an explicit
  inverse of the phylogenetic VCV matrix before estimating the Brownian Motion
  residual variance. The optimized path uses Cholesky factorization and
  triangular solves for positive-definite VCV matrices, while retaining the
  original inverse-based calculation as a fallback. A repeated effect-size pass
  now solves `[1, x]` as one Cholesky right-hand side and derives the residual
  solve from those columns, avoiding two duplicate triangular solves per
  positive-definite R2 evaluation. The white-noise variance term now uses the
  trait vector's ndarray `var()` method, avoiding generic NumPy dispatch in both
  the Cholesky path and inverse fallback.
- `PhylogeneticSignal._log_likelihood` baseline time formed an explicit inverse
  of every candidate lambda covariance matrix. The optimized path uses Cholesky
  factorization and triangular solves for positive-definite matrices, while
  retaining the original inverse-based calculation as a fallback. `_pagels_lambda`
  benefits directly because each bounded optimizer evaluation calls this
  likelihood kernel. A later startup pass replaced eager `scipy.linalg` and
  `scipy.optimize` imports with same-name lazy wrappers, preserving the
  Cholesky and scalar-optimization behavior while avoiding linalg/optimize
  startup on module import. A repeated-call pass now caches those resolved
  SciPy Cholesky and optimizer callables after first use, avoiding import
  wrapper overhead inside repeated likelihood and lambda-optimization calls.
  A follow-up startup pass defers direct NumPy imports and localizes the PGLS
  `max_lambda` import to Pagel lambda runs, so Blomberg's K command discovery
  avoids PGLS/NumPy startup. Cached read-only
  `PhylogeneticSignal.run` now uses the explicit unmodified tree read helper to
  avoid copying the cached parsed tree before VCV and statistic calculation.
  A later startup pass defers JSON output and multi-trait parsing behind local
  forwarding wrappers, preserving patch points while avoiding those helper
  imports for import-only callers. Another startup pass converts annotation-only
  `typing` aliases to built-in postponed annotations, so command discovery no
  longer loads `typing`. A repeated-likelihood pass now solves `[1, x]` as one
  combined Cholesky right-hand side and derives the residual solve from those
  columns, avoiding two extra triangular solves per positive-definite likelihood
  evaluation. Lambda matrix transforms inside `_pagels_lambda` now copy the
  cached source diagonal through ndarray access and restore it via a flat
  diagonal stride, avoiding `np.fill_diagonal` dispatch for each bounded-search
  likelihood evaluation. The two-column trait parser now streams rows directly
  from the file handle, uses a single tab partition on valid rows, and builds
  taxa sets directly from parsed dictionaries while preserving comment/blank
  filtering, mismatch warnings, and detailed validation errors. A later parser
  pass returns immediately for exact tree/trait taxon matches after all rows are
  validated and at least three taxa are shared, preserving too-few-shared-taxa
  errors while avoiding shared/warning set construction.
- `PhyloImpute._estimate_complete_case_stats` baseline time formed an explicit
  inverse or pseudoinverse of the complete-case VCV matrix before estimating
  GLS trait means and residual trait covariance. The optimized path computes the
  same quantities from one Cholesky factorization for positive-definite VCV
  matrices. The inverse/pseudoinverse fallback now computes `C^-1 Y` once,
  derives all GLS trait means from that matrix product, and reuses centered
  right-hand sides when forming the residual trait covariance. A repeated-stats
  pass now applies the same reuse to the Cholesky path by solving `[1, Y]` as one
  right-hand side and deriving the centered weighted trait matrix from those
  columns, avoiding an extra `n x p` triangular solve.
- `PhyloImpute._impute_taxon` baseline time formed explicit inverses for each
  phylogenetic conditional prediction and each trait-correlation adjustment.
  The optimized path uses Cholesky solves for positive-definite conditional
  systems and retains the inverse/pseudoinverse path as a fallback. A later pass
  caches the per-taxon observed-trait covariance factor/inverse and observed
  trait delta across all missing traits for that taxon, avoiding repeated setup
  when many traits are missing for the same taxon. A repeated Cholesky pass now
  solves each missing trait's phylogenetic residual and variance-weight vectors
  as one two-column RHS, avoiding a duplicate triangular solve for the same
  factor. A subsequent startup pass
  replaced eager `scipy.linalg` imports with same-name lazy wrappers, so
  import-only callers avoid linalg startup while imputation still calls the same
  Cholesky implementations. A follow-up startup pass defers NumPy behind a
  module-level proxy and postpones annotations while preserving the module-level
  `np` access pattern used by tests. Cached read-only `PhyloImpute.run` now
  avoids copying the cached parsed tree because validation, tip matching, VCV
  construction, imputation, and output do not mutate it. A later startup pass
  keeps JSON output behind a lazy forwarding wrapper so import-only callers avoid
  the shared JSON helper. Another startup pass converts annotation-only `typing`
  aliases to postponed built-in annotations, keeping command discovery free of
  `typing`. Trait parsing with NA markers now finds the first non-comment
  header and parses subsequent rows in one file pass instead of materializing
  both `readlines()` and a stripped `data_lines` list, preserving missing-value
  metadata, shared-taxa filtering, warnings, and logical data-row error
  numbering. A later parser pass validates off-tree rows for column counts and
  non-numeric values, but skips retaining their numeric values and missing-value
  metadata because those rows are discarded after shared-taxa filtering. A
  follow-up parser pass converts all-numeric rows with a list comprehension and
  falls back to the detailed NA/error loop only for rows containing missing
  markers or invalid values, preserving missing-value metadata and exact
  validation messages.
  Exact all-shared trait files now return the already retained trait mapping
  after the same row validation and missing-value metadata collection, avoiding
  mismatch set differences and the filtered dictionary copy when no warnings can
  be emitted. Off-tree rows now validate numeric values without allocating
  discarded trait-value lists, preserving off-tree nonnumeric errors and
  missing-marker handling.
  Discordance-VCV shared-taxa subsetting now builds the shared-taxon
  membership set once before row filtering instead of rebuilding it for every
  ordered taxon. The replacement trait-table writer now formats each NumPy row
  directly, preserving taxon order and six-decimal TSV output while avoiding
  repeated two-dimensional indexing. The run loop now reuses its initial
  missing-value mask for complete-case detection and per-taxon trait index
  selection, iterating only taxa with missing values and passing NumPy index
  arrays through to imputation instead of rebuilding scans and lists per row.
  The complete-case total now uses `np.count_nonzero` on that boolean mask
  instead of reducing booleans through `np.sum`. A
  follow-up writer pass formats bounded mixed string/numeric chunks through
  `np.savetxt`, preserving the exact TSV text while reducing per-float Python
  formatting overhead without materializing one object array for the full table.
  Text output now batches the imputation summary and per-missing-value rows into
  one newline-joined print while preserving exact stdout text. A follow-up
  text-output pass calculates taxon and trait column widths in one loop instead
  of two generator scans, preserving byte-identical table formatting.
- `TraitCorrelation._compute_correlation_matrices` baseline time formed a full
  VCV inverse before GLS-centering all traits and computing the phylogenetic
  covariance matrix. The optimized path solves the required right-hand sides
  from one Cholesky factorization. The inverse fallback now computes `C^-1 Y`
  once, derives all GLS trait means from that matrix product, and reuses the
  centered right-hand sides when forming the covariance matrix. P-value matrix
  assembly now computes off-diagonal t-statistics and survival probabilities in
  one vectorized pass, while preserving diagonal p-values of `1.0` and the
  scalar `0.0` p-value for perfect off-diagonal correlations. A later pass
  removed the eager `scipy.stats` import by evaluating two-tailed Student-t
  p-values through lazy `scipy.special.stdtr`. A subsequent startup pass
  replaced eager `scipy.linalg` imports with same-name lazy wrappers, so
  import-only callers avoid linalg startup while correlation computation still
  uses the same Cholesky implementations. A later p-value pass evaluates only
  upper-triangle trait pairs and mirrors the results into the lower triangle,
  halving the Student-t survival-probability work while preserving symmetric
  p-values. The lazy `stdtr` wrapper now caches the resolved SciPy special
  function after first use, avoiding repeated import dispatch on subsequent
  p-value batches. A follow-up startup pass defers NumPy behind a module-level
  proxy and postpones annotations while preserving the module-level `np` access
  pattern. A later startup pass keeps JSON output and trait parsing behind
  local forwarding wrappers and localizes `PlotConfig` to argument processing,
  avoiding helper imports for import-only callers. Cached read-only
  `TraitCorrelation.run` now uses the unmodified tree
  reader because validation, tip extraction, VCV construction, correlation
  matrix calculation, plotting, and output do not mutate the parsed species tree.
  Heatmap significance stars now render through grouped scatter collections for
  `*`, `**`, and `***`, preserving the existing p-value thresholds while avoiding
  one Matplotlib text artist per significant trait pair. Text output now builds
  matrix rows with list joins, emits the report with one print, and counts
  significant upper-triangle pairs with vectorized indexing while preserving
  exact stdout text. Heatmaps with no significant cells now skip off-diagonal
  mask construction and marker grouping before returning. A follow-up
  text-output pass inlines the star-threshold
  checks while formatting each matrix row, avoiding one method call per
  off-diagonal cell and preserving byte-identical output. A subsequent text
  pass joins local row cells without suffix duplication in every append,
  preserving byte-identical output while reducing per-cell string work. JSON payload
  construction now formats correlation and p-value matrix rows directly and
  reuses those row slices while collecting significant pairs, preserving rounded
  payload values and pair order. A later JSON pass rounds dense matrices through
  NumPy and selects significant pairs from the upper triangle in vectorized form,
  preserving rounded payload values and upper-triangle pair order. Sparse
  significant-pair output now scans upper-triangle rows directly and keeps the
  vectorized gather for dense outputs, avoiding full triangular-mask allocation
  when few pairs pass `alpha`. Another startup pass removes the annotation-only
  `typing` import by using postponed built-in collection annotations. A
  repeated-correlation pass now solves `[1, Y]` as one
  combined Cholesky right-hand side and derives the centered weighted trait
  matrix from those columns, avoiding an extra `n x p` triangular solve per
  positive-definite correlation computation.
- `NetworkSignal._blombergs_k` baseline time recomputed each permutation's
  GLS-centered residual ratio one shuffle at a time. The optimized path
  batches permutations, reuses precomputed inverse weights, and computes
  permutation numerators and denominators with vectorized row operations. A
  later permutation-core pass uses invariant permutation sums to compute the
  centered numerator algebraically and derives the denominator from the
  uncentered quadratic form, avoiding the centered residual matrix allocation
  while preserving the fixed-seed permutation stream. K p-values now count
  extreme permutation statistics with `np.count_nonzero` before dividing by the
  permutation count, avoiding boolean mean reductions. K setup now builds the
  required covariance inverse from a Cholesky identity solve for
  positive-definite network VCVs while preserving the explicit inverse fallback
  for non-Cholesky cases. The algebraic permutation setup now computes the
  invariant trait sum through the ndarray `sum()` method, avoiding generic NumPy
  dispatch while preserving the same centered numerator formula.
- `NetworkSignal._validate_tree` baseline time materialized terminal clades,
  traversed the tree again to check branch lengths, and then scanned all clades
  a third time to count unresolved polytomies. The optimized path gathers tip
  count, missing branch lengths, and polytomy count in one direct stack
  traversal while preserving the previous error priority and warning behavior.
  A later validation pass pushes child lists directly because the traversal only
  accumulates tip count, missing-branch existence, and polytomy count.
- `NetworkSignal._compute_network_vcv` baseline time copied every node's
  covariance row one cell at a time and then extracted the terminal VCV
  submatrix with nested Python loops. The optimized path uses slice assignment
  for tree and hybrid-node covariance rows and `np.ix_` indexing for the
  ordered tip submatrix.
- Cached read-only `NetworkSignal.run` now uses the unmodified tree reader
  because validation, tip matching, DAG construction, hybrid-edge handling,
  VCV construction, and output operate on read-only tree data or separate DAG
  dictionaries. Text output now also builds a parent-to-child index once for
  hybrid-edge donor reporting and emits the report in a single newline-joined
  print, avoiding a full tip-map scan for every reported hybrid recipient.
  Quartet JSON hybrid-edge inference now uses a stable top-two topology-index
  helper, avoiding per-hybrid-quartet `max()+index()` and three-item `sorted()`
  calls while preserving first-index tie behavior. A follow-up inference pass
  caches repeated topology-string validation, reuses the three derived quartet
  topology tuples for repeated taxa sets, and accumulates gamma sums directly
  instead of storing one gamma list per inferred edge. Another pass replaces
  the pair-evidence `Counter` with a plain dictionary and stable descending
  sort, preserving `most_common()` tie ordering while reducing per-quartet
  count-update overhead.
- `NetworkSignal._log_likelihood` baseline time formed an explicit inverse of
  each candidate network VCV matrix. The optimized path uses Cholesky solves
  and the Cholesky log determinant for positive-definite matrices, retaining
  the inverse implementation as a fallback. `_pagels_lambda` benefits directly
  because each bounded optimizer evaluation calls this likelihood kernel. A
  later startup pass replaced eager `scipy.linalg` and `scipy.optimize` imports
  with same-name lazy wrappers, preserving Cholesky and bounded-optimizer
  behavior while avoiding linalg/optimize startup on module import. A follow-up
  startup pass defers direct NumPy imports behind a lazy proxy, so command
  discovery avoids array startup until network VCV construction or signal
  calculations run. A later startup pass keeps JSON output and optional quartet
  JSON parsing behind lazy wrappers so plain imports avoid stdlib JSON and the
  shared JSON output helper. A follow-up startup pass removes the remaining
  annotation-only `typing` import by using a built-in postponed annotation. A
  repeated-likelihood pass now solves `[1, x]` as one combined Cholesky right-hand
  side and derives the residual solve from those columns, avoiding two extra
  triangular solves per positive-definite network-likelihood evaluation. A
  later repeated-call pass caches the resolved SciPy Cholesky and optimizer
  callables after first use, preserving lazy imports while avoiding import
  wrapper overhead inside repeated network likelihood and lambda-optimization
  calls. Lambda matrix transforms inside `_pagels_lambda` now copy the cached
  source diagonal through ndarray access and restore it via a flat diagonal
  stride, avoiding `np.fill_diagonal` dispatch for each bounded-search
  likelihood evaluation.
  The two-column trait parser now streams rows directly from the file handle,
  uses a single tab partition on valid rows, and builds taxa sets directly from
  parsed dictionaries while preserving comment/blank filtering, mismatch
  warnings, and detailed validation errors. A later parser pass returns
  immediately for exact tree/trait taxon matches after all rows are validated
  and at least three taxa are shared, preserving too-few-shared-taxa errors
  while avoiding shared/warning set construction. A follow-up strict-split
  pass parses valid two-column rows with one bounded split instead of
  `partition` plus a separate tab scan, preserving missing/extra-column counts
  and nonnumeric error text. The regime parser now applies the same strict split
  in binary mode, decoding only retained taxon/regime fields while preserving
  missing/extra-column counts and mismatch warnings.
- `RateHeterogeneity` setup baseline time used repeated Bio.Phylo traversal
  APIs while building the parent map, assigning Fitch branch regimes, and
  mapping tips to root-to-tip paths before per-regime VCV construction. The
  optimized setup helpers use direct stack traversals and allow the run path to
  share preorder/postorder clade lists. A later traversal-helper pass keeps
  preorder output identical while using a one-/two-child fast path instead of
  `reversed(children)` for direct preorder materialization, and replaces the
  standalone postorder helper's visited-tuple stack with reverse-preorder
  materialization. The parent-map
  helper now uses a direct
  map-only stack traversal when no preorder list is supplied, while retaining
  the supplied-preorder path for callers that already materialized traversal
  state. The branch-regime Fitch downpass now merges binary child state sets
  directly; the non-binary path now scans child state sets once with direct
  dictionary lookups instead of building a temporary child-set list and slicing
  it for the intersection loop. Copied-tree pruning setup now uses the shared terminal-name traversal
  instead of materializing terminal clade objects. Cached read-only
  `RateHeterogeneity.run` now avoids an extra cached tree copy before creating
  its explicit protective copy for pruning and
  optional ladderized plotting. A later all-shared run-setup pass skips that
  protective copy entirely when no taxa are pruned and ladderized plotting is
  disabled, while preserving copy isolation for pruning and ladderized plots.
  Text output now batches the summary and per-regime sigma rows into one
  newline-joined print while preserving exact stdout text. Trait-file parsing
  now streams directly over the file handle instead of materializing all input
  lines first, preserving comment/blank filtering and validation messages. A
  follow-up parser pass uses a single tab partition for valid two-column rows
  and builds taxa sets directly from parsed dictionaries, avoiding temporary
  split lists and keys views on the common path. Regime-file parsing now uses
  the same streaming valid-row parser and returns immediately for exact
  tree/regime taxon matches after all rows are validated, avoiding shared and
  warning set construction while retaining mismatch warnings for nonidentical
  taxon sets. Trait-file parsing now uses the same exact-match shortcut after
  all rows are validated, retaining mismatch warnings and shared-taxa errors
  for partial overlaps and nonidentical tree/trait taxon sets. Regime
  effect-size weighting now counts tip regimes in
  one `Counter` pass instead of rescanning all assignments once per regime,
  preserving zero counts for regimes with no assigned tips.
- `RateHeterogeneity._build_per_regime_vcv` baseline time compared every tip
  pair's root-to-tip path prefixes to assign shared path lengths by regime. The
  optimized path accumulates each branch length into the descendant-tip
  submatrix for that branch's regime, preserving the same per-regime covariance
  decomposition. A later pass writes terminal-only branches directly to the
  diagonal, avoiding one-element index arrays and `np.ix_` tuples for single-tip
  branch groups.
- `RateHeterogeneity._fit_single_rate` and `_fit_multi_rate` baseline time
  formed explicit covariance inverses and computed log determinants through
  `slogdet`. The optimized paths use Cholesky solves and Cholesky log
  determinants for positive-definite covariance matrices, retaining inverse
  likelihood helpers as fallbacks. `_fit_multi_rate` benefits across every
  optimizer likelihood evaluation. A repeated-likelihood pass now solves
  `[1, y]` as one Cholesky right-hand side for single-rate and multi-rate
  likelihood helpers and derives the residual solve from those columns, avoiding
  two duplicate triangular solves per positive-definite covariance evaluation.
  A later startup pass replaced eager
  `scipy.linalg` and `scipy.optimize` imports with same-name lazy wrappers, so
  command import avoids loading SciPy linalg/optimize until rate fitting runs.
  A repeated-call pass now caches the resolved SciPy Cholesky and optimizer
  callables after first use, avoiding import-wrapper overhead across repeated
  single-rate, multi-rate, and optimizer likelihood evaluations.
  A follow-up startup pass also defers NumPy behind a lazy proxy, leaving array
  startup until VCV construction, fitting, bootstrap, or plotting code first
  needs numeric arrays. A later startup pass keeps JSON output behind a local
  forwarding wrapper and localizes pickle, `PlotConfig`, circular-layout, and
  color-annotation helpers to the run or plot paths that need them. Another
  startup pass converts annotation-only `typing` aliases to built-in postponed
  annotations so command discovery no longer loads `typing`. A follow-up
  p-value pass handles the common three-regime LRT case with the
  closed-form chi-square survival function for two degrees of freedom, avoiding
  `scipy.stats` dispatch while preserving the existing SciPy fallback for
	  higher degrees of freedom.
- `OUShiftDetection` setup baseline time repeatedly used Bio.Phylo traversal
  APIs while building parent maps, tip lineages, candidate edge lists, and
  descendant counts. The optimized setup helpers use direct stack traversals,
  allow the run path to share one preorder list, and compute descendant counts
  for all candidate edges in one postorder pass. The parent-map helper now uses
  a direct map-only stack traversal when no preorder list is supplied, preserving
  the shared-preorder path for run setup. A later traversal-helper pass keeps
  preorder output identical while using a one-/two-child fast path instead of
  `reversed(children)` for the shared preorder materialization. Copied-tree
  pruning setup now uses the shared terminal-name traversal instead of
  materializing terminal clade objects.
  Cached read-only `OUShiftDetection.run` now avoids the extra cached tree copy
  before creating its protective working copy for pruning. A
  later all-shared run-setup pass skips that protective copy entirely when all
  tree tips have trait data, while still copying before pruning missing taxa.
  Trait parsing now streams directly over the file handle, validates
  two-column rows with a single tab partition, and builds the trait taxa set
  directly from the parsed dictionary. A later parser pass returns immediately
  for exact tree/trait taxon matches after all rows are validated, avoiding
  shared/warning set construction while retaining mismatch warnings for
  nonidentical taxon sets. Shift-edge descriptions now also use the shared
  direct terminal-name traversal for parsed clades, preserving existing label
  text while avoiding terminal object materialization for each reported edge. A
  later descendant-count pass builds the direct postorder iterator by reversing
  root-right-left preorder, preserving Biopython postorder while avoiding
  per-node visited stack tuples. A follow-up descendant-count pass handles
  binary clades with direct child-count addition and keeps the generic loop for
  multifurcations. A later
  startup pass keeps JSON output behind a lazy forwarding wrapper so command
  discovery avoids the shared JSON helper. Text output now batches the summary
  and detected-shift rows into one newline-joined print while preserving exact
  stdout text. A follow-up startup pass defers pickle until copied-tree pruning
  is required.
- `OUShiftDetection._fit_ou1_for_alpha` and `_fit_and_score_config` baseline
  time formed explicit OU VCV inverses for each candidate alpha likelihood
  evaluation. The optimized paths route profiled GLS likelihood calculations
  through a shared Cholesky-solve helper with inverse fallback. A repeated-GLS
  pass now solves `[W, x]` as one Cholesky right-hand side and derives the
  residual solve from those columns, avoiding two duplicate triangular solves
  per positive-definite GLS likelihood evaluation. The LASSO
  whitening transform also avoids explicitly inverting the Cholesky factor and
  instead applies triangular solves. Shift-weight matrix construction now uses
  the row ndarray reduction for the baseline-regime weight, avoiding generic
  `np.sum` dispatch on every taxon row. Indicator design matrix construction
  now caches the tip row indices under each lineage clade id and fills selected
  shift columns from those rows, avoiding a full root-to-tip path scan for every
  LASSO candidate configuration. A later pass deferred scikit-learn's
  `lars_path` import until LASSO path extraction, preserving the fitting path
  while avoiding `sklearn` on normal module import. A subsequent startup pass
  replaced eager `scipy.linalg` and `scipy.optimize` imports with same-name lazy
  wrappers, so import-only callers avoid linalg/optimize startup while fitting
  still calls the same Cholesky, triangular-solve, and optimizer
  implementations. Those lazy wrappers now cache the resolved SciPy callables
  after first use, avoiding repeated import dispatch inside GLS, whitening, pBIC
  scoring, and bounded alpha-optimization loops. LASSO path setup now computes
  residualized shift-column L2 norms with an `einsum` reduction instead of
  `np.linalg.norm`, preserving the same scaling while avoiding linalg dispatch.
  The pBIC Fisher-information
  correction now computes the log determinant of the scaled inverse directly
  from `slogdet(info)`, avoiding the coefficient-covariance inverse when only
  its determinant is needed. A follow-up startup pass also
  defers NumPy behind a lazy proxy, leaving array startup until trait parsing,
  VCV construction, LASSO setup, or scoring first needs numeric arrays. The
  cache-inclusive first measured batch was
  still faster (`0.005829s` vs. `0.161885s`), with the table reporting the
  steady-state repeated-candidate cost. pBIC scoring now solves only for
  `V^-1 W` and forms the small `W' V^-1 W` information matrix directly instead
  of materializing the full `V^-1` covariance inverse for each final candidate
  score. A subsequent startup pass removes annotation-only `typing` aliases
  under postponed annotations, so command discovery no longer imports
  `typing`.
- `PhyloAnova._cholesky_transform` baseline time explicitly inverted the
  Cholesky factor before whitening responses and design matrices. The optimized
  path applies triangular solves. `_run_anova` and `_run_manova` baseline time
  refit least-squares models for every RRPP permutation even though the design
  matrices are fixed; the optimized paths compute QR projection bases once and
  reuse them across observed and permuted fits. A later ANOVA pass computes
  each permuted full-model residual sum of squares from `y'y - ||Q'y||^2`,
  avoiding per-permutation fitted-value and residual-vector materialization
  while preserving the seeded RRPP permutation sequence. A matching MANOVA pass
  computes each residual SSCP matrix as `Y'Y - (Q'Y)'(Q'Y)`, avoiding
  per-permutation fitted-value and residual-matrix materialization before the
  Pillai trace calculation. Pillai trace reductions now use the ndarray
  reduction method directly for the observed and permuted eigenvalue vectors,
  avoiding generic `np.sum` dispatch in the MANOVA permutation loop.
  Permutation p-value/z-score summaries now compute permutation means and
  standard deviations through ndarray reductions, avoiding lazy NumPy proxy
  dispatch in the shared ANOVA/MANOVA and pairwise summary helper.
  Phylomorphospace PCA axis-label variance setup now uses a dot product for the
  singular-value total variance, avoiding a temporary squared array reduction
  after SVD while preserving the reported PC percentages.
  `_run_pairwise` now
  preserves the legacy `RandomState.permutation` sequence while computing
  permuted group residual means and distances in bulk for each group pair. A
  later pairwise pass computes permuted mean differences with a contrast-weight
  matrix and one residual matrix multiply, avoiding the larger
  `permutations x pair_size x traits` residual tensor for each group pair.
  Pairwise observed and permuted L2 distances now use direct dot/`einsum`
  reductions instead of `np.linalg.norm` dispatch. A design-matrix setup pass
  removes an unused all-ones allocation and caches
  group-name-to-column indices once, avoiding a linear `unique_groups.index()`
  search for every taxon while preserving treatment coding. ANOVA and MANOVA
  summary calculation now share a helper that computes permutation p-values,
  means, and standard deviations once, avoiding the duplicate `std` pass in
  z-score calculation. Permutation p-values now count extreme simulations with
  `np.count_nonzero` before dividing by the permutation count, avoiding boolean
  mean reductions in ANOVA/MANOVA summaries and pairwise tests.
  A subsequent startup pass replaced the eager `scipy.linalg.solve_triangular`
  import with a same-name lazy wrapper, so import-only callers avoid linalg
  startup while whitening still calls SciPy's triangular solve. A follow-up
  startup pass defers NumPy behind a module-level proxy and postpones
  annotations while preserving the module-level `np` access pattern used by
  tests. A later startup pass keeps JSON output behind a local forwarding
  wrapper and localizes `PlotConfig` to argument processing, avoiding helper
  imports for import-only callers. A follow-up startup pass converts
  annotation-only `typing` aliases to postponed built-in annotations.
  Phylomorphospace plotting now batches
  gray parent-child branch segments into one `LineCollection`, preserving branch
  styling while leaving group-colored tip scatter points and labels unchanged.
  `_print_results` now batches ANOVA/MANOVA summary and pairwise comparison
  rows into one newline-joined print while preserving exact text formatting.
  A later run-setup pass reads the cached tree without copying when branch
  lengths are already present and copies before validation when default branch
  lengths would be assigned. The branch-length preflight now pushes child lists
  directly because it only returns whether any branch length is missing.
  Group and response-matrix setup now precomputes response column indices and
  writes into a preallocated numeric matrix, avoiding the per-row conditional
  response-column scan while preserving single-response two-dimensional shape.
  Trait parsing now streams rows directly from the file handle instead of
  building both `readlines()` and a stripped `data_lines` list, preserving
  mixed categorical/numeric parsing and logical data-row error numbering.
  A follow-up parser pass precomputes the categorical group column and numeric
  columns once, binds `float`, and returns the validated trait dictionary
  directly when all parsed taxa are present on the tree, preserving mismatch
  warnings for partial-overlap inputs.
- `CharacterMap.run` taxon-set setup baseline time materialized terminal clade
  objects while only terminal names were needed for matching the character
  matrix. The optimized path uses the shared direct terminal-name traversal.
  Cached read-only `CharacterMap.run` avoids the extra cached tree copy while
  retaining its explicit working copy for pruning, polytomy resolution,
  branch-length filling, and optional ladderizing. A later all-clean run-setup
  pass skips that protective copy entirely when all taxa are shared, the tree
  is already bifurcating, all branch lengths are present, and ladderizing is
  disabled; the mutating cases still copy before changing the tree. Character
  matrix parsing now streams nonblank rows directly from the file handle instead
  of materializing a full filtered line list before splitting, preserving
  logical nonblank row numbering while reducing parse time and peak memory for
  large matrices. Node-label assignment now reuses the direct preorder helper
  and checks child lists instead of `find_clades(order="preorder")` plus
  `is_terminal()`, preserving tip labels and internal `node_N` numbering.
- `CharacterMap._summarize_character_states` baseline time scanned every
  character column once to count distinct states, then scanned the transposed
  character matrix again to count parsimony-informative characters. The
  optimized helper transposes once with `zip(*tip_states.values())` and reuses
  one wildcard-filtered `Counter` for both summaries. A later pass counts each
  full column with `Counter`, drops wildcard keys, and short-circuits the
  informative-state scan once two repeated states are found, avoiding filtered
  generator overhead for every character column. Subsequent byte-matrix passes
  compute counts for single-character ASCII state columns in bulk, and the full
  summary now derives returned per-character state lists from the same byte
  matrix transpose instead of a second Python `zip` transpose. The generic
  Counter path remains the fallback for multi-character, ragged, or non-ASCII
  states. Large observed ASCII alphabets now count per-character states with one
  column-offset `bincount` pass instead of one equality reduction per symbol,
  while small state alphabets keep the previous equality-count path to avoid
  `bincount` setup overhead.
- `retention_index` baseline time rebuilt a filtered non-wildcard list before
  counting each character column. The optimized path counts the full column
  once, removes wildcard keys from the `Counter`, and derives the observed
  taxon count from the removed wildcard count. A later shared-helper pass uses a
  byte matrix for single-character ASCII state columns and computes per-symbol
  column counts in bulk, with the Counter path retained for multi-character,
  ragged, or non-ASCII states.
- `CharacterMap.run` now keeps the wildcard-filtered per-character counters
  from state summarization and feeds them into a count-aware RI helper. This
  preserves the public `_summarize_character_states` return shape while avoiding
  a second `Counter` pass over each transposed character column during normal
  character-map runs. A later pass lets `run()` request counts only; common
  single-character ASCII matrices use a byte matrix and per-symbol column counts,
  while multi-character or non-ASCII states retain the Counter fallback. A
  follow-up counts-only pass streams rows directly into byte chunks while
  validating shape and ASCII compatibility, avoiding the temporary row list and
  giant intermediate text string. The
  public full-summary helper now also reuses those byte-derived counts when
  available and constructs only the required per-character state lists, avoiding
  one `Counter` per column without changing the returned state-list shape.
- `CharacterMap._print_verbose` baseline time scanned every classified branch
  change for every character while printing one line at a time. The optimized
  path groups changes by character once, then emits the same verbose report with
  one newline-joined print.
- `CharacterMap._print_summary` now batches the fixed summary report into one
  newline-joined print while preserving exact stdout text and broken-pipe
  handling.
- `CharacterMap._plot_character_map` layout setup baseline time rebuilt tip
  order, rectangular coordinates, and branch drawing traversal with
  `get_terminals()` plus separate preorder/postorder `find_clades()` scans.
  The optimized path uses the shared direct `compute_node_positions` helper and
  reuses one direct preorder list for rectangular/circular branch drawing and
  clade-color overlays. A later setup pass passes that preorder list into
  `compute_node_positions` and passes the same preorder list plus the prepared
  tip list into `compute_circular_coords`, avoiding additional direct tree
  walks while preserving circular tip order. A later traversal-helper pass keeps
  preorder output identical while using a one-/two-child fast path instead of
  `reversed(children)` for direct preorder materialization. A later rectangular
  plotting pass routes base branch
  drawing through the shared batched `LineCollection` helper, avoiding one
  Matplotlib `Line2D` artist per horizontal and vertical branch. A later
  color-overlay rendering pass batches rectangular clade horizontal/vertical
  branches and circular clade radial branches into `LineCollection`s, preserving
  color-file styling while avoiding one artist per highlighted branch segment. A
  later marker-rendering pass batches classified change circles into one
  `scatter` call per plot while preserving per-change colors, black marker
  edges, and individual annotations. A later import-time pass benefits from the
  shared circular-layout lazy NumPy arc cache, so command discovery avoids NumPy
  startup while circular drawing still materializes the cached arc array on
  demand. Later startup work localizes the protective-copy `pickle` import,
  `PlotConfig`, layout helpers, circular/color helpers, and keeps JSON and
  parsimony helper access behind command-module forwarding wrappers. Another
  startup pass converts annotation-only `typing` aliases to built-in postponed
  annotations so command discovery no longer imports `typing`.
- `RateHeterogeneity._plot_regime_tree` setup baseline time built y positions
  from terminal clades, x positions with a separate preorder scan, internal
  y positions with a postorder scan, and then scanned branches again for
  rectangular/circular drawing. The optimized path uses the shared direct
  `compute_node_positions` helper and reuses one direct preorder list for plot
  loops. A later setup pass passes that preorder list into
  `compute_node_positions` and passes the same preorder list plus the prepared
  tip list into `compute_circular_coords`, avoiding additional direct tree
  walks while preserving circular tip order. A later rectangular plotting pass batches gray vertical connectors and
  regime-colored horizontal branches into two `LineCollection`s, preserving the
  existing branch colors and z-order while avoiding one Matplotlib `Line2D`
  artist per segment. A later circular plotting pass batches regime-colored
  radial branches and internal arcs into two `LineCollection`s while preserving
  branch and arc colors.
- `StochasticCharacterMap._prune_tree_to_tip_states` baseline time collected tip
  names and pruned by name before SIMMAP setup. The optimized shared helper
  prunes terminal clade objects directly and is reused by stochastic character
  mapping, density-map, and SIMMAP summary workflows. A later pass uses the
  standard-tree batch pruning path when at least one terminal is retained,
  collapsing all missing-tip removals into one postorder traversal while
  retaining the object-prune fallback for all-tip and nonstandard trees. A
  later `run` setup pass reads the cached tree without copying, counts missing
  tip states with direct traversal, and only creates a working copy when pruning
  or ladderizing would mutate the tree. The missing-tip-state count now pushes
  child lists directly because the helper only returns a count. The prune-target
  scan now preserves terminal order while pushing binary children right-then-left
  and indexing multifurcations backward, avoiding `reversed(children)` iterator
  setup before batch pruning.
- `SimmapSummary` branch label generation baseline time called
  `clade.get_terminals()` for every internal branch label in text, JSON, and CSV
  output. The optimized path computes descendant tip-name tuples once in
  postorder and reuses them for branch and node labels. Text output now batches
  the summary headers, Q matrix, branch rows, and node posterior rows into one
  newline-joined print, and reuses one preorder list for the branch and node
  sections while preserving exact stdout text. A later cache-construction pass
  builds postorder by reversing root-right-left preorder, removing visited-flag
  stack entries while preserving the same descendant tip-name tuples.
- `SimmapSummary._plot_posterior_pie` setup baseline time materialized terminal
  clades and ran separate preorder scans for branch coloring and node pies. The
  optimized path reuses one direct preorder list for tip counting, branch loops,
  and posterior pie placement. A later setup pass passes that preorder list into
  `compute_node_positions`, avoiding another direct tree walk while preserving
  the same rectangular coordinates. A later rendering pass batches dominant-state
  horizontal branches and gray vertical connectors into two `LineCollection`s,
  preserving branch colors and alpha while avoiding one `Line2D` artist per
  segment. A later startup pass keeps JSON output behind a
  forwarding wrapper and localizes `PlotConfig` plus plot helper imports to
  argument processing and plotting paths, avoiding those helpers during command
  discovery. A subsequent `run` setup pass starts from the cached read-only
  tree, counts missing tip states directly, and only copies when pruning or
  ladderizing would mutate the tree.
- `QuartetPie` JSON/CSV output used the same repeated descendant-tip extraction
  for every reported internal branch. The optimized output path computes
  descendant tip-name tuples once in postorder and reuses them for JSON nodes,
  CSV branch labels, and the JSON taxon count. A later pass moves that postorder
  cache and the JSON/CSV branch loops onto direct traversal helpers, avoiding
  Bio.Phylo traversal overhead while preserving output ordering. A later
  cache-construction pass builds that postorder iterator by reversing
  root-right-left preorder, preserving Biopython order while avoiding per-node
  visited stack tuples.
- `QuartetPie.run` native mode previously scanned the species tree from the
  root once for every branch while attaching confidence/support values to
  branch info. The optimized path builds a clade-id-to-confidence map in one
  direct preorder traversal and reuses it for every branch. A follow-up startup
  pass converts annotation-only `typing` aliases to postponed built-in
  annotations, keeping command discovery free of `typing`. A later traversal
  pass pushes binary children right-then-left and indexes multifurcations
  backward, preserving Biopython preorder while avoiding per-node reversed
  iterator setup in the shared preorder helper.
- `QuartetPie._plot_quartet_pie` setup baseline time materialized terminal
  clades and ran separate preorder scans for pie placement, branch labels, and
  color overlays. The optimized path reuses one direct preorder list across
  rectangular and circular plot setup. A later setup pass passes that preorder
  list into `compute_node_positions` and passes the same preorder list plus the
  prepared tip list into `compute_circular_coords`, avoiding additional direct
  tree walks while preserving circular tip order. A later color-overlay rendering pass
  batches rectangular clade horizontal/vertical branches and circular clade
  radial branches into `LineCollection`s, preserving color-file styling while
  avoiding one artist per highlighted branch segment. A later startup pass
  removed an unused NumPy import, leaving numeric and plotting libraries out of
  plain command discovery imports. A follow-up startup pass keeps JSON output
  behind a forwarding wrapper and localizes plot, quartet, circular-layout, and
  color-annotation helper imports to the code paths that use them.
- `InternalBranchStats.calculate_internal_branch_stats` non-verbose baseline
  time still built descendant terminal-name rows for every internal branch even
  though summary output only reports aggregate statistics. The optimized run
  path requests those rows only for verbose text or verbose JSON output, while
  keeping the default helper API unchanged. A later pass replaced the remaining
  `find_clades()`/`get_nonterminals()` setup in `get_internal_branch_lengths`
  with direct standard Bio.Phylo traversals, retaining the original generic
  traversal fallback for nonstandard tree objects. Verbose `run` now batches
  internal length/descendant-tip rows into
  one newline-joined print while preserving the same space-separated stdout
  text, JSON payloads, summary output, and broken-pipe handling. Later startup
  passes defer the annotation-only Bio.Phylo import and the shared statistics
  helper's NumPy import until summary calculation. The non-verbose standard-tree
  stats path now collects internal branch lengths directly into a NumPy array
  before summary calculation, avoiding an intermediate Python list conversion.
  Cached read-only `InternalBranchStats.run` now also uses the explicit
  unmodified tree read helper to avoid copying the cached parsed tree before
  collecting internal lengths. A later verbose helper pass reuses the direct
  preorder output list and computes descendant terminal names in reverse
  preorder, removing visited-flag stack entries while preserving row order.
  A follow-up verbose helper pass uses order-preserving indexed child pushes
  instead of `reversed(children)` setup and specializes the common binary
  descendant-name merge, retaining the generic path for multifurcations.
  A later startup pass keeps JSON output behind a
  lazy forwarding wrapper so import-only callers avoid the JSON helper module.
  The latest startup pass applies the same forwarding-wrapper pattern to
  summary-statistics helpers, avoiding `phykit.helpers.stats_summary` during
  command module import. A follow-up startup pass converts annotation-only
  `typing` aliases to built-in postponed annotations, so command discovery no
  longer loads `typing`. The direct array helper and non-verbose list helper now
  push binary children explicitly and index larger child lists backward,
  preserving internal branch length order while avoiding the per-node
  `reversed(children)` iterator.
- `TerminalBranchStats.calculate_terminal_branch_stats` non-verbose used the
  same pattern at a smaller scale, building terminal name rows even for summary
  output. The optimized run path skips those rows unless verbose output needs
  them, while preserving the default helper return shape. A later pass replaced
  the remaining `get_terminals()` setup in `get_terminal_branch_lengths` with
  direct standard Bio.Phylo terminal traversal, retaining the original generic
  terminal fallback for nonstandard tree objects. Verbose `run` now batches
  terminal rows into one newline-joined print while preserving stdout text,
  JSON payloads, summary output, and broken-pipe handling. The non-verbose
  standard-tree stats path now collects terminal lengths directly into a NumPy
  array for summary calculation, then converts back to a list for the public
  helper return shape. Later startup passes defer the annotation-only
  Bio.Phylo import and the shared statistics helper's NumPy import until
  summary calculation. Cached read-only `TerminalBranchStats.run` now also uses
  the explicit unmodified tree read helper to avoid copying the cached parsed
  tree before collecting terminal lengths. A later startup pass keeps JSON
  output behind a lazy forwarding wrapper so import-only callers avoid the JSON
  helper module. The latest startup pass applies the same forwarding-wrapper
  pattern to summary-statistics helpers, avoiding `phykit.helpers.stats_summary`
  during command module import. A follow-up startup pass removes annotation-only
  `typing` aliases by using built-in postponed annotations. The direct array
  helper now pushes binary children explicitly and indexes larger child lists
  backward, preserving terminal-length order while avoiding the per-node
  `reversed(children)` iterator. The list-returning direct helper now uses the
  same order-preserving child push for verbose and non-verbose API callers.
  Command execution can now request `return_lengths=False`, so non-verbose
  `run()` keeps the default API return shape available to direct callers while
  avoiding an unused NumPy-array-to-list conversion before printing summaries.
- `calculate_summary_statistics_from_arr` and
  `calculate_summary_statistics_from_dict` baseline time mixed multiple Python
  `statistics` passes with repeated NumPy conversions. The optimized shared
  helper converts values once, computes all summary values from that array, and
  retains the existing fewer-than-two-values message/`None` behavior. A later
  dictionary-specific pass builds numeric arrays directly from dictionary values
  with `np.fromiter()`, retaining the generic list conversion fallback when
  one-pass numeric conversion is unavailable. The helper now also requests the
  25th percentile, median, and 75th percentile in one NumPy call instead of
  making separate percentile and median passes over the same array. A later pass
  derives sample standard deviation from the already computed sample variance,
  avoiding a second full-array reduction while preserving the reported values.
  Constant arrays now return after a guarded equality check, preserving identical
  summary values while skipping percentile, variance, and standard-deviation
  reductions for common all-equal branch-length/support summaries. Nonconstant
  summaries now compute the mean through the ndarray method, avoiding generic
  `np.mean` dispatch while preserving exact integer-mean formatting. The same
  direct ndarray reduction path now covers minimum, maximum, and sample
  variance, avoiding lazy proxy dispatch for the remaining non-quantile summary
  reductions.
  A startup pass wraps the NumPy module in a lazy proxy so importing commands
  that only reference summary helpers does not import NumPy until a summary is
  calculated. Summary output now emits the existing eight-line report through
  one newline-joined print call, reducing stdout overhead in commands that
  report shared statistics while preserving exact text. A later formatter pass
  keeps that single print call but builds the report with a static `%s`
  template fed by the existing `round(..., 4)` values, preserving integer and
  scientific-notation text while reducing per-report formatting overhead. The
  no-values diagnostic now also emits the existing three-line message with one
  print call, preserving exact stdout text while reducing repeated print
  overhead for empty summary requests. Sized empty and single-value array/dict
  inputs now return that same diagnostic before constructing NumPy arrays,
  preserving `None` results and stdout text while avoiding unnecessary NumPy
  startup and conversion work.
- `ThresholdModel._run_mcmc` baseline time recomputed full bivariate BM
  quadratic forms for every Metropolis-Hastings proposal and used scalar
  `truncnorm.rvs` calls for every discrete liability Gibbs update. The
  optimized path caches C-inverse sufficient statistics whenever liabilities
  are unchanged and evaluates proposal likelihoods from scalar arithmetic. The
  Gibbs sampler now draws threshold-truncated normal liabilities by inverse CDF
  with the existing SciPy sampler retained for degenerate tail intervals.
  One-sided liability bounds now use a specialized inverse-CDF path that avoids
  repeated generic bound normalization in the per-tip Gibbs loop while matching
  the previous inverse-CDF random stream. Liability Gibbs sweeps also keep an
  incremental residual vector during sequential tip updates instead of
  rebuilding the full residual array for every tip. A later MCMC pass caches
  invariant Gibbs diagonal/conditional-scale data once per discrete trait and
  refreshes only the changed liability vector's `C_inv @ x` product after each
  Gibbs sweep, preserving exact sampled chains for the benchmark seed. The
  truncated-normal sampler now also reuses cached float clipping bounds instead
  of calling `np.finfo(float)` for every liability draw. A later pass deferred
  the remaining `scipy.stats.truncnorm` import to the degenerate-interval
  fallback path, so normal threshold-model imports no longer load `scipy.stats`.
  A subsequent scalar-path pass moved hot likelihood/proposal arithmetic from
  the lazy NumPy proxy to `math`, calls cached SciPy special functions directly
  inside inverse-CDF liability draws, and trusts prevalidated Gibbs context
  standard deviations instead of repeating per-tip NumPy finite checks. A later
  Gibbs-sweep pass loads the cached inverse-CDF helpers once per sweep and
  inlines the one-sided threshold draw math, avoiding the generic sampler
  dispatcher and per-draw global helper checks in the discrete-liability loop.
  A subsequent startup pass deferred `scipy.special.ndtr`/`ndtri` and the shared
  VCV helper's `scipy.linalg.eigvalsh` import, avoiding SciPy special/linalg
  startup until inverse-CDF sampling or PSD correction actually needs it.
  The lazy `ndtr`/`ndtri` wrappers now cache the resolved SciPy special
  functions after first use, preserving cold-import behavior while avoiding a
  repeated import dispatch on every scalar inverse-CDF draw. A follow-up
  truncated-normal pass replaces scalar `ndtr` CDF calls in the direct
  one-sided sampler with `math.erfc`, preserving the same inverse-CDF samples
  while keeping SciPy special cached only for `ndtri` in that path. A later
  startup pass inherits the shared VCV helper's lazy Bio.Phylo parser, further
  reducing threshold-model import cost before tree parsing is needed. A follow-up
  startup pass defers NumPy behind a module-level proxy, postpones annotations,
  and uses `sys.float_info` for cached floating-point bounds so import-only
  callers avoid loading NumPy while tests retain the module-level `np` patch
  point. A later startup pass keeps JSON output behind a local forwarding
  wrapper, localizes `PlotConfig` to argument processing, and imports the VCV
  helper only when `_build_vcv_matrix` runs, avoiding helper and pickle startup
  for import-only callers. A later startup pass converts annotation-only
  `typing` aliases to built-in postponed annotations, so command discovery no
  longer imports `typing`.
  Discrete liability initialization now calls SciPy's truncated-normal sampler
  once with vectorized per-taxon bounds instead of once per taxon, preserving
  the same fixed-seed scalar draw stream while removing Python call overhead.
  The MCMC loop now also carries cached proposal standard deviations and
  accepted log-sigma values, refreshing proposal scales only when adaptive
  tuning changes the underlying proposal variance. MCMC setup now derives both
  `C^-1` and `logdet(C)` from one Cholesky factor for positive-definite VCV
  matrices, while preserving the explicit inverse plus `slogdet` fallback for
  non-Cholesky cases. Cached bivariate quadratic statistics now sum the
  precomputed `C^-1x` product vectors through each ndarray's `sum()` method,
  avoiding lazy NumPy proxy dispatch in repeated MCMC setup/proposal paths.
  Continuous-trait initialization now computes each liability variance once and
  reuses it for the positive-variance fallback, then uses each liability
  vector's ndarray `var()` and `mean()` methods for initial continuous-trait
  scaling. The covariance diagonal mean now comes from `trace(C) / n`,
  preserving the same scalar while avoiding a copied diagonal array. Posterior
  summaries and plots now use each small NumPy trace's `mean()` method while
  keeping the generic `np.mean` path for larger traces, avoiding lazy proxy
  dispatch on common sampled-chain summaries.
  `ThresholdModel.run` now
  reads the cached tree directly for validation and copies only when parsed
  trait taxa omit one or more tree tips before pruning, so all-shared trait
  matrices avoid the protective copy and no-op prune traversal. Text output now
  batches the fixed threshold-model summary report into one newline-joined print
  while preserving exact stdout text and broken-pipe handling. A follow-up
  output pass formats the same fixed report as one f-string, preserving exact
  text while avoiding the temporary list construction. A later output pass uses
  a static percent-format report template for the same fixed fields, preserving
  byte-identical stdout while reducing repeated summary-report formatting cost.
  Trait-table
  parsing now streams rows instead of materializing `readlines()` and collects
  discrete value sets during parsing, preserving header validation, extra-column
  handling, and filtered shared-taxon output while lowering peak memory. An
  exact tree/trait taxon match now returns after constructing the same sorted
  shared-name list and filtered trait dictionaries, avoiding mismatch set work
  when no warnings can be emitted.
  A later all-shared parser pass bounds row splitting to the header-width
  contract so ignored trailing columns are not fully tokenized, then returns
  the parsed trait dictionaries directly with the same sorted ordered-name list
  when every parsed taxon is present in the tree.
- `OUwie._concentrated_ll_bm` and `_fit_bms` baseline time formed explicit VCV
  inverses and log determinants for BM likelihoods. The optimized paths use
  Cholesky factorization, Cholesky log determinants, and triangular solves for
  positive-definite covariance matrices while retaining inverse-based helpers
  as fallbacks. `_fit_bms` benefits across every optimizer likelihood
  evaluation. A later repeated-likelihood pass solves `[1, x]` as one Cholesky
  right-hand side for BM1 and BMS likelihood helpers and derives the residual
  solve from those columns, avoiding two duplicate triangular solves per
  positive-definite covariance evaluation. `OUwie.run` copied-tree pruning
  setup now uses the shared
  terminal-name traversal instead of materializing terminal clade objects.
  Cached read-only `OUwie.run` now avoids the extra cached tree copy before
  creating its protective copy for pruning. A later all-shared run-setup pass
  skips that protective copy entirely when every tree tip has trait and regime
  data, while still copying before pruning missing taxa. A follow-up regime
  setup pass builds the object parent map with a direct stack traversal and
  computes branch regime assignments and the root regime from one direct
  state-set traversal instead of running the same Fitch-style pass twice.
  A later state-set pass merges the common two-child case directly; a follow-up
  model-comparison pass computes the small AICc weight vector with scalar
  `math.exp` calls instead of allocating a NumPy array. A follow-up
  shares the same direct non-binary child state-set merge as
  `RateHeterogeneity`, avoiding a temporary child-set list and slice for
  multifurcations.
  A later parent-map pass localizes stack operations and pushes children
  directly because the result is a map rather than an ordered traversal.
  Root-to-tip path setup now filters terminal clades with a one-time set of
  ordered names instead of scanning the full ordered-name list for every tip.
  Lineage-info setup applies the same set-backed terminal filter before OU
  weight and covariance construction, preserving lineage order and contents
  while avoiding another O(n²) ordered-name membership scan. The run path now
  passes the already-built lineage info into per-regime VCV construction,
  avoiding a second root-to-tip path rebuild before model fitting.
  Model-comparison R2 setup now averages unweighted multi-regime sigma2
  dictionaries directly over their values view, avoiding a temporary list when
  no regime-tip weights are available.
  A later startup pass replaced eager `scipy.linalg` and `scipy.optimize`
  imports with same-name lazy wrappers, so command import avoids loading SciPy
  linalg/optimize until likelihood evaluation or fitting actually runs. A
  follow-up startup pass defers NumPy behind a module-level proxy and postpones
  annotations so import-only callers avoid loading NumPy while preserving the
  module-level `np` access pattern used by tests. A later startup pass keeps
  JSON output behind a lazy forwarding wrapper so plain imports avoid the shared
  JSON helper. A follow-up startup pass defers pickle until copied-tree pruning
  is needed. Another startup pass converts annotation-only `typing` aliases to
  built-in postponed annotations so command discovery no longer imports
  `typing`. Regime-file parsing now mirrors the trait parser's all-shared fast
  return, preserving validation and warning behavior while avoiding mismatch
  set construction and shared-taxon filtering when every parsed regime taxon is
  present in the tree.
  `typing`. Text output now batches the model table and best-model parameter
  rows into one newline-joined print while preserving exact stdout text.
  Trait and regime parsing now stream directly over file handles, validate
  two-column rows with a single tab partition, and build taxa sets directly from
  parsed dictionaries, avoiding full-file materialization and temporary split
  lists on the valid-row path. A later trait parser pass returns immediately
  for exact tree/trait taxon matches after all rows are validated, avoiding
  shared/warning set construction while retaining mismatch warnings for
  nonidentical taxon sets.
- `OUwie._fit_ou1` and `_fit_oum` baseline time formed explicit OU VCV
  inverses and log determinants for each candidate alpha likelihood
  evaluation. The optimized paths share a profiled GLS likelihood helper that
  uses Cholesky solves and Cholesky log determinants for positive-definite OU
  covariance matrices, retaining the inverse implementation as a fallback.
  A later Cholesky-helper pass derives `V^-1 e` from the already solved
  `V^-1 W` and `V^-1 x` columns, removing one residual triangular solve from
  each profiled and GLS likelihood evaluation.
  The single-alpha OU VCV builder now also caches tip heights, shared-path
  lengths, and unique path lengths once per fit. Cache preparation now groups
  descendant tip indices by branch and adds branch lengths to the shared-path
  matrix in bulk, and each candidate alpha covariance matrix is rebuilt with
  NumPy array operations instead of walking root-to-tip path prefixes for every
  tip pair. A follow-up cache-building pass writes terminal branches directly
  to the shared-path diagonal instead of allocating one-element index arrays and
  `np.ix_` tuples for single-tip branch groups. The OUM/OUMV single-alpha weight
  matrix now uses the same pattern:
  per-tip branch rows, regime columns, branch lengths, and tail distances are
  cached once, then each optimizer alpha rebuilds the weight matrix with
  vectorized exponentials and indexed accumulation.
- `OUwie._build_per_regime_vcv` baseline time compared every tip-pair path
  prefix to find shared branches. The optimized path groups descendant tip
  indices by branch once and adds each branch length to the appropriate
  per-regime covariance submatrix in bulk. A follow-up pass writes one-tip
  terminal branches directly to the diagonal instead of allocating one-element
  index arrays and `np.ix_` tuples for every terminal branch.
- `OUwie._build_ou_H_matrices` used the same pairwise shared-prefix scan for
  OUMV per-regime OU covariance components. The optimized path groups
  descendant tip indices by branch and adds each branch's OU-decayed
  contribution with a vectorized outer product over descendant tip decays. A
  follow-up pass writes one-tip branches directly to the diagonal instead of
  constructing one-element index arrays and 1x1 outer products.
- `OUwie._build_ou_vcv_multi_alpha` had the same shared-prefix and repeated
  downstream-decay work for OUMA/OUMVA covariance matrices. It now computes
  branch-end-to-tip decays once per tip path and accumulates each branch's
  regime-specific alpha/sigma contribution with descendant-tip outer products.
  A later pass applies the same direct diagonal update for one-tip branches.
