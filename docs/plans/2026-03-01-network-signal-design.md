# Network Phylogenetic Signal: Design Document

## Problem

All phylogenetic signal methods (Blomberg's K, Pagel's lambda) assume a tree. When the evolutionary history includes hybridization or introgression, the tree assumption is violated and signal estimates may be biased. No Python tool computes phylogenetic signal on networks. Blomberg's K on networks has never been implemented anywhere. Pagel's lambda on networks exists only in Julia (PhyloNetworks.jl).

## Goal

Add a `network_signal` command to PhyKIT that computes Blomberg's K and Pagel's lambda on phylogenetic networks, using the Bastide et al. (2018, Systematic Biology) network VCV algorithm.

## Architecture

The existing `phylogenetic_signal.py` K and lambda implementations only need a VCV matrix. The core new work is:

1. **Network representation**: Convert a species tree + hybrid edge specifications into an internal DAG
2. **Network VCV computation**: Bastide et al. 2018 recursive algorithm
3. **Reuse existing K/lambda math**: Same formulas, different VCV

## Input Specification

### Required arguments

- `-t / --tree`: Rooted species tree with branch lengths (Newick format)
- `-d / --trait-data`: Tab-separated trait file (taxon<TAB>value, same format as `phylogenetic_signal`)

### Network specification (one required)

**Option A: Explicit hybrid edges**
```bash
--hybrid donor:recipient:gamma [donor2:recipient2:gamma2 ...]
```
- `donor`: taxon name whose lineage is the source of gene flow
- `recipient`: taxon name whose lineage receives gene flow
- `gamma`: inheritance probability from the donor lineage (0 < gamma < 0.5)
- Multiple hybrid edges supported

**Option B: From quartet_network JSON output**
```bash
--quartet-json quartets.json
```
- Reads the JSON output from `phykit quartet_network --json`
- Auto-infers hybrid edges from hybrid quartet classifications
- Estimates gamma from concordance factor ratios

### Optional arguments

- `--method`: `both` (default), `blombergs_k`, or `lambda`
- `--permutations`: Number of permutations for K p-value (default 1000)
- `--json`: JSON output mode
- `-v / --verbose`: Print network VCV matrix and hybrid edge details

## Algorithm

### 1. Network Construction

**From tree + explicit hybrids:**
1. Parse rooted Biopython tree
2. Convert to internal DAG: assign integer IDs, build parent-edge list
3. For each hybrid spec `donor:recipient:gamma`:
   - Recipient node becomes a hybrid node with two parents:
     - Original tree parent (weight 1-gamma, original edge length)
     - Donor's tree parent (weight gamma, donor's edge length)
   - Validate: no cycles created (donor's parent must not be a descendant of recipient)

**From tree + quartet_network JSON:**
1. Parse JSON, extract hybrid quartets
2. For each hybrid quartet, identify the "swap pair":
   - Compare dominant and elevated minor topologies
   - Taxa that change sides between topologies form the swap pair
3. Aggregate swap counts per taxon pair
4. Select top pair(s) exceeding evidence threshold
5. Estimate gamma from CF ratios: `gamma = mean(CF_elevated_minor / (CF_elevated_minor + CF_lower_minor))`
6. Construct network as in Option A

### 2. Network VCV Computation (Bastide et al. 2018)

Process nodes in topological order (root first, tips last). Maintain a full node-by-node covariance matrix V.

**Root**: `V[root, root] = 0`

**Tree node** c with parent p, edge length l:
```
V[c, c] = V[p, p] + l
V[c, j] = V[p, j]  for all previously processed nodes j
```

**Hybrid node** h with parents p1 (weight gamma, length l1) and p2 (weight 1-gamma, length l2):
```
V[h, h] = gamma^2 * (V[p1,p1] + l1) + (1-gamma)^2 * (V[p2,p2] + l2) + 2*gamma*(1-gamma)*V[p1,p2]
V[h, j] = gamma * V[p1, j] + (1-gamma) * V[p2, j]  for all previously processed nodes j
```

Extract the tip-by-tip submatrix as the final VCV.

### 3. Blomberg's K (same as phylogenetic_signal.py)

```
C_inv = inv(VCV)
a_hat = (1' C_inv x) / (1' C_inv 1)     # GLS mean
e = x - a_hat                             # residuals
observed_ratio = (e'e) / (e' C_inv e)
expected_ratio = (trace(C) - n/sum(C_inv)) / (n-1)
K = observed_ratio / expected_ratio
```

P-value: permutation test (shuffle trait values, recompute K).

### 4. Pagel's Lambda (same as phylogenetic_signal.py)

```
V(lambda) = lambda * VCV + (1-lambda) * diag(VCV)
LL(lambda) = concentrated log-likelihood with sigma^2 profiled out
lambda_hat = argmax LL(lambda) over [0, 1]
LRT: 2*(LL(lambda_hat) - LL(0)) ~ chi2(df=1)
```

Multi-interval bounded optimization (10 sub-intervals, matching phytools).

## Output Format

**Text mode (default):**
```
Hybrid edges: B -> C (gamma=0.300)
Network taxa: 6
---
Blomberg's K: 0.8234    p-value: 0.0320
Pagel's lambda: 0.7651    log-likelihood: -12.3456    p-value: 0.0012
```

**JSON mode:**
```json
{
  "hybrid_edges": [{"donor": "B", "recipient": "C", "gamma": 0.3}],
  "n_taxa": 6,
  "blombergs_k": {"K": 0.8234, "p_value": 0.032, "permutations": 1000},
  "pagels_lambda": {"lambda": 0.7651, "log_likelihood": -12.3456, "p_value": 0.0012}
}
```

**Verbose mode** additionally prints the network VCV matrix.

## Key Design Decisions

1. **No Extended Newick parser** in v1: Adds significant complexity for limited immediate value. The tree + hybrid edges approach is simpler and sufficient. Extended Newick can be added later.

2. **Reuse existing K/lambda code**: The formulas are identical; only the VCV changes. Rather than duplicating code, import or refactor shared methods.

3. **Hybrid edge placement**: The recipient inherits from its tree parent (weight 1-gamma) and the donor's tree parent (weight gamma). This is a level-1 network approximation suitable for most biological scenarios.

4. **Gamma convention**: gamma is the inheritance weight from the *secondary* (donor) parent. Range: 0 < gamma < 0.5 (the minority contribution). The tree parent always contributes the majority.

## Validation

- **Tree-only case**: When gamma=0 or no hybrid edges, the network VCV should exactly equal the tree VCV from `phylogenetic_signal`
- **Known network**: Verify against PhyloNetworks.jl VCV computation on a simple 4-taxon network
- **Biological reasonableness**: K and lambda on networks should be between K/lambda on the tree and K/lambda on a star tree

## References

- Bastide, P., Solís-Lemus, C., Kriebel, R., Sparks, K. W., & Ané, C. (2018). Phylogenetic comparative methods on phylogenetic networks with reticulations. Systematic Biology, 67(5), 800-820.
- Blomberg, S. P., Garland, T., & Ives, A. R. (2003). Testing for phylogenetic signal in comparative data. Evolution, 57(4), 717-745.
- Pagel, M. (1999). Inferring the historical patterns of biological evolution. Nature, 401(6756), 877-884.
