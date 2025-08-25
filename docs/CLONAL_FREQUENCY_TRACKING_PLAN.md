# Clonal Frequency Tracking & Visualization – Implementation Plan

_This document describes how to extend the longitudinal pipeline to compute and visualise **clone-level frequencies** over time.  Copy-number adjustments are **explicitly ignored** as per requirement._

---
## 1. Current Longitudinal Workflow (high-level recap)
1. **`data_loader.py`** – loads per-sample VAF tables & metadata into tidy `pandas` frames.
2. **Optimisation (`optimize*.py`, `iterative_pipeline_real.py`)** – searches for the best phylogenetic tree; each node represents a clone; leaves contain mutation lists.
3. **`tree_updater.py` / `longitudinal_update*.py`** – refines trees across timepoints and finalises a _converged_ tree plus a mutation ➜ clone mapping.
4. **Visualisation (`visualize.py`, `longitudinal_visualizer.py`)** – currently plots **mutation-level** VAF trajectories only.

> Result: After convergence we still lack clone-level frequency curves.

---
## 2. Definition: Clone Frequency (no CN adjustment)
For a sample _t_ and a clone _C_:
* Let **M<sub>C</sub>** be all mutations assigned to clone _C_.
* For each mutation _m ∈ M<sub>C</sub>_ we already have VAF<sub>m,t</sub>.
* **Raw clone frequency**
  \[ F_{raw}(C,t) = \text{mean}\_{m∈M_C}\big(\text{VAF}_{m,t}\big) \]

### Root-correction demanded by PI
1. Compute **F<sub>subtree</sub>(t)** – frequency of the _entire unrooted subtree_ (all nodes below the root).
2. Compute **F<sub>leafsum</sub>(t) = Σ F_{raw}(L_i,t)** for all leaf clones _L_ under that subtree.
3. Derive **remainder**
   \[ F_{remainder}(t) = F_{subtree}(t) - F_{leafsum}(t) \]
4. Treat the remainder as a pseudo-clone so total frequency still sums to 1.

---
## 3. New Module: `clone_frequency.py` (<150 LOC)
| Function | Purpose |
| --- | --- |
| `compute_clone_frequencies(tree, vaf_df)` | Return tidy DF `[sample, time, clone_id, freq]` using mapping mutation ➜ clone. |
| `compute_subtree_remainder(freq_df, tree)` | Append pseudo-clone `"unrooted_remainder"` per sample with values from step 2. |

Intermediate result saved as **`clone_frequencies.csv`** via `output_manager.py`.

---
## 4. Visualisation
* Add **`clone_visualizer.py`** (or extend existing visualiser if size permits).
* Function `plot_clone_trajectories(freq_df, palette=None)` mirrors mutation-level plot: one line per clone, remainder plotted last in semi-transparent grey.
* Interactive and static outputs placed beside existing VAF plots.

---
## 5. Pipeline Integration
1. After convergence in `iterative_pipeline_real.py`:
   ```python
   import clone_frequency, output_manager, clone_visualizer
   freq_df = clone_frequency.compute_clone_frequencies(best_tree, vaf_df)
   freq_df = clone_frequency.compute_subtree_remainder(freq_df, best_tree)
   output_manager.write_clone_frequencies(freq_df)
   clone_visualizer.plot_clone_trajectories(freq_df)
   ```
2. Introduce config flag `track_clone_freq: true` (YAML) & CLI option `--plot_clonal_freq` handled by `config_handler.py`.

---
## 6. Testing & Validation
* **Unit tests** (pytest):
  * Correct aggregation on toy tree/VAF matrix.
  * Edge cases: single leaf, no leaves.
* **Plot regression test**: deterministic mini dataset → hash PNG to check for visual regressions.

---
## 7. Incremental Checklist
1. Create `clone_frequency.py` with both compute functions.
2. Add `write_clone_frequencies` to `output_manager.py`.
3. Propagate config flag through `longitudinal_main.py`.
4. Implement `clone_visualizer.py` (≤120 LOC).  Update docs.
5. Write unit & plot tests under `tests/longitudinal/`.
6. Update CI and README.

---
### File-touch Summary
```
src/longitudinal/clone_frequency.py    # new (core logic)
src/longitudinal/output_manager.py     # +write_clone_frequencies
src/longitudinal/clone_visualizer.py   # new (or extend existing visualizer)
src/longitudinal/config_handler.py     # parse new flag
docs/CLONAL_FREQUENCY_TRACKING_PLAN.md # this file
``` 