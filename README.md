# Benders and C&CG for two-stage robust optimization

Python implementations of Benders decomposition and column-and-constraint generation (C&CG) for a two-stage robust location–transportation example. The models use IBM CPLEX through DOcplex, with synthetic instances from the included generator.

## Files

- `benders.py` — Benders implementation
- `ccg.py` — C&CG implementation
- `data_generator_fixed.py` — synthetic instance generator

## Run

Install IBM CPLEX and the Python packages in `requirements.txt`. Then run from this folder:

```bash
python -m pip install -r requirements.txt
python ccg.py
python benders.py
```

The instance settings are near the top of each algorithm script. As supplied, the scripts use different instance sizes, so their printed results are not a direct comparison. For a fair comparison, use the same generator settings in both scripts, including `p`, `q`, `density`, `seed`, and `total_cap_slack`.

## References and attribution

The model and decomposition methods follow [1]. Parts of the implementation were adapted from the Gurobi-based example [2] for CPLEX/DOcplex; that example also follows [1].

1. Zeng B, Zhao L. Solving two-stage robust optimization problems using a column-and-constraint generation method[J]. *Operations Research Letters*, 2013, **41**(5): 457–461. [doi:10.1016/j.orl.2013.05.003](https://doi.org/10.1016/j.orl.2013.05.003).
2. sometimesstudy. *Two-Stage-Robust-Optimization*: Python/Gurobi implementations of C&CG and Benders decomposition [source code]. GitHub. [Repository link](https://github.com/sometimesstudy/two-stage-robust-optimization).

**AI assistance:** The synthetic data generator was created using AI tools. I developed and adapted the Benders and C&CG implementations myself, with AI help on some coding and debugging steps.