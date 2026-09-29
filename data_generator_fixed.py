import numpy as np
from scipy.optimize import linprog

def generate_instance(p, q, density=1.0, seed=None,
                       cap_unit_range=(8,20), open_cost_range=(8,15),
                       cap_cost_range=(2,4), route_cost_range=(3,8),
                       nominal_demand_range=(5,12), delta_range=(3,8),
                       total_cap_slack=1.0, budget_frac=0.6, sub_budget_frac=0.4,
                       sparsify=('d',)):
    """
    Generalizes the robust facility-location structure:
      y: p binaries (open facility i)
      z: p continuous (capacity level at facility i)
      x: p*q continuous (route flow facility i -> customer j), facility-major order
      g: q uncertainty vars in [0,1] (extra demand fraction per customer)

    Structural matrices (A, G, E, M) are built directly from p, q -- their zero
    pattern is dictated by the model's logic, not randomized:
      A (p+1 x 2p): rows 0..p-1: cap_unit_i * y_i - z_i >= 0  (z_i <= cap_unit_i * y_i)
                    row p:       sum_i z_i >= total_cap_req
      G ((p+q) x p*q): rows 0..p-1: -sum_j x_ij >= -z_i  (capacity: sum_j x_ij <= z_i)
                       rows p..p+q-1: sum_i x_ij >= demand_j  (demand satisfaction)
      E ((p+q) x 2p): only rows 0..p-1 have a 1 at column p+i (ties z_i into capacity RHS)
      M ((p+q) x q): only rows p..p+q-1 have -delta_j on the diagonal (uncertain demand)
      h (p+q): rows 0..p-1 = 0; rows p..p+q-1 = nominal_demand_j

    density: fraction of entries in the cost/parameter vectors listed in `sparsify`
             that are actual nonzero random values -- the rest are hard zeros.
             density=1.0 -> every entry is a real (nonzero) value ("dense" data).
             density=0.3 -> only ~30% of entries are nonzero, rest are 0 ("sparse" data).
    sparsify: which parameter vectors density applies to. Options: 'd' (routing costs,
              the natural default -- this is the largest vector), 'c' (opening costs),
              'a' (capacity costs), 'cap_unit', 'nominal_demand', 'delta'.
              Pass a tuple of any subset, e.g. ('d','c','a').
    """
    rng = np.random.default_rng(seed)

    def sparse_values(size, low, high, name):
        vals = rng.uniform(low, high, size=size)
        if name in sparsify and density < 1.0:
            mask = rng.random(size) < density
            if not mask.any():          # guarantee at least one nonzero
                mask[rng.integers(size)] = True
            vals = np.where(mask, vals, 0.0)
        return vals

    nominal_demand = sparse_values(q, *nominal_demand_range, 'nominal_demand')
    delta          = sparse_values(q, *delta_range, 'delta')
    d              = sparse_values(p*q, *route_cost_range, 'd')

    # Scale cap_unit_range by the customer-to-facility ratio (q/p). Without this,
    # a fixed cap_unit_range means max achievable TOTAL capacity only scales with
    # p, while total demand scales with q -- so whenever q >> p (e.g. p=2, q=10),
    # even every facility at max capacity can't reach the required total, making
    # the master problem structurally infeasible regardless of any algorithm bug.
    ratio = q / p
    scaled_cap_unit_range = (cap_unit_range[0] * ratio, cap_unit_range[1] * ratio)
    cap_unit = sparse_values(p, *scaled_cap_unit_range, 'cap_unit')
    c        = sparse_values(p, *open_cost_range, 'c')
    a        = sparse_values(p, *cap_cost_range, 'a')

    # if cap_unit ended up 0 for some facility (only happens if 'cap_unit' is
    # in sparsify and density<1), that facility can never open usefully --
    # guard against a degenerate all-zero row by giving it a tiny nonzero value
    cap_unit = np.where(cap_unit == 0.0, 1e-3, cap_unit)

    # --- A, b ---
    A = np.zeros((p+1, 2*p))
    for i in range(p):
        A[i, i] = cap_unit[i]
        A[i, p+i] = -1.0
    A[p, p:2*p] = 1.0

    # uncertainty budget (sum of all g <= budget, plus a partial sub-budget on first half)
    total_budget = budget_frac * q
    sub_size = max(1, q // 2)
    sub_budget = sub_budget_frac * sub_size

    # --- Complete-recourse capacity sizing ---
    # For this fully-connected bipartite transportation structure (every
    # facility can route to every customer -- density/sparsify only ever
    # zeroes out COST values, never removes a route), a classical
    # transportation-problem feasibility result applies: a feasible routing
    # exists for ANY demand realization if and only if total installed
    # capacity >= total demand. So sizing total capacity to cover the
    # worst-case TOTAL demand the uncertainty set can produce is both
    # necessary and sufficient to guarantee complete recourse -- no
    # per-facility or per-customer case analysis needed.
    #
    # worst-case total demand = sum(nominal_demand) + max_{g in U} sum(delta_j * g_j)
    # solved exactly via a small LP (not a heuristic), matching the actual
    # uncertainty set (g in [0,1]^q, sub-budget on first half, total budget).
    A_ub_worst = np.vstack([
        np.concatenate([np.ones(sub_size), np.zeros(q - sub_size)]),
        np.ones(q)
    ])
    b_ub_worst = [sub_budget, total_budget]
    worst_case_lp = linprog(
        c=-delta, A_ub=A_ub_worst, b_ub=b_ub_worst,
        bounds=[(0, 1)] * q, method='highs'
    )
    worst_case_demand_increase = -worst_case_lp.fun

    total_cap_req = total_cap_slack * (nominal_demand.sum() + worst_case_demand_increase)
    b = np.zeros(p+1)
    b[p] = total_cap_req

    # --- G (fixed structural pattern, does not depend on density) ---
    G = np.zeros((p+q, p*q))
    for i in range(p):
        for j in range(q):
            col = i*q + j
            G[i, col] = -1.0        # capacity row i
            G[p+j, col] = 1.0       # demand row j

    # --- E (fixed structural pattern) ---
    E = np.zeros((p+q, 2*p))
    for i in range(p):
        E[i, p+i] = 1.0

    # --- M (fixed structural pattern; -delta_j values can themselves be sparse
    # if 'delta' is in sparsify) ---
    M = np.zeros((p+q, q))
    for j in range(q):
        M[p+j, j] = -delta[j]

    # --- h ---
    h = np.zeros(p+q)
    h[p:p+q] = nominal_demand


    return {
        'p': p, 'q': q, 'A': A, 'b': b, 'G': G, 'E': E, 'M': M, 'h': h,
        'c': c, 'a': a, 'd': d,
        'total_budget': total_budget, 'sub_size': sub_size, 'sub_budget': sub_budget,
        'bigM': 10 * max(d.max(), nominal_demand.max()*2, cap_unit.max(), 1.0),
    }

if __name__ == '__main__':
    print("Dense (density=1.0):")
    inst = generate_instance(p=3, q=3, density=1.0, seed=0)
    print("  d:", inst['d'])

    print("\nSparse (density=0.3), sparsifying d:")
    inst = generate_instance(p=3, q=3, density=0.3, seed=0, sparsify=('d',))
    print("  d:", inst['d'])

    print("\nSparse (density=0.3), sparsifying d, c, and a:")
    inst = generate_instance(p=3, q=3, density=0.3, seed=0, sparsify=('d','c','a'))
    print("  d:", inst['d'])
    print("  c:", inst['c'])
    print("  a:", inst['a'])

    for k in ['A','b','G','E','M','h','c','a','d']:
        print(k, np.shape(inst[k]))
