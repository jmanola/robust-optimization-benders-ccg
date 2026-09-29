# Example results

These runs use the same generated instance for both algorithms. Values below are copied from the scripts' printed iteration output.

## Parameters used

| Generator parameter | Value used |
| --- | --- |
| Facilities (`p`) | [`3`] |
| Customers (`q`) | [`100`] |
| `density` | [`1.0`] |
| `seed` | [`0`] |
| `total_cap_slack` | [`1.0`] |
| `budget_frac` | [default `0.6`] |
| `sub_budget_frac` | [default `0.4`] |
| `sparsify` | [default `('d',)`] |
| `cap_unit_range` | [default `(8, 20)`] |
| `open_cost_range` | [default `(8, 15)`] |
| `cap_cost_range` | [default `(2, 4)`] |
| `route_cost_range` | [default `(3, 8)`] |
| `nominal_demand_range` | [default `(5, 12)`] |
| `delta_range` | [default `(3, 8)`] |

| Solver/run setting | Value used |
| --- | --- |
| Stopping tolerance (`epsilon`) | `1e-4` |
| Maximum iterations | `100` |
| CPLEX version | [`22.2.0.0`] |

## Final results

| Algorithm | Final lower bound | Final upper bound | Iterations completed |
| --- | ---: | ---: | ---: |
| C&CG | [10144.286175930918] | [10144.286175930936] | [2] |
| Benders | [10144.286175906509] | [10144.286175962206] | [6] |

## C&CG iteration history

| Iteration | Lower bound | Upper bound |
| ---: | ---: | ---: |
| 1 | [10133.374242990489] | [10146.814507898296] |
| 2 | [10144.286175930918] | [10144.286175930936] |

## Benders iteration history

| Iteration | Lower bound | Upper bound |
| ---: | ---: | ---: |
| 1 | [9967.458883345995] | [10185.662725396016] |
| 2 | [10089.748385145665] | [10146.159395154678] |
| 3 | [10139.885161936018] | [10146.159395154678] |
| 4 | [10143.567977671284] | [10144.44388062215] |
| 5 | [10144.172877921219] | [10144.431213491072] |
| 6 | [10144.286175906509] | [10144.286175962206] |
