from docplex.mp.model import Model
from data_generator_fixed import generate_instance
import numpy as np

#Parameters generated with "generator.py"

inst = generate_instance(p=3, q=100, density=1.0, seed=0, total_cap_slack = 1.0)
A, b, G, E, M, h = inst['A'], inst['b'], inst['G'], inst['E'], inst['M'], inst['h']
c, a, d = inst['c'], inst['a'], inst['d']
bigM = inst['bigM']

#Create MP model

MP = Model(name = "MP")

#Ensure mipgap and integrality gap are tighter than epsilon

MP.parameters.mip.tolerances.mipgap = 1e-8
MP.parameters.mip.tolerances.integrality = 1e-9

#Continue with algorithm

LB = -1*float('inf')

UB = float('inf')

epsilon = 1e-4


k = 1

y = MP.binary_var_list(c.shape[0], name = 'y')
z = MP.continuous_var_list(a.shape[0], lb = 0, name = 'z')
eta = MP.continuous_var(name = 'eta', lb = 0)

#Make helper methods

def matvec(mat, variables):
    return [
        sum(mat[i,j]*variables[j] for j in range(mat.shape[1])) for i in range(mat.shape[0])
    ]

def dotvec(vec1, vec2):
    if len(vec1) != len(vec2):
        raise ValueError("They're dimensions, don't match")
    else:
        return sum(vec1[i]*vec2[i] for i in range(len(vec1)))

#Add constraints

Ay = matvec(A[:,:len(y)], y)
Az = matvec(A[:,len(y):], z)
MP.add_constraints(Ay[i] + Az[i] >= b[i] for i in range(A.shape[0]))

MP.minimize(dotvec(c, y) + dotvec(a, z) + eta)
mp_soln = MP.solve()

if mp_soln is None:
    raise RuntimeError("The master is infeasible or did not solve properly")

LB = max(LB, mp_soln.objective_value)


#Create SP

SP = Model(name = 'SP')

#Define mipgap and integrality tolerances

SP.parameters.mip.tolerances.mipgap = 1e-8
SP.parameters.mip.tolerances.integrality = 1e-9

#Create variables

x = SP.continuous_var_list(G.shape[1], name = 'x')
pi = SP.continuous_var_list(G.shape[0], name = 'pi')
g = SP.continuous_var_list(M.shape[1], lb = 0, ub = 1, name = 'g')
v = SP.binary_var_list(G.shape[0], name = 'v')
w = SP.binary_var_list(G.shape[1], name = 'w')

#Add constraints

GT = G.T

Gx = matvec(G,x)
yz = np.concatenate([
    [var.solution_value for var in y],
    [var.solution_value for var in z]
    ])
Eyz = matvec(E,yz)
Mg = matvec(M,g)
GTpi = matvec(GT, pi)

#Primal feasiblity

G1 = SP.add_constraints(Gx[i] >= h[i] - Eyz[i] - Mg[i] for i in range(G.shape[0]))

#Dual feasibility

SP.add_constraints(GTpi[i] <= d[i] for i in range(G.shape[1]))

#First complementarity pair

bigM_pi = 2*d.max()
SP.add_constraints(pi[j] <= bigM_pi * v[j] for j in range(len(pi)))
G2 = SP.add_constraints(Gx[i] - h[i] + Eyz[i] + Mg[i] <= bigM*(1-v[i]) for i in range(G.shape[0]))

#Second complementarity pair

SP.add_constraints(x[j] <= bigM*w[j] for j in range(len(x)))
SP.add_constraints((d-GTpi)[j] <= bigM*(1-w[j]) for j in range(G.shape[1]))

#Uncertainty set (Bertsimas-Sim style "budget of uncertianty")

sub_size, sub_budget, total_budget = inst['sub_size'], inst['sub_budget'], inst['total_budget']
SP.add_constraint(sum(g[:sub_size]) <= sub_budget)
SP.add_constraint(sum(g) <= total_budget)

dx = dotvec(d,x)
SP.maximize(dx)
sp_soln = SP.solve()

if sp_soln is None:
    raise RuntimeError("The subproblem is not feasible or did not properly solve")

cy_value = sum(c[i]*y[i].solution_value for i in range(len(c)))
az_value = sum(a[i]*z[i].solution_value for i in range(len(a)))

UB = min(UB, cy_value + az_value + sp_soln.objective_value)

pi_value = [var.solution_value for var in pi]
g_value = [var.solution_value for var in g]

Ey_sym = matvec(E[:,:len(y)], y)
Ez_sym = matvec(E[:,len(y):], z)

Mg_value = matvec(M, g_value)

cut_rhs = sum(
    (h[i] - Ey_sym[i] - Ez_sym[i] - Mg_value[i]) * pi_value[i]
    for i in range(len(h))
)

MP.add_constraint(eta >= cut_rhs)

#Start iterative loop

max_iterations = 100

while abs(UB-LB) >= epsilon and k <= max_iterations:

    mp_soln = MP.solve()

    if mp_soln is None:
        raise RuntimeError("The master is infeasible")
        
    LB = max(LB, mp_soln.objective_value)

    SP.remove_constraints(G1)
    SP.remove_constraints(G2)

    yz = np.concatenate([
                    [var.solution_value for var in y],
                    [var.solution_value for var in z]
                        ])
    Eyz = matvec(E, yz)

    G1 = SP.add_constraints(Gx[i] >= h[i] - Eyz[i] - Mg[i] for i in range(G.shape[0]))
    G2 = SP.add_constraints(Gx[i] - h[i] + Eyz[i] + Mg[i] <= bigM*(1-v[i]) for i in range(G.shape[0]))

    sp_soln = SP.solve()

    if sp_soln is None:
        raise RuntimeError("The subproblem is infeasible")

    cy_value = sum(c[i]*y[i].solution_value for i in range(len(c)))
    az_value = sum(a[i]*z[i].solution_value for i in range(len(z)))

    UB = min(UB, cy_value + az_value + sp_soln.objective_value)

    pi_value = [var.solution_value for var in pi]
    g_value = [var.solution_value for var in g]

    print(
        f"iter {k} max |pi|:",
        max(abs(v) for v in pi_value),
        "bigM_pi:",
        bigM_pi
    )

    Ey_sym = matvec(E[:,:len(y)], y)
    Ez_sym = matvec(E[:,len(y):], z)

    Mg_value = matvec(M, g_value)

    cut_rhs = sum(
        (h[i] - Ey_sym[i] - Ez_sym[i] - Mg_value[i]) * pi_value[i]
        for i in range(len(h))
    )

    MP.add_constraint(eta >= cut_rhs)


    print(f"UB: {UB}")
    print(f"LB: {LB}")
    print(f"iteration: {k}")

    k+=1

