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

#Cplex doesn't support Numpy function, '@', so make a helper method

def matvec(mat, variables):
    return [
        sum(mat[i,j]*variables[j] for j in range(mat.shape[1])) for i in range(mat.shape[0])
    ]

def dotvec(vec1, vec2):
    if vec1.shape[0] != len(vec2):
        raise ValueError("They're dimensions, don't match")
    else:
        return sum(vec1[i]*vec2[i] for i in range(vec1.shape[0]))
    

#Add constraint and minimize MP

Ay = matvec(A[:,:len(y)], y)
Az = matvec(A[:,len(y):], z)
MP.add_constraints(Ay[i] + Az[i] >= b[i] for i in range(A.shape[0]))

MP.minimize(dotvec(c, y) + dotvec(a, z) + eta)
mp_soln = MP.solve()

if mp_soln is None:
    raise RuntimeError("Master problem is infeasible or could not be solved")

#MP_obj = MP_soln.objective_value
LB = max(LB, mp_soln.objective_value)

#create SP

SP = Model(name = "SP")

SP.parameters.mip.tolerances.mipgap = 1e-8
SP.parameters.mip.tolerances.integrality = 1e-9

x = SP.continuous_var_list(G.shape[1], lb = 0, name = 'x')
pi = SP.continuous_var_list(G.shape[0], lb = 0, name = 'pi')
g = SP.continuous_var_list(M.shape[1], lb = 0, ub = 1, name = 'g')
v = SP.binary_var_list(G.shape[0], name = 'v')
w = SP.binary_var_list(G.shape[1], name = 'w')

yz = np.concatenate([
    [var.solution_value for var in y],
    [var.solution_value for var in z]
])

Gx = matvec(G,x)
Eyz = matvec(E,yz)
Mg = matvec(M,g)
GTpi = matvec(G.T,pi)

#Primal feasibility

G1 = SP.add_constraints(Gx[i] >= h[i] - Eyz[i] - Mg[i] for i in range(G.shape[0]))

#Dual feasibility

SP.add_constraints(GTpi[i] <= d[i] for i in range(G.shape[1]))

#First complementarity pair

SP.add_constraints(pi[j] <= bigM*v[j] for j in range(len(pi)))
G2 = SP.add_constraints([Gx[i] - h[i] + Eyz[i] + Mg[i] <= bigM*(1-v[i]) for i in range(G.shape[0])])

#Second complementarity pair

SP.add_constraints(x[j] <= bigM*w[j] for j in range(len(x)))
SP.add_constraints(d[j] - GTpi[j] <= bigM*(1-w[j]) for j in range(len(d)))

#Uncertainty set (Bertsimas-Sim style "budget of uncertianty")

sub_size, sub_budget, total_budget = inst['sub_size'], inst['sub_budget'], inst['total_budget']
SP.add_constraint(sum(g[:sub_size]) <= sub_budget)
SP.add_constraint(sum(g) <= total_budget)

dx = dotvec(d,x)
SP.maximize(dx)

sp_soln = SP.solve()

if sp_soln is None:
    SP_obj = float('inf')
else:
    SP_obj = sp_soln.objective_value
    cy_value = sum(c[i]*y[i].solution_value for i in range(len(c)))
    az_value = sum(a[i]*z[i].solution_value for i in range(len(z)))

    UB = min(UB, cy_value + az_value + SP_obj)

#Start Iterative loop

max_iterations =  100

Ey = matvec(E[:, :len(y)], y)
Ez = matvec(E[:, len(y):], z)

while abs(UB-LB) >= epsilon and k <= max_iterations:
    if SP_obj < float('inf'):
        #Create x^{k+1}
        x_new = MP.continuous_var_list(G.shape[1], lb = 0, name = f'x_new_{k}')
        #Add constraint (12)
        MP.add_constraint(eta >= dotvec(d, x_new))

        #Add constraint (13)
        g_values = [var.solution_value for var in g]

        Gx_new = matvec(G, x_new)
        Mg_values = matvec(M, g_values)

        MP.add_constraints(Ey[i] + Ez[i] + Gx_new[i] >= h[i] - Mg_values[i] for i in range(G.shape[0]))

        MP_soln = MP.solve()
        if MP_soln is None:
            raise RuntimeError("Updated master problem did not solve")

        MP_obj = MP_soln.objective_value
        LB = max(LB, MP_obj)
    
    else:
         #Create x^{k+1}
        x_new = MP.continuous_var_list(G.shape[1], lb = 0, name = f'x_new_{k}')
        
        #Add constraint (14)
        g_values = [var.solution_value for var in g]

        Gx_new = matvec(G, x_new)
        Mg_values = matvec(M, g_values)

        MP.add_constraints(Ey[i] + Ez[i] + Gx_new[i] >= h[i] - Mg_values[i] for i in range(G.shape[0]))

        MP_soln = MP.solve()
        if MP_soln is None:
            raise RuntimeError("Updated master problem did not solve")

        MP_obj = MP_soln.objective_value
        LB = max(LB, MP_obj)
    
    #update the SP constraints in response to MP soln
    SP.remove_constraints(G1)
    SP.remove_constraints(G2)

    yz = np.concatenate([
        [var.solution_value for var in y],
        [var.solution_value for var in z]
    ])

    Eyz = matvec(E, yz)

    G1 = SP.add_constraints(Gx[i] >= h[i] - Eyz[i] - Mg[i] for i in range(G.shape[0]))
    G2 = SP.add_constraints([Gx[i] - h[i] + Eyz[i] + Mg[i] <= bigM*(1-v[i]) for i in range(G.shape[0])])

    SP_soln = SP.solve()
    if SP_soln is None:
        SP_obj = float('inf')
    else:
        SP_obj = SP_soln.objective_value

        #print("g:", [var.solution_value for var in g])

        cy_value = sum(c[i]*y[i].solution_value for i in range(len(c)))
        az_value = sum(a[i]*z[i].solution_value for i in range(len(z)))

#        print("Opening cost:", cy_value)
#        print("Capacity cost:", az_value)
#        print("Worst-case routing cost:", SP_obj)
#        print("UB candidate:", cy_value + az_value + SP_obj)
#        print("MP objective:", MP_obj)
#        print("Eta:", eta.solution_value)
        # print("y: ", [var.solution_value for var in y])
        # print("x: ", [var.solution_value for var in x])
        # print("z: ", [var.solution_value for var in z])

        UB = min(UB, cy_value + az_value + SP_obj)


    print(f"Number of Iterations: {k}")
    print(f"Upper Bound: {UB}")
    print(f"Lower Bound: {LB}")

    k += 1
