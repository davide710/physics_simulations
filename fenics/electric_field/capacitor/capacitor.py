import dolfinx
from dolfinx import default_scalar_type
import numpy as np
from ufl import Measure, dot, lhs, grad, rhs, as_vector, sqrt, TrialFunction, TestFunction, curl
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.fem import Function, FunctionSpace, Constant, Expression, locate_dofs_topological, locate_dofs_geometrical, dirichletbc
from dolfinx.mesh import compute_midpoints, locate_entities_boundary


def f_di(f, x0, mesh):
    tree = dolfinx.geometry.bb_tree(mesh, mesh.geometry.dim)
    cell_candidates = dolfinx.geometry.compute_collisions_points(tree, x0)
    cells = dolfinx.geometry.compute_colliding_cells(mesh, cell_candidates, x0)
    return f.eval(x0, cells[0])

domain, cell_tags, facet_tags = gmshio.read_from_msh("electric_field/capacitor/extrude.msh", MPI.COMM_WORLD, gdim=3)

C_POS_SUP = 1
C_NEG_SUP = 2
C_POS_VOL = 3
C_NEG_VOL = 4
EXT_SUP = 5
VACUUM = 6

dx = Measure('dx', domain=domain, subdomain_data=cell_tags)
ds = Measure('ds', domain=domain, subdomain_data=facet_tags)

Function_space = FunctionSpace(domain, ("Lagrange", 1))
Vector_space = FunctionSpace(domain, ("DG", 0, (domain.geometry.dim,)))

#tdim = domain.topology.dim
#facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
#dofs = locate_dofs_topological(Function_space, tdim - 1, facets)
#def on_boundary(x):
#    return np.isclose(x[0], -5) + np.isclose(x[0], 15) + np.isclose(x[1], -5) + np.isclose(x[1], 15) + np.isclose(x[2], -5) + np.isclose(x[2], 5.4)
#
#dofs = locate_dofs_geometrical(Function_space, on_boundary)
#bc = dirichletbc(default_scalar_type(0), dofs, Function_space)


u = TrialFunction(Function_space)
v = TestFunction(Function_space)

"""
C = eps * 10*5 / 0.2 = 7.38e-10 F
e.g. DV := 30 V
=> Q = 2.21e-8 C, rho = 2.21 nC / m^3
"""
eps = 8.85e-12
rho = 2.21e-9

F = dot(grad(u), grad(v))*dx(4) + dot(grad(u), grad(v))*dx(5) + dot(grad(u), grad(v))*dx(6) - rho/eps*v*dx(4) + rho/eps*v*dx(5)

a, L = lhs(F), rhs(F)
problem = LinearProblem(a, L)#, bcs=[bc])

V = Function(Function_space)
V = problem.solve()

import pandas as pd

df = pd.DataFrame(columns=['x', 'y', 'z', 'V'])
df.x = domain.geometry.x[:, 0]
df.y = domain.geometry.x[:, 1]
df.z = domain.geometry.x[:, 2]
df.V = V.x.array               
df.to_csv('results.txt', index=False) # as expected it's linear, and V(10, 5.8, 2.5) - V(10, 5.2, 2.5) = 28.9, close enough to the expected 30