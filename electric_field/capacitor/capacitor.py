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

domain, cell_tags, facet_tags = gmshio.read_from_msh("electric_field/capacitor/capacitor.msh", MPI.COMM_WORLD, gdim=3)

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
def on_boundary(x):
    return np.isclose(x[0], -5) + np.isclose(x[0], 15) + np.isclose(x[1], -5) + np.isclose(x[1], 15) + np.isclose(x[2], -5) + np.isclose(x[2], 5.4)

dofs = locate_dofs_geometrical(Function_space, on_boundary)
bc = dirichletbc(default_scalar_type(0), dofs, Function_space)


u = TrialFunction(Function_space)
v = TestFunction(Function_space)

"""
C = eps * 10^2 / 0.2 = 4.43 nF
e.g. DV := 5 V
=> Q = 22.1 nC, rho = 2.21 nC / m^3 oppure sigma = 0.22 nC / m^2
"""
eps = 8.85e-12
rho = 2.21e-9
sigma = 0.22e-9

F = dot(grad(u), grad(v))*dx(C_POS_VOL) + dot(grad(u), grad(v))*dx(C_NEG_VOL) + dot(grad(u), grad(v))*dx(VACUUM) - rho/eps*v*dx(C_POS_VOL) + rho/eps*v*dx(C_NEG_VOL)
#F = dot(grad(u), grad(v))*dx(C_POS_VOL) + dot(grad(u), grad(v))*dx(C_NEG_VOL) + dot(grad(u), grad(v))*dx(VACUUM) - sigma/eps*v*dx(C_POS_SUP) + sigma/eps*v*dx(C_NEG_SUP)

a, L = lhs(F), rhs(F)
problem = LinearProblem(a, L, bcs=[bc])

V = Function(Function_space)
V = problem.solve()

E = Function(Vector_space)
E.interpolate(Expression(-grad(V), Vector_space.element.interpolation_points()))


vector_space_x = Vector_space.tabulate_dof_coordinates()[:, 0]
vector_space_y = Vector_space.tabulate_dof_coordinates()[:, 1]

import pandas as pd

df = pd.DataFrame(columns=['x', 'y', 'z', 'V'])
df.x = domain.geometry.x[:, 0]
df.y = domain.geometry.x[:, 1]
df.z = domain.geometry.x[:, 2]
df.V = V.x.array               
df.to_csv('results.txt', index=False)

vector_space_x = Vector_space.tabulate_dof_coordinates()[:, 0]
vector_space_y = Vector_space.tabulate_dof_coordinates()[:, 1]
vector_space_z = Vector_space.tabulate_dof_coordinates()[:, 2]

df = pd.DataFrame(columns=['x', 'y', 'z', 'E', 'EzsuE'])
df.x = vector_space_x
df.y = vector_space_y
df.z = vector_space_z
df.E = np.sqrt(E.x.array[0::3]**2 + E.x.array[1::3]**2 + E.x.array[2::3]**2)
df.EzsuE = 0
df.EzsuE = df.EzsuE + E.x.array[2::3] / np.sqrt(E.x.array[0::3]**2 + E.x.array[1::3]**2 + E.x.array[2::3]**2)
df.to_csv('E.txt', index=False)