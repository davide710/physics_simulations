import dolfinx
from dolfinx import default_scalar_type
import numpy as np
from ufl import Measure, dot, lhs, grad, rhs, as_vector, sqrt, TrialFunction, TestFunction, curl, div
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.fem import Function, FunctionSpace, Constant, Expression, locate_dofs_topological, dirichletbc, locate_dofs_geometrical
from dolfinx.mesh import compute_midpoints, locate_entities_boundary


BOX = 1
VACUUM = 2

def f_di(f, x0, mesh):
    tree = dolfinx.geometry.bb_tree(mesh, mesh.geometry.dim)
    cell_candidates = dolfinx.geometry.compute_collisions_points(tree, x0)
    cells = dolfinx.geometry.compute_colliding_cells(mesh, cell_candidates, x0)
    return f.eval(x0, cells[0])

domain, cell_tags, facet_tags = gmshio.read_from_msh("electric_field/sphere/boxes.msh", MPI.COMM_WORLD, gdim=3)

dx = Measure('dx', domain=domain, subdomain_data=cell_tags)
ds = Measure('ds', domain=domain, subdomain_data=facet_tags)

Function_space = FunctionSpace(domain, ("Lagrange", 1))
Vector_space = FunctionSpace(domain, ("DG", 0, (domain.geometry.dim,)))

"""
(lap V) * v = - rho / eps * v
div (v * grad V) - dot(grad V, grad v) = - rho / eps * v
dot(grad V, grad v) - bc = rho / eps * v
"""

eps = 8.85e-12

#tdim = domain.topology.dim
#facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
#dofs = locate_dofs_topological(Function_space, tdim - 1, facets)
#bc = dirichletbc(default_scalar_type(450000), dofs, Function_space)

def on_boundary(x):
    condition1 = (x[0] > 5.95)
    condition2 = (x[0] < -4.95)
    condition3 = (x[1] > 5.95)
    condition4 = (x[1] < -4.95)
    condition5 = (x[2] > 5.95)
    condition6 = (x[2] < -4.95)
    return np.any(np.vstack([condition1, condition2, condition3, condition4, condition5, condition6]), axis=0)

dofs = locate_dofs_geometrical(Function_space, on_boundary)
bc = dirichletbc(default_scalar_type(0), dofs, Function_space)


rho = 1e-9
zero = Constant(domain, default_scalar_type(0))

u = TrialFunction(Function_space)
v = TestFunction(Function_space)

F = dot(grad(u), grad(v))*dx(BOX) + dot(grad(u), grad(v))*dx(VACUUM) - rho/eps*v*dx(BOX) + zero/eps*v*dx(VACUUM)

a, L = lhs(F), rhs(F)
problem = LinearProblem(a, L, bcs=[bc])

V = Function(Function_space)
V = problem.solve()

with open('box.csv', 'w') as f:
    #f.write(f'r,V\n')
    for r in np.arange(1, 4, 0.1):
        f.write(f'{r} {f_di(V, np.array([0.5, r, 0.5]), domain)[0]}\n')

