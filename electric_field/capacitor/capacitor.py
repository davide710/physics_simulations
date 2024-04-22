import dolfinx
from dolfinx import default_scalar_type
import numpy as np
from ufl import Measure, dot, lhs, grad, rhs, as_vector, sqrt, TrialFunction, TestFunction, curl
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.fem import Function, FunctionSpace, Constant, Expression, locate_dofs_topological, dirichletbc
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
VACUUM = 6

dx = Measure('dx', domain=domain, subdomain_data=cell_tags)
ds = Measure('ds', domain=domain, subdomain_data=facet_tags)

Function_space = FunctionSpace(domain, ("Lagrange", 1))
Vector_space = FunctionSpace(domain, ("DG", 0, (domain.geometry.dim,)))

tdim = domain.topology.dim
facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
dofs = locate_dofs_topological(Function_space, tdim - 1, facets)
bc = dirichletbc(default_scalar_type(0), dofs, Function_space)


u = TrialFunction(Function_space)
v = TestFunction(Function_space)

"""
C = eps * 10^2 / 0.2 = 4.43 nF
e.g. DV := 5 V
=> Q = 22.1 nC, rho = 2.21 nC / m^3 oppure sigma = 0.22 nC / m^3
"""
eps = 8.85e-12
rho = 2.21e-9
sigma = 0.22e-9

F = dot(grad(u), grad(v))*dx(C_POS_VOL) + dot(grad(u), grad(v))*dx(C_NEG_VOL) + dot(grad(u), grad(v))*dx(VACUUM) - rho/eps*v*dx(C_POS_VOL) + rho/eps*v*dx(C_NEG_VOL)

a, L = lhs(F), rhs(F)
problem = LinearProblem(a, L, bcs=[bc])

V = Function(Function_space)
V = problem.solve()

E = Function(Vector_space)
E.interpolate(Expression(grad(V), Vector_space.element.interpolation_points()))

#print(np.max(V.x.array))

for i in range(10):
    z = 0.1 + 0.2/10*i
    print(f'z = {z} --> V = {f_di(V, np.array([5, 5, z]), domain)[0]}')