import dolfinx
from dolfinx import default_scalar_type
import numpy as np
from ufl import Measure, dot, lhs, grad, rhs, as_vector, sqrt, TrialFunction, TestFunction, curl
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.fem import Function, FunctionSpace, Constant, Expression, locate_dofs_topological, dirichletbc
from plotting import plot_scalar_function, plot_vector_function
from dolfinx.mesh import compute_midpoints, locate_entities_boundary


def f_di(f, x0, mesh):
    tree = dolfinx.geometry.bb_tree(mesh, mesh.geometry.dim)
    cell_candidates = dolfinx.geometry.compute_collisions_points(tree, x0)
    cells = dolfinx.geometry.compute_colliding_cells(mesh, cell_candidates, x0)
    return f.eval(x0, cells[0])

domain, cell_tags, facet_tags = gmshio.read_from_msh("spira.msh", MPI.COMM_WORLD, gdim=2)

dx = Measure('dx', domain=domain, subdomain_data=cell_tags)
ds = Measure('ds', domain=domain, subdomain_data=facet_tags)

Function_space = FunctionSpace(domain, ("Lagrange", 1))
Vector_space = FunctionSpace(domain, ("DG", 0, (domain.geometry.dim,)))

def Jx(x):
    y = x[1, :]
    x = x[0, :]
    jx = -y / (np.sqrt(x**2 + y**2) + 1e-9)
    jx[np.abs(np.sqrt(x**2+y**2) - 6) > 0.1] = 0
    return jx

def Jy(x):
    y = x[1, :]
    x = x[0, :]
    jy = x / (np.sqrt(x**2 + y**2) + 1e-9)
    jy[np.abs(np.sqrt(x**2+y**2) - 6) > 0.1] = 0
    return jy

tdim = domain.topology.dim
facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
dofs = locate_dofs_topological(Function_space, tdim - 1, facets)
bc = dirichletbc(default_scalar_type(0), dofs, Function_space)

"""
(dB/dy, -dB/dx) = J
B = dAy/dx - dAx/dy, A=(Ax, Ay)
J = (d2/dxdy Ay - d2/dy2 Ax, - d2/dx2 Ay + d2/dxdy Ax)

I:    || Jx*v1 dxdy = - || (d/dy Ax)*(d/dy v1) dxdy
II:    || Jy*v2 dxdy = - || (d/dx Ay)*(d/dx v2) dxdy
"""