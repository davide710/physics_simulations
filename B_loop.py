import dolfinx
import numpy as np
from ufl import Measure, dot, lhs, grad, rhs, as_vector, sqrt, TrialFunction, TestFunction
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.fem import Function, FunctionSpace, Constant, Expression
from plotting import plot_scalar_function, plot_vector_function


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

# rot B = J
# 

r = Function(Function_space)
r.interpolate(lambda x: np.sqrt(x[0]**2 + x[1]**2))

def Jx(x):
    if x[0]**2 + x[1]**2 >= 6 and x[0]**2 + x[1]**2 >= 6.5:
        return -x[1] / np.sqrt(x[0]**2 + x[1]**2)
    else:
        return 0

def Jy(x):
    if x[0]**2 + x[1]**2 >= 6 and x[0]**2 + x[1]**2 >= 6.5:
        return x[0] / np.sqrt(x[0]**2 + x[1]**2)
    else:
        return 0

B_x = TrialFunction(Function_space)
B_y = TrialFunction(Function_space)
v_x = TestFunction(Function_space)
v_y = TestFunction(Function_space)

jx = Function(Function_space)
jy = Function(Function_space)
jx.interpolate(lambda x: Jx(x))
jy.interpolate(lambda x: Jy(x))

plot_scalar_function(domain, jx, Function_space, 'jx.png')
plot_scalar_function(domain, jy, Function_space, 'jy.png')
