import dolfinx
import numpy as np
from ufl import Measure, dot, lhs, grad, rhs, FacetNormal, FiniteElement, TrialFunctions, TestFunctions, MixedElement, as_vector, sqrt, atan2, conditional
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

domain = dolfinx.mesh.create_rectangle(MPI.COMM_WORLD, [np.array([-2, -2]), np.array([2, 2])],[100, 100], dolfinx.mesh.CellType.triangle)
Function_space = FunctionSpace(domain, ("Lagrange", 1))
Vector_space = FunctionSpace(domain, ("DG", 0, (domain.geometry.dim,)))

r = Function(Function_space)
r.interpolate(lambda x: np.sqrt(x[0]**2 + x[1]**2))

u_theta = Function(Vector_space)
u_theta.interpolate(lambda x: np.array([-x[1], x[0]]) / np.sqrt(x[0]**2 + x[1]**2))

mu0 = 4*3.14e-7
B = Function(Vector_space)
B_expr = u_theta / r
B.interpolate(Expression(B_expr, Vector_space.element.interpolation_points()))

plot_scalar_function(domain, r, Function_space, 'r.png')
plot_vector_function(domain, u_theta, Vector_space, 'u_theta.png')
plot_vector_function(domain, B, Vector_space, 'B.png')

