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

# rot B = J
# 

r = Function(Function_space)
r.interpolate(lambda x: np.sqrt(x[0]**2 + x[1]**2))

def Jx(x):
    y = x[1, :]
    x = x[0, :]
    jx = -y / np.sqrt(x**2 + y**2)
    jx[np.abs(np.sqrt(x**2+y**2) - 6) > 0.5] = 0
    return jx

def Jy(x):
    y = x[1, :]
    x = x[0, :]
    jy = x / np.sqrt(x**2 + y**2)
    jy[np.abs(np.sqrt(x**2+y**2) - 6) > 0.5] = 0
    return jy

tdim = domain.topology.dim
facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
dofs = locate_dofs_topological(Function_space, tdim - 1, facets)
bc = dirichletbc(default_scalar_type(0), dofs, Function_space)

B_z1 = TrialFunction(Function_space)
B_z2 = TrialFunction(Function_space)
v1 = TestFunction(Function_space)
v2 = TestFunction(Function_space)

jx = Function(Function_space)
jy = Function(Function_space)
jx.interpolate(lambda x: Jx(x))
jy.interpolate(lambda x: Jy(x))

J = Function(Vector_space)
J.interpolate(Expression(as_vector([jx, jy]), Vector_space.element.interpolation_points()))

plot_scalar_function(domain, jx, Function_space, 'jx.png')
plot_scalar_function(domain, jy, Function_space, 'jy.png')
plot_vector_function(domain, J, Vector_space, 'j.png')

F1 = dot(grad(B_z1), grad(v1))/12*10e7*dx(2)-jx*v1*dx(2)

a, L = lhs(F1), rhs(F1)
problem = LinearProblem(a, L, bcs=[bc])

sol1 = Function(Function_space)
sol1 = problem.solve()

plot_scalar_function(domain, sol1, Function_space, 'sol1.png')

F2 = dot(grad(B_z2), grad(v2))/12*10e7*dx-jy*v2*dx

a, L = lhs(F2), rhs(F2)
problem = LinearProblem(a, L, bcs=[bc])

sol2 = Function(Function_space)
sol2 = problem.solve()

plot_scalar_function(domain, sol2, Function_space, 'sol2.png')

#u1 = dot(as_vector([1, 0]), grad(solnUC))
#u2 = dot(as_vector([1, 0]), grad(solnUS))

B = Function(Function_space)
B.interpolate(Expression(sol2.dx(0) - sol1.dx(1), Function_space.element.interpolation_points()))

plot_scalar_function(domain, B, Function_space, 'sol.png')
print(f_di(B, np.array([0,0,0]), domain))
print(f_di(B, np.array([3,3,0]), domain))
