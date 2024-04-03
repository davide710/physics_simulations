import dolfinx
from dolfinx import default_scalar_type
import numpy as np
from ufl import Measure, dot, lhs, grad, rhs, as_vector, sqrt, TrialFunction, TestFunction, curl
from dolfinx.io import gmshio
from mpi4py import MPI
from dolfinx.fem.petsc import LinearProblem
from dolfinx.fem import Function, FunctionSpace, Constant, Expression, locate_dofs_topological, dirichletbc
#from plotting import plot_scalar_function, plot_vector_function
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

def jx(x):
    y = x[1, :]
    x = x[0, :]
    Jx = -y / (np.sqrt(x**2 + y**2) + 1e-9)
    Jx[np.abs(np.sqrt(x**2+y**2) - 6) > 0.1] = 0
    return Jx

def jy(x):
    y = x[1, :]
    x = x[0, :]
    Jy = x / (np.sqrt(x**2 + y**2) + 1e-9)
    Jy[np.abs(np.sqrt(x**2+y**2) - 6) > 0.1] = 0
    return Jy

tdim = domain.topology.dim
facets = locate_entities_boundary(domain, tdim - 1, lambda x: np.full(x.shape[1], True))
dofs = locate_dofs_topological(Function_space, tdim - 1, facets)
bc = dirichletbc(default_scalar_type(0), dofs, Function_space)

"""
r = Function(Function_space)
r.interpolate(lambda x: np.sqrt(x[0]**2 + x[1]**2))
x_coord = Function(Function_space)
x_coord.interpolate(lambda x: x[0])
y_coord = Function(Function_space)
y_coord.interpolate(lambda x: x[1])
"""

Jx = Function(Function_space)
Jy = Function(Function_space)

Jx.interpolate(lambda x: jx(x))
Jy.interpolate(lambda x: jy(x))

J = Function(Vector_space)
J.interpolate(Expression(as_vector([Jx, Jy]), Vector_space.element.interpolation_points()))

#plot_scalar_function(domain, Jx, Function_space, 'jx.png')
#plot_scalar_function(domain, Jy, Function_space, 'jy.png')
#plot_vector_function(domain, J, Vector_space, 'j.png')

"""
(dB/dy, -dB/dx) = J
B = dAy/dx - dAx/dy, A=(Ax, Ay)
J = (d2/dxdy Ay - d2/dy2 Ax, - d2/dx2 Ay + d2/dxdy Ax)

I:    || Jx*v1 dxdy = || (d/dy Ax)*(d/dy v1) dxdy
II:    || Jy*v2 dxdy = || (d/dx Ay)*(d/dx v2) dxdy
"""

ux = TrialFunction(Function_space)
uy = TrialFunction(Function_space)
v1 = TestFunction(Function_space)
v2 = TestFunction(Function_space)

F1 = ux.dx(1)*v1.dx(1)*dx(2)-Jx*v1*dx(2)
a1, L1 = lhs(F1), rhs(F1)
problem1 = LinearProblem(a1, L1, bcs=[bc])

Ax = Function(Function_space)
Ax = problem1.solve()

F2 = uy.dx(0)*v2.dx(0)*dx(2)-Jy*v2*dx(2)
a2, L2 = lhs(F2), rhs(F2)
problem2 = LinearProblem(a2, L2, bcs=[bc])

Ay = Function(Function_space)
Ay = problem2.solve()

B = Function(Function_space)
B.interpolate(Expression(Ay.dx(0) - Ax.dx(1), Function_space.element.interpolation_points()))

#plot_scalar_function(domain, B, Function_space, 'B.png')


x = np.linspace(0, 5.5, 1000)
y = 0*x
r = np.sqrt(x**2 + y**2)
P = [[x[i], y[i], 0] for i in range(x.shape[0])]
z = [f_di(B, P[i], domain) for i in range(x.shape[0])]

with open('results.csv', 'w') as f:
    for i in range(x.shape[0]):
        f.write(f'{r[i]},{z[i]}\n')