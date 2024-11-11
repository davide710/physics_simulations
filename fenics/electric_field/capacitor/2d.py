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

domain, cell_tags, facet_tags = gmshio.read_from_msh("electric_field/capacitor/2d.msh", MPI.COMM_WORLD, gdim=2)

C_POS = 1
C_NEG = 2
VACUUM = 3

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

eps = 8.85e-12
#f = Constant(domain, dolfinx.default_scalar_type(0))

F = dot(grad(u), grad(v))*dx(VACUUM) + dot(grad(u), grad(v))*dx(C_POS) + dot(grad(u), grad(v))*dx(C_NEG) - 5e-6*v*dx(C_POS) + 5e-6*v*dx(C_NEG)

a, L = lhs(F), rhs(F)
problem = LinearProblem(a, L, bcs=[bc])

V = Function(Function_space)
V = problem.solve()

x = np.linspace(5.2, 5.8, 10)
y = [f_di(V, np.array([10, i, 0]), domain)[0] for i in x]
with open('graph.txt', 'w') as f:
    for i in range(len(x)):
        f.write(f'{x[i]} {y[i]}\n') # perfectly linear as it should be


import dolfinx
import numpy as np
import pyvista
from dolfinx import plot

def plot_scalar_function(mesh, function, space, filename):
    pyvista.start_xvfb()
    topology, cell_types, geometry = plot.vtk_mesh(mesh, mesh.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    u_topology, u_cell_types, u_geometry = plot.vtk_mesh(space)
    u_grid = pyvista.UnstructuredGrid(u_topology, u_cell_types, u_geometry)
    u_grid.point_data["u"] = function.x.array.real
    u_grid.set_active_scalars("u")
    u_plotter = pyvista.Plotter(off_screen=True)
    u_plotter.add_mesh(u_grid, show_edges=False)
    u_plotter.view_xy()
    figure = u_plotter.screenshot(filename)

plot_scalar_function(domain, V, Function_space, 'V.png')



