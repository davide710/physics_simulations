from dolfinx import default_scalar_type
from dolfinx.fem import (dirichletbc, Expression, Function, FunctionSpace, VectorFunctionSpace, locate_dofs_topological)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io.gmshio import model_to_mesh
from dolfinx.mesh import compute_midpoints, locate_entities_boundary
from dolfinx.plot import vtk_mesh
from ufl import TestFunction, TrialFunction, as_vector, dot, dx, grad, inner
from mpi4py import MPI
import gmsh
import numpy as np
import pyvista
from plotting import plot_scalar_function, plot_vector_function


rank = MPI.COMM_WORLD.rank

gmsh.initialize()
r = 0.1   # Radius of copper wires
R = 2     # Radius of domain
N = 1
c_1 = 0.8  # Radius of inner copper wires
c_2 = 1.4  # Radius of outer copper wires
gdim = 2  # Geometric dimension of the mesh
model_rank = 0
mesh_comm = MPI.COMM_WORLD
if mesh_comm.rank == model_rank:

    # Define geometry for background
    background = gmsh.model.occ.addDisk(0, 0, 0, R, R)
    gmsh.model.occ.synchronize()

    # Define the copper-wires inside iron cylinder
    angles_N = [i * 2 * np.pi / N for i in range(N)]
    wires_N = [(2, gmsh.model.occ.addDisk(c_1 * np.cos(v), c_1 * np.sin(v), 0, r, r)) for v in angles_N]

    # Resolve all boundaries of the different wires in the background domain
    all_surfaces = [(2, wires_N)]
    whole_domain = gmsh.model.occ.fragment([(2, background)], all_surfaces)
    gmsh.model.occ.synchronize()
    # Create physical markers for the different wires.
    # We use the following markers:
    # - Vacuum: 0
    # - Iron cylinder: 1
    # - Inner copper wires: $[2,3,\dots,N+1]$
    # - Outer copper wires: $[N+2,\dots, 2\cdot N+1]
    inner_tag = 2
    outer_tag = 2 + N
    background_surfaces = []
    other_surfaces = []
    for domain in whole_domain[0]:
        com = gmsh.model.occ.getCenterOfMass(domain[0], domain[1])
        mass = gmsh.model.occ.getMass(domain[0], domain[1])
        if np.allclose(com, [0, 0, 0]):
            background_surfaces.append(domain[1])

        # Identify the inner circles by their center of mass
        elif np.isclose(np.linalg.norm(com), c_1):
            gmsh.model.addPhysicalGroup(domain[0], [domain[1]], inner_tag)
            inner_tag += 1
            other_surfaces.append(domain)
    # Add marker for the vacuum
    gmsh.model.addPhysicalGroup(2, background_surfaces, tag=0)
    # Create mesh resolution that is fine around the wires and
    # iron cylinder, coarser the further away you get
    gmsh.model.mesh.field.add("Distance", 1)
    edges = gmsh.model.getBoundary(other_surfaces, oriented=False)
    gmsh.model.mesh.field.setNumbers(1, "EdgesList", [e[1] for e in edges])
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "IField", 1)
    gmsh.model.mesh.field.setNumber(2, "LcMin", r / 3)
    gmsh.model.mesh.field.setNumber(2, "LcMax", 6 * r)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 4 * r)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 10 * r)
    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(5)
    # Generate mesh
    gmsh.option.setNumber("Mesh.Algorithm", 7)
    gmsh.model.mesh.generate(gdim)
    gmsh.model.mesh.optimize("Netgen")
    mesh, ct, _ = model_to_mesh(gmsh.model, mesh_comm, model_rank, gdim=2)
    gmsh.finalize()

Q = FunctionSpace(mesh, ("DG", 0))
material_tags = np.unique(ct.values)
mu = Function(Q)
J = Function(Q)
# As we only set some values in J, initialize all as 0
J.x.array[:] = 0
for tag in material_tags:
    cells = ct.find(tag)
    # Set values for mu
    if tag == 0:
        mu_ = 4 * np.pi * 1e-7  # Vacuum
    elif tag == 1:
        mu_ = 1e-5  # Iron (This should really be 6.3e-3)
    else:
        mu_ = 1.26e-6  # Copper
    mu.x.array[cells] = np.full_like(cells, mu_, dtype=default_scalar_type)
    if tag in range(2, 2 + N):
        J.x.array[cells] = np.full_like(cells, 1, dtype=default_scalar_type)
    elif tag in range(2 + N, 2 * N + 2):
        J.x.array[cells] = np.full_like(cells, -1, dtype=default_scalar_type)

V = FunctionSpace(mesh, ("Lagrange", 1))
tdim = mesh.topology.dim
facets = locate_entities_boundary(mesh, tdim - 1, lambda x: np.full(x.shape[1], True))
dofs = locate_dofs_topological(V, tdim - 1, facets)
bc = dirichletbc(default_scalar_type(0), dofs, V)

u = TrialFunction(V)
v = TestFunction(V)
a = (1 / mu) * dot(grad(u), grad(v)) * dx
L = J * v * dx

A_z = Function(V)
problem = LinearProblem(a, L, u=A_z, bcs=[bc])
problem.solve()

W = VectorFunctionSpace(mesh, ("DG", 0))
B = Function(W)
B_expr = Expression(as_vector((A_z.dx(1), -A_z.dx(0))), W.element.interpolation_points())
B.interpolate(B_expr)

plotter = pyvista.Plotter(off_screen=True)
plotter.set_position([0, 0, 5])

grid = pyvista.UnstructuredGrid(*vtk_mesh(mesh, mesh.topology.dim))
num_local_cells = mesh.topology.index_map(mesh.topology.dim).size_local
grid.cell_data["Marker"] = ct.values[ct.indices < num_local_cells]
grid.set_active_scalars("Marker")

# We include ghosts cells as we access all degrees of freedom (including ghosts) on each process
top_imap = mesh.topology.index_map(mesh.topology.dim)
num_cells = top_imap.size_local + top_imap.num_ghosts
midpoints = compute_midpoints(mesh, mesh.topology.dim, range(num_cells))

num_dofs = W.dofmap.index_map.size_local + W.dofmap.index_map.num_ghosts
assert (num_cells == num_dofs)
values = np.zeros((num_dofs, 3), dtype=np.float64)
values[:, :mesh.geometry.dim] = B.x.array.real.reshape(num_dofs, W.dofmap.index_map_bs)
valuesa=values
for q in [0,1]:
    valuesa[:,q]=values[:,q]/np.sqrt(np.power(values[:,0],2)+np.power(values[:,1],2))

cloud = pyvista.PolyData(midpoints[0::100, :])
cloud["B"] = valuesa[0::100, :]
glyphs = cloud.glyph("B", factor=.5)
actor = plotter.add_mesh(grid, style="wireframe", color="k")
actor2 = plotter.add_mesh(glyphs)
#plotter.camera.zoom('tight')
B_fig = plotter.screenshot("B.png")

