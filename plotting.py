import dolfinx
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pyvista
from dolfinx import plot


def plot_vector_function(mesh, grad_u, Vector_space, fname, show_grid=False):
    top_imap = mesh.topology.index_map(mesh.topology.dim)
    num_cells = top_imap.size_local + top_imap.num_ghosts
    midpoints = dolfinx.mesh.compute_midpoints(mesh, mesh.topology.dim, range(num_cells))
    num_dofs = Vector_space.dofmap.index_map.size_local + Vector_space.dofmap.index_map.num_ghosts
    #assert (num_cells == num_dofs)
    values = np.zeros((num_dofs, 3), dtype=np.float64)
    values[:, :mesh.geometry.dim] = grad_u.x.array.real.reshape(num_dofs, Vector_space.dofmap.index_map_bs)
    valuesa=values
    for q in [0,1]:
        valuesa[:,q]=values[:,q]/np.sqrt(np.power(values[:,0],2)+np.power(values[:,1],2))

    pyvista.start_xvfb()
    topology, cell_types, geometry = plot.vtk_mesh(mesh, mesh.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)

    cloud = pyvista.PolyData(midpoints[0::100, :])
    cloud["grad_u"] = valuesa[0::10, :]
    plotter = pyvista.Plotter(off_screen=True)
    #plotter.camera_position = [(-400, -300, 1000), (-400, -300, 0), (0, 1, 0)]
    glyphs = cloud.glyph("grad_u", factor=1)
    actor = plotter.add_mesh(grid)
    actor2 = plotter.add_mesh(glyphs)
    plotter.camera.zoom(2)
    B_fig = plotter.screenshot(fname)


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


def create_xdmf(mesh, function, filename):
    with dolfinx.io.XDMFFile(mesh.comm, filename, "w") as xdmf:
        xdmf.write_mesh(mesh)
        xdmf.write_function(function)

