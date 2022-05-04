"""
Render a decent representaion of a point cloud.
"""
import open3d as o3d
import numpy as np
import pyvista as pv


def render(pcd, w=2000, h=2000, off_screen=False):
    """
    Parameters
    ----------
    pcd:
        Open3D point cloud

    Returns
    -------
    p:
        PyVista plotter
    """
    pl = pv.Plotter(
        window_size=(w, h),
        off_screen=off_screen,
        multi_samples=8,
        line_smoothing=True,
        polygon_smoothing=True,
    )

    # Get data from point cloud
    centers = np.asarray(pcd.points)
    colors = np.asarray(pcd.colors)

    # Add colored spheres to represent the points
    for i in range(centers.shape[0]):
        sphere = pv.Sphere(
            radius=0.1, center=centers[i, :], phi_resolution=100, theta_resolution=100
        )
        pl.add_mesh(sphere, color=colors[i, :])

    pl.set_background("white")
    pl.enable_anti_aliasing()
    pl.reset_camera()

    return pl


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="Optimise and score point clouds")
    p.add_argument("input", type=str, help="Input point cloud")
    p.add_argument("output", type=str, help="Output image")
    p.add_argument("--height", type=float, default=2000, help="Image height")
    p.add_argument("--width", type=float, default=2000, help="Image width")
    args = p.parse_args()

    pcd = o3d.io.read_point_cloud(args.input)
    plotter = render(pcd, args.height, args.width, off_screen=True)
    plotter.save_graphic(args.output, raster=False)
