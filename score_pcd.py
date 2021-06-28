import argparse
import open3d as o3d

p = argparse.ArgumentParser(description="Optimise and score point clouds")
p.add_argument("pcd1", type=str, help="Point cloud 1")
p.add_argument("pcd2", type=str, help="Point cloud 2")
p.add_argument("-r", "--resolution", type=float, default=0.5, help="Grid resolution")
p.add_argument("-t", "--threshold", type=float, default=0.5, help="Threshold")

args = p.parse_args()

# Voxel size
voxel_size = args.resolution

# Point clouds
pcds = [o3d.io.read_point_cloud(args.pcd1), o3d.io.read_point_cloud(args.pcd2)]

# + ---------------- +
# | Compute features |
# + ---------------- +

radius_normal = voxel_size * 2
radius_feature = voxel_size * 5

fpfhs = []
for pcd in pcds:
    pcd.estimate_normals(
        o3d.geometry.KDTreeSearchParamHybrid(radius=radius_normal, max_nn=30)
    )

    fpfh = o3d.pipelines.registration.compute_fpfh_feature(
        pcd, o3d.geometry.KDTreeSearchParamHybrid(radius=radius_feature, max_nn=100),
    )

    fpfhs.append(fpfh)

# + ------------------- +
# | Global registration |
# + ------------------- +

distance_threshold = voxel_size * 1.5

# RANSAC is non-deterministic; impossible to compare with the original SENSAAS implementation
#   open3d/#288     | https://github.com/intel-isl/Open3D/issues/288
#   open3d/#1263    | https://github.com/intel-isl/Open3D/issues/1263
gresult = o3d.pipelines.registration.registration_ransac_based_on_feature_matching(
    pcds[0],
    pcds[1],
    fpfhs[0],
    fpfhs[1],
    mutual_filter=True,
    max_correspondence_distance=distance_threshold,
    estimation_method=o3d.pipelines.registration.TransformationEstimationPointToPoint(
        False
    ),
    ransac_n=4,
    checkers=[
        o3d.pipelines.registration.CorrespondenceCheckerBasedOnEdgeLength(0.9),
        o3d.pipelines.registration.CorrespondenceCheckerBasedOnDistance(
            distance_threshold
        ),
    ],
    criteria=o3d.pipelines.registration.RANSACConvergenceCriteria(400000, 1000),
)

gfit = o3d.pipelines.registration.evaluate_registration(
    pcds[0], pcds[1], args.threshold, gresult.transformation
)

print(f"GLOBAL | Fitness: {gfit.fitness:.3f} | RMSE: {gfit.inlier_rmse:.3f}")

# + -------------------- +
# | Colored registration |
# + -------------------- +

cresult = o3d.pipelines.registration.registration_colored_icp(
    pcds[0],
    pcds[1],
    max_correspondence_distance=radius_normal,
    init=gresult.transformation,
    estimation_method=o3d.pipelines.registration.TransformationEstimationForColoredICP(
        lambda_geometric=0.8
    ),
    criteria=o3d.pipelines.registration.ICPConvergenceCriteria(
        relative_fitness=1e-6, relative_rmse=1e-6, max_iteration=100
    ),
)

cfit = o3d.pipelines.registration.evaluate_registration(
    pcds[0], pcds[1], args.threshold, cresult.transformation
)

print(f"COLORED | Fitness: {cfit.fitness:.3f} | RMSE: {cfit.inlier_rmse:.3f}")
