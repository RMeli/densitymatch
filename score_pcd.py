import open3d as o3d


def fit_and_score(pcds, voxel_size, threshold):
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
            pcd,
            o3d.geometry.KDTreeSearchParamHybrid(radius=radius_feature, max_nn=100),
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
        pcds[0], pcds[1], threshold, gresult.transformation
    )

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
        pcds[0], pcds[1], threshold, cresult.transformation
    )

    # The 4x4 transformation matrix represents and affine transformation:
    #   r00 r01 r02 t0
    #   r10 r11 r12 t1
    #   r20 r21 r22 t2
    #   0   0   0   1

    return gfit, cfit, cresult.transformation


if __name__ == "__main__":
    import argparse
    import os

    p = argparse.ArgumentParser(description="Optimise and score point clouds")
    p.add_argument("pcd1", type=str, help="Point cloud 1")
    p.add_argument("pcd2", type=str, help="Point cloud 2")
    p.add_argument(
        "-r", "--resolution", type=float, default=0.5, help="Grid resolution"
    )
    p.add_argument("-t", "--threshold", type=float, default=0.5, help="Threshold")
    p.add_argument("-v", "--verbose", action="store_true", help="Verbose")

    args = p.parse_args()

    # Check files exist
    for f in [args.pcd1, args.pcd2]:
        if not os.path.exists(f):
            raise FileNotFoundError(f"File {f} not found.")

    # Point clouds
    pcds = [o3d.io.read_point_cloud(args.pcd1), o3d.io.read_point_cloud(args.pcd2)]

    gfit, cfit, transformation = fit_and_score(
        pcds, voxel_size=args.resolution, threshold=args.threshold
    )

    if args.verbose:
        print(f"Transformation:\n{transformation}")

    print(f"GLOBAL | Fitness: {gfit.fitness:.3f} | RMSE: {gfit.inlier_rmse:.3f}")
    print(f"COLORED | Fitness: {cfit.fitness:.3f} | RMSE: {cfit.inlier_rmse:.3f}")