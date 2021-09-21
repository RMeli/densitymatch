"""
Fit and score point clouds using SENSAAS-based algorithms.

Changes to the original algorithm:
* Works directly with point clouds
* Does not perform additional scoring

---

BSD 3-Clause License

Copyright (c) 2018-2021, CNRS, Inserm, Université Côte d'Azur, Dominique Douguet and Frédéric Payan.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

import open3d as o3d
import numpy as np


def fit_and_score(pcds, voxel_size, threshold, fast=False):
    """
    Parameters
    ----------
    pcds:
        Tuple (pair) of point clouds
    voxel_size:
        PCD voxel_size
    threshold:
        Distance threshold (maximum correspondence distance)
    """

    assert len(pcds) == 2

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

    if not fast:
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
    else:
        # This is a faster version of global point cloud registration
        gresult = o3d.pipelines.registration.registration_fast_based_on_feature_matching(
            pcds[0],
            pcds[1],
            fpfhs[0],
            fpfhs[1],
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
            lambda_geometric=0.8 # Weight of the geometric objective; (1-lambda) for the photometric objective
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
