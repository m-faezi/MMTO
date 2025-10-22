"""We improved parameters in this version"""

import os
import glob
import higra as hg
import numpy as np
import mmto
import helper
import background

"""Example program - using original settings"""
trees, latitudes, longitudes, mu_list, graphs, depths, ids, areas, reffs = [], [], [], [], [], [], [], [], []

for fit_file in glob.glob(os.path.join('./data', '*.fits')):

    image, header = helper.read_image_data(fit_file)  # use fit_file here
    image = helper.image_value_check(image)
    image_processed = helper.smooth_filter(image, 2)

    bg_mean, bg_var, bg_gain = background.estimate_background(image_processed)
    image_calibrated = image_processed - bg_mean

    graph_structure, tree_structure, altitudes = helper.image_to_hierarchical_structure(image_calibrated)

    x, y = helper.centroid(tree_structure, image.shape[:2])
    distances = np.sqrt((x[tree_structure.parents()] - x) ** 2 + (y[tree_structure.parents()] - y) ** 2)
    mean, variance = hg.attribute_gaussian_region_weights_model(tree_structure, image_calibrated)
    area = hg.attribute_area(tree_structure)
    parent_area = area[tree_structure.parents()]

    gaussian_intensities = helper.compute_gaussian_profile(
        mean,
        variance,
        distances,
        altitudes / area
    )

    volume = hg.attribute_volume(tree_structure, altitudes)
    depth = hg.attribute_depth(tree_structure)
    first_hu_moment = hg.attribute_moment_of_inertia(tree_structure, graph_structure)
    parent_altitude = altitudes[tree_structure.parents()]
    gamma = hg.attribute_topological_height(tree_structure)
    parent_gamma = gamma[tree_structure.parents()]

    significant_nodes = helper.attribute_statistical_significance(
        tree_structure, altitudes, volume, area, bg_var, bg_gain
    )
    objects = helper.select_objects(tree_structure, significant_nodes)

    modified_isophote = helper.move_up(
        tree_structure, altitudes, area, parent_area, distances, objects, bg_var, bg_gain,
        parent_gamma - gamma, gaussian_intensities, 0, None, .0
    )

    tree_of_segments, n_map_segments = hg.simplify_tree(tree_structure, np.logical_not(modified_isophote))

    segment_ids = np.arange(tree_of_segments.num_leaves(), tree_of_segments.num_vertices())
    label_data = np.full(tree_of_segments.num_vertices(), -1, dtype=np.int32)
    label_data[segment_ids] = np.arange(len(segment_ids))
    seg_array = hg.reconstruct_leaf_data(tree_of_segments, label_data)
    unique_segment_ids = np.arange(tree_of_segments.num_vertices())[::-1]

    coords_per_segment = [[] for _ in range(len(segment_ids))]

    for y_ in range(seg_array.shape[0]):
        for x_ in range(seg_array.shape[1]):
            label = seg_array[y_, x_]
            if label >= 0:
                coords_per_segment[label].append((y_, x_))

    r_eff = [helper.half_light_radius(image, coords) for coords in coords_per_segment]
    r_fwhm = [helper.compute_r_fwhm(image, coords) for coords in coords_per_segment]

    centroids = [helper.weighted_centroid_coords_from_segments(image, coords) for coords in coords_per_segment]

    y = [cen[0] for cen in centroids]
    x = [cen[1] for cen in centroids]
    ra, dec = helper.sky_coordinates(y, x, header)

    a, b, theta = helper.second_order_moments(tree_of_segments, image.shape[:2], image)

    trees.append(tree_of_segments)
    latitudes.append(x[::-1])
    longitudes.append(y[::-1])
    mu_list.append(np.log10(volume[tree_of_segments.num_leaves():][::-1]))
    graphs.append(graph_structure)
    depths.append(depth[tree_of_segments.num_leaves():][::-1])
    ids.append(unique_segment_ids[tree_of_segments.num_leaves():][::-1])
    areas.append(area[n_map_segments][tree_of_segments.num_leaves():][::-1])
    reffs.append(r_eff[tree_of_segments.num_leaves():][::-1])

tree_map = mmto.tree_map(*trees, *mu_list, *latitudes, *longitudes, *depths, *ids, *areas, *reffs)


