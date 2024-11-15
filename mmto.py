"""We improved parameters in this version"""

import mtolib.main as mto
import os
import glob
import higra as hg
import numpy as np
from PIL import Image
import mdmto
import operator
from helpers import parameter_helper as ph


"""Example program - using original settings"""
trees, latitudes, longitudes, mu_list, graphs, depths, ids, areas, moments = [], [], [], [], [], [], [], [], []

for fit_file in glob.glob(os.path.join('./data', '*.fits')):

    # Get the input image and parameters
    image, params = mto.setup(fit_file)

    # Pre-process the image
    processed_image = mto.preprocess_image(image, params, n=2)

    # Build a max tree
    mt = mto.build_max_tree(processed_image, params)

    # Filter the tree and find objects
    id_map, sig_ancs = mto.filter_tree(mt, processed_image, params)

    # Relabel objects for clearer visualisation
    id_map = mto.relabel_segments(id_map, shuffle_labels=False)

    '''Revaluing id_map to construct max-tree'''
    old_id_map = id_map
    value, counts = np.unique(old_id_map, return_counts=True)
    tuple_value_count = list(zip(value.tolist(), counts.tolist()))
    tuple_value_count.sort(key=operator.itemgetter(1))
    tuple_value_count.reverse()
    id_map = np.ones(image.shape)
    id_map = -id_map

    for i in tuple_value_count:
        if i[0] == -1:
            continue
        else:
            _id = (old_id_map == i[0])
            id_map[_id] = tuple_value_count.index(i)
    id_map = id_map.astype(int)
    # end of revaluing

    # Generate output files
    mto.generate_image(image, id_map, params)
    para = mto.generate_parameters(image, id_map, sig_ancs, params)

    graph = hg.get_8_adjacency_graph(image.shape)
    t, a = hg.component_tree_max_tree(graph, id_map)
    depth = hg.attribute_depth(t)
    area = hg.attribute_area(t)
    first_hu_moment = hg.attribute_moment_of_inertia(t, graph)

    # divide area to mean(area) of leaves
    area[t.num_leaves(): -1] = area[t.num_leaves(): -1] / np.mean(area[t.num_leaves(): -1])

    new_mu = ph.total_flux(t, image.shape, processed_image)
    mu = [-1] * (len(a) - len(new_mu))
    new_mu = mu + new_mu.tolist()

    new_x, new_y = ph.centroid(t, image.shape)
    x = [-1] * (len(a) - len(new_x))
    new_x = x + new_x.tolist()
    y = [-1] * (len(a) - len(new_y))
    new_y = y + new_y.tolist()

    trees.append(t)
    latitudes.append(new_x)
    longitudes.append(new_y)
    mu_list.append(new_mu)
    graphs.append(graph)
    depths.append(depth)
    ids.append(a)
    areas.append(area)
    moments.append(first_hu_moment)

tree_map = mdmto.tree_map(*trees, *mu_list, *latitudes, *longitudes, *depths, *ids, *areas, *moments)

fusion_tree, fusion_altitude = hg.component_tree_max_tree(
    graphs[0],
    tree_map
)

colors = np.random.randint(
    0,
    256,
    (fusion_tree.num_vertices(), 3),
    dtype=np.uint8,
)
colors[fusion_tree.root(), :] = 0
seg = hg.reconstruct_leaf_data(fusion_tree, colors)
Image.fromarray(np.flipud(seg)).save("mdmto-out.png")


