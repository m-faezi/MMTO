import higra as hg
import numpy as np
from PIL import Image, ImageOps
from mto2lib.utils import io_utils as io_uts
import os
from random import randint

def get_segmentation_map(tree_structure, modified_isophote, header, arguments):

    tree_of_segments, n_map_segments = hg.simplify_tree(tree_structure, np.logical_not(modified_isophote))
    colors = np.random.randint(0, 254, (tree_of_segments.num_vertices(), 3), dtype=np.uint8)
    colors[tree_of_segments.root(), :] = 0

    unique_segment_ids = np.arange(tree_of_segments.num_vertices())[::-1]
    seg_with_ids = hg.reconstruct_leaf_data(tree_of_segments, unique_segment_ids)

    results_dir = os.path.join("./results", arguments.time_stamp)
    output_fits = os.path.join(results_dir, str(randint(4000, 5000))+"_segmentation_map.fits")


    io_uts.save_fits_with_header(seg_with_ids, header, output_fits)



    return tree_of_segments, n_map_segments, unique_segment_ids

