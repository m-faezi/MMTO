from image import Image
from dark_frame import DarkFrame

from run import Run
from max_tree import MaxTree
from extractor import Extractor
from mto2lib.utils import base_utils as uts

import os
import glob
import higra as hg
import numpy as np
import mmto


"""Example program - using original settings"""
trees, latitudes, longitudes, fluxes, volumes, gammas, areas, ids = [], [], [], [], [], [], [], []


def mmto_run():

        run = Run()

        run.setup_args()

        for fit_file in glob.glob(os.path.join('./data', '*.fits')):

            run.arguments.file_path = fit_file
            image = Image()
            dark_frame = DarkFrame()

            image.get_image(run)
            image.preprocess_image(run.arguments.s_sigma)

            try:

                dark_frame.estimate_const_bg(image.smooth_image)
                dark_frame.create_reduced_image(image, run.results_dir)

                maxtree = MaxTree()
                maxtree.construct_max_tree(image.reduced_image)
                maxtree.compute_attributes(run, image)

            except Exception as e:

                run.arguments.background_mode = 'morph'

                print(f"Note: Background mode switched from 'const' to '{run.arguments.background_mode}'!")

                maxtree = MaxTree()
                maxtree.construct_max_tree(image.smooth_image)
                maxtree.compute_attributes(run, image)
                dark_frame.estimate_morph_bg(image, maxtree)

            maxtree.detect_significant_objects(dark_frame)
            maxtree.move_up(dark_frame, run)

            extractor = Extractor()
            extractor.create_segmentation(maxtree, image, run)

            tree_of_segments = extractor.maxtree_of_segment
            segment_ids = np.arange(tree_of_segments.num_leaves(), tree_of_segments.num_vertices())
            label_data = np.full(tree_of_segments.num_vertices(), -1, dtype=np.int32)
            #label_data[segment_ids] = np.arange(len(segment_ids))
            #seg_array = hg.reconstruct_leaf_data(tree_of_segments, label_data)
#
            coords_per_segment = [[] for _ in range(len(segment_ids))]
#
            #for y_ in range(seg_array.shape[0]):
#
            #    for x_ in range(seg_array.shape[1]):
#
            #        label = seg_array[y_, x_]
#
            #        if label >= 0:
            #            coords_per_segment[label].append((y_, x_))
#
            #r_eff = [uts.half_light_radius(image.image, coords) for coords in coords_per_segment]

            centroids = [uts.weighted_centroid_coords_from_segments(image.image, coords) for coords in coords_per_segment]

            y = [cen[0] for cen in centroids]
            x = [cen[1] for cen in centroids]

            flux = hg.accumulate_sequential(tree_of_segments, image.image, hg.Accumulators.sum)


            trees.append(extractor.maxtree_of_segment)
            latitudes.append(maxtree.x[extractor.segment_node_map])
            longitudes.append(maxtree.y[extractor.segment_node_map])
            fluxes.append(flux)
            gammas.append(maxtree.gamma[extractor.segment_node_map])
            ids.append(label_data)
            areas.append(maxtree.area[extractor.segment_node_map])
            volumes.append(maxtree.volume[extractor.segment_node_map])


        mmto.tree_map(trees, latitudes, longitudes, fluxes, gammas, areas, volumes, ids)


if __name__ == "__main__":

    mmto_run()

