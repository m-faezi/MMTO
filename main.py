from mto2lib.image import Image
from mto2lib.dark_frame import DarkFrame

from mto2lib.run import Run
from mto2lib.max_tree import MaxTree
from mto2lib.extractor import Extractor

import os
import glob
import higra as hg
import numpy as np
import mmto
import uuid


trees, latitudes, longitudes, fluxes, volumes, gammas, areas, ids, tree_ids = [], [], [], [], [], [], [], [], []


def mmto_run():

    run = Run()

    run.setup_args()

    for fit_file in glob.glob(os.path.join('data', '*.fits')):

        tree_id = str(uuid.uuid4())[:8]
        tree_ids.append(tree_id)

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
        extractor.create_segmentation(maxtree, image, run, tree_id)

        tree_of_segments = extractor.maxtree_of_segment
        segment_ids = np.arange(tree_of_segments.num_leaves(), tree_of_segments.num_vertices())
        label_data = np.full(tree_of_segments.num_vertices(), -1, dtype=np.int32)
        label_data[segment_ids] = np.arange(len(segment_ids))

        flux = hg.accumulate_sequential(tree_of_segments, image.image, hg.Accumulators.sum)

        trees.append(extractor.maxtree_of_segment)
        latitudes.append(maxtree.x[extractor.segment_node_map])
        longitudes.append(maxtree.y[extractor.segment_node_map])
        fluxes.append(flux)
        gammas.append(maxtree.gamma[extractor.segment_node_map])
        ids.append(label_data)
        areas.append(maxtree.area[extractor.segment_node_map])
        volumes.append(maxtree.volume[extractor.segment_node_map])

    mmto.tree_map(trees, latitudes, longitudes, fluxes, gammas, areas, volumes, ids, tree_ids)


if __name__ == "__main__":

    mmto_run()

