"""We improved parameters in this version"""
from image import Image
from dark_frame import DarkFrame
from run import Run
from max_tree import MaxTree
from extractor import Extractor

from mto2lib.utils import io_utils
import sys

import os
import glob
import higra as hg
import numpy as np
import mmto
import helper
import background

"""Example program - using original settings"""
trees, latitudes, longitudes, mu_list, graphs, depths, ids, areas, reffs = [], [], [], [], [], [], [], [], []


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


            trees.append(maxtree.tree_structure)
            latitudes.append(maxtree.x)
            longitudes.append(maxtree.y)
            mu_list.append(maxtree.volume)
            graphs.append(maxtree.graph)
            depths.append(maxtree.gamma)
            ids.append(extractor.ids)
            areas.append(maxtree.area)
            reffs.append(maxtree.area) #change the value to reff

        tree_map = mmto.tree_map(*trees, *mu_list, *latitudes, *longitudes, *depths, *ids, *areas, *reffs)



if __name__ == "__main__":

    mmto_run()

