from mto2lib.image import Image
from mto2lib.dark_frame import DarkFrame
from mto2lib.run import Run
from mto2lib.max_tree import MaxTree
from mto2lib.extractor import Extractor
import higra as hg
import numpy as np
import mmto
import uuid
from mto2lib.utils import io_utils, base_utils
import sys


trees, latitudes, longitudes, fluxes, volumes, gammas, areas, ids, tree_ids = [], [], [], [], [], [], [], [], []


def mmto_run():

    run = Run()
    run.setup_args()

    for band_name, band_args in run.config['bands'].items():

        print(f"Processing band: {band_name}")
        _id = str(uuid.uuid4())[:8]
        fit_file = band_args['file_path']

        image = Image()

        dark_frame = DarkFrame()

        try:

            image.get_image(fit_file)
            image.preprocess_image(band_args['s_sigma'])

            if band_args['background_mode'] == 'const':

                try:

                    dark_frame.estimate_const_bg(image.smooth_image)


                    if band_args['skip_reduction']:

                        dark_frame.bg_map = np.full_like(image.image, 0, dtype=np.float32)

                    dark_frame.create_reduced_image(image)

                    maxtree = MaxTree()
                    maxtree.construct_max_tree(image.smooth_reduced_image)
                    maxtree.compute_attributes(band_args, image)

                except Exception as e:

                    maxtree = MaxTree()
                    maxtree.construct_max_tree(image.smooth_image)
                    maxtree.compute_attributes(band_args, image)
                    dark_frame.estimate_morph_bg(image, maxtree)

                    if band_args['skip_reduction']:

                        dark_frame.bg_map = np.full_like(image.image, 0, dtype=np.float32)

                    else:

                        dark_frame.bg_map = np.full_like(image.image, dark_frame.bg_mean, dtype=np.float32)

                    dark_frame.create_reduced_image(image)

            else:

                maxtree = MaxTree()
                maxtree.construct_max_tree(image.smooth_image)
                maxtree.compute_attributes(band_args, image)
                dark_frame.estimate_morph_bg(image, maxtree)

                if band_args['skip_reduction']:

                    dark_frame.bg_map = np.full_like(image.image, 0, dtype=np.float32)

                dark_frame.create_reduced_image(image)

            tree_ids.append(_id)

            maxtree.detect_significant_objects(dark_frame, band_args)
            maxtree.move_up(dark_frame, band_args)

            extractor = Extractor()
            extractor.create_segmentation(maxtree, image, run, _id)
            extractor.extract_parameters(extractor, maxtree, run, image, _id)

            tree_of_segments = extractor.maxtree_of_segment
            unique_segment_ids = extractor.ids

            label_data = np.full(tree_of_segments.num_vertices(), -1, dtype=np.int32)
            segment_ids = np.arange(tree_of_segments.num_leaves(), tree_of_segments.num_vertices())
            label_data[segment_ids] = unique_segment_ids[segment_ids]

            flux = hg.accumulate_sequential(tree_of_segments, image.image, hg.Accumulators.sum)

            trees.append(extractor.maxtree_of_segment)
            latitudes.append(maxtree.x[extractor.segment_node_map])
            longitudes.append(maxtree.y[extractor.segment_node_map])
            fluxes.append(base_utils.normalize(flux))
            gammas.append(maxtree.gamma[extractor.segment_node_map])
            ids.append(label_data)
            areas.append(base_utils.normalize(maxtree.area[extractor.segment_node_map]))
            volumes.append(base_utils.normalize(maxtree.volume[extractor.segment_node_map]))

            run.status = "Completed"
            io_utils.save_run_metadata(run, band_args, _id)

        except KeyboardInterrupt:

            run.status = "Interrupted"
            io_utils.save_run_metadata(run, band_args, _id)

            print("\nMMTO run interrupted by user!")
            sys.exit(1)

        except Exception as e:

            run.status = "Terminated"
            io_utils.save_run_metadata(run, band_args, _id)

            print(f"\nMMTO run terminated with error: {e}")

            sys.exit(1)

    mmto.tree_map(
        trees,
        latitudes,
        longitudes,
        fluxes,
        gammas,
        areas,
        volumes,
        ids,
        tree_ids,
        run.time_stamp,
        run.arguments.co_sim,
        run.arguments.pix_dist,
    )


if __name__ == "__main__":

    mmto_run()

