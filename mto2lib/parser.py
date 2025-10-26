import argparse
import mto2lib.validators as validators


def make_parser():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--pix_dist',
        type=validators.restricted_non_negative,
        default=3,
        help='cross-band central distance (default = 0)'
    )

    parser.add_argument(
        '--co_sim',
        type=validators.restricted_normal,
        default=0.90,
        help='cross_band cosine similarity (default = .90)'
    )

    return parser

