from mto2lib import io_utils
from mto2lib.parser import make_parser
from mto2lib.validators import validate_crop_coordinates
from datetime import datetime
import os


def setup_args():

    arguments = make_parser().parse_args()

    arguments.time_stamp = datetime.now().isoformat()

    print("Run ID: " + arguments.time_stamp)

    results_dir = os.path.join("./results", arguments.time_stamp)

    os.makedirs(results_dir, exist_ok=True)

    return arguments, results_dir


def get_image(arguments, results_dir):

    image, header = io_utils.read_image_data(arguments.file_path)

    return image, header

