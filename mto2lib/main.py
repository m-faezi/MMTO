from mto2lib import io_utils


def get_image(arguments, results_dir):

    image, header = io_utils.read_image_data(arguments.file_path)

    return image, header

