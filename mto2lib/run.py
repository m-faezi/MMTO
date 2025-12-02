from datetime import datetime
import os
import yaml
from mto2lib.parser import make_parser


class Run:

    def __init__(self):

        self.arguments = None
        self.results_dir = None
        self.status = None
        self.time_stamp = None
        self.config = None
        self.trees = []
        self.latitudes = []
        self.longitudes = []
        self.fluxes = []
        self.volumes = []
        self.gammas = []
        self.areas = []
        self.ids = []
        self.tree_ids = []

    def setup_args(self):

        self.status = "Running"
        self.arguments = make_parser().parse_args()
        self.time_stamp = datetime.now().isoformat()
        self.results_dir = os.path.join("./results", self.time_stamp)

        os.makedirs(self.results_dir, exist_ok=True)

        with open('config.yaml', 'r') as f:

            self.config = yaml.safe_load(f)

        return self

