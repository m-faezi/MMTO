from datetime import datetime
import os
import yaml


class Run:

    def __init__(self):

        self.arguments = None
        self.results_dir = None
        self.status = None
        self.time_stamp = None
        self.config = None

    def setup_args(self):

        self.status = "Running"
        self.time_stamp = datetime.now().isoformat()
        self.results_dir = os.path.join("./results", self.time_stamp)

        os.makedirs(self.results_dir, exist_ok=True)

        with open('config.yaml', 'r') as f:

            self.config = yaml.safe_load(f)

        return self

