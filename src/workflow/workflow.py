from .methods import *
import os

class Workflow:
    """
    This class is meant to house the computational workflow involved with CAVaLRi
    """

    def __init__(self, cohort):
        self.cohort = cohort

    def run(self):
        print([case.case_id for case in self.cohort.cases])
