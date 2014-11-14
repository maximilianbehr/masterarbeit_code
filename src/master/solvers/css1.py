# -*- coding: utf-8 -*-
from css import Solver as CSS

class Solver(CSS):
    "First-order consistent splitting scheme by Guermond and Shen."

    def __init__(self, options):
        CSS.__init__(self, options, 1)

    def __str__(self):
        return "CSS1"
