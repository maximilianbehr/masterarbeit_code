# -*- coding: utf-8 -*-

from src.solvers.css import Solver as CSS


class Solver(CSS):
    "Second-order consistent splitting scheme by Guermond and Shen."

    def __init__(self, options):
        CSS.__init__(self, options, 2)

    def __str__(self):
        return "CSS2"
