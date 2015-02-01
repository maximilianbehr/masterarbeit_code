# -*- coding: utf-8 -*-
from src.bws.mesh.bws import MeshBuilder


if __name__ == "__main__":

        # generate a mesh builder
        meshbuilder = MeshBuilder()

        for ref in range(0, 5):
            # refine and save mesh
            meshbuilder.refine()
            meshbuilder.save()





