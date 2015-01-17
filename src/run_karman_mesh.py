# -*- coding: utf-8 -*-
from src.mesh.karman import MeshBuilder


if __name__ == "__main__":

        # generate a mesh builder
        meshbuilder = MeshBuilder()

        for ref in range(0, 6):
            # refine and save mesh
            meshbuilder.refine()
            meshbuilder.save()




