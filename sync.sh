#!/bin/sh

rsync -avz  --exclude=*.dmg \
            --exclude=*.pyc \
            --exclude=*.pvd \
            --exclude=*.gz  \
            --exclude=*.xml \
            --exclude=*.vtu \
            --exclude=*.h5  \
            --exclude=*.xdmf  . behr@ssh.mpi-magdeburg.mpg.de:/home/behr
