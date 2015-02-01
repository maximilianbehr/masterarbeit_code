import time
import os

def createdir(file):
    if not os.path.exists(os.path.dirname(file)):
        os.makedirs(os.path.dirname(file))


def gettime():
    return time.strftime("%d/%m/%Y %H:%M:%S")


