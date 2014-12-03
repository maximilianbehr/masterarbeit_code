import time
import sys
import os


class Tee(object):
    def __init__(self, *files):
        self.files = files

    def write(self, obj):
        for f in self.files:
            f.write(obj)

    def flush(self):
        for f in self.files:
            f.flush()


class TeeHandler():
    def __init__(self, file):
        dir = os.path.dirname(file)
        if not os.path.exists(dir):
            os.makedirs(dir)
        self.log = open(file, "w")

    def start(self):
        self.original = sys.stdout
        sys.stdout = Tee(sys.stdout, self.log)

    def stop(self):
        sys.stdout = self.original
        self.log.close()


def deletedir(dir):
    fileList = os.listdir(dir)
    if os.path.exists(dir):
        for fileName in fileList:
            os.remove(os.path.join(dir, fileName))
        os.rmdir(dir)


def gettime():
    return time.strftime("%d/%m/%Y %H:%M:%S")

