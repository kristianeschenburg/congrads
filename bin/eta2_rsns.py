import argparse
import nibabel as nb
import numpy as np
from niio import loaded, write
import os

from congrads import conmap

parser = argparse.ArgumentParser()
parser.add_argument('-f', '--features', help='Features to compute similarity matrices from.',
    required=True, type=str)
parser.add_argument('-r', '--rsn', help='Resting-state network mask file.',
    required=True, type=str)
parser.add_argument('-d', '--dir', help='Output directory.',
    required=True, type=str)
parser.add_argument('-ob', '--out_base', help='Output base name.',
    required=True, type=str)

args = parser.parse_args()

nb_r = nb.load(args.rsn)

