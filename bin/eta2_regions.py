import argparse
import numpy as np
import scipy.io as sio

from niio import loaded
import os
from fragmenter import RegionExtractor as re

from congrads import conmap

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--subject', 
                        help='Input subject name.',
                        required=True, type=str)
parser.add_argument('-f', '--features', 
                        help='Feature data file.',
                        required=True, type=str)
parser.add_argument('-l', '--label', 
                        help='Label file.', required=True,
                        type=str)
parser.add_argument('-sr', '--sroi', 
                        help='Source rois.', required=True,
                        type=str, nargs='+')
parser.add_argument('-tr', '--troi', 
                        help='Target rois.', required=False,
                        type=str, default=None, nargs='+')
parser.add_argument('-hemi', '--hemisphere', 
                        help='Hemisphere to process.',
                        required=False, default='L', 
                        choices=['L', 'R'], type=str)
parser.add_argument('-d', '--dir', 
                        help='Output directory.', 
                        required=True, type=str)
parser.add_argument('-bo', '--base_out', 
                        help='Base output name, without extension.',
                        required=True, type=str)

parser.add_argument('-pf', '--full', 
                        help='Process full.', default=True, 
                        required=False, type=bool, choices=[True, False])
parser.add_argument('-pi', '--iters', 
                        help='Process iterations.', default=True,
                        required=False, type=bool, choices=[True, False])

args = parser.parse_args()

def eta2(data, sinds, tinds):

    """
    Sub-method for generating eta2 and correlation matrices.

    Parameters:
    - - - - -
    data: float, array
        input data array
    sinds, tinds: list
        source and target indices
    """

    data = (data-data.mean(1)[:, None]) / (data.std(1)[:, None])

    A = data[sinds, :]

    A[np.isnan(A)] = 0
    A[np.isinf(A)] = 0
    A = A.T

    # get target region data matrix
    B = data[tinds, :]

    B[np.isnan(B)] = 0
    B[np.isinf(B)] = 0

    zeros = (np.abs(B).sum(1) == 0)
    B = B[~zeros, :]
    B = B.T

    print('Computing voxel-wise connectivity fingerprints...')
    [evecs, Bhat, evals] = conmap.pca(B)
    R = conmap.corr(A, Bhat)

    print('Computing similarity matrix.')
    E2 = conmap.eta2(R)

    return [E2, R]

if not os.path.exists(args.dir):
    print('Output directory does not exist -- creating now.')
    os.mkdir(args.dir)

# Load region of interest
try:
    assert os.path.isfile(args.label)
except:
    raise('Label file does not exist.')
else:
    label = loaded.load(args.label)
    R = re.Extractor(args.label)
    region_map = R.map_regions()

# get source and target indices
sinds = R.indices(region_map, args.sroi)
if args.troi:
    tinds = R.indices(region_map, args.troi)
else:
    tinds = list(set(np.arange(label.shape[0])).difference(set(sinds)))

# Load feature matrix
try:
    assert os.path.isfile(args.features)
except:
    raise('Feature file for {:} does not exist.'.format(args.subject))
else:
    print('Loading feature data.')
    F = loaded.load(args.features)
    n, p = F.shape

    # n_samples should be greater than n_features
    if n < p:
        F = F.T

F[np.isnan(F)] = 0
F[np.isinf(F)] = 0

if args.full:

    print('Processing full.')

    fext_eta = '%s%s.%s.Eta2.%s.Full.mat' % (
        args.dir, args.subject, args.hemisphere, args.base_out)

    fext_cor = '%s%s.%s.Corr.%s.Full.mat' % (
        args.dir, args.subject, args.hemisphere, args.base_out)

    if not os.path.exists(fext_eta) and not os.path.exists(fext_cor):

        [E, R] = eta2(F, sinds, tinds)
        r = {'r2': R}
        e = {'eta2': E}

        sio.savemat(file_name=fext_cor, mdict=r)
        sio.savemat(file_name=fext_eta, mdict=e)

if args.iters:

    print('Processing iterations.')
    ranges = [(0, 1200), (1200, 2400), (2400, 3600), (3600, 4800)]

    for itx, inds in enumerate(ranges):

        r_data = F[:, inds[0]:inds[1]]
        [E, R] = eta2(r_data, sinds, tinds)

        r = {'r2': R}
        e = {'eta2': E}

        fext_eta = '%s%s.%s.Eta2.%s.Iter.%i.mat' % (
            args.dir, args.subject, args.hemisphere, args.base_out, itx)

        fext_cor = '%s%s.%s.Corr.%s.Iter.%i.mat' % (
            args.dir, args.subject, args.hemisphere, args.base_out, itx)

        sio.savemat(file_name=fext_cor, mdict=r)
        sio.savemat(file_name=fext_eta, mdict=e)
