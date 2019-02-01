import h5py
import argparse
import numpy as np
from niio import loaded
import os

from congrads import conmap

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--subject',
                    help='Input subject name.', required=True, type=str)
parser.add_argument('-f', '--features',
                    help='Feature data file.', required=True, type=str)
parser.add_argument('-c', '--cluster', help='Cluster map -- generate eta2 matrix for each label in a map.', 
                    required=True, type=str)
parser.add_argument('-d', '--dir', help='Output directory.',
                    required=True, type=str)
parser.add_argument('-bo', '--base_out',
                    help='Base output name, without extension.', required=True, type=str)

args = parser.parse_args()

if not os.path.exists(args.dir):
    print('Output directory does not exist -- creating now.')
    os.mkdir(args.dir)

# Load region of interest
try:
    assert os.path.isfile(args.cluster)
except:
    raise('Mask ROI does not exist.')
else:
    print('Loading ROI.')
    cluster = loaded.load(args.cluster)

# Load feature matrix
try:
    assert os.path.isfile(args.features)
except:
    raise('Feature file for {:} does not exist.'.format(args.subject))
else:
    print('Loading feature data.')
    F = loaded.load(args.features)
    n, p = F.shape

    if n < p:
        F = F.T

F = (F-F.mean(1)[:, None]) / (F.std(1)[:, None])

labels = list(set(np.unique(cluster)).difference({0, -1}))

fext = '{:}{:}.L.{:}'.format(args.dir, args.subject, args.base_out)

sims = h5py.File(name='{:}.ConnectopySims.h5'.format(fext),
                 mode='a')

for lab in labels:

    print('Processing label {:}'.format(lab))

    bool_inds = (cluster == lab)
    indices = np.where(cluster == lab)[0]

    sims.create_group(name=str(lab))
    sims[str(lab)].create_dataset(name='indices', data=indices)

    A = F[bool_inds, :]

    print('ROI shape: {:}'.format(A.shape))

    print('Transpose to generate time X samples matrix.')
    A = A.T

    B = F[~bool_inds, :]
    zeros = np.isnan(np.abs(B).sum(1))
    B = B[~zeros, :]

    print('Target shape: {:}'.format(B.shape))
    print('Transpose target to generate time X samples matrix.')
    B = B.T

    print('Computing voxel-wise connectivity fingerprints...')
    [evecs, Bhat, evals] = conmap.pca(B)
    R = conmap.corr(A, Bhat)

    print('Computing similarity matrix.')
    E2 = conmap.eta2(R)

    sims[str(lab)].create_dataset(name='eta2', data=E2)
    sims[str(lab)].create_dataset(name='r2', data=R)

sims.close()
