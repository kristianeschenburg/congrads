import argparse
import numpy as np
from niio import loaded
import os
from fragmenter import RegionExtractor as re

from congrads import conmap

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--subject', help='Input subject name.', required=True, type=str)
parser.add_argument('-f', '--features', help='Feature data file.', required=True, type=str)
parser.add_argument('-l', '--label', help='Label file.', required=True, type=str)
parser.add_argument('-sr', '--sroi', help='Source rois.', required=True, type=str, nargs='+')
parser.add_argument('-tr', '--troi', help='Target rois.', required=False, type=str, default=None, nargs='+')
parser.add_argument('-d', '--dir', help='Output directory.', required=True, type=str)
parser.add_argument('-bo', '--base_out', help='Base output name, without extension.', required=True, type=str)

args = parser.parse_args()

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

sindices=[]
tindices=[]

# get indices of source region
for sr in args.sroi:
    sindices.append(region_map[sr])
sindices = np.concatenate(sindices)

# if targets are supplied, get indices of targets
# otherwise set to False
if args.troi:
   target_exists=True
   for tr in args.troi:
      tindices.append(region_map[tr])
   tindices = np.concatenate(tindices)
else:
   target_exists=False

# Load feature matrix
try:
    assert os.path.isfile(args.features)
except:
    raise('Feature file for {:} does not exist.'.format(args.subject))
else:
    print ('Loading feature data.')
    F = loaded.load(args.features)
    n, p = F.shape

    # n_samples should be greater than n_features
    if n < p:
        F = F.T

# standardize
F = (F-F.mean(1)[:, None]) / (F.std(1)[:, None])

# get source region data matrix and transpose
A = F[sindices, :]
print('ROI shape: {:}'.format(A.shape))

print('Transpose to generate time X samples matrix.')
A = A.T

# if target regions supplied, get target region data matrix
if target_exists:
    B = F[tindices, :]
else:
    B = F[sindices, :]

zeros = np.isnan(np.abs(B).sum(1))
print('Zero index target region shape: {:}'.format(zeros.shape))
B = B[~zeros, :]

print('Target shape: {:}'.format(B.shape))
print('Transpose target to generate time X samples matrix.')
B = B.T

print('Computing voxel-wise connectivity fingerprints...')
[evecs, Bhat, evals] = conmap.pca(B)
R = conmap.corr(A, Bhat)

print('Computing similarity matrix.')
E2 = conmap.eta2(R)

print('Saving correlation and eta^2 matrix.')
import scipy.io as sio
r = {'r2': R}
e = {'eta2': E2}

fext_eta = '{:}{:}.L.Eta2.{:}.mat'.format(args.dir, args.subject, args.base_out)
fext_cor = '{:}{:}.L.Corr.{:}.mat'.format(args.dir, args.subject, args.base_out)

sio.savemat(file_name=fext_cor, mdict=r)
sio.savemat(file_name=fext_eta, mdict=e)
