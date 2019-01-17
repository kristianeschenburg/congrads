import argparse
import numpy as np
from niio import loaded
import os

from congrads import conmap

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--subject', help='Input subject name.', required=True, type=str)
parser.add_argument('-f', '--features', help='Feature data file.', required=True, type=str)
parser.add_argument('-r', '--roi', help='ROI mask.', required=True, type=str)
parser.add_argument('-t', '--target', help='Target ROI mask.', required=False, type=str, default=None)
parser.add_argument('-d', '--dir', help='Output directory.', required=True, type=str)
parser.add_argument('-bo', '--base_out', help='Base output name, without extension.', required=True, type=str)

args = parser.parse_args()

if not os.path.exists(args.dir):
    print('Output directory does not exist -- creating now.')
    os.mkdir(args.dir)

# Load region of interest
try:
    assert os.path.isfile(args.roi)
except:
    raise('Mask ROI does not exist.')
else:
    print('Loading ROI.')
    roi = loaded.load(args.roi)
    roi = (roi>0)

# Load target region
if args.target:
    try:
        assert os.path.isfile(args.target)
    except:
        raise('Target mask does not exist.')
    else:
        print('Loading target.')
        target_exists = True
        target = loaded.load(args.target)
        target = (target>0)
else:
    target_exists = False

# Load feature matrix
try:
    assert os.path.isfile(args.features)
except:
    raise('Feature file for {:} does not exist.'.format(args.subject))
else:
    print ('Loading feature data.')
    F = loaded.load(args.features)
    n,p = F.shape

    if n < p:
        F = F.T

F = (F-F.mean(1)[:, None]) / (F.std(1)[:, None])

A = F[roi, :]
print('ROI shape: {:}'.format(A.shape))

print('Transpose to generate time X samples matrix.')
A = A.T

if target_exists:
    B = F[target, :]
else:
    B = F[~roi, :]

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

print('Saving correlation and eta^2 matrix.')
import scipy.io as sio
r = {'r2': R}
e = {'eta2': E2}

fext = '{:}{:}.L.{:}'.format(args.dir, args.subject, args.base_out)

sio.savemat(file_name='{:}.Corr'.format(fext), mdict=r)
sio.savemat(file_name='{:}.Eta2'.format(fext), mdict=e)
