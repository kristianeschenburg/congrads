import argparse
import numpy as np
import nibabel as nb
import os
from niio import write, loaded
from scipy.linalg import eigh

from congrads import conmap

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--subject', help='Subject to process.', required=True, type=str)
parser.add_argument('-r', '--roi', help='ROI to process.', required=True, type=str)
parser.add_argument('-sim', '--similarity', help='Similarity matrix.', required=True, type=str)
parser.add_argument('-d', '--dir', help='Output directory.', required=True, type=str)
parser.add_argument('-o', '--output', help='Output name, without extension.', required=True, type=str)
parser.add_argument('-e', '--evecs', help='Number of eigenvalues to compute.', required=False, default=5, type=int)
parser.add_argument('-n', '--normalize', help='Normalize eigenvectors to [0,1].', required=False, default=True, type=bool)

args = parser.parse_args()

try:
    assert os.path.isfile(args.roi)
except:
    raise('Mask ROI does not exist.')
else:
    print('Loading ROI.')
    roi = loaded.load(args.roi)
    roi = (roi>0)

try:
    assert os.path.isfile(args.similarity)
except:
    raise('Similarity matrix does not exist.')
else:
    sim = loaded.load(args.similarity)

print('Computing distance matrix.')
distance = conmap.norm(sim)

print('Computing adjacency matrix.')
adj = conmap.adjacency(distance)

print('Computing laplacian.')
W = np.multiply(adj, sim)
D = np.diag(np.sum(W, 0))
L = np.subtract(D, W)

print('Computing the dominant ' + str(args.evecs) + ' connectopic maps...')
l,y = eigh(L,D, eigvals=(0, args.evecs))

print('Normalizing eigenvectors and correcting sign.')
corr_vec = np.arange(roi.sum())

for evec in range(1, y.shape[1]):		
    y[:,evec] = np.multiply(y[:, evec], np.sign(np.corrcoef(y[:, evec], corr_vec)[0, 1]))
    if args.normalize: 
        tmp = y[:, evec] - min(y[:, evec])
        y[:,evec] = np.divide(tmp, max(tmp))

print('Writing connectopic maps to file...')
z = np.zeros((roi.shape[0], y.shape[1]-1))
for evec in range(0, y.shape[1]-1):
    z[roi, evec] = y[:, evec+1]

outEvecs = ''.join([args.dir, args.output])
write.save(z, outEvecs, 'CortexLeft')