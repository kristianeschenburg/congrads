import argparse
import numpy as np
import nibabel as nb
import os
from niio import write, loaded
from scipy.linalg import eigh

from fragmenter import RegionExtractor as re
from congrads import conmap

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--subject', help='Subject to process.',
    required=True, type=str)
parser.add_argument('-l', '--label', help='Label file.',
    required=True, type=str)
parser.add_argument('-sr', '--sroi', help='Source rois.', required=True,
    type=str, nargs='+')
parser.add_argument('-sim', '--similarity', help='Similarity matrix.',
    required=True, type=str)
parser.add_argument('-d', '--dir', help='Output directory.',
    required=True, type=str)
parser.add_argument('-o', '--output', help='Output name, without extension.',
    required=True, type=str)
parser.add_argument('-e', '--evecs', help='Number of eigenvalues to compute.',
    required=False, default=5, type=int)
parser.add_argument('-n', '--normalize',
    help='Normalize eigenvectors to [0,1].',
    required=False, default=True, type=bool)
parser.add_argument('-hemi', '--hemisphere', help='Hemisphere to process.',
    required=False, default='L', choices=['L', 'R'], type=str)

args = parser.parse_args()

outEvecs = ''.join([args.dir, args.output])

hemi_map = {'L': 'CortexLeft',
            'R': 'CortexRight'}

try:
    assert os.path.isfile(args.label)
except:
    raise('Label file does not exist.')
else:
    R = re.Extractor(args.label)
    regions = R.map_regions()
    inds = R.indices(regions, args.sroi)

try:
    assert os.path.isfile(args.similarity)
except:
    raise('Similarity matrix does not exist.')
else:
    sim = loaded.load(args.similarity)

sim[np.isnan(sim)] = 0
sim[np.isinf(sim)] = 0

row_sums = np.where(np.abs(sim).sum(1) != 0)[0]
col_sums = np.where(np.abs(sim).sum(0) != 0)[0]

assert np.all(row_sums == col_sums)
sim = sim[:, row_sums][col_sums, :]
inds = inds[row_sums]

if not os.path.exists(outEvecs):

    print('Computing distance matrix.')
    distance = conmap.norm(sim)

    print('Computing adjacency matrix.')
    adj = conmap.adjacency(distance)

    print('Computing laplacian.')
    W = np.multiply(adj, sim)
    D = np.diag(np.sum(W, 0))
    L = np.subtract(D, W)

    print('Computing the dominant ' + str(args.evecs) + ' connectopic maps...')
    l, y = eigh(L, D, eigvals=(0, args.evecs))

    print('Normalizing eigenvectors and correcting sign.')
    corr_vec = np.arange(len(inds))

    for evec in range(1, y.shape[1]):		
        y[:, evec] = np.multiply(y[:, evec],
                                np.sign(np.corrcoef(y[:, evec], corr_vec)[0, 1]))


    print('Writing sign-flipped eigenvectors.')
    z = np.zeros((32492,y.shape[1]-1))
    for evec in range(0, y.shape[1]-1):
        z[inds, evec] = y[:, evec+1]
    write.save(z, '%s.Signed' % (outEvecs), hemi_map[args.hemisphere])

    if args.normalize:
        for evec in range(0, y.shape[1] - 1):
            # normalize to range 0-1
            tmp = y[:, evec] - min(y[:, evec])
            y[:, evec] = np.divide(tmp, (max(y[:, evec]) - min(y[:, evec])))

    print('Writing connectopic maps to file...')
    z = np.zeros((32492, y.shape[1]-1))
    for evec in range(0, y.shape[1]-1):
        z[inds, evec] = y[:, evec+1]

    write.save(z, outEvecs, hemi_map[args.hemisphere])