import argparse
import numpy as np
import nibabel as nb
import os
import h5py
from niio import write, loaded
from scipy.linalg import eigh

from congrads import conmap

parser = argparse.ArgumentParser()

parser.add_argument('-s', '--subject', help='Subject to process.', 
    required=True, type=str)
parser.add_argument('-sim', '--similarity', help='Similarity matrix file.', 
    required=True, type=str)
parser.add_argument('-c', '--cluster', help='Cluster map.',
    required=True, type=str)
parser.add_argument('-d', '--dir', help='Output directory.', 
    required=True, type=str)
parser.add_argument('-o', '--outbase', help='Output name, without extension.', 
    required=True, type=str)
parser.add_argument('-e', '--evecs', help='Number of eigenvalues to compute.', 
    required=False, default=5, type=int)
parser.add_argument('-n', '--normalize', help='Normalize eigenvectors to [0,1].', 
    required=False, default=True, type=bool)

args = parser.parse_args()

clusters = loaded.load(args.cluster)

try:
    assert os.path.isfile(args.similarity)
except:
    raise('Similarity file does not exist.')
else:
    print('Loading similarity file.')
    sims = h5py.File(name=args.similarity, mode='r')

labels = list(sims.keys())

outEvecs = ''.join([args.dir, args.outbase, '.Evecs.h5'])
evecs = h5py.File(name=outEvecs, mode='a')

z = np.zeros((clusters.shape[0], args.evecs-1))

for lab in labels:

    print('Processing label {:}'.format(lab))

    evecs.create_group(name=str(lab))
    indices = np.asarray(sims[lab]['indices'])
    evecs[lab].create_dataset(name='indices', data=indices)
    
    simmat = np.asarray(sims[lab]['eta2'])

    print('Computing distance matrix.')
    distance = conmap.norm(simmat)

    print('Computing adjacency matrix.')
    adj = conmap.adjacency(distance)

    print('Computing laplacian.')
    W = np.multiply(adj, simmat)
    D = np.diag(np.sum(W, 0))
    L = np.subtract(D, W)

    print('Computing the dominant ' + str(args.evecs) + ' connectopic maps...')
    l, y = eigh(L,D, eigvals=(0, args.evecs))

    print('Normalizing eigenvectors and correcting sign.')
    corr_vec = np.arange(len(indices))

    for evec in range(1, y.shape[1]):		
        y[:,evec] = np.multiply(y[:, evec], np.sign(np.corrcoef(y[:, evec], corr_vec)[0, 1]))
        if args.normalize: 
            tmp = y[:, evec] - min(y[:, evec])
            y[:,evec] = np.divide(tmp, max(tmp))

    print('Eigenvector shape: {:}'.format(y.shape))
    print('Indices shape: {:}'.format(indices.shape))

    print('Writing connectopic maps to file...')
    for evec in range(0, y.shape[1]-1):
        z[indices, evec] = y[:, evec+1]

    evecs[lab].create_dataset(name='evecs', data=y)

outEvecs = ''.join([args.dir, args.output, '.func.gii'])
write.save(z, outEvecs, 'CortexLeft')
