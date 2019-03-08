import argparse
import nibabel as nb
import numpy as np
import scipy.io as sio
import os

from sklearn.cluster import SpectralClustering as spect

from fragmenter import RegionExtractor as re
import networkx as nx
from surface_utilities import adjacency as adj

from niio import loaded, write
from congrads import conmap

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--subject', help='Subject ID.', required=True, type=str)
parser.add_argument('-l', '--label', help='Label file.', required=True, type=str)
parser.add_argument('-surf', '--surface', help='Surface file.', required=True, type=str)
parser.add_argument('-e', '--eta', help='Eta2 file.', required=True, type=str)
parser.add_argument('-c', '--clusters', help='Cluster range.', required=True,
    type=int, nargs='+')
parser.add_argument('-r', '--rois', help='ROI names.', required=True, type=str,
    nargs='+')

parser.add_argument('-hops', '--hop_distance', help='Maximum hop distance.',
    required=True, type=int)
parser.add_argument('-dm', '--distmat', help='Precomputed distance matrix.',
    required=False, type=str, default=None)

parser.add_argument('-d', '--dir', help='Output directory.', required=True, type=str)
parser.add_argument('-bo', '--outbase', help='Output base name.', required=True, type=str)

print('Generating region map...')
args = parser.parse_args()
label = loaded.load(args.label)
R = re.Extractor(args.label)
indices = R.indices(R.map_regions(), args.rois)
sort_inds = np.argsort(indices)

print('Loading eta matrix...')
eta = loaded.load(args.eta)
eta[np.isinf(eta)] = 0
eta[np.isnan(eta)] = 0

cmin = args.clusters[0]
cmax = args.clusters[1]

print('Generating adjacency matrix...')
surf = nb.load(args.surface)
vertices = surf.darrays[0].data
faces = surf.darrays[1].data
S = adj.SurfaceAdjacency(vertices=vertices, faces=faces)
S.generate(indices=indices)

adjmat = np.zeros((len(indices), len(indices)))
index2coords = dict(zip(list(S.adj.keys()), np.arange(len(indices))))
coords2index = dict(zip(index2coords.values(), index2coords.keys()))

for k, v in S.adj.items():
    for n in v:
        adjmat[index2coords[k], index2coords[n]] = 1

if not args.distmat:
    print('Generating distance matrix...')
    G = nx.from_numpy_array(adjmat)
    apsp = nx.floyd_warshall_numpy(G)
else:
    print('Loading distance matrix...')
    apsp = loaded.load(args.distmat)

apsp = np.asarray(apsp)
A = {'apsp': apsp}

# save distance matrix if it doesn't exist yet
distfile = '{:}{:}.L.IPL.DistanceMatrix.mat'.format(args.dir, args.subject)
if not os.path.isfile(distfile):
    print('Saving distance matrix.')
    sio.savemat(file_name=distfile, mdict=A)

hops = list(np.arange(10, args.hop_distance, 10))
hops = hops + ['Full']

# loop over distances
for hp in hops:

    print('Clustering at distance: {:}'.format(hp))

    # sort eta matrix
    sorted_eta = eta[:, sort_inds]
    sorted_eta = sorted_eta[sort_inds, :]

    try:
        float(hp)
    except ValueError:
        print('Processing whole similarity matrix.')
    else:
        # threshold distance matrix
        print('Processing distance-thresholded similarity matrix.')
        distmat = np.asarray(apsp <= hp)
        sorted_eta = sorted_eta*distmat
    finally:
        pass
        
    print('{:} non-zero entries in eta matrix.'.format((sorted_eta != 0).sum()))

    # loop over cluster counts
    for clust_count in np.arange(cmin, cmax+1):

        print('Clusters: {:}'.format(clust_count))

        S = spect(n_clusters=clust_count, affinity='precomputed')
        S.fit(sorted_eta)

        labs = S.labels_
        labs[labs==0] = (labs.max() + 1)

        z = np.zeros((label.shape))
        z[indices[sort_inds]] = labs

        out_path = '{:}{:}.L.{:}.Cluster.{:}.Distance.{:}.func.gii'.format(
            args.dir, args.subject, args.outbase, clust_count, hp)
        write.save(z, out_path, 'CortexLeft')
