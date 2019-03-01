import argparse
import nibabel as nb
import numpy as np

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

parser.add_argument('-d', '--dir', help='Output directory.', required=True, type=str)
parser.add_argument('-bo', '--outbase', help='Output base name.', required=True, type=str)

print('Generating region map...')
args = parser.parse_args()
label = loaded.load(args.label)
R = re.Extractor(args.label)
region_map = R.map_regions()

print('Loading eta matrix...')
eta = loaded.load(args.eta)
eta[np.isinf(eta)] = 0
eta[np.isnan(eta)] = 0

cmin = args.clusters[0]
cmax = args.clusters[1]

indices = []
for r in args.rois:
    indices.append(region_map[r])
indices = np.concatenate(indices)
sort_inds = np.argsort(indices)

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

print('Generating distance matrix...')
G = nx.from_numpy_array(adjmat)
apsp = nx.floyd_warshall_numpy(G)

for dist in np.arange(1, args.hop_distance+1):

    print('Clustering at distance: {:}'.format(dist))

    distmat = (apsp<=dist).astype(np.float32)
    print('{:} non-zero entries in distance matrix.'.format(distmat.sum()))

    sorted_eta = eta[:, sort_inds]
    sorted_eta = sorted_eta[sort_inds, :]
    sorted_eta = sorted_eta*distmat
    print('{:} non-zero entries in eta matrix.'.format((sorted_eta != 0).sum()))

    for c in np.arange(cmin, cmax+1):

        print('Clusters: {:}'.format(c))

        S = spect(n_clusters=c, affinity='precomputed')
        S.fit(sorted_eta)

        labs = S.labels_
        labs[labs==0] = (labs.max() + 1)

        z = np.zeros((label.shape))
        z[indices[sort_inds]] = labs

        out_path = '{:}{:}.L.{:}.Cluster.{:}.Distance.{:}.func.gii'.format(args.dir, args.subject, args.outbase, c, dist)
        write.save(z, out_path, 'CortexLeft')