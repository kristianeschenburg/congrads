import argparse
import numpy as np

from sklearn.cluster import SpectralClustering as spect

from fragmenter import RegionExtractor as re
from niio import loaded, write

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--subject', help='Subject ID.', required=True, type=str)
parser.add_argument('-l', '--label', help='Label file.', required=True, type=str)
parser.add_argument('-e', '--eta', help='Eta2 file.', required=True, type=str)
parser.add_argument('-c', '--clusters', help='Cluster range.', required=True,
    type=int, nargs='+')
parser.add_argument('-r', '--rois', help='ROI names.', required=True, type=str,
    nargs='+')
parser.add_argument('-d', '--dir', help='Output directory.', required=True, type=str)
parser.add_argument('-bo', '--outbase', help='Output base name.', required=True, type=str)

args = parser.parse_args()

label = loaded.load(args.label)
R = re.Extractor(args.label)
region_map = R.map_regions()

indices = []
for r in args.rois:
    indices.append(region_map[r])
indices = np.concatenate(indices)

eta = loaded.load(args.eta)
eta[np.isinf(eta)] = 0
eta[np.isnan(eta)] = 0

cmin = args.clusters[0]
cmax = args.clusters[1]

for c in np.arange(cmin, cmax+1):

    print('Clusters: {:}'.format(c))

    S = spect(n_clusters=c, affinity='precomputed')
    S.fit(eta)

    labs = S.labels_
    labs[labs==0] = (labs.max() + 1)

    z = np.zeros((label.shape))
    z[indices] = labs

    out_path = '{:}{:}.L.{:}.func.gii'.format(args.dir, args.subject, args.outbase)
    write.save(z, out_path, 'CortextLeft')