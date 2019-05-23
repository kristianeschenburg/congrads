import argparse
import numpy as np
import scipy.io as sio

from niio import loaded
import os
from fragmenter import RegionExtractor as re

from congrads import conmap

parser = argparse.ArgumentParser()
parser.add_argument('-s', '--subject', help='Input subject name.',
    required=True, type=str)
parser.add_argument('-f', '--features', help='Feature data file.',
    required=True, type=str)
parser.add_argument('-rsns', '--networks', help='Resting state network map.',
    required=True, type=str)
parser.add_argument('-lt', '--label_table', help='RSN label table file.',
    required=True, type=str)
parser.add_argument('-n2r', '--network2rsn', help='Mapping of label name to RSN name.',
    required=True, type=str)
parser.add_argument('-d', '--dir', help='Output directory.', required=True,
    type=str)
parser.add_argument('-hemi', '--hemisphere', help='Hemisphere to process.',
    required=False, type=str, choices=['L', 'R'], default='L')

args = parser.parse_args()

with open(args.label_table, 'r') as inlabtab:
    lab_tab = inlabtab.readlines()

lab_dict = {}
for i in np.arange(0, len(lab_tab), 2):
    parts = lab_tab[i+1].strip()
    parts = parts.split(' ')
    
    temp = {'alpha': np.int(parts[0]),
            'rgb': parts[1:4]}
            
    lab_dict[lab_tab[i].strip()] = temp

networks = {k: lab_dict[k] for k in lab_dict.keys() if k[0] == '7'}
labels = loaded.load(args.networks)[:, 0]

with open(args.network2rsn, 'r') as innet:
    net2rsn = innet.readlines()
net2rsn = [k.strip().split(' ') for k in net2rsn]
net2rsn = {k[0]: k[1] for k in net2rsn}

if not os.path.exists(args.dir):
    print('Output directory does not exist -- creating now.')
    os.mkdir(args.dir)

features = loaded.load(args.features)
if features.shape[0] < features.shape[1]:
    features = features.T

data_inds = (np.abs(features).sum(1) != 0)

for net in networks.keys():

    print('Processing: {:}'.format(net2rsn[net]))
    alpha = networks[net]['alpha']
    temp_idx = (labels == alpha)

    gInds = (temp_idx * data_inds)
    eInds = (~temp_idx * data_inds)

    region_data = features[gInds, :]
    extern_data = features[eInds, :]

    print('Network data shape: {:}'.format(region_data.shape))
    print('External data shape: {:}'.format(extern_data.shape))

    region_data = region_data.T
    extern_data = extern_data.T

    print('Computing voxel-wise connectivity fingerprints...')
    [evecs, ehat, evals] = conmap.pca(extern_data)
    R = conmap.corr(region_data, ehat)

    print('Computing similarity matrix.')
    E2 = conmap.eta2(R)

    print('Saving correlation and eta^2 matrix.')
    r = {'r2': R}
    e = {'eta2': E2}

    fext_eta = '{:}{:}.{:}.Eta2.{:}.{:}.mat'.format(args.dir, args.subject,
        args.hemisphere, net2rsn[net])
    fext_cor = '{:}{:}.{:}.Corr.{:}.{:}.mat'.format(args.dir, args.subject,
        args.hemisphere, net2rsn[net])

    sio.savemat(file_name=fext_cor, mdict=r)
    sio.savemat(file_name=fext_eta, mdict=e)
