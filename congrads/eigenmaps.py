import numpy as np
from scipy.linalg import eigh
from congrads import conmap

def eigenmap(sim, inds, normalize=True, evecs=6):
    
    """
    Compute eigenmaps of eta2 matrix.
    
    Parameters:
    - - - - -
    sim: float, array
        eta2 matrix
    inds: list
        list of indices in cortical map to which eigenmap projects onto
    """
    
    sim[np.isnan(sim)] = 0
    sim[np.isinf(sim)] = 0
    
    row_sums = np.where(np.abs(sim).sum(1) != 0)[0]
    col_sums = np.where(np.abs(sim).sum(0) != 0)[0]

    assert np.all(row_sums == col_sums)
    sim = sim[:, row_sums][col_sums, :]
    inds = inds[row_sums]
    
    # compute distance and adjacency matrix
    distance = conmap.norm(sim)
    adj = conmap.adjacency(distance)
    
    # compute graph laplacian
    W = np.multiply(adj, sim)
    D = np.diag(np.sum(W, 0))
    L = np.subtract(D, W)
    
    # compute eigendecomposition of laplacian
    l, y = eigh(L, D, eigvals=(0, evecs))
    
    # correct vectors for sign
    corr_vec = np.arange(len(inds))
    for evec in range(1, y.shape[1]):
        y[:, evec] = np.multiply(y[:, evec],
                                np.sign(np.corrcoef(y[:, evec], corr_vec)[0, 1]))
    
    sign_flipped = np.zeros((32492, y.shape[1]-1))
    for evec in range(0, y.shape[1]-1):
        sign_flipped[inds, evec] = y[:, evec+1]
        
    # normalized to range [0, 1]
    if normalize:
        for evec in range(0, y.shape[1] - 1):
            # normalize to range 0-1
            tmp = y[:, evec] - min(y[:, evec])
            y[:, evec] = np.divide(tmp, (max(y[:, evec]) - min(y[:, evec])))
            
    normalized = np.zeros((32492, y.shape[1]-1))
    for evec in range(0, y.shape[1]-1):
        normalized[inds, evec] = y[:, evec+1]
        
    return [sign_flipped, normalized]