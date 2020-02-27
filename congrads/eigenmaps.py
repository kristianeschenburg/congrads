import numpy as np
from scipy.linalg import eigh
from congrads import conmap

def eigenmaps(similarity, sinds, evecs):
    """
    Compute the eigenmaps of a given similarity matrix.

    Parameters:
    - - - - -
    similarity: float, array
        symmetric similarity matrix
    evecs: int
        number of eigenvectors to compute
    normalize: bool
        generate [0, 1] ranging eigenvectors

    Returns:
    - - - -
    vectors: float, array
        eigenvectors of similarity matrix
    """

    similarity[np.isnan(similarity)] = 0
    similarity[np.isinf(similarity)] = 0

    row_sums = (np.abs(similarity).sum(1) != 0)
    col_sums = (np.abs(similarity).sum(0) != 0)

    assert np.all(row_sums == col_sums)
    S = similarity[:, row_sums][col_sums, :]
    inds = sinds[row_sums]

    print('Computing distance matrix.')
    distance = conmap.norm(S)

    print('Computing adjacency matrix.')
    adj = conmap.adjacency(distance)

    print('Computing laplacian.')
    W = np.multiply(adj, S)
    D = np.diag(np.sum(W, 0))
    L = np.subtract(D, W)

    print('Computing the dominant ' + str(evecs) + ' connectopic maps...')
    l, y = eigh(L, D, eigvals=(0, evecs))

    corr_vec = np.arange(len(inds))

    # compute sign-flipped eigenvectors
    sign_flipped = np.zeros((y.shape))
    for evec in range(1, y.shape[1]):
        temp = np.multiply(y[:, evec], 
                        np.sign(np.corrcoef(y[:, evec], corr_vec)[0, 1]))

        sign_flipped[:, evec] = temp

    sign_flipped[:, 0] = y[:, 0]

    signed = np.zeros((32492, sign_flipped.shape[1]-1))
    for evec in range(0, sign_flipped.shape[1]-1):
        signed[inds, evec] = sign_flipped[:, evec+1]

    return signed
