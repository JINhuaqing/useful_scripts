import numpy as np

def mat_power(mat, pv=-1):
    """mat is a symmetric positive-semi-definite matrix
    """
    eps = 1e-10
    
    # p-inv
    S, U = np.linalg.eig(mat)
    assert np.sum(S[np.abs(S)>eps] <0) == 0, "mat should be PSD mat"
    
    S = np.abs(S) # S >=0 for PSD
    S1 = S.copy()
    S1[S<=eps] = 0
    S1[S>eps] = (S[S>eps])**(pv)
    
    mat_r = U @ np.diag(S1) @ U.T
    return mat_r

def geodesic_dist(Q1, Q2):
    """Calculate the geodesic distance between two sys-PSD matrices. 
       Strictly, Q1 and Q2 should be invertible, but I let it can be PSD
       Follows https://github.com/makto-toruk/FC_geodesic/blob/master/utils/distance_FC/distance_FC.py
    """
    eps = 1e-10
    Q1_neg_half = mat_power(Q1, -1/2)
    Q = Q1_neg_half @ Q2 @ Q1_neg_half
    eigvs, _ = np.linalg.eig(Q)
    assert np.sum(eigvs[np.abs(eigvs)>eps] <0) == 0, "mat should be PSD mat"
    
    eigvs = np.abs(eigvs) # Q is PSD
    eigvs_part = eigvs[eigvs>eps]
    dist = np.sqrt(np.sum(np.log(eigvs_part)**2))
    return dist
