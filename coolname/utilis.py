from typing import Literal, get_args, Union, Annotated

import numpy as np
import pandas as pd
import os
import logging

from scipy import optimize, integrate, interpolate
from scipy.stats import bootstrap


def ParseArgs(args):
    """
    Parse arguments from the commandline.

    Parameters
    ----------
    args : :class:`list`
        List of arguments

    Returns
    -------
    fad, fimg, args.diskClosing, args.diskOpening, outdir, sample

    """

    if args.sample is None:
        logging.info('sample name is None. Setting to default.')
        sample = 'BudoidObject'
    else:
        sample = str(args.sample)

    if (args.dir is not None) and os.path.isdir(args.dir):
        logging.info('Directory is provided. Looking for dapi_segmentation_mask.tif and spatial_anndata_spotipy.h5ad in it.')
        fimg = os.path.join(args.dir, 'dapi_segmentation_mask.tif')
        fad = os.path.join(args.dir, 'spatial_anndata_spotipy.h5ad')
        indir = args.dir

    elif args.dir is None:
        if os.path.isfile(args.fimg):
            fimg = args.fimg
        else:
            raise FileNotFoundError('Invalid image path.')

        if os.path.isfile(args.fad):
            fad = args.fad
        else:
            raise FileNotFoundError('Invalid anndata path.')
        indir = os.path.split(fimg)[0]


    outdir = indir if args.outdir is None else args.outdir
    os.makedirs(outdir, exist_ok = True)

    return fad, fimg, args.diskClosing, args.diskOpening, outdir, sample



def ScaleMinMax(x):
    """
    Scale the input vector into the range between zero and one.

    Parameters
    ----------
    x : :class:`array_like`

    Returns
    -------
    :class:`array_like`
        Scaled x

    """

    if not isinstance(x, np.ndarray):
        x = np.array(x)

    return (x - min(x)) / (max(x) - min(x))



def EuclideanDist(pt, pts):
    """
    Calculate the Euclidean distance between one point and other point(s).

    Parameters
    ----------
    pt : :class:`array_like` of size one
    pts : :class:`array_like`

    Returns
    -------
    dist : :class:`float` or :class:`numpy.ndarray`
        Euclidean distance
    """

    if not instances(pts, np.ndarray):
        pts = np.array(pts)

    if len(pt) != 1:
        raise ValueError("pt must be of size one.")

    if len(pts) == 1:
        dist = np.linalg.norm(pts - pt)
    else:
        dist = np.linalg.norm(pts - pt, axis=1)

    return dist



# get average exp by condition
def grouped_obs(adata, groupby, method, layer=None, gene_symbols=None):
    """
    Get average exp by condition.

    Parameters
    ----------
    adata : :class:`anndata.AnnData`
        Annotated data matrix.
    groupby : :class:`str`
        The key of the observations grouping to consider.
    method : :class:`str`
        Method used to aggregate the expression. Must be one of ``['sum','mean']``
    layer : :class:`str` (default: `None`)
        Key from `adata.layers` whose value will be used to. If None, adata.X will be used.
    gene_symbols : :class:`list` | :class:`None` (default: `None`)
        Genes to aggregate. If None, calculation will be done for all genes


    Returns
    -------
    :class:`pandas.DataFrame`
        A gene by group dataframe

    """

    if method not in ['sum', 'mean']:
        raise ValueError(f"Method must be one of {['sum', 'mean']}.")

    if layer is not None:
        try:
            getX = lambda x: x.layers[layer]
        except:
            raise ValueError(f"layer must be one of {adata.layers.keys()}")
    else:
        getX = lambda x: x.X


    if gene_symbols is not None:
        excl = set(gene_symbols) - set(adata.var_names)
        genes = set(gene_symbols) & set(adata.var_names)

        if len(excl) > 0:
            logging.warning('%s is not in the adata, continue without them', excl)
        adata = adata[:,list(genes)]

    grouped = adata.obs.groupby(groupby)
    out = pd.DataFrame(
        np.zeros((adata.n_vars, len(grouped)), dtype=np.float64),
        columns=list(grouped.groups.keys()),
        index=adata.var_names
    )

    if method == 'sum':
        for group, idx in grouped.indices.items():
            X = getX(adata[idx, :])
            out[group] = np.ravel(X.sum(axis=0, dtype=np.float64)).tolist()

    else: # mean
        for group, idx in grouped.indices.items():
            X = getX(adata[idx, :])
            out[group] = np.ravel(X.mean(axis=0, dtype=np.float64)).tolist()

    return out



## ================== plot on path ====================


def get_gene_expr(adata, gene, norm, coords):

    g_idx = np.where(adata.var.index == gene)[0][0]
    major_coor = adata.obs[coords]

    aggr = pd.DataFrame({'coords': major_coor, 'expr': adata.raw.X[:, g_idx]})

    if norm:
        aggr = aggr.groupby('coords', as_index=False).mean()
    return aggr



def bin_aggr(x, idx, method, sf = 10, q = 0.95):
    if method not in ['sum', 'mean', 'median', 'max', 'q', 'LogMean']:
        sys.exit('method must be one of [sum, mean, median, max, q, LogMean]')
    if method == 'LogMean' and isinstance(sf, int):
        return np.log(1 + sf * np.mean(np.take(x, idx, axis = 0)))
    elif method == 'sum':
        return np.sum(np.take(x, idx))
    elif method == 'mean':
        return np.mean(np.take(x, idx))
    elif method == 'median':
        return np.median(np.take(x, idx))
    elif method == 'max':
        return np.max(np.take(x, idx))
    elif method == 'q' and q <=1 :
        return np.quantile(np.take(x, idx), q)



def get_bins(x, nbin):
    cut_index, cut_bins = pd.cut(x, nbin, labels = False, duplicates='drop', retbins=True)
    mid_pt = [(a+b) / 2 for a, b in zip(cut_bins[1:], cut_bins[:-1])]

    return cut_index, cut_bins, mid_pt



def expr_binning(x, bins, bin_method = 'LogMean', bootstrapping = False):

    aggr_bin = []
    for n in np.sort(np.unique(bins)):
        idx = np.where(bins == n)[0]
        if bootstrapping:
            idx = np.random.choice(idx, size = len(idx)*5, replace = True)
        aggr_bin.append(bin_aggr(x, idx, method = bin_method))

    return aggr_bin



def minimize_linear(X, y):

    def linear_reg(params):
        return params[0] * X + params[1]

    def loss(params):
        return np.sum((linear_reg(params) - y)**2)

    def con(params):
        return sum(params)

    bnds = ((None, None), (0, None))
    cons = ({'type': 'ineq', 'fun': con})
    res = optimize.minimize(loss, [1, 1], method='SLSQP', bounds=bnds, constraints = cons)
    return np.poly1d(res.x)



def minimize_sigmoid(X, y):

    def sigmoid_reg(params):
        c1, c2, c3 = params
        return c3 / (1 + np.exp(-c1 * (X - c2)))

    def loss(params):
        return np.sum((sigmoid_reg(params) - y)**2)

    bnds = ((0, None), (0, 1), (0, None))
    res = optimize.minimize(loss, [1, 1, 1], method='SLSQP', bounds=bnds)

    c1, c2, c3 = res.x
    # print(c1, c2, c3)
    return lambda x: c3 / (1 + np.exp(-c1 * (x - c2)))


from scipy.special import gamma
def minimize_gamma(X, y):

    def gamma_reg(params):
        k, th = params
        return (X**(k-1) * np.exp(-X/th)) / (gamma(k) * th**k)

    def loss(params):
        return np.sum((gamma_reg(params) - y)**2)

    bnds = ((1, None), (1e-5, None))
    res = optimize.minimize(loss, [1, 0], method='SLSQP', bounds = bnds)

    a, b = res.x
    return lambda x: (x**(a-1) * np.exp(-b * x) * b**a)/gamma(a)



def minimize_splqua(X, y):

    def splrep_qua(s):
        return interpolate.splrep(X, y, k = 2, s = s)

    def loss(s):
        tck = splrep_qua(s)
        y_ = interpolate.splev(X, tck, der=0)
        der = interpolate.splev(X, tck, der=1)
        return 0.65 * np.sum((y_ - y)**2) + 0.35 * np.quantile(np.abs(der), 0.95)

    bnds = [(0.1, None)] # s being too small causes warnings
    res = optimize.minimize(loss, [1], method='SLSQP', bounds=bnds)

    return interpolate.splrep(X, y, k = 2, s = res.x), res.x



from scipy.signal import peak_prominences, argrelextrema

def find_peaks_scipy(y):
    y2 = np.array([0] + y.tolist() + [0])
    peaks = argrelextrema(y2, np.greater)[0]
    vally = argrelextrema(y2, np.less)[0]
    breaks = np.sort(np.concatenate([peaks, vally]))

    proms = peak_prominences(y2, peaks = peaks)[0]
    return len(peaks), np.std(proms)
