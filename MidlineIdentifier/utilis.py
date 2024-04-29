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



def EuclideanDist(pts, pt):
    """
    Calculate the Euclidean distance between one point and other point(s).

    Parameters
    ----------
    pts : :class:`array_like`
    pt : :class:`array_like` of size one

    Returns
    -------
    dist : :class:`float` or :class:`numpy.ndarray`
        Euclidean distance
    """

    if not isinstance(pts, np.ndarray):
        pts = np.array(pts)

    if not isinstance(pt, np.ndarray):
        pt = np.array(pt)

    if pt.shape[0] != 1:
        raise ValueError("pt must be of size one.")

    if pts.shape[0] == 1:
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
