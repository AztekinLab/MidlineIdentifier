import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import scanpy as sc
import seaborn as sns



# plot along coords
def trend_plot(budoid, features, groupby, coords = 'major_coor_used', save = False, **kwargs):
    """
    Makes a *trend plot* of the expression values of `var_names` as a function of `coords`

    For each var_name and each `groupby` category a dot is plotted.
    Each dot represents two values: mean expression within each category
    (visualized by color) and fraction of cells expressing the `var_name` in the
    category (visualized by the size of the dot). If `groupby` is not given,
    the dotplot assumes that all data belongs to a single category.

    This function use :func:`seaborn.lmplot`. If you need more flexibility, you should use :func:`seaborn.lmplot` directly.


    Parameters
    ----------
    feature : :class:`str` | :class:`list`
        Column name in `.var` DataFrame that stores gene symbols. By default `var_names` refer to the index column of the `.var` DataFrame.
    groupby : :class:`str`
        The key of the observation grouping to consider. Must be one of `obs.columns`
    coords : :class:`str` (default: `'major_coor_used'`)
        To which the gene expression should be consider to. Must be one of `obs.columns`.
    save : :class:`bool` (default: `False`)
        If `True` or a `str`, save the figure. A string is appended to the default filename. Infer the filetype if ending on {`'.pdf'`, `'.png'`, `'.svg'`}.
    kwargs
        Additonal arguments to pass to :func:`seaborn.lmplot`


    Returns
    -------
    :meth:`seaborn.lmplot` object.


    Examples
    --------

    Create a trend plot using the given markers using an example dataset grouped by the category 'batch'.

    .. highlight:: python
    .. code-block:: python

        import PSUils as ps

        budoid1 = ps.io.ReadObj('testdata/Budoid_1A/Budoids.pkl')
        budoid2 = ps.io.ReadObj('testdata/Budoid_3H/Budoids.pkl')
        budoid1.Concat(budoid2)

        markers = ['Col9a2','Col3a1']
        sc.pl.dotplot(budoid1, markers, groupby='batch')

    """

    adata = budoid.data.adata

    if groupby is None or groupby not in adata.obs.columns:
        raise ValueError(f"groupby = {groupby} must be one of {adata.obs.columns}")

    if isinstance(features, str):
        features = [features]

    kwargs_final = {
        'order': 2,
        'line_kws' :{'lw':5},
        'x_bins' : np.linspace(0, 1, 8)[1:-1],
        'x_estimator': lambda x: np.log(np.mean(x) + 1),
        'truncate' : True
        }

    kwargs_final.update(kwargs)

    genes = list(set(features) & set(adata.var_names))
    if len(genes) == 0:
        print('None of the requested genes is in the data.')
    if len(genes) < len(features):
        print("%s is not in the data, continue without them." % ', '.join(set(features) - set(genes)))


    df = sc.get.obs_df(adata, keys = genes + [coords, groupby])
    df = pd.melt(df,id_vars = [coords, groupby], var_name='genes', value_name='exp')

    lm = sns.lmplot(df, x = coords, y = 'exp', hue = groupby, col = 'genes', **kwargs_final)

    if isinstance(save, bool):
        fn = 'trendplot.pdf'
        plt.savefig(fn)
    elif isinstance(save , str):
        fn = save if save.endswith((".svg", ".pdf", ".png")) else save + '.pdf'
        plt.savefig(fn)

    return lm
