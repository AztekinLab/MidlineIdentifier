import scanpy as sc
from anndata import AnnData

import numpy as np
import pandas as pd


import matplotlib.pyplot as plt
from adjustText import adjust_text

import os
import logging

from . import utilis

class Adata:
    """
    Adata object to wrap anndata

    Attributes
    ----------
    adata : :class:`anndata.AnnData`
        anndata object to store the single cell data. Compatible with all scanpy functions.
    outdir : :class:`str`
        output directory to save files


    Methods
    -------
    Preprocessing
    EnrichBins
    FindOrientation
    FindDEG
    FindSVG

    """


    def __init__(self, fad, sample, outdir):

        """
        Constructs all the necessary attributes for the :class:`~coolname.adata_class.Adata`.

        Parameters
        ----------

        fad : :class:`str`
            Path to single cell object. Will be load into an anndata object.

        sample : :class:`str`
            Sample identifier. Will be store in ``adata.obs['sample']``. Useful when concanating multiple :class:`~coolname.adata_class.Adata` instances.
        outdir : :class:`str`
            Output directory to save files

        """

        self.adata =  sc.read_h5ad(fad)
        self.adata.obs['sample'] = sample
        self.outdir = outdir

    def Preprocessing(self):
        """
        Preprocessing of single cell dataset. This function calls :func:`scanpy.pp.normalize_total` and :func:`scanpy.pp.log1p`.

        Raw data, normalized data and log data will be stored into ``.layers['counts']``, ``.layers["norm_counts"]`` and ``.layers["lognorm_counts"]`` respectively.
        """

        adata = self.adata

        if adata.raw is None and isinstance(adata.X.max(), (int, np.integer)) :
            adata = adata[:, adata.var_names != 'nan']

        if adata.n_vars == 0:
            logging.warning('No detected genes in the adata.')

        adata.raw = adata
        adata.layers["counts"] = adata.X.copy() # storing raw counts in layers

        sc.pp.normalize_total(adata)
        adata.layers["norm_counts"] = adata.X.copy()

        sc.pp.log1p(adata)
        adata.layers["lognorm_counts"] = adata.X.copy()

        self.adata = adata


    def EnrichBins(self, genes, coords, nbin = 4, score_name = 'score', **kwargs):

        """
        Perform gene set enrichment in bins. This function calls :func:`scanpy.tl.score_genes`. Result will be stored into ``.uns[score_name]``.

        Parameters
        ----------
        genes : :class:`list` | :class:`str`
            The list of gene names used for score calculation
        coords : :class:`str`
            The key of the observations to consider. Must be one of the ``.obs.columns``
        nbin : :class:`int`
            The number of bins
        score_name : :class:`str`
            Name of the field to be added in ``.uns``
        kwargs
            Additonal arguments to pass to :func:`scanpy.tl.score_genes`
        """

        adata = self.adata

        genes = set(genes) & set(adata.var_names)
        logging.info(', '.join(genes) + 'are used to calculate enrichment.')

        cut_index = pd.cut(adata.obs[coords], nbin, labels = False, duplicates='drop', retbins=False)
        adata.obs['bins'] = cut_index

        out = utilis.grouped_obs(adata, 'bins', method = 'mean')
        adata_bins = sc.AnnData(out.transpose(), out.columns.to_frame(), out.index.to_frame())

        try:
            sc.tl.score_genes(adata_bins, genes, n_bins = 25, **kwargs)
            score = adata_bins.obs['score'].tolist()
        except IndexError:
            sc.tl.score_genes(adata_bins, genes, n_bins = 10, **kwargs)
            score = adata_bins.obs['score'].tolist()
        except ValueError:
            score = None

        adata.uns[score_name] = score
        self.adata = adata


    def FindOrientation(self, coords = 'major_coor_scaled', start_genes = ['Sox9','Acan','Col2a1','Col9a1','Col9a2','Col11a1'], end_genes = ['Col1a1', 'Col3a1'], plot = True, **kwargs):
        """
        Orient the coords based on the provided genelists. This allows cross-dataset/structure comparisons. Result will be stored into ``.uns[start_score]`` and ``.uns[end_score]``.

        Parameters
        ----------
        coords : :class:`str`
            The key of the observations to consider. Must be one of the ``.obs.columns``
        start_genes : :class:`list`
            The list of gene names used to calculate the start
        end_genes : :class:`list`
            The list of gene names used to calculate the end
        plot : :class:`bool` (default: `True`)
            Set to ``True`` by default. If True, save teh plot into ``'Orientation_score.pdf'`` in the output directory (`.outdir`)
        kwargs
            Additonal arguments to pass to :func:`~coolname.adata_class.Adata.EnrichBins`


        """

        self.EnrichBins(start_genes, coords, score_name = 'start_score', **kwargs)
        self.EnrichBins(end_genes, coords, score_name = 'end_score', **kwargs)

        adata = self.adata
        ss, es = adata.uns['start_score'], adata.uns['end_score']
        if ss[0] < ss[-1]:
            # reverse the path
            # self.start, self.end = self.end, self.start
            adata.obs['major_coor_used'] = 1 - adata.obs[coords]
        else:
            adata.obs['major_coor_used'] = adata.obs[coords]


        # max_s, max_e = max(ss), max(es)
        # if max_s > Thre and max_e > Thre:
        #     adata.obs['loc'] = 'P'
        #     adata.obs.loc[adata.obs['major_coor_used'] > 0.5, 'loc'] = 'D'
        # else:
        #     adata.obs['loc'] = 'A'

        self.adata = adata


        if plot:
            fig, ax = plt.subplots(1, 2, figsize = (8, 4))
            if ss is not None:
                ax[0].scatter(range(len(ss)), ss)
                ax[0].axhline(0, color = 'r', linestyle = 'dashed')
                ax[0].set_title('Start score')
            if es is not None:
                ax[1].scatter(range(len(es)), es)
                ax[1].axhline(0, color = 'r', linestyle = 'dashed')
                ax[1].set_title('End score')

            plt.savefig(os.path.join(self.outdir, 'Orientation_score.pdf'))




    def FindDEG(self, groupby, condition, method, **kwargs):
        """
        Finds differentially expressed genes (DEGs) for each of the identity classes in a dataset.

        Parameters
        ----------
        adata : :class:`anndata.AnnData`
            Annotated data matrix.
        groupby : :class:`str`
            The key of the observations grouping to consider.
        method : :class:`str` (`default: 'DESeq2'`)
            Method used to calcualte DEGs.
            ``'DESeq2'`` and ``'DESeq2_pb'`` use `pydeseq2 <https://github.com/owkin/PyDESeq2?tab=readme-ov-file>`_, the python implementation of the `DESeq2 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8>`_ method. ``DESeq2`` calculates DEGs on single cell level while ``DESeq2_pb`` generate pseudobulk expression based on ``condition``.
            ``'t-test'``, ``'t-test_overestim_var'`` overestimates variance of each group, ``'wilcoxon'`` uses Wilcoxon rank-sum, ``'logreg'`` uses logistic regression.

            If method is one of ``['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']``, This function directly calls :func:`scanpy.tl.rank_genes_groups`.
        kwargs
            Additonal arguments to pass to :func:`scanpy.tl.rank_genes_groups`


        Returns
        -------
        :class:`pandas.DataFrame`

        """

        # groups = [ref], reference = g,

        if method.startswith("DESeq2"):
            if method == "DESeq2_pb":
                groupby_ = condition
            else:
                groupby_ = groupby
            all_ = self.adata.obs[groupby_].unique().tolist()

            # check groups
            if 'groups' not in kwargs.keys():
                groups = all_
            else:
                groups = kwargs['groups']
            if isinstance(groups, (str, int)):
                groups = [groups]


            # check reference
            if 'reference' not in kwargs.keys():
                reference = all_
            else:
                reference = kwargs['reference']
            if isinstance(reference, str):
                reference = [reference]

                diff = set(reference) - set(all_)
                if len(diff) > 0:
                    raise ValueError(
                        f"reference = {diff} is not one of groupby = {all_}.")

            logging.info("Accessing adata.layers['counts'] to get the raw counts")
            adata_comp = self.adata[self.adata.obs[groupby_].isin(set(groups + reference))].copy()

            if method == "DESeq2_pb":
                df = utilis.grouped_obs(adata_comp, groupby = groupby, method = 'sum', layer = 'counts')

                g2c = dict(zip(adata_comp.obs[groupby],adata_comp.obs[condition]))
                meta = np.array([g2c.get(g) for g in df.columns])

                res_list = []
                for group_i in groups:
                    reference_i = list(set(reference)-set(group_i))

                    if len(reference_i)  == 0:
                        logging.warning('skipping comparison')
                        continue

                    col1_idx = np.where(meta == group_i)[0]
                    col2_idx = np.where(meta == reference_i)[0]

                    res_tmp = _deseq2(df, col1_idx, col2_idx)
                    res_tmp.insert(0, 'group', [group_i] * len(res_tmp))
                    res_list.append(res_tmp)

                res = pd.concat(res_list, axis = 0)

            else: # TODO: DESeq2 on single cell, very slow
                res_list = []
                for group_i in groups:
                    reference_i = list(set(reference)-set(group_i))

                    if len(reference_i) == 0:
                        logging.warning('skipping comparison')
                        continue

                    cells1_idx = np.where(adata_comp.obs[groupby] == group_i)[0]
                    cells2_idx = np.where(adata_comp.obs[groupby].isin(reference_i))[0]

                    res_list.append(_deseq2(adata_comp, cells1_idx, cells2_idx))
                res = pd.concat(res_list, axis = 1)

        else:
            sc.tl.rank_genes_groups(self.adata, groupby = groupby, method=method, **kwargs)
            # All groups are returned if groups is None
            res = sc.get.rank_genes_groups_df(self.adata, group = None)


            # add group column to res if len(groups) == 1
            if 'group' not in res.columns:
                g = list(self.adata.uns['rank_genes_groups']["names"].dtype.names)
                res.insert(0, 'group', g * len(res))

        return res


    def FindSVG(self, coords, sample, layer = 'counts', min_exp_gene = 0, min_exp_cell = 0):
        """
        Finds spatially variable genes (SVGs) for each of the identity classes in a dataset. This function incoporate :func:`SpatialDE`. Raw counts should be used.

        Parameters
        ----------
        sample : :class:`str` (default: `'sample'`)
            Sample identifier. Must be one of the `.obs.columns`
        coords : :class:`str`
            Spatial coordinates for each cell. Can be one of the ``.obs.columns`` or a :class:`pandas.DataFrame` with rows as cells and columns as spatial dimensions.
        layer : :class:`str` (default: `'counts'`)
            Key from `adata.layers` whose value will be used to. If None, .`adata.layers['counts']` will be used.
        min_exp_gene : :class:`int` (default: `'0'`)
            Filter genes whose expression lower than this
        min_exp_cell : :class:`int` (default: `'0'`)
            Filter cells whose total expression lower than this
        """

        if isinstance(coords, str) and (coords in self.adata.obs.columns):
            coords = self.adata.obs[coords]
        elif isinstance(coords, pd.DataFrame):
            if len(coords) != self.adata.n_obs:
                raise ValueError(f"The provided coords must have the same size as the number of testing cells {self.adata.n_obs}.")
        else:
            raise ValueError(f"The provided coords should be either one columns in adata.obs or a pandas.DataFrame.")

        # print(coords[counts.index])
        sample_unqiue = self.adata.obs[sample].unique()

        for s in sample_unqiue:
            adata = self.adata[self.adata.obs[sample] == s]
            counts = pd.DataFrame(adata.layers[layer], columns=adata.var_names, index=adata.obs_names)

            # Filter 'nan' genes
            counts = counts.loc[:,counts.columns != 'nan']
            # Filter practically unobserved genes
            counts = counts.loc[:,counts.sum(0) >= min_exp_gene]
            counts = counts.loc[counts.sum(1) >= min_exp_cell,: ]

            sample_info = pd.DataFrame({
                'coords': coords.loc[counts.index],
                'total_counts':counts.sum(1)},
                index = counts.index)

            _spatialde(counts, sample_info, self.outdir, s)





    def PolarizationScoring(self, genes, norm = False, coords = 'major_coor_used', bootstrapping = False, n_bs = 1000, random_state = 1234):

        """
        Calculate the polarization score for genes.

        Parameters
        ----------
        genes : :class:`list` | :class:`str`
            genes to calculate
        bootstrapping : :class:`bool` (default : `False`)
            Whether to bootstrapping or not
        n_bs : :class:`int` (default : `100`)
            The number of bootstrapping to perform
        random_state : :class:`int` (default : `1234`)
            Set random seed for reproducibility

        """

        adata = self.adata.adata

        if g not in adata.var_names:
            logging.warning("The requested gene %s is not in adata! Returning None", g)
            return None

        aggr_ = utilis.get_gene_expr(adata, g, norm = norm, coords = coords)

        if bootstrapping:
            np.random.seed(random_state)
            score = []
            for _ in range(n_bs):
                idx = np.random.choice(aggr_.index, size = len(aggr_), replace = True)
                df = aggr_.iloc[idx,:]
                score.append(_PolarizationScoring(df))

            # score = np.array(score) / (10 * np.std(score))
            bs = bootstrap((score, ), np.median, confidence_level = 0.95,random_state = random_state, method = 'percentile')
            CI = bs.confidence_interval
            return (score, CI)

        else:
            return _PolarizationScoring(aggr_)



from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats


# gene_mat = dr_utilis.grouped_obs_single(adata, groupby = 'sample', use_raw = True, method = 'sum') # raw counts
# gene_mat_PD = dr_utilis.grouped_obs_single(adata_PD, groupby = 'comp', use_raw = True, method = 'sum') # raw counts
# gene_mat_P = gene_mat_PD.iloc[:,gene_mat_PD.columns.str.endswith('P')]
# gene_mat_D = gene_mat_PD.iloc[:,gene_mat_PD.columns.str.endswith('D')]




def _deseq2(data, cells1, cells2, **kwagrs):
    """
    Finds differentially expressed genes (DEGs) using DESeq2 for one comparison.

    Parameters
    ----------
    data : :class:`anndata.AnnData` | :class:`pandas.DataFrame`
        Raw expression counts. If :class:`anndata.AnnData`, `.layers["counts"]` will be used. If :class:`pandas.DataFrame`, genes as rows and cells/samples as columns
    cells1 : :class:`list`
        Cells index for condition 1
    cells2 : :class:`list`
        Cells index for condition 2 (reference)
    **kwagrs
        Additonal arguments to pass to :func:`pydeseq2.dds.DeseqDataSet`

    Returns
    -------
    :class:`pandas.DataFrame`

    """

    if isinstance(data, AnnData):
        df = np.concatenate([data.layers["counts"][cells1], data.layers["counts"][cells2]], axis = 0)
        meta = pd.DataFrame({'comp': ['G1']* len(cells1) +  ['G2']* len(cells2)}, index = df.index)
    elif isinstance(data, pd.DataFrame):
        df = pd.concat([data.iloc[:,cells1], data.iloc[:,cells2]], axis = 1).transpose()
        meta = pd.DataFrame({'comp': ['G1']* len(cells1) +  ['G2']* len(cells2)}, index = df.index)
    else:
        raise ValueError('Unrecognized input format.')

    # print(df)
    # print(meta)
    dds = DeseqDataSet(
        counts=df,
        metadata=meta,
        design_factors="comp",
        ref_level = ['comp', 'G2'],
        refit_cooks=False,
        inference=DefaultInference(n_cpus=1),
        **kwagrs
    )

    dds.deseq2()
    stat_res = DeseqStats(dds)
    stat_res.summary()
    res = stat_res.results_df.sort_values('padj')

    return res



import NaiveDE
import SpatialDE
import SpatialDE.plot

def _spatialde(counts, sample_info, outdir, sample = 'sample', cal_pattern = True, random_state = 1234):
    """
    Finds spatially variable genes (SVGs) for one sample.

    THis function will save two / four files in the specified folder
    - '1.FSV_sig_' + s + '.pdf': Plot of Fraction Spatial Variance vs Q-value
    - '2.SVG_filter_' + s + '.csv': Statistics of significant SVG. See `here <https://github.com/Teichlab/SpatialDE>`_ for more detailed documentatiton.

    If `cal_pattern = True`:
    - '3.Pattern_' + s + '.pdf': Spatial patterns of gene expression
    - '4.Pattern_members_' + s + '.csv': Gene membership of spatial patterns


    Parameters
    ----------
    counts : :class:`pandas.DataFrame`
        Raw expression counts. Genes as columns and cells/samples as rows.
    outdir : :class:`str`
        Output directory
    sample: :class:`str` (default: `'sample'`)
        Sample identifier. Will be used as suffix of the output file.
    cal_pattern: :class:`bool` (default: `'True'`)
        If True, to detect the spatial patterns of gene expression and assign genes to patterns. This callls the :func:`SpatialDE.spatial_patterns`
    random_state : :class:`int` (default : `'1234'`)
        Set random seed for reproducibility

    """

    np.random.seed(random_state)

    counts_n = NaiveDE.stabilize(counts.T).T
    counts_r = NaiveDE.regress_out(sample_info, counts_n.T, 'np.log(total_counts)').T

    X = sample_info['coords'].values[:, None] # do not change. it's here for a reason
    results = SpatialDE.run(X, counts_r)
    results['pval'] = results['pval'].clip(lower = results.query('pval > 0')['pval'].min() / 2)
    results['qval'] = results['qval'].clip(lower = results.query('qval > 0')['qval'].min() / 2)


    plt.clf()
    plt.figure(figsize=(10,6))
    SpatialDE.plot.FSV_sig(results)

    text_list = []
    for i in range(results.shape[0]):
        if results.qval[i] < 0.05:
            text_list.append(plt.text(results.FSV[i], results.pval[i], results.g[i], fontsize = 12))

    if len(text_list) >0:
        adjust_text(text_list, arrowprops=dict(arrowstyle="-", color='k', lw=0.5))
    plt.title(sample)
    plt.tight_layout()
    fn = '1.FSV_sig_' + sample + '.pdf'
    plt.savefig(os.path.join(outdir, fn))
    plt.clf()

    sig_res = results.query('qval < 0.05').copy()
    sig_res['sample'] = sample

    fn = '2.SVG_filter_' + sample + '.csv'
    if not cal_pattern:
        sig_res.to_csv(os.path.join(outdir, fn),index = False)

    else:
        logging.info('Pattern calculation...')
        pattern_results, patterns = SpatialDE.spatial_patterns(X, counts_n, sig_res, C = 4, l = 0.2, verbosity=1)
        comb_results = sig_res.join(pattern_results.set_index('g'), on='g')
        comb_results.to_csv(os.path.join(outdir, fn),index = False)


        fig, ax = plt.subplots(2, 2, figsize = (6, 6))
        axs = ax.ravel()
        for i, p in enumerate(patterns.columns):
            C = patterns[p]

            axs[i].scatter(X, C, s = 10)
            axs[i].vlines(0.5, min(C), max(C), linestyles = 'dashed', colors = 'r')
            axs[i].set_yticks([])
            axs[i].set_title('Pattern {} - {} genes'.format(i, pattern_results.query('pattern == @i').shape[0]))

        plt.tight_layout()
        fn = '3.Pattern_' + sample + '.pdf'
        plt.savefig(os.path.join(outdir, fn))

        fn = '4.Pattern_members_' + sample + '.csv'
        pattern_results.sort_values('pattern').to_csv(os.path.join(outdir, fn), index = False)



def _PolarizationScoring(aggr_):

    m = np.mean(aggr_['expr'])
    s = np.sum(aggr_['expr'] > 0)
    score = np.sum((aggr_['coords']-0.5) * aggr_['expr'])
    scale1 = np.sum((aggr_['expr'] - m)**2)

    if scale1 == 0:
        scale1 = 1

    scale =  scale1 * np.sum((aggr_['coords'] - 0.5)**2)
    return s * score / scale
