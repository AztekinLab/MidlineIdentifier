from typing import Literal, get_args

import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
import os, sys
import logging

from .image_class import Img
from .adata_class import Adata
from . import utilis
from . import io

logging.basicConfig(stream=sys.stdout, level=logging.INFO, format="%(name)s %(asctime)s %(levelname)s %(message)s")


_Method = Literal["logreg", "t-test", "wilcoxon", "t-test_overestim_var", "DESeq2","DESeq2_pb"]
_CorrMethod = Literal["benjamini-hochberg", "bonferroni"]


class Budoid:
    """
    Budoid object class

    Attributes
    ----------
    img : :class:`~coolname.image_class.Img`
        Image object
    data : :class:`~coolname.adata_class.Adata`
        Single cell object
    sample : :class:`str`
        Sample identifier. Will be store in `adata.obs['sample']`. Useful when concanating multiple Adata objects.
    outdir : :class:`str`
        Output directory where files will be saved


    Methods
    -------
    FindPath
    ADProcess
    RMOutliers
    ProjectCells
    FindOrientation
    run_wrapper
    FindDEG
    FindSVG
    Concat

    """


    def __init__(self, args):
        """
        Constructs all the necessary attributes for the Adata object.

        Parameters
        ----------
        fad : :class:`str`
            Path to single cell object. Will be load into an anndata object.

        sample : :class:`str`
            Sample identifier. Will be store in ``adata.obs['sample']``. Useful when concanating multiple Adata objects.
        outdir : :class:`str`
            Output directory to save files

        """
        fad, fimg, dc, do, outdir, sample = utilis.ParseArgs(args)
        logging.info("Initializing the object using %s" % args.sample)

        self.img = Img(fimg, dc, do, outdir)
        self.data = Adata(fad, sample, outdir)
        self.sample = sample
        self.outdir = outdir

    def FindPath(self, **kwagrs):
        """
        Identify the morphological midline of the structure. A wrapper function of :func:`~coolname.image_class.Img.FindPath`. See :func:`~coolname.image_class.Img.FindPath` for more detail.

        Parameters
        ----------
        kwagrs
            Additonal arguments to pass to :func:`~coolname.image_class.Img.FindPath`
        """

        self.img.FindPath(**kwagrs)

    def ADProcess(self):
        """
        Preprocessing of single cell dataset. A wrapper function of :func:`~coolname.adata_class.Aata.Preprocessing`. See :func:`~coolname.adata_class.Aata.Preprocessing` for more detail.

        """

        self.data.Preprocessing()


    # load image to filter cells outside of the boundary
    # possibly remove addtional cells on the boundary, future fix
    def RMOutliers(self, plot=True):

        """
        Remove cells that fall out of the structure segmentation.

        Parameters
        ----------
        plot : :class:`bool`
            If True, save the plot into ``'Cells_remove.pdf'`` in the output directory (``.outdir``)

        """

        adata = self.data.adata

        x = adata.obsm['centers'][:,0]
        y = adata.obsm['centers'][:,1]
        keep_cells = self.img.seg[x, y].nonzero()[0]
        self.data.adata = adata[keep_cells,:]

        if plot:
            fig, ax = plt.subplots()
            c = ['b' if i in keep_cells else 'r' for i in range(len(adata)) ]
            ax.scatter(x, y, c = c)
            plt.savefig(os.path.join(self.outdir, 'Cells_remove.pdf'))



    def ProjectCells(self, alpha = 0.01, plot = True):
        """
        Project cells onto the nearest coordinate on the morphological midline. We developed a scoring scheme which takes into account the distance between coordinates and cells and the number of cells associated with the coordinates. The score of coordinate-cell pair :math:`(i,c)` is defined as

        .. math::
            S_{ic}  = D_{ic}  e^{αN_{i}}

        where :math:`D_{ic}` represents the Euclidian distance, :math:`N_i` is the number of cells associated with :math:`i` and :math:`α` is the scaling factor. Each cell was then projected to the coordinate with the highest score.

        Parameters
        ----------
        alpha : :class:`float` (default: `0.01`)
            alpha (:math:`α`) that control the level of penalty.

        plot : :class:`bool` (default: `True`)
            If True, save teh plot into ``'Cells_remove.pdf'`` in the output directory (``.outdir``)


        """

        adata, path = self.data.adata, self.img.path

        cells = adata.obsm['centers']
        prob_take, num_take = np.ones(len(path)), np.zeros(len(path))


        np_list, pairs = [], []
        # find match loc on path for each cell
        for cell in tqdm(cells, desc = 'Projecting each cell to the path'):

            dist = utilis.EuclideanDist(path, cell)
            dist = dist / prob_take

            idx = np.argmin(dist)
            np_ = path[idx]
            np_list.append(np_)
            pairs.append([cell, np_])

            # update prop
            # lower the prop if chosen multiple times
            num_take[idx] += 1
            prob_take[idx] = np.exp(-alpha * num_take[idx])


        x, y = path[:,0], path[:,1]
        np_list = np.array(np_list)
        if (max(x)- min(x)) > (max(y)- min(y)):
            logging.info('Major axis is more vertical')
            major_coor = np_list[:,0]
            self.tag = 'vertical'
            # start = path[np.where(x == start)[0][0]]
            # end = path[np.where(x == end)[0][0]]
        else:
            logging.info('Major axis is more horizontal')
            major_coor = np_list[:,1]
            self.tag = 'horizontal'
            # start = path[np.where(y == start)[0][0]]
            # end = path[np.where(y == end)[0][0]]

        adata.obs['major_coor'] = major_coor
        adata.obs['major_coor_scaled'] = utilis.ScaleMinMax(major_coor)

        self.data.adata = adata
        # self.start = start
        # self.end = end


        if plot:
            fig, ax = plt.subplots(1, 2, figsize = (8, 4), width_ratios=[1, 2])

            ax[0].hist(num_take)
            ax[0].set_xlabel('#Cells matched to each position')
            ax[0].set_ylabel('Frequency')

            ax[1].scatter(*zip(*cells))

            for pt in pairs:
                ax[1].plot(*zip(*pt), c = 'm', linewidth = 0.5)
            ax[1].plot(*zip(*path), c = 'y')

            plt.tight_layout()
            plt.savefig(os.path.join(self.outdir, 'Cell_projection.pdf'))

    def FindOrientation(self, **kwargs):
        """
        Orient the coords based on the provided genelists. A wrapper function of :func:`~coolname.adata_class.Adata.FindOrientation`. See :func:`~coolname.adata_class.Adata.FindOrientation` for more detail.

        Parameters
        ----------
        kwagrs
            Additonal arguments to pass to :func:`~coolname.adata_class.Adata.FindOrientation`

        Examples
        --------
        To define the proximal (start) and distal (end) ends of the midline using an example datset. If the dataset with both proximal score and distal score greater than self-defined threshold (Thre = 0.01), it will be considered as polarized; Otherwise, it will be considered as non-polarized.

        .. highlight:: python
        .. code-block:: python

            >>> import PSUils as ps
            >>> budoid = ps.io.ReadObj('testdata/Budoid_1A/Budoids.pkl')

            >>> start_genes = ['Sox9','Acan','Col2a1','Col9a1','Col9a2','Col11a1']
            >>> end_genes = ['Col1a1', 'Col3a1']
            >>> coords = 'major_coor_scaled' # previsouly stored midline coordinates

            >>> budoid.FindOrientation(coords, start_genes, end_genes)

            >>> adata = budoid.data.adata
            >>> Thre = 0.01 # self-defined threshold
            >>> max_s, max_e = adata.uns['start_score'], adata.uns['end_score']

            >>> if max_s > Thre and max_e > Thre:
            ... idx = adata.obs['major_coor_used'] > 0.5
            ... adata.obs.loc[idx, 'loc'] = 'Proximal'
            ... adata.obs.loc[idx, 'loc'] = 'Distal'
            >>> else:
                adata.obs['loc'] = 'Round'

        """

        self.data.FindOrientation(**kwargs)


    def FindSVG(self, coords, sample = 'sample', **kwargs):
        """
        Finds spatially variable genes (SVGs) for each of the identity classes in a dataset. This should be done on the sample level. A wrapper function of :func:`~coolname.adata_class.Adata.FindSVG`.

        Parameters
        ----------
        sample : :class:`str` (Default: `sample`)
            Sample identifier. Must be one of the `.obs.columns`
        kwagrs
            Additonal arguments to pass to :func:`~coolname.adata_class.Adata.FindSVG`


        """
        ALL_ = self.data.adata.obs.columns

        if (sample is None) or (sample not in ALL_):
            raise ValueError(f"sample={sample} must be one of {ALL_}.")

        self.data.FindSVG(coords, sample, **kwargs)



    def FindDEG(
        self,
        groupby,
        condition, # condition for each sample in groupby, required for DESeq2 pseudobulk method
        method: _Method = "DESeq2",
        corr_method: _CorrMethod = "benjamini-hochberg",
        **kwagrs):

        """
        Finds differentially expressed genes (DEGs) for each of the identity classes in a dataset. A wrapper function of :func:`~coolname.adata_class.Adata.FindDEG`

        Parameters
        ----------
        groupby : :class:`str`
            The key of the observations grouping to consider.
        method : :class:`str` (default: `'DESeq2'`)
            Method used to calcualte DEGs.
            ``DESeq2`` and ``DESeq2_pb`` apply `pydeseq2 <https://github.com/owkin/PyDESeq2?tab=readme-ov-file>`_, the python implementation of the `DESeq2 <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8>`_ method. ``DESeq2`` calculates DEGs on single cell level while ``DESeq2_pb`` generate pseudobulk expression based on ``condition``.
            ``t-test``, ``'t-test_overestim_var'`` overestimates variance of each group, ``'wilcoxon'`` uses Wilcoxon rank-sum, ``'logreg'`` uses logistic regression.

            If method is one of ``['logreg', 't-test', 'wilcoxon', 't-test_overestim_var']``, This function directly calls :func:`scanpy.tl.rank_genes_groups`.
        condition : :class:`str`
            Required for ``DESeq2_pb`` method.
        kwagrs
            Additonal arguments to pass to :func:`scanpy.tl.rank_genes_groups`


        Examples
        --------
        .. highlight:: python
        .. code-block:: python

            >>> import PSUils as ps
            >>> budoid = ps.io.ReadObj('testdata/Budoid_1A/Budoids.pkl')

            >>> groupby = 'condition'
            >>> cond = 'loc'
            >>> budoid.data.adata.obs

            >>> # test DESeq2_pb method
            >>> budoid.FindDEG(groupby, cond, method = 'DESeq2_pb', groups = 'P', reference = 'D')

            >>> # test wilcoxon method
            >>> budoid.FindDEG(groupby, method = 'wilcoxon')
        """

        ALL_ = self.data.adata.obs.columns

        if (groupby is None) or (groupby not in ALL_):
            raise ValueError(f"groupby={groupby} must be one of {ALL_}.")

        if method is None:
            method = "DESeq2"

        avail_methods = set(get_args(_Method))
        if method not in avail_methods:
            raise ValueError(f"Method must be one of {avail_methods}.")

        avail_corr = set(get_args(_CorrMethod))
        if corr_method not in avail_corr:
            raise ValueError(f"Correction method must be one of {avail_corr}.")


        if method == "DESeq2_pb":
            if condition is None:
                raise ValueError(f"Condition information is required in DESeq2 pseudobulk method.")
            if not isinstance(condition, str) or (condition not in ALL_):
                raise ValueError(f"condition = {condition} must be one of {ALL_}.")

        return self.data.FindDEG(groupby, condition, method, corr_method, **kwagrs)



    def run_wrapper(self, save = True, **kwagrs):
        """
        A wrapper function to process the datset.

        Parameters
        ----------
        save : :class:`bool` (default: `True`)
            If `True`, save the processed data into pickle file.
        kwagrs
            Additonal arguments to pass to :func:`~coolname.io.SaveObj`

        """
        logging.info("Running wrapper function...")
        logging.info("Finding path...")
        self.FindPath()
        logging.info("Processing single cell dataset...")
        self.ADProcess()
        logging.info("Removing cells fall out of the structure...")
        self.RMOutliers()
        logging.info("Projecting cells on the path...")
        self.ProjectCells()
        logging.info("Identitying the orientation of the path...")
        self.FindOrientation()

        if save:
            logging.info("Saving the object...")
            io.SaveObj(self, **kwagrs)
        logging.info("Done!")


    def Concat(self, object_list):
        """
        Merge multiple objects.

        Parameters
        ----------
        object_list : :class:`list` of :class:`~coolname.Budoids_class.Budoid`
            A list of :class:`~coolname.Budoids_class.Budoid` to merge

        """
        for attribute, value in self.__dict__.items():
            if attribute == 'outdir':
                pass
            elif attribute == 'data':
                for obj in object_list:
                    del getattr(obj.data, 'adata').obsm

                merge_adata = value.adata.concatenate([getattr(obj.data, 'adata') for obj in object_list])
                setattr(self.data, 'adata', merge_adata)
            else:
                setattr(self, attribute, [value] + [getattr(obj, attribute) for obj in object_list])
