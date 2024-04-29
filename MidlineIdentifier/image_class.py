from skimage.io import imread
from skimage.morphology import binary_closing, binary_opening, skeletonize, disk
from skimage.measure import label, regionprops, regionprops_table

from skimage.filters import  meijering
from skimage.feature import corner_harris, corner_peaks
from skimage.graph import route_through_array

from scipy.ndimage import binary_fill_holes, distance_transform_edt

import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from .utilis import EuclideanDist

class Img:

    """
    Image object

    Attributes
    ----------
    img : :class:`np.ndarray`
        Original image
    imgb : :class:`np.ndarray`
        Binary image
    dc : :class:`int`
        Disk size for image closing
    do : :class:`int`
        Disk size for image openning
    paddings : :class:`numpy.ndarray`
        Record if paddings has been performed for the four borders of the image.
    outdir : :class:`str`
        output directory to save files

    ar : :class:`float`
        The aspect ratio of the major structure
    measure : :class:`pandas.DataFrame`
    seg : :class:`numpy.ndarray`
        The segmentation of the major structure
    ridge : :class:`numpy.ndarray`
        The segmentation of the major structure
    path : :class:`numpy.ndarray`
        The morphological midline of the major structure in the image
    start : :class:`tuple`
        The start of the morphological midline
    end : :class:`tuple`
        The end of the morphological midline


    Methods
    -------
    Padding
    Segmentation
    Image_measurement
    RMSmallRegion
    GetStartEnd
    GetAspectRatio
    FindRidge
    FindPath

    """


    def __init__(self, fimg, dc, do, outdir):

        """
        Constructs all the necessary attributes for the Img object.

        Parameters
        ----------

        fimg : :class:`str`
            Path to image
        dc : :class:`int` (default: `50`)
            Disk size for image closing
        do : :class:`int` (default: `20`)
            Disk size for image openning
        outdir : :class:`str`
            Output directory to save files

        """

        self.img = imread(fimg)
        self.imgb = self.img > 1
        self.dc = dc
        self.do = do
        self.paddings = [False, False, False, False] # t, b, l, r
        self.outdir = outdir


    def Padding(self, tolerance = 50, num_pads = 100):
        """
        Image padding. Paddings will be performed if the target structure locates too close to the image borders. This will allow better structure segmentation.

        If padding is performed, the `paddings` attribute will be modified accordingly.

        Parameters
        ----------
        tolerance : :class:`int` (default: `50`)
            The number of pixel to tolerant.
        num_pads : :class:`int` (default: `100`)
            The number of pixel to pad


        """

        all_trues = self.imgb.nonzero()

        edge_t, edge_b = all_trues[0].min(), all_trues[0].max()
        edge_l, edge_r = all_trues[1].min(), all_trues[1].max()

        self.npads = num_pads
        self.imgp = self.imgb
        # top padding
        if edge_t < tolerance:
            pads_x = np.full((num_pads, self.imgp.shape[1]), False) # padding size may change during padding
            self.imgp = np.vstack((pads_x, self.imgp))
            self.paddings[0] = True

        # bottom padding
        if (self.imgp.shape[0] - edge_b) < tolerance:
            pads_x = np.full((num_pads, self.imgp.shape[1]), False)
            self.imgp = np.vstack((self.imgp, pads_x))
            self.paddings[1] = True

        # left padding
        if edge_l < tolerance:
            pads_y = np.full((self.imgp.shape[0], num_pads), False)
            self.imgp = np.hstack((pads_y, self.imgp))
            self.paddings[2] = True

        # right padding
        if  (self.imgp.shape[1] - edge_r) < tolerance:
            pads_y = np.full((self.imgp.shape[0], num_pads), False)
            self.imgp = np.hstack((self.imgp, pads_y))
            self.paddings[3] = True



    def Segmentation(self):
        """
        Get a closed segmentation of the major structure

        """

        self.Padding()
        seg = binary_closing(self.imgp, disk(self.dc))
        seg = binary_fill_holes(seg).astype(int)
        seg = binary_opening(seg, disk(self.do))

        # to prevent removing too many peripheral cells
        if self.do > 20:
            seg = dilation(seg, disk(5))

        # remove padding from segmentation
        if self.paddings[0]:
            seg = seg[self.npads:,:]
        if self.paddings[1]:
            seg = seg[:-self.npads,:]
        if self.paddings[2]:
            seg = seg[:,self.npads:]
        if self.paddings[3]:
            seg = seg[:,:-self.npads]


        self.seg = seg

        self.Image_measurement()
        self.RMSmallRegion()



    def Image_measurement(self, plot = True):
        """
        Measure properties of all image regions. Measurements will be stored in ``.measure``

        Parameters
        ----------
        plot : :class:`bool` (default: `True`)
            If True, save teh plot into ``'region_measure.pdf'`` in the output directory (``.outdir``)
        """

        label_img = label(self.seg)
        df = pd.DataFrame(regionprops_table(label_img, properties=('centroid', 'orientation', 'axis_major_length', 'axis_minor_length', 'bbox')))
        df.sort_values(by = ['axis_major_length'], ascending = False, inplace = True, ignore_index=True)
        self.measure = df

        if plot:
            fig, ax = plt.subplots()
            ax.imshow(self.seg, cmap=plt.cm.gray)

            for props in regionprops(label_img):
                y0, x0 = props.centroid
                orientation = props.orientation
                x1 = x0 + math.cos(orientation) * 0.5 * props.axis_minor_length
                y1 = y0 - math.sin(orientation) * 0.5 * props.axis_minor_length
                x2 = x0 - math.sin(orientation) * 0.5 * props.axis_major_length
                y2 = y0 - math.cos(orientation) * 0.5 * props.axis_major_length

                ax.plot((x0, x1), (y0, y1), '--r', linewidth=2.5)
                ax.plot((x0, x2), (y0, y2), '-r', linewidth=2.5)
                ax.plot(x0, y0, '.g', markersize=15)

                minr, minc, maxr, maxc = props.bbox
                bx = (minc, maxc, maxc, minc, minc)
                by = (minr, minr, maxr, maxr, minr)
                ax.plot(bx, by, '-b', linewidth=2.5)

            plt.savefig(os.path.join(self.outdir, 'region_measure.pdf'))


    def RMSmallRegion(self):
        """
        Remove small regions (noise) by setting the corresponding pixels to false. Only the largest segment was kept for structural segmentation.

        """

        if self.measure is None:
            self.Image_measurement()

        df, seg = self.measure, self.seg

        if len(df) > 1:
            for _, df_sub in df[1:].iterrows():
                x1 = int(df_sub['bbox-0'])
                y1 = int(df_sub['bbox-1'])
                x2 = int(df_sub['bbox-2'])
                y2 = int(df_sub['bbox-3'])
                seg[x1:x2, y1:y2] = False

        self.seg = seg


    def GetStartEnd(self):
        """
        Return the start and end of the mophological midline.

        Returns
        -------
        start : :class:`tuple` of `( :class:`int`, :class:`int`)`
            Start pixel coordiantes of the morphological midline
        end : :class:`tuple` of `( :class:`int`, :class:`int`)`
            End pixel coordiantes of the morphological midline
        """

        df = self.measure

        ori = df['orientation'][0]
        maj_ax = df['axis_major_length'][0]
        x_start = df['centroid-0'][0]
        y_start = df['centroid-1'][0]

        x1 = int(x_start - math.cos(ori) * maj_ax )
        y1 = int(y_start - math.sin(ori) * maj_ax )
        x2 = int(x_start + math.cos(ori) * maj_ax )
        y2 = int(y_start + math.sin(ori) * maj_ax )

        return (x1, y1), (x2, y2)


    def GetAspectRatio(self):
        """
        Return the aspect ratio of the structure.

        Returns
        -------
        AspectRatio : :class:`float`
            Aspect ratio of the major structure

        """
        if self.ar is None:
            long = self.df['axis_major_length'][0]
            short =  self.df['axis_minor_length'][0]
            self.ar = round(long/short, 2)

        return self.ar


    # def Skeletonization(self, pruning = True):

    #     # skeletonization on ridge
    #     edt = distance_transform_edt(self.seg)
    #     ridge_f = filters.meijering(edt, black_ridges = False)
    #     sk = skeletonize(ridge_f > 0.1)

    #     if pruning:
    #         sk = skel_pruning_DSE(sk, edt, 100)

    #     self.sk = sk

    def FindRidge(self, **kwargs):
        """
        Filter the Euclidean distance transform of the image with the Meijering neuriteness filter. This function calls :func:`scipy.ndimage.distance_transform_edt` and :func:`skimage.filters.meijering`.

        Parameters
        ----------
        kwargs
            Additonal arguments to pass to :func:`skimage.filters.meijering`
        """

        edt = distance_transform_edt(self.seg)
        ridge = meijering(edt, black_ridges = False, **kwargs)

        self.ridge = ridge


    # get all possible path of the pruned skeleton
    def FindPath(self, plot = True):
        """
        Identify the morphological midline of the major structure in the image.The midline will be store in ``.path``.

        Parameters
        ----------
        plot : :class:`bool` (default: `True`)
            If True, save teh plot into ``'Major_axis_sk_on_ridge.pdf'`` in the output directory (``.outdir``)

        kwargs
            Additonal arguments to pass to :func:`skimage.filters.meijering`
        """

        if not hasattr(self, 'seg'):
            self.Segmentation()

        if not hasattr(self, 'ridge'):
            self.FindRidge()

        corners = _CornerDetector(self.seg, min_distance = 30, correct = False)
        (x1, y1), (x2, y2) = self.GetStartEnd()

        dist1 = EuclideanDist(corners, [x1, y1])
        dist2 = EuclideanDist(corners, [x2, y2])

        pt1 = corners[np.argmin(dist1)]
        pt2 = corners[np.argmin(dist2)]

        path, _ = route_through_array(1 - self.ridge, pt1, pt2)

        self.path = np.array(path)
        self.start, self.end = pt1, pt2

        if plot:
            path_r = [tup[::-1] for tup in path] # reverse to plot
            corners_i = [[y, x] for x, y in corners]

            fig, axes = plt.subplots(2, 2, figsize=(8, 8))
            ax = axes.ravel()

            ax[0].set_title('Original picture')
            ax[0].imshow(self.imgb, cmap=plt.cm.gray)
            ax[0].set_axis_off()

            ax[1].set_title('Transformed picture')
            ax[1].imshow(self.seg, cmap=plt.cm.gray)
            ax[1].plot(*zip(*path_r))
            ax[1].set_axis_off()

            ax[2].set_title('Ridge')
            ax[2].imshow(self.ridge, cmap=plt.cm.gray)
            ax[2].scatter(*zip(*corners_i), s = 1, c = 'red')

            ax[3].set_title('Ridge with path')
            ax[3].imshow(self.ridge, cmap=plt.cm.gray)
            ax[3].scatter(y1, x1, c = 'r')
            ax[3].scatter(y2, x2, c = 'y')
            ax[3].scatter(pt1[1], pt1[0], c = 'r', marker = '*')
            ax[3].scatter(pt2[1], pt2[0],c = 'y', marker = '*')
            ax[3].plot(*zip(*path_r))
            ax[3].set_axis_off()

            plt.tight_layout()
            plt.savefig(os.path.join(self.outdir, 'Major_axis_sk_on_ridge.pdf'))
            plt.clf()



def _CornerDetector(mat, correct = False, threshold_abs = 0, min_distance = 15, **kwargs):

    corners = corner_peaks(corner_harris(mat), threshold_abs = threshold_abs,min_distance = min_distance, **kwargs)

    if correct:
        true_idx = np.argwhere(mat)
        corrected = []

        for corner in tqdm(corners, total = len(corners), desc = 'Correcting corners'):
            if mat[corner[0], corner[1]] == False:
                dist = np.linalg.norm(true_idx - corner, ord=2, axis=1.)
                corrected.append(true_idx[np.argmin(dist)])
            else:
                corrected.append(corner)

        corners = np.array(corrected)

    return corners
