'''

This script aims to provide examples of the functionalities of coolname

'''

import MidlineIdentifier as MIL
import argparse
import sys, os
from anndata import ImplicitModificationWarning

# Necessary!
# print(sys.getrecursionlimit()) # 1000
sys.setrecursionlimit(5000)

import warnings

warnings.filterwarnings("default", category=FutureWarning)
warnings.filterwarnings("default", category=ImplicitModificationWarning)


parser = argparse.ArgumentParser()
parser.add_argument('--dir', default = None, help = 'Full path to the directory containing image and anndata files. If provided, separate path to files will be ignored')
parser.add_argument('--img', default = None, help = 'Full path to the segamentation image')
parser.add_argument('--ad', default = None, help = 'Full path to the anndata object')

# Find path
parser.add_argument('--diskClosing', default = 50, type=int, help = 'Disk size for closing')
parser.add_argument('--diskOpening', default = 20, type=int, help = 'Disk size for opening')

parser.add_argument('--sample', default = None, help = 'Name of the object')
parser.add_argument('--outdir', default = None, help = 'Dictrectory to the output. If not provided, output files will be stored in the same folder as the input file')


args = parser.parse_args()



# budoid = MIL.Budoids_class.Budoid(args)

##========= run wrapper
# budoid.run_wrapper()


##========= run step by step
# budoid.img.Segmentation()
# budoid.img.FindRidge()
# budoid.img.FindPath()

# budoid.adata.Preprocessing()

# budoid.RMOutliers()
# budoid.ProjectCells()
# budoid.DefineDirection()




# budoid1 = MIL.io.ReadObj('data/Budoid_1A/Budoids.pkl')
# budoid2 = MIL.io.ReadObj('data/Budoid_3H/Budoids.pkl')

##========= DEG test
# groupby = 'condition'
# cond = 'loc' # same len of groupby.unique()
# method = 'DESeq2_pb'


# adata_comp =  budoid.data.adata
# adata_comp.obs['condition'] = adata_comp.obs['loc']
# adata_comp.obs['condition'][:400] = adata_comp.obs['condition'][:400] + '_1'
# adata_comp.obs['condition'][400:900] = adata_comp.obs['condition'][400:900] + '_2'
# adata_comp.obs['condition'][900:] = adata_comp.obs['condition'][900:] + '_3'


# print(budoid.data.adata.obs)
# print(budoid.FindDEG(groupby, cond, method, groups = ['P'], reference = 'D'))



##========= SVG test
# budoid1.FindSVG(coords = 'major_coor_used', sample = 'loc', min_exp_gene = 10, min_exp_cell = 10)



##========= concat and plotting
# budoid1.Concat([budoid2])
# MIL.plotting.trend_plot(budoid1, features = 'Acan', groupby = 'batch', save = True)