import pandas as pd
from anndata import AnnData
from scipy.sparse import csr_matrix
import numpy as np
from itertools import product
from statsmodels.stats.multitest import fdrcorrection

def non_mirrored_product(list1, list2):
    pairs = list(product(list1, list2))
    non_mirrored_pairs = [(x, y) for x, y in pairs if ((x <= y) & (x != y))]
    return non_mirrored_pairs

def onehot_groupby(adata, groupby='cell_type'):
    cts = pd.DataFrame(pd.get_dummies(adata.obs[groupby]).values.astype(np.float32),
                   index = adata.obs.index,
                   columns = pd.get_dummies(adata.obs[groupby]).columns
                   )
    ctdata = AnnData(X=csr_matrix(cts.values), obs=adata.obs, var=pd.DataFrame(index=cts.columns))
    ctdata.obsm = adata.obsm
    ctdata.obsp = adata.obsp
    
    return ctdata

def format_truth(truth, x_name='source', y_name='target'):
    truth['truth'] = fdrcorrection(truth['morans_pvals'], alpha=0.05)[0]
    truth = truth.loc[truth['truth'], [x_name, y_name, 'truth']]
    
    return truth