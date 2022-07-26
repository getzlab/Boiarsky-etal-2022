import numpy as np
import pandas as pd 
import scanpy as sc
import scipy
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from adjustText import adjust_text

#for scanpy functions
import collections.abc as cabc
from itertools import product
from typing import Optional, Union, Mapping  # Special
from pandas.api.types import is_categorical_dtype
    
# function to convert to subscript
def get_sub(x):
        normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
        sub_s = "ₐ₈CDₑբGₕᵢⱼₖₗₘₙₒₚQᵣₛₜᵤᵥwₓᵧZₐ♭꜀ᑯₑբ₉ₕᵢⱼₖₗₘₙₒₚ૧ᵣₛₜᵤᵥwₓᵧ₂₀₁₂₃₄₅₆₇₈₉₊₋₌₍₎"
        res = x.maketrans(''.join(normal), ''.join(sub_s))
        return x.translate(res)
        
        
#for plotting top genes on NMF sigs based on W matrix
def marker_scatter(R, W, FE, sig,nlabel=10,ax=None,fontsize= 'x-small', figheight=5, figwidth=5):
    if ax is None:
        f,ax = plt.subplots(1)
        f.set_figheight(figheight)
        f.set_figwidth(figwidth)
    
    marker_genes=R[sig].sort_values(ascending=False).index[0:nlabel]

    ax.scatter(W[sig],FE[sig],4)

    t=[ax.text(W.loc[m,sig],FE.loc[m,sig],m,weight="bold", fontsize=fontsize) for m in marker_genes]
    adjust_text(t,arrowprops=dict(arrowstyle='->', color='red'),ax=ax)

    ax.set_xlabel('Gene Loading on Signature (Fraction of Sig)',fontsize=14)
    ax.set_ylabel('Specificity of Gene to Signature',fontsize=14)
    ax.set_ylim(0,1)

# functions to calculate and plot DEGs (some are original, some are edited versions of scanpy fxns, and some are verbatim internal scanpy fxns that are required to run the edited fxns)
def prepare_rankgenesgroupsdf(df):
    #create columns with neglogp, absolute score, and "my_score" which is product of neglogp and log fold change
    df['negative_log_padj'] = np.minimum(300, -np.log10(df.pvals_adj))
    df['abs_score'] = np.absolute(df.scores)
    df['my_score'] = np.absolute(df['logfoldchanges']*df['negative_log_padj'])
    df['direction'] = np.where(df.logfoldchanges>0, 'up', 'down')
    return(df)


def plot_all_markers_volcano(marker_genes, n_annot=5, ncol=4, pos_only=False, fontsize=11, groupname="unknown", rankby_col="my_score", logfc_thresh=np.log2(1.5), pval_adj_thresh=0.1, height=4, arrows=True):
    #ncol is number of columns in the facet grid
    #n_annot is the number of genes to annotate with text. chosen by rankby_col (default "my_score" column)
    #groupname is the name to assign to the "cluster" in the event that there is only one group in the marker_genes df, and thus no columns named 'group' (so it needs to be added)
    #height=height of each facet, passed to sns.FacetGrid
    marker_genes = prepare_rankgenesgroupsdf(marker_genes)
    marker_genes['passthresh'] = (marker_genes.pvals_adj<pval_adj_thresh) & (np.abs(marker_genes.logfoldchanges)>logfc_thresh) #add column denoting whether gene passes thresholds, to use for coloring
    if ~np.isin('group', marker_genes.columns):
        marker_genes['group'] = groupname
    g = sns.FacetGrid(marker_genes, col="group",  col_wrap=ncol, sharex=False, sharey=False, col_order = marker_genes.group.sort_values().unique(), height=height, 
                      hue="passthresh", palette={True:'orange', False:'blue'})
    g = g.map(plt.scatter, "logfoldchanges", "negative_log_padj", edgecolor="w")
    g.set_ylabels("-log{}(q-value)".format(get_sub('10'))); 
    g.set_xlabels("log{}(fold change)".format(get_sub('2')));
    axes_map = dict(zip(marker_genes.group.sort_values().unique(), np.arange(len(marker_genes.group.sort_values().unique()))))#which cluster is on which axis
    ## label top DEGs (MIGHT BREAK IF THERE ARE NO SIGNIF DEG FOR A FACET, NOT SURE)
    top_marker_genes=marker_genes[marker_genes.passthresh].sort_values(rankby_col, ascending=False).groupby('group').head(n_annot) #only label if pass thresholds
    #if you want to force it to plot a specific set of genes: top_marker_genes = pd.concat([top_marker_genes,marker_genes[marker_genes.names=="LILRB4"]])
    try: #handle the case where the groups you are comparing have names that cant be cast to int
        clusts = np.unique(top_marker_genes.group.astype('int').sort_values())
    except:
        clusts = top_marker_genes.group.sort_values().unique()
    for i in np.arange(len(clusts)): 
        this_clust_top_marker_genes = top_marker_genes.loc[top_marker_genes.group == str(clusts[i]),:]
        ax = g.axes[axes_map[str(clusts[i])]] 
        ax.axvline(c="grey", ls='--') #vertical line at 0
        if pos_only == False: #if positive and negative markers, center plot
            xmax = np.max(np.absolute(marker_genes.loc[marker_genes.group == str(clusts[i])].logfoldchanges)) # calculate xlim for subplot
            ax.set_xlim(-xmax, xmax)
        xs = list(this_clust_top_marker_genes.logfoldchanges)
        ys = list(this_clust_top_marker_genes.negative_log_padj)
        text = list(this_clust_top_marker_genes.names)
        texts=[]
        for x, y, s in zip(xs, ys, text):
            texts.append(ax.text(x, y, s, size=fontsize))
        if arrows:
            adjust_text(texts, only_move={'points':'y', 'texts':'y'}, arrowprops=dict(arrowstyle="->", color='r', lw=0.5), ax=ax)
        else:
            adjust_text(texts, only_move={'points':'y', 'texts':'y'}, ax=ax)
    
    
def _prepare_dataframe(
    adata,
    var_names,
    groupby: Optional[str] = None,
    use_raw: Optional[bool] = None,
    log: bool = False,
    num_categories: int = 7,
    layer=None,
    gene_symbols: Optional[str] = None,
):
    """
    Given the anndata object, prepares a data frame in which the row index are the categories
    defined by group by and the columns correspond to var_names.
    Parameters
    ----------
    adata
        Annotated data matrix.
    var_names
        `var_names` should be a valid subset of  `adata.var_names`.
    groupby
        The key of the observation grouping to consider. It is expected that
        groupby is a categorical. If groupby is not a categorical observation,
        it would be subdivided into `num_categories`.
    use_raw
        Use `raw` attribute of `adata` if present.
    log
        Use the log of the values
    num_categories
        Only used if groupby observation is not categorical. This value
        determines the number of groups into which the groupby observation
        should be subdivided.
    gene_symbols
        Key for field in .var that stores gene symbols.
    Returns
    -------
    Tuple of `pandas.DataFrame` and list of categories.
    """
    from scipy.sparse import issparse

    adata._sanitize()
    if use_raw is None and adata.raw is not None:
        use_raw = True
    if isinstance(var_names, str):
        var_names = [var_names]

    if groupby is not None:
        if groupby not in adata.obs_keys():
            raise ValueError(
                'groupby has to be a valid observation. '
                f'Given {groupby}, valid observations: {adata.obs_keys()}'
            )

    if gene_symbols is not None and gene_symbols in adata.var.columns:
        # translate gene_symbols to var_names
        # slow method but gives a meaningful error if no gene symbol is found:
        translated_var_names = []
        for symbol in var_names:
            if symbol not in adata.var[gene_symbols].values:
                print(
                    f"Gene symbol {symbol!r} not found in given "
                    f"gene_symbols column: {gene_symbols!r}"
                )
                return
            translated_var_names.append(
                adata.var[adata.var[gene_symbols] == symbol].index[0]
            )
        symbols = var_names
        var_names = translated_var_names
    if layer is not None:
        if layer not in adata.layers.keys():
            raise KeyError(
                f'Selected layer: {layer} is not in the layers list. '
                f'The list of valid layers is: {adata.layers.keys()}'
            )
        matrix = adata[:, var_names].layers[layer]
    elif use_raw:
        matrix = adata.raw[:, var_names].X
    else:
        matrix = adata[:, var_names].X

    if issparse(matrix):
        matrix = matrix.toarray()
    if log:
        matrix = np.log1p(matrix)

    obs_tidy = pd.DataFrame(matrix, columns=var_names)
    if groupby is None:
        groupby = ''
        categorical = pd.Series(np.repeat('', len(obs_tidy))).astype('category')
    else:
        if not is_categorical_dtype(adata.obs[groupby]):
            # if the groupby column is not categorical, turn it into one
            # by subdividing into  `num_categories` categories
            categorical = pd.cut(adata.obs[groupby], num_categories)
        else:
            categorical = adata.obs[groupby]

    obs_tidy.set_index(categorical, groupby, inplace=True)
    if gene_symbols is not None:
        # translate the column names to the symbol names
        obs_tidy.rename(
            columns=dict([(var_names[x], symbols[x]) for x in range(len(var_names))]),
            inplace=True,
        )
    categories = obs_tidy.index.categories

    return categories, obs_tidy
    
    
def my_filter_rank_genes_groups(
    #edited from scanpy version to handle both up and downregulated genes correctly
    adata,
    key=None,
    groupby=None,
    use_raw=True,
    log=True,
    key_added='rank_genes_groups_filtered',
    min_in_group_fraction=0.25,
    min_fold_change=2,
    max_out_group_fraction=0.5,
) -> None:
    """\
    Filters out genes based on fold change and fraction of genes expressing the
    gene within and outside the `groupby` categories.
    See :func:`~scanpy.tl.rank_genes_groups`.
    Results are stored in `adata.uns[key_added]`
    (default: 'rank_genes_groups_filtered').
    To preserve the original structure of adata.uns['rank_genes_groups'],
    filtered genes are set to `NaN`.
    Parameters
    ----------
    adata
    key
    groupby
    use_raw
    log
        If true, it means that the values to work with are in log scale
    key_added
    min_in_group_fraction
    min_fold_change
    max_out_group_fraction
    Returns
    -------
    Same output as :func:`scanpy.tl.rank_genes_groups` but with filtered genes names set to
    `nan`
    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, 'bulk_labels', method='wilcoxon')
    >>> sc.tl.filter_rank_genes_groups(adata, min_fold_change=3)
    >>> # visualize results
    >>> sc.pl.rank_genes_groups(adata, key='rank_genes_groups_filtered')
    >>> # visualize results using dotplot
    >>> sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups_filtered')
    """
    if key is None:
        key = 'rank_genes_groups'

    if groupby is None:
        groupby = str(adata.uns[key]['params']['groupby'])

    # convert structured numpy array into DataFrame
    gene_names = pd.DataFrame(adata.uns[key]['names'])

    fraction_in_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )
    fold_change_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )
    fraction_out_cluster_matrix = pd.DataFrame(
        np.zeros(gene_names.shape),
        columns=gene_names.columns,
        index=gene_names.index,
    )
    print(
        f"Filtering genes using: "
        f"min_in_group_fraction: {min_in_group_fraction} "
        f"min_fold_change: {min_fold_change}, "
        f"max_out_group_fraction: {max_out_group_fraction}"
    )

    for cluster in gene_names.columns:
        # iterate per column
        var_names = gene_names[cluster].values

        # add column to adata as __is_in_cluster__. This facilitates to measure
        # fold change of each gene with respect to all other clusters
        adata.obs['__is_in_cluster__'] = pd.Categorical(adata.obs[groupby] == cluster)

        # obs_tidy has rows=groupby, columns=var_names
        categories, obs_tidy = _prepare_dataframe(
            adata,
            var_names,
            groupby='__is_in_cluster__',
            use_raw=use_raw,
        )

        # for if category defined by groupby (if any) compute for each var_name
        # 1. the mean value over the category
        # 2. the fraction of cells in the category having a value > 0

        # 1. compute mean value
        mean_obs = obs_tidy.groupby(level=0).mean()

        # 2. compute fraction of cells having value >0
        # transform obs_tidy into boolean matrix
        obs_bool = obs_tidy.astype(bool)

        # compute the sum per group which in the boolean matrix this is the number
        # of values >0, and divide the result by the total number of values in the group
        # (given by `count()`)
        fraction_obs = obs_bool.groupby(level=0).sum() / obs_bool.groupby(level=0).count()

        # Because the dataframe groupby is based on the '__is_in_cluster__' column,
        # in this context, [True] means __is_in_cluster__.
        # Also, in this context, fraction_obs.loc[True].values is the row of values
        # that is assigned *as column* to fraction_in_cluster_matrix to follow the
        # structure of the gene_names dataFrame
        fraction_in_cluster_matrix.loc[:, cluster] = fraction_obs.loc[True].values
        fraction_out_cluster_matrix.loc[:, cluster] = fraction_obs.loc[False].values

        # compute fold change.
        if log:
            fold_change_matrix.loc[:, cluster] = (np.exp(mean_obs.loc[True]) / np.exp(mean_obs.loc[False])).values
        else:
            fold_change_matrix.loc[:, cluster] = (mean_obs.loc[True] / mean_obs.loc[False]).values

    # remove temporary columns
    adata.obs.drop(columns='__is_in_cluster__')
    # filter original_matrix
    gene_names = gene_names[
        (fraction_in_cluster_matrix > min_in_group_fraction) &
        (fraction_out_cluster_matrix < max_out_group_fraction) &
        (np.absolute(np.log2(fold_change_matrix)) > np.log2(min_fold_change)) #Rebecca: needs to be in log scale so that up and downregulated magnitudes are comparable
    ]
    # create new structured array using 'key_added'.
    adata.uns[key_added] = adata.uns[key].copy()
    adata.uns[key_added]['names'] = gene_names.to_records(index=False)

def my_rank_genes_groups(
    adata,
    groupby: str,
    zero_offset, #added by Rebecca; should be set by user to half of lowest non-zero value in dataset, instead of using scanpy's 1e-9
    groups = "all",
    use_raw: bool = True,
    reference: str = 'rest',
    n_genes: Optional[int] = None,
    rankby_abs: bool = False,
    pts: bool = False,
    key_added: Optional[str] = None,
    copy: bool = False,
    method = None,
    corr_method = 'benjamini-hochberg',
    tie_correct: bool = False,
    layer: Optional[str] = None,
    **kwds,
):
    """\
    Rank genes for characterizing groups.
    Expects logarithmized data.
    Parameters
    ----------
    adata
        Annotated data matrix.
    groupby
        The key of the observations grouping to consider.
    use_raw
        Use `raw` attribute of `adata` if present.
    layer
        Key from `adata.layers` whose value will be used to perform tests on.
    groups
        Subset of groups, e.g. [`'g1'`, `'g2'`, `'g3'`], to which comparison
        shall be restricted, or `'all'` (default), for all groups.
    reference
        If `'rest'`, compare each group to the union of the rest of the group.
        If a group identifier, compare with respect to this group.
    n_genes
        The number of genes that appear in the returned tables.
        Defaults to all genes.
    method
        The default method is `'t-test'`,
        `'t-test_overestim_var'` overestimates variance of each group,
        `'wilcoxon'` uses Wilcoxon rank-sum,
        `'logreg'` uses logistic regression. See [Ntranos18]_,
        `here <https://github.com/theislab/scanpy/issues/95>`__ and `here
        <http://www.nxn.se/valent/2018/3/5/actionable-scrna-seq-clusters>`__,
        for why this is meaningful.
    corr_method
        p-value correction method.
        Used only for `'t-test'`, `'t-test_overestim_var'`, and `'wilcoxon'`.
    tie_correct
        Use tie correction for `'wilcoxon'` scores.
        Used only for `'wilcoxon'`.
    rankby_abs
        Rank genes by the absolute value of the score, not by the
        score. The returned scores are never the absolute values.
    pts
        Compute the fraction of cells expressing the genes.
    key_added
        The key in `adata.uns` information is saved to.
    **kwds
        Are passed to test methods. Currently this affects only parameters that
        are passed to :class:`sklearn.linear_model.LogisticRegression`.
        For instance, you can pass `penalty='l1'` to try to come up with a
        minimal set of genes that are good predictors (sparse solution meaning
        few non-zero fitted coefficients).
    Returns
    -------
    **names** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the gene
        names. Ordered according to scores.
    **scores** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the z-score
        underlying the computation of a p-value for each gene for each
        group. Ordered according to scores.
    **logfoldchanges** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Structured array to be indexed by group id storing the log2
        fold change for each gene for each group. Ordered according to
        scores. Only provided if method is 't-test' like.
        Note: this is an approximation calculated from mean-log values.
    **pvals** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        p-values.
    **pvals_adj** : structured `np.ndarray` (`.uns['rank_genes_groups']`)
        Corrected p-values.
    **pts** : `pandas.DataFrame` (`.uns['rank_genes_groups']`)
        Fraction of cells expressing the genes for each group.
    **pts_rest** : `pandas.DataFrame` (`.uns['rank_genes_groups']`)
        Only if `reference` is set to `'rest'`.
        Fraction of cells from the union of the rest of each group
        expressing the genes.
    Notes
    -----
    There are slight inconsistencies depending on whether sparse
    or dense data are passed. See `here <https://github.com/theislab/scanpy/blob/master/scanpy/tests/test_rank_genes_groups.py>`__.
    Examples
    --------
    >>> import scanpy as sc
    >>> adata = sc.datasets.pbmc68k_reduced()
    >>> sc.tl.rank_genes_groups(adata, 'bulk_labels', method='wilcoxon')
    >>> # to visualize the results
    >>> sc.pl.rank_genes_groups(adata)
    """
    if method is None:
        print(
            "Default of the method has been changed to 't-test' from 't-test_overestim_var'"
        )
        method = 't-test'

    if 'only_positive' in kwds:
        rankby_abs = not kwds.pop('only_positive')  # backwards compat

    avail_methods = {'t-test', 't-test_overestim_var', 'wilcoxon', 'logreg'}
    if method not in avail_methods:
        raise ValueError(f'Method must be one of {avail_methods}.')

    avail_corr = {'benjamini-hochberg', 'bonferroni'}
    if corr_method not in avail_corr:
        raise ValueError(f'Correction method must be one of {avail_corr}.')

    adata = adata.copy() if copy else adata
    adata._sanitize() #_utils.sanitize_anndata(adata)
    # for clarity, rename variable
    if groups == 'all':
        groups_order = 'all'
    elif isinstance(groups, (str, int)):
        raise ValueError('Specify a sequence of groups')
    else:
        groups_order = list(groups)
        if isinstance(groups_order[0], int):
            groups_order = [str(n) for n in groups_order]
        if reference != 'rest' and reference not in set(groups_order):
            groups_order += [reference]
    if reference != 'rest' and reference not in adata.obs[groupby].cat.categories:
        cats = adata.obs[groupby].cat.categories.tolist()
        raise ValueError(
            f'reference = {reference} needs to be one of groupby = {cats}.'
        )

    if key_added is None:
        key_added = 'rank_genes_groups'
    adata.uns[key_added] = {}
    adata.uns[key_added]['params'] = dict(
        groupby=groupby,
        reference=reference,
        method=method,
        use_raw=use_raw,
        layer=layer,
        corr_method=corr_method,
    )

    test_obj = _RankGenes(adata, groups_order, groupby, reference, use_raw, layer, pts)

    # for clarity, rename variable
    n_genes_user = n_genes
    # make sure indices are not OoB in case there are less genes than n_genes
    # defaults to all genes
    if n_genes_user is None or n_genes_user > test_obj.X.shape[1]:
        n_genes_user = test_obj.X.shape[1]

    test_obj.compute_statistics(
        method, zero_offset, corr_method, n_genes_user, rankby_abs, tie_correct, **kwds
    )

    if test_obj.pts is not None:
        groups_names = [str(name) for name in test_obj.groups_order]
        adata.uns[key_added]['pts'] = pd.DataFrame(
            test_obj.pts.T, index=test_obj.var_names, columns=groups_names
        )
    if test_obj.pts_rest is not None:
        adata.uns[key_added]['pts_rest'] = pd.DataFrame(
            test_obj.pts_rest.T, index=test_obj.var_names, columns=groups_names
        )

    test_obj.stats.columns = test_obj.stats.columns.swaplevel()

    dtypes = {
        'names': 'O',
        'scores': 'float32',
        'logfoldchanges': 'float32',
        'pvals': 'float64',
        'pvals_adj': 'float64',
    }

    for col in test_obj.stats.columns.levels[0]:
        adata.uns[key_added][col] = test_obj.stats[col].to_records(
            index=False, column_dtypes=dtypes[col]
        )

    return adata if copy else None


class _RankGenes:
    def __init__(
        self,
        adata,
        groups,
        groupby,
        reference='rest',
        use_raw=True,
        layer=None,
        comp_pts=False,
    ):

        if 'log1p' in adata.uns_keys() and adata.uns['log1p']['base'] is not None:
            self.expm1_func = lambda x: np.expm1(x * np.log(adata.uns['log1p']['base']))
        else:
            self.expm1_func = np.expm1

        self.groups_order, self.groups_masks = sc._utils.select_groups(
            adata, groups, groupby
        )

        # Singlet groups cause division by zero errors
        invalid_groups_selected = set(self.groups_order) & set(
            adata.obs[groupby].value_counts().loc[lambda x: x < 2].index
        )

        if len(invalid_groups_selected) > 0:
            raise ValueError(
                "Could not calculate statistics for groups {} since they only "
                "contain one sample.".format(', '.join(invalid_groups_selected))
            )

        adata_comp = adata
        if layer is not None:
            if use_raw:
                raise ValueError("Cannot specify `layer` and have `use_raw=True`.")
            X = adata_comp.layers[layer]
        else:
            if use_raw and adata.raw is not None:
                adata_comp = adata.raw
            X = adata_comp.X

        # for correct getnnz calculation
        if scipy.sparse.issparse(X):
            X.eliminate_zeros()

        self.X = X
        self.var_names = adata_comp.var_names

        self.ireference = None
        if reference != 'rest':
            self.ireference = np.where(self.groups_order == reference)[0][0]

        self.means = None
        self.vars = None

        self.means_rest = None
        self.vars_rest = None

        self.comp_pts = comp_pts
        self.pts = None
        self.pts_rest = None

        self.stats = None

        # for logreg only
        self.grouping_mask = adata.obs[groupby].isin(self.groups_order)
        self.grouping = adata.obs.loc[self.grouping_mask, groupby]

    def _basic_stats(self):
        n_genes = self.X.shape[1]
        n_groups = self.groups_masks.shape[0]

        self.means = np.zeros((n_groups, n_genes))
        self.vars = np.zeros((n_groups, n_genes))
        self.pts = np.zeros((n_groups, n_genes)) if self.comp_pts else None

        if self.ireference is None:
            self.means_rest = np.zeros((n_groups, n_genes))
            self.vars_rest = np.zeros((n_groups, n_genes))
            self.pts_rest = np.zeros((n_groups, n_genes)) if self.comp_pts else None
        else:
            mask_rest = self.groups_masks[self.ireference]
            X_rest = self.X[mask_rest]
            self.means[self.ireference], self.vars[self.ireference] = sc.pp._simple._get_mean_var(
                X_rest
            )
            # deleting the next line causes a memory leak for some reason
            del X_rest

        if scipy.sparse.issparse(self.X):
            get_nonzeros = lambda X: X.getnnz(axis=0)
        else:
            get_nonzeros = lambda X: np.count_nonzero(X, axis=0)

        for imask, mask in enumerate(self.groups_masks):
            X_mask = self.X[mask]

            if self.comp_pts:
                self.pts[imask] = get_nonzeros(X_mask) / X_mask.shape[0]

            if self.ireference is not None and imask == self.ireference:
                continue

            self.means[imask], self.vars[imask] = sc.pp._simple._get_mean_var(X_mask)

            if self.ireference is None:
                mask_rest = ~mask
                X_rest = self.X[mask_rest]
                self.means_rest[imask], self.vars_rest[imask] = sc.pp._simple._get_mean_var(X_rest)
                # this can be costly for sparse data
                if self.comp_pts:
                    self.pts_rest[imask] = get_nonzeros(X_rest) / X_rest.shape[0]
                # deleting the next line causes a memory leak for some reason
                del X_rest

    def t_test(self, method):
        from scipy import stats

        self._basic_stats()

        for group_index, mask in enumerate(self.groups_masks):
            if self.ireference is not None and group_index == self.ireference:
                continue

            mean_group = self.means[group_index]
            var_group = self.vars[group_index]
            ns_group = np.count_nonzero(mask)

            if self.ireference is not None:
                mean_rest = self.means[self.ireference]
                var_rest = self.vars[self.ireference]
                ns_other = np.count_nonzero(self.groups_masks[self.ireference])
            else:
                mean_rest = self.means_rest[group_index]
                var_rest = self.vars_rest[group_index]
                ns_other = self.X.shape[0] - ns_group

            if method == 't-test':
                ns_rest = ns_other
            elif method == 't-test_overestim_var':
                # hack for overestimating the variance for small groups
                ns_rest = ns_group
            else:
                raise ValueError('Method does not exist.')

            # TODO: Come up with better solution. Mask unexpressed genes?
            # See https://github.com/scipy/scipy/issues/10269
            with np.errstate(invalid="ignore"):
                scores, pvals = stats.ttest_ind_from_stats(
                    mean1=mean_group,
                    std1=np.sqrt(var_group),
                    nobs1=ns_group,
                    mean2=mean_rest,
                    std2=np.sqrt(var_rest),
                    nobs2=ns_rest,
                    equal_var=False,  # Welch's
                )

            # I think it's only nan when means are the same and vars are 0
            scores[np.isnan(scores)] = 0
            # This also has to happen for Benjamini Hochberg
            pvals[np.isnan(pvals)] = 1

            yield group_index, scores, pvals

    def wilcoxon(self, tie_correct):
        from scipy import stats

        self._basic_stats()

        n_genes = self.X.shape[1]
        # First loop: Loop over all genes
        if self.ireference is not None:
            # initialize space for z-scores
            scores = np.zeros(n_genes)
            # initialize space for tie correction coefficients
            if tie_correct:
                T = np.zeros(n_genes)
            else:
                T = 1

            for group_index, mask in enumerate(self.groups_masks):
                if group_index == self.ireference:
                    continue

                mask_rest = self.groups_masks[self.ireference]

                n_active = np.count_nonzero(mask)
                m_active = np.count_nonzero(mask_rest)

                if n_active <= 25 or m_active <= 25:
                    print(
                        'Few observations in a group for '
                        'normal approximation (<=25). Lower test accuracy.'
                    )

                # Calculate rank sums for each chunk for the current mask
                for ranks, left, right in _ranks(self.X, mask, mask_rest):
                    scores[left:right] = np.sum(ranks.iloc[0:n_active, :])
                    if tie_correct:
                        T[left:right] = _tiecorrect(ranks)

                std_dev = np.sqrt(
                    T * n_active * m_active * (n_active + m_active + 1) / 12.0
                )

                scores = (
                    scores - (n_active * ((n_active + m_active + 1) / 2.0))
                ) / std_dev
                scores[np.isnan(scores)] = 0
                pvals = 2 * stats.distributions.norm.sf(np.abs(scores))

                yield group_index, scores, pvals
        # If no reference group exists,
        # ranking needs only to be done once (full mask)
        else:
            n_groups = self.groups_masks.shape[0]
            scores = np.zeros((n_groups, n_genes))
            n_cells = self.X.shape[0]

            if tie_correct:
                T = np.zeros((n_groups, n_genes))

            for ranks, left, right in _ranks(self.X):
                # sum up adjusted_ranks to calculate W_m,n
                for imask, mask in enumerate(self.groups_masks):
                    scores[imask, left:right] = np.sum(ranks.iloc[mask, :])
                    if tie_correct:
                        T[imask, left:right] = _tiecorrect(ranks)

            for group_index, mask in enumerate(self.groups_masks):
                n_active = np.count_nonzero(mask)

                if tie_correct:
                    T_i = T[group_index]
                else:
                    T_i = 1

                std_dev = np.sqrt(
                    T_i * n_active * (n_cells - n_active) * (n_cells + 1) / 12.0
                )

                scores[group_index, :] = (
                    scores[group_index, :] - (n_active * (n_cells + 1) / 2.0)
                ) / std_dev
                scores[np.isnan(scores)] = 0
                pvals = 2 * stats.distributions.norm.sf(np.abs(scores[group_index, :]))

                yield group_index, scores[group_index], pvals

    def logreg(self, **kwds):
        # if reference is not set, then the groups listed will be compared to the rest
        # if reference is set, then the groups listed will be compared only to the other groups listed
        from sklearn.linear_model import LogisticRegression

        # Indexing with a series causes issues, possibly segfault
        X = self.X[self.grouping_mask.values, :]

        if len(self.groups_order) == 1:
            raise ValueError('Cannot perform logistic regression on a single cluster.')

        clf = LogisticRegression(**kwds)
        clf.fit(X, self.grouping.cat.codes)
        scores_all = clf.coef_
        for igroup, _ in enumerate(self.groups_order):
            if len(self.groups_order) <= 2:  # binary logistic regression
                scores = scores_all[0]
            else:
                scores = scores_all[igroup]

            yield igroup, scores, None

            if len(self.groups_order) <= 2:
                break

    def compute_statistics(
        self,
        method,
        zero_offset, #added by Rebecca; defaults to 1e-9 but should be set by user to half of lowest non-zero value in dataset
        corr_method='benjamini-hochberg',
        n_genes_user=None,
        rankby_abs=False,
        tie_correct=False,
        **kwds,
    ):

        if method in {'t-test', 't-test_overestim_var'}:
            generate_test_results = self.t_test(method)
        elif method == 'wilcoxon':
            generate_test_results = self.wilcoxon(tie_correct)
        elif method == 'logreg':
            generate_test_results = self.logreg(**kwds)

        self.stats = None

        n_genes = self.X.shape[1]

        for group_index, scores, pvals in generate_test_results:
            group_name = str(self.groups_order[group_index])

            if n_genes_user is not None:
                scores_sort = np.abs(scores) if rankby_abs else scores
                global_indices = _select_top_n(scores_sort, n_genes_user)
                first_col = 'names'
            else:
                global_indices = slice(None)
                first_col = 'scores'

            if self.stats is None:
                idx = pd.MultiIndex.from_tuples([(group_name, first_col)])
                self.stats = pd.DataFrame(columns=idx)

            if n_genes_user is not None:
                self.stats[group_name, 'names'] = self.var_names[global_indices]

            self.stats[group_name, 'scores'] = scores[global_indices]

            if pvals is not None:
                self.stats[group_name, 'pvals'] = pvals[global_indices]
                if corr_method == 'benjamini-hochberg':
                    from statsmodels.stats.multitest import multipletests

                    pvals[np.isnan(pvals)] = 1
                    _, pvals_adj, _, _ = multipletests(
                        pvals, alpha=0.05, method='fdr_bh'
                    )
                elif corr_method == 'bonferroni':
                    pvals_adj = np.minimum(pvals * n_genes, 1.0)
                self.stats[group_name, 'pvals_adj'] = pvals_adj[global_indices]

            if self.means is not None:
                mean_group = self.means[group_index]
                if self.ireference is None:
                    mean_rest = self.means_rest[group_index]
                else:
                    mean_rest = self.means[self.ireference]
                #add offset to numerator and denominator prior to calculating fold change
                foldchanges = (self.expm1_func(mean_group) + zero_offset) / ( #1e-9
                    self.expm1_func(mean_rest) + zero_offset ) #1e-9 # add small value to remove 0's
                # REBECCA EDIT - use half of lowest value to replace num&denom when zeros encountered in fold change calculation (current code replaces zeros, but that could flip sign of fold change if numerator is very small, so should update to replace num & denom)
                #numerators = np.where(self.expm1_func(mean_group) == 0, zero_offset, self.expm1_func(mean_group))
                #denominators = np.where(self.expm1_func(mean_rest) == 0, zero_offset, self.expm1_func(mean_rest))
                #foldchanges = numerators/ denominators
                self.stats[group_name, 'logfoldchanges'] = np.log2(
                    foldchanges[global_indices]
                )

        if n_genes_user is None:
            self.stats.index = self.var_names
        
def _ranks(X, mask=None, mask_rest=None):
    from math import floor
    
    CONST_MAX_SIZE = 10000000

    n_genes = X.shape[1]

    if scipy.sparse.issparse(X):
        merge = lambda tpl: scipy.sparse.vstack(tpl).toarray()
        adapt = lambda X: X.toarray()
    else:
        merge = np.vstack
        adapt = lambda X: X

    masked = mask is not None and mask_rest is not None

    if masked:
        n_cells = np.count_nonzero(mask) + np.count_nonzero(mask_rest)
        get_chunk = lambda X, left, right: merge(
            (X[mask, left:right], X[mask_rest, left:right])
        )
    else:
        n_cells = X.shape[0]
        get_chunk = lambda X, left, right: adapt(X[:, left:right])

    # Calculate chunk frames
    max_chunk = floor(CONST_MAX_SIZE / n_cells)

    for left in range(0, n_genes, max_chunk):
        right = min(left + max_chunk, n_genes)

        df = pd.DataFrame(data=get_chunk(X, left, right))
        ranks = df.rank()
        yield ranks, left, right


def _select_top_n(scores, n_top):
    n_from = scores.shape[0]
    reference_indices = np.arange(n_from, dtype=int)
    partition = np.argpartition(scores, -n_top)[-n_top:]
    partial_indices = np.argsort(scores[partition])[::-1]
    global_indices = reference_indices[partition][partial_indices]

    return global_indices

#Bayesian purity estimation    
import scipy.integrate as integrate

def _integrand_kappa_n(kappa_n, mean_kappa_normal, std_kappa_normal, rho, kappa_t, num_kappa_cells, num_total_cells):
    prior_kappa_n = scipy.stats.norm.pdf(kappa_n, mean_kappa_normal, std_kappa_normal)
    p_l = rho*kappa_t + (1-rho)*kappa_n
    likelihood_K = scipy.stats.binom.pmf(num_kappa_cells, num_total_cells, p_l)
    integrand_kappa_n = likelihood_K*prior_kappa_n
    return integrand_kappa_n

def estimate_purity_lc_only(cell_lightchain_type_df, sample_col, this_sample,
                            mean_kappa_normal, std_kappa_normal, rho_granularity=.01, a=1, b=1):
    #rho_granularity is used for calculating the cdf for different slices of the distribution over rho
    num_kappa_cells = cell_lightchain_type_df[cell_lightchain_type_df[sample_col]==this_sample].kappa_bool.sum()
    num_total_cells = len(cell_lightchain_type_df[cell_lightchain_type_df[sample_col]==this_sample])

    rhos=[]
    marginals_rho=[]
    for rho in np.arange(rho_granularity,1+rho_granularity,rho_granularity): #will calculate the posterior cdf for rho in [0,0.01], [0.01,0.02],...[0.99, 1.00]
        prior_rho = scipy.stats.beta.cdf(rho, a, b) - scipy.stats.beta.cdf(rho-rho_granularity, a, b) #default a, b assumes uniform prior on rho
        
        marginal_rho = 0
        for kappa_t in [0,1]: #to get the marginal P(rho|K,N), we sum P(rho, kappa_t|K,N) over possible values of kappa_t 
            prior_kappa_t=0.5 #could update this to match the literature for how likely an MM tumor is to be kappa; for now, uniform     
            integral_over_kappa_n = integrate.quad(_integrand_kappa_n, 0, 1, args=(mean_kappa_normal, std_kappa_normal, rho, kappa_t, num_kappa_cells, num_total_cells))[0]         
            marginal_rho += prior_kappa_t*prior_rho*integral_over_kappa_n
        
        rhos.append(rho)
        marginals_rho.append(marginal_rho)
        
    marginals_rho = marginals_rho/np.sum(marginals_rho) #normalize prob to sum to 1
    return rhos, marginals_rho

#plots of gene expression and signature activity per healthy/malig cells per patient
def plot_per_patient_healthymalig(summary_df, things_to_plot, facet_col='signature', value_col='mean_activity', sample_col='person', label_col="ground_truth", normalcells="healthy plasma", ds_col="disease_stage", ds_colors=['cornflowerblue', '#F1CE46', '#DC7209', '#880B0B'], ylabel='mean activity', ncols=1, ptsize=6, custom_ylim=None, sharey=False, sharex=False, figsize=None, legend_fontsize="small", filename=None, dpi=150):
    #sigs_to_plot: a list of signatures or genes to plot. must contain a subset of the values in the facet_col column of your summary_df
    #facet_col: column name in summary_df which defines the facets
    #value_col: column name in summary_df for column that contained values that will be plotted
    
    assert isinstance(things_to_plot,list), "things_to_plot must be a list (even if it only contains one element)"

    #initialize figure+subplots
    if figsize is None:
        figsize = (10,5*len(things_to_plot))
    nrows = int(np.ceil(len(things_to_plot)/ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharey=False, sharex=False, figsize=figsize, squeeze=False)
    plt.xticks(np.arange(0, len(summary_df[sample_col].drop_duplicates())+3, 1.0), rotation=90) #additional 4 (or 3 if now distinction b.w SMMl/h) blank spaces b/w disease stages
    plt.subplots_adjust(hspace=0.3) #increase vertical space between plots otherwise labels may overlap
    if custom_ylim is not None:
        plt.ylim(custom_ylim)
    
    if 'stderror' not in summary_df.columns: #enables user to pass df without calculating error bars
        summary_df['stderror'] = None

    #plot each signature on a subplot
    sig_counter = 0
    for i in np.arange(nrows):
        for j in np.arange(ncols):
            if sig_counter > len(things_to_plot)-1:
                break;
            ax = axes[i,j]

            ## plot aesthetics
            #despine
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_ylabel(ylabel)
            ax.set_xlabel('patient ID')
            ax.grid(False)

            blue_patch = mpatches.Patch(color=ds_colors[0], label='normal plasma cells')
            mgus_patch = mpatches.Patch(color=ds_colors[1], label='MGUS neoplastic cells')
            smm_patch = mpatches.Patch(color=ds_colors[2], label='SMM neoplastic cells')
            mm_patch = mpatches.Patch(color=ds_colors[3], label='MM neoplastic cells')

            plt.legend(handles=[blue_patch, mgus_patch, smm_patch, mm_patch], loc="best", fontsize=legend_fontsize)

            ## signature specific settings
            #get signature
            sig = things_to_plot[sig_counter] #signature for this facet

            #patients in descending order
            this_sig_df = summary_df[summary_df[facet_col]==sig].copy()
            #order by disease stage + descending
            this_sig_df[ds_col] = pd.Categorical(this_sig_df[ds_col], categories=[
                'MM','SMM','MGUS','NBM'], ordered=True)
            this_sig_df.sort_values([ds_col,value_col], ascending=False, inplace=True)
            
            #set title, xticklabels
            ax.set_title(sig) #IDK why i cant get this to work!! this_sig_df.title[0])
            myxlabels = []
            for d in this_sig_df[ds_col].drop_duplicates():
                myxlabels+=list(this_sig_df[this_sig_df[ds_col]==d][sample_col].drop_duplicates())
                myxlabels+=['']
            myxlabels=myxlabels[:-1]
            ax.set_xticklabels(myxlabels)        
            #x axis locations per person, for this sig
            tmp_person_xaxis_loc_dict = dict(zip(myxlabels, np.arange(len(myxlabels))))

            ## plot values with error bars
            # 1. plot the points for each person, with error bars
            # order first by disease stage
            # within each disease stage, order by decreasing values (malignant) (perhaps should be able to choose increasing or decreasing)
            # color by malignant vs. healthy
            # plot malignant points in shades of orange, increasing in intensity with disease stage
            
            if np.all(pd.isna(this_sig_df.stderror)):
                print("WARNING: plotting {} without error bars".format(sig))
                
                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="MGUS")][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[ds_col]=="MGUS"),value_col], 
                            yerr=None,
                            fmt='o', color=ds_colors[1], ecolor=ds_colors[1], elinewidth=2, capsize=0);

                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="SMM")][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[ds_col]=="SMM"),value_col], 
                            yerr=None,
                            fmt='o', color=ds_colors[2], ecolor=ds_colors[2], elinewidth=2, capsize=0);

                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="MM")][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[ds_col]=="MM"),value_col], 
                            yerr=None,
                            fmt='o', color=ds_colors[3], ecolor=ds_colors[3], elinewidth=2, capsize=0);

                #plot healthy points in blue
                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[label_col]==normalcells)][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[label_col]==normalcells),value_col], 
                            yerr=None,
                            fmt='o', color=ds_colors[0], ecolor=ds_colors[0], elinewidth=2, capsize=0);
            
            else: #plot with error bars
                
                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="MGUS")][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[ds_col]=="MGUS"),value_col], 
                            yerr=this_sig_df[(this_sig_df[ds_col]=="MGUS")].stderror, 
                            fmt='o', color=ds_colors[1], ecolor=ds_colors[1], elinewidth=2, capsize=0);

                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="SMM")][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[ds_col]=="SMM"),value_col], 
                            yerr=this_sig_df[(this_sig_df[ds_col]=="SMM")].stderror, 
                            fmt='o', color=ds_colors[2], ecolor=ds_colors[2], elinewidth=2, capsize=0);

                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="MM")][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[ds_col]=="MM"),value_col], 
                            yerr=this_sig_df[(this_sig_df[ds_col]=="MM")].stderror, 
                            fmt='o', color=ds_colors[3], ecolor=ds_colors[3], elinewidth=2, capsize=0);

                #plot healthy points in blue
                ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[label_col]==normalcells)][sample_col]], 
                            y=this_sig_df.loc[(this_sig_df[label_col]==normalcells),value_col], 
                            yerr=this_sig_df[(this_sig_df[label_col]==normalcells)].stderror, 
                            fmt='o', color=ds_colors[0], ecolor=ds_colors[0], elinewidth=2, capsize=0);

            # 2. draw line connecting the points
            #find ppl who have both healthy and malignant cells (same for every facet)
            ppl_who_need_lines = this_sig_df[[sample_col,label_col]].drop_duplicates().groupby(sample_col).size()
            ppl_who_need_lines = ppl_who_need_lines.index[ppl_who_need_lines>1]

            for p in ppl_who_need_lines:
                ax.plot([tmp_person_xaxis_loc_dict[p],tmp_person_xaxis_loc_dict[p]], this_sig_df[this_sig_df[sample_col]==p].sort_values("group")[value_col], linestyle="--", color='lightgray', zorder=1)

            ## add horizontal bars over each disease stage
            barheight=ax.get_ylim()[1]
            ax.axhline(y=barheight, 
                       xmin=0+1/(max(tmp_person_xaxis_loc_dict.values())+2),
                       xmax=(1+max([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='NBM'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       color=ds_colors[0]) #over NBM
            ax.axhline(y=barheight, 
                       xmin=(1+min([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='MGUS'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       xmax=(1+max([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='MGUS'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       color=ds_colors[1]) #
            ax.axhline(y=barheight, 
                       xmin=(1+min([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='SMM'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       xmax=(1+max([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='SMM'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       color=ds_colors[2]) #
            ax.axhline(y=barheight, 
                       xmin=(1+min([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='MM'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       xmax=1-1/(max(tmp_person_xaxis_loc_dict.values())+2),
                       color=ds_colors[3]) #


            sig_counter += 1

    if filename is not None:
        plt.savefig(filename, bbox_inches="tight", dpi=dpi)

#adjusted function to work on Amit data (eg. "Normal" instead of "NBM")
def plot_per_patient_healthymalig_Amit(summary_df, things_to_plot, facet_col='signature', value_col='mean_activity', sample_col='person', ds_col="disease_stage", ds_colors=['cornflowerblue', '#F1CE46', '#DC7209', '#880B0B'], ylabel='mean activity', ncols=1, ptsize=6, sharey=False, sharex=False, figsize=None, filename=None):
    #sigs_to_plot: a list of signatures or genes to plot. must contain a subset of the values in the facet_col column of your summary_df
    #facet_col: column name in summary_df which defines the facets
    #value_col: column name in summary_df for column that contained values that will be plotted
    
    assert isinstance(things_to_plot,list), "things_to_plot must be a list (even if it only contains one element)"

    #initialize figure+subplots
    if figsize is None:
        figsize = (10,5*len(things_to_plot))
    nrows = int(np.ceil(len(things_to_plot)/ncols))
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, sharey=False, sharex=False, figsize=figsize, squeeze=False)
    plt.xticks(np.arange(0, len(summary_df[sample_col].drop_duplicates())+3, 1.0), rotation=90) #additional 4 (or 3 if now distinction b.w SMMl/h) blank spaces b/w disease stages
    plt.subplots_adjust(hspace=0.3) #increase vertical space between plots otherwise labels may overlap

    #plot each signature on a subplot
    sig_counter = 0
    for i in np.arange(nrows):
        for j in np.arange(ncols):
            if sig_counter > len(things_to_plot)-1:
                break;
            ax = axes[i,j]

            ## plot aesthetics
            #despine
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.set_ylabel(ylabel)
            ax.set_xlabel('patient ID')
            ax.grid(False)

            blue_patch = mpatches.Patch(color=ds_colors[0], label='normal plasma cells')
            mgus_patch = mpatches.Patch(color=ds_colors[1], label='MGUS neoplastic cells')
            smm_patch = mpatches.Patch(color=ds_colors[2], label='SMM neoplastic cells')
            mm_patch = mpatches.Patch(color=ds_colors[3], label='MM neoplastic cells')

            plt.legend(handles=[blue_patch, mgus_patch, smm_patch, mm_patch], loc="best")

            ## signature specific settings
            #get signature
            sig = things_to_plot[sig_counter] #signature for this facet

            #patients in descending order
            this_sig_df = summary_df[summary_df[facet_col]==sig].copy()
            #order by disease stage + descending
            this_sig_df[ds_col] = pd.Categorical(this_sig_df[ds_col], categories=[
                'MM','SMM','MGUS','Normal'], ordered=True)
            this_sig_df.sort_values([ds_col,value_col], ascending=False, inplace=True)
            
            #set title, xticklabels
            ax.set_title(sig) #IDK why i cant get this to work!! this_sig_df.title[0])
            myxlabels = []
            for d in this_sig_df[ds_col].drop_duplicates():
                myxlabels+=list(this_sig_df[this_sig_df[ds_col]==d][sample_col].drop_duplicates())
                myxlabels+=['']
            myxlabels=myxlabels[:-1]
            ax.set_xticklabels(myxlabels)        
            #x axis locations per person, for this sig
            tmp_person_xaxis_loc_dict = dict(zip(myxlabels, np.arange(len(myxlabels))))

            ## plot values with error bars
            # 1. plot the points for each person, with error bars
            # order first by disease stage
            # within each disease stage, order by decreasing values (malignant) (perhaps should be able to choose increasing or decreasing)
            # color by malignant vs. healthy
            # plot malignant points in shades of orange, increasing in intensity with disease stage
            ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="MGUS")][sample_col]], 
                        y=this_sig_df.loc[(this_sig_df[ds_col]=="MGUS"),value_col], 
                        yerr=this_sig_df[(this_sig_df[ds_col]=="MGUS")].stderror, 
                        fmt='o', color=ds_colors[1], ecolor=ds_colors[1], elinewidth=2, capsize=0);

            ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="SMM")][sample_col]], 
                        y=this_sig_df.loc[(this_sig_df[ds_col]=="SMM"),value_col], 
                        yerr=this_sig_df[(this_sig_df[ds_col]=="SMM")].stderror, 
                        fmt='o', color=ds_colors[2], ecolor=ds_colors[2], elinewidth=2, capsize=0);

            ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df[ds_col]=="MM")][sample_col]], 
                        y=this_sig_df.loc[(this_sig_df[ds_col]=="MM"),value_col], 
                        yerr=this_sig_df[(this_sig_df[ds_col]=="MM")].stderror, 
                        fmt='o', color=ds_colors[3], ecolor=ds_colors[3], elinewidth=2, capsize=0);

            #plot healthy points in blue
            ax.errorbar(x=[tmp_person_xaxis_loc_dict[p] for p in this_sig_df[(this_sig_df.ground_truth=="healthy plasma")][sample_col]], 
                        y=this_sig_df.loc[(this_sig_df.ground_truth=="healthy plasma"),value_col], 
                        yerr=this_sig_df[(this_sig_df.ground_truth=="healthy plasma")].stderror, 
                        fmt='o', color=ds_colors[0], ecolor=ds_colors[0], elinewidth=2, capsize=0);

            # 2. draw line connecting the points
            #find ppl who have both healthy and malignant cells (same for every facet)
            ppl_who_need_lines = this_sig_df[[sample_col,'ground_truth']].drop_duplicates().groupby(sample_col).size()
            ppl_who_need_lines = ppl_who_need_lines.index[ppl_who_need_lines>1]

            for p in ppl_who_need_lines:
                ax.plot([tmp_person_xaxis_loc_dict[p],tmp_person_xaxis_loc_dict[p]], this_sig_df[this_sig_df[sample_col]==p].sort_values("group")[value_col], linestyle="--", color='lightgray', zorder=1)

            ## add horizontal bars over each disease stage
            barheight=ax.get_ylim()[1]
            ax.axhline(y=barheight, 
                       xmin=0+1/(max(tmp_person_xaxis_loc_dict.values())+2),
                       xmax=(1+max([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='Normal'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       color=ds_colors[0]) #over NBM
            ax.axhline(y=barheight, 
                       xmin=(1+min([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='MGUS'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       xmax=(1+max([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='MGUS'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       color=ds_colors[1]) 
            ax.axhline(y=barheight, 
                       xmin=(1+min([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='SMM'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       xmax=(1+max([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='SMM'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       color=ds_colors[2]) #
            ax.axhline(y=barheight, 
                       xmin=(1+min([tmp_person_xaxis_loc_dict[p] for p in this_sig_df[this_sig_df[ds_col]=='MM'][sample_col]]))/(max(tmp_person_xaxis_loc_dict.values())+2), 
                       xmax=1-1/(max(tmp_person_xaxis_loc_dict.values())+2),
                       color=ds_colors[3]) #

            sig_counter += 1

    if filename is not None:
        plt.savefig(filename, bbox_inches="tight")

#functions to calculate standard error to plot on the 'per_patient_healthymalig' plots above

#bootstrap standard error - group by person + ground_truth
#computation time does NOT scale well with the number of genes
def bootstrap_means_genes(adata, n_bootstrap=10000, groupby_vars=['person','ground_truth'], layer='lognorm', genes=None):
    #if genes =None: defaults to highly variable gens from adata
    if genes==None:
        genes = adata.var.index[adata.var.highly_variable]
                          
    tmp = pd.concat([pd.DataFrame(adata[:,genes].layers[layer].todense(), index = adata.obs.index, columns = genes), 
                           adata.obs[groupby_vars]], axis=1) #expression + metadata
    all_means = pd.DataFrame()
    
    for i in np.arange(n_bootstrap):
        #print(i)
        #resample cells per person/ground_truth pair
        #print("sampling...")
        smpld_data = tmp.groupby(groupby_vars).sample(frac=1, replace=True)
        

        #calculate mean on bootstrapped samples
        #print("calculating mean...")
        boot_means = smpld_data.groupby(groupby_vars).mean().dropna(axis=0)
        #print("melting...")
        boot_means = boot_means.melt(var_name='gene', value_name='mean', ignore_index=False).reset_index()
        
        #record bootstrapped sample mean
        #print("concatenating...")
        all_means = pd.concat([all_means, boot_means], axis=0)
    return(all_means)

def bootstrap_means_sigs(adata, sig_colnames, n_bootstrap=10000, groupby_vars=['person','ground_truth']):
    #print("calculating tmp...")
    
    #check that arguments were given as lists
    if not isinstance(sig_colnames, list):
        sig_colnames = [sig_colnames]
    if not isinstance(groupby_vars, list):
        groupby_vars = [groupby_vars]
        
    tmp = adata.obs.loc[:,groupby_vars+sig_colnames] #signature activity + metadata
    all_means = pd.DataFrame()
    
    for i in np.arange(n_bootstrap):
        #print(i)
        #resample cells per person/ground_truth pair
        #print("sampling...")
        smpld_data = tmp.groupby(groupby_vars).sample(frac=1, replace=True)
        

        #calculate mean on bootstrapped samples
        #print("calculating mean...")
        boot_means = smpld_data.groupby(groupby_vars).mean().dropna(axis=0)
        #print("melting...")
        boot_means = boot_means.melt(var_name='signature', value_name='mean', ignore_index=False).reset_index()
        
        #record bootstrapped sample mean
        #print("concatenating...")
        all_means = pd.concat([all_means, boot_means], axis=0)
    return(all_means)

def calc_boot_SE(bootstrapped_means, groupby_vars=['person','ground_truth']):
    #handle both signature and gene inputs
    if "gene" in bootstrapped_means.columns:
        bootstrapped_SE = bootstrapped_means.groupby(groupby_vars+['gene']).std().reset_index().rename(columns={'mean':'stderror'})
    elif "signature" in bootstrapped_means.columns:
        bootstrapped_SE = bootstrapped_means.groupby(groupby_vars+['signature']).std().reset_index().rename(columns={'mean':'stderror'})
    return(bootstrapped_SE)