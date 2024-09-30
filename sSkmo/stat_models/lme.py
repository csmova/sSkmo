import pandas as pd
import numpy as np
import warnings
import anndata as ad
import statsmodels.api as sm
import statsmodels.stats.multitest as smm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

def lme_differential_expression(adata, leiden_group, group1, group2):
    '''builds a linear mixed model comparing young and old groups for given Leiden cluster(s)'''
    # Suppress warnings
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")

        # Subset the data to include only the cells from the specified Leiden cluster
        adata_subset = adata[adata.obs['leiden'].isin(leiden_group)].copy()

        # Prepare lists to store results
        genes = []
        coefs = []
        p_values = []

        for gene in adata.var_names:
            # Prepare the data for the mixed model
            gene_data = pd.DataFrame({'expression': adata_subset[:, gene].X.flatten(),
                                      'group_rep': adata_subset.obs['group_rep'],
                                      'group': adata_subset.obs['group']})

            # Check if there are at least three individuals in at least one of the groups using adata_subset.raw
            raw_gene_data = adata_subset.raw[:, gene].X.flatten()
            raw_gene_data_df = pd.DataFrame({'expression': raw_gene_data,
                                             'group_rep': adata_subset.obs['group_rep'],
                                             'group': adata_subset.obs['group']})

            raw_group_counts = raw_gene_data_df[raw_gene_data_df['expression'] != 0].groupby('group')[
                'group_rep'].nunique()

            if (raw_group_counts.get(group1, 0) < 3) and (raw_group_counts.get(group2, 0) < 3):
                # If the condition is not met, skip this gene
                genes.append(gene)
                coefs.append(np.nan)
                p_values.append(np.nan)
                continue

            # Fit the mixed model
            try:
                mixed_lm = smf.mixedlm("expression ~ group", gene_data, groups=gene_data["group_rep"])
                mixed_lm_result = mixed_lm.fit()

                # Extract the coefficient for the group effect and its p-value
                coef = mixed_lm_result.params[
                    'group[T.' + group2 + ']'] if 'group[T.' + group2 + ']' in mixed_lm_result.params else None
                p_value = mixed_lm_result.pvalues[
                    'group[T.' + group2 + ']'] if 'group[T.' + group2 + ']' in mixed_lm_result.pvalues else None

                # Append results to lists
                genes.append(gene)
                coefs.append(coef)
                p_values.append(p_value)

            except Exception as e:
                print(f"Warning: Model for gene '{gene}' did not converge.")

        # Create DataFrame from lists
        results_df = pd.DataFrame({'gene': genes, 'coef': coefs, 'p_value': p_values})

        # Exclude genes with None p-values
        valid_results_df = results_df[results_df['p_value'].notna()]

        # Adjust p-values for multiple testing using BH correction
        valid_results_df['FDR'] = multipletests(valid_results_df['p_value'], method='fdr_bh')[1]

        # Return the sorted results DataFrame
        return valid_results_df.sort_values('FDR')