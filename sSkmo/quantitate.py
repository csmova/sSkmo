import pandas as pd
import matplotlib.pyplot as plt

# METHOD=quantitate; type = iBAQ
# TODO: add argument/additional function for normalization by user input protein
def ibaq(df, mydict):
    '''
    iBAQ: loop through each unique protein, sum all the peptides that correspond to that protein,
    add those summed values to a new df, divide by # possible peptides for that protein(len_dict)

    Parameters:
        df (DataFrame): data to quantitate
    Returns:
        df_ibaq (DataFrame): data normalized by iBAQ
    '''
    proteins = df['Protein.Ids'].unique()
    tmplist = []

    for x in proteins:
        if ';' not in x:
            tmpdf = df[df['Protein.Ids'] == x]
            tmpser = tmpdf.loc[:, tmpdf.columns[3:]].sum()
            tmplist.append(tmpser.to_frame(name=x).T / mydict[x])

    df_ibaq = pd.concat(tmplist)
    return df_ibaq

def normalize_by_protein(df, protein):
    '''Normalizes data by given protein. Caution: drops samples without any protein present'''
    df_norm = [col for col in df if df[col].loc[protein]>0]
    df_norm = df[df_norm]
    df_norm = df_norm.div(df_norm.loc[protein])
    assert(all(df_norm.loc[protein]==1))
    print('Dropped ', str(df.shape[1]-df_norm.shape[1]), ' samples containing 0 ', protein)
    return df_norm

def count_zeros_by_plate_columns(df):
    """
    This function returns the number of zeros and the percentage of zeros in each sample grouped by plate number when plate numbers are in the column headers.

    Parameters:
    df (pandas.DataFrame): The DataFrame with plate numbers in the column headers.

    Returns:
    pandas.DataFrame: A DataFrame containing the count and percentage of zeros for each plate number.
    """
    # Count the number of zeros in each column
    zero_counts = (df == 0).sum()

    # Extract the plate number from the column headers
    zero_counts.index = zero_counts.index.str.extract(r'(P\d+)')[0]

    # Group by plate number and sum the zero counts
    zero_counts_by_plate = zero_counts.groupby(zero_counts.index).sum()

    # Calculate the total counts for each plate number
    total_counts_by_plate = df.columns.str.extract(r'(P\d+)')[0].value_counts() * len(df)

    # Calculate the percentage of zeros
    percentage_zeros_by_plate = (zero_counts_by_plate / total_counts_by_plate) * 100

    # Create a DataFrame with the counts and percentages
    result = pd.DataFrame({
        'ZeroCount': zero_counts_by_plate,
        'Percentage': percentage_zeros_by_plate
    })

    return result


def get_most_abundant(df):
    '''prints list of the number of samples in which a protein was the most abundant'''
    tmplist = []
    for x in df.columns:
        tmplist.append(df[x].sort_values(ascending=False).iloc[:1, ])
    tmplist

    # print number of samples with most abundant protein
    tmpdict = {}
    for x in df.index:
        tmpdict[x] = sum([1 if x in i else 0 for i in tmplist])

    for key, value in sorted(tmpdict.items(), key=lambda x: x[1], reverse=True):
        if value > 0:
            print("{} : {}".format(key, value))
    print("All remaining proteins were never most abundant in provided samples.")


def get_redundant_pairs(df):
    '''Get diagonal and lower triangular pairs of correlation matrix'''
    pairs_to_drop = set()
    cols = df.columns
    for i in range(0, df.shape[1]):
        for j in range(0, i+1):
            pairs_to_drop.add((cols[i], cols[j]))
    return pairs_to_drop

def get_top_abs_correlations(df, n=5):
    au_corr = df.corr().abs().unstack()
    labels_to_drop = get_redundant_pairs(df)
    au_corr = au_corr.drop(labels=labels_to_drop).sort_values(ascending=False)
    return au_corr[0:n]

def get_highest_correlations(df,top_n=5,mthd='spearman'):
    df_corr=df.T.corr(method=mthd)
    plt.matshow(df_corr)
    print("Top Absolute Correlations")
    top_abs_correlations=pd.DataFrame(get_top_abs_correlations(df_corr,top_n))
    print(top_abs_correlations)
    #TODO: add ability to search by analyte
    return top_abs_correlations


def get_pct_isoform(df, protein_code,exclude=''):
    # extract all MYH rows
    for i in protein_code:
        tmpdf=0
        tmpdf=df.filter(like=i, axis=0).astype('float')
        #if exclude:#TODO: optional exlclusion
        tmpdf=tmpdf.div(tmpdf.sum())
        return tmpdf

def get_pct_isoform_by_regex(df, protein_family):
    res=pd.DataFrame(columns=df.columns)
    for i in protein_family.keys():
        tmpdf=0
        tmpdf=df.filter(items=protein_family[i], axis=0).astype('float')
        tmpdf=tmpdf.div(tmpdf.sum())
        res=pd.concat([res,tmpdf])
    return res

def old_young_ratio(scdata, leiden_group):
    '''returns the ratio of old:(young+old) fibers in each Leiden cluster'''
    tmp_obs = scdata.obs[scdata.obs['leiden'] == leiden_group]
    count_young = (tmp_obs['group'] == 'young').sum()
    count_old = (tmp_obs['group'] == 'old').sum()

    return count_old / (count_young + count_old)


def get_myh_hybrid_counts_by_age(df, ages, fiber_types):
    ##find a way to get how many and what type mixed for all young and old

    res = pd.DataFrame(columns=ages, index=df.loc[fiber_types].unique())
    res.fillna(0, inplace=True)
    for col in df.columns:
        for ft in df.loc[fiber_types].unique():
            if df[col].loc[fiber_types] == ft:
                for age in ages:
                    if df[col].loc['group'] == age:
                        res[age].loc[ft] += 1

    return res