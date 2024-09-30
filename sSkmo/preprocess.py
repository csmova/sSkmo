import pandas as pd
import string

def drop_qcs(df):
    '''Takes a dataframe, removes all columns where column name contains 'QC', returns a new dataframe.


    Parameters:
        df (DataFrame): The DataFrame to remove columns from.

    Returns:
        df_no_qcs (DataFrame): A new dataframe with columns containing 'QC' removed.
    '''
    # todo: extend to any phrase other than qc?
    qc_cols = [col for col in df.columns if 'QC' in col]
    df_no_qcs = df.drop(qc_cols, axis=1)
    print("Dropped " + str(len(qc_cols)) + " samples containing 'QC'.")
    return df_no_qcs


def assign_fiber_age(df):
    '''adds young/old labels to dataframe'''

    df_yo=df.copy()
    df_yo=df_yo.T

    # young 1
    letters = list(string.ascii_uppercase[:8])
    numbers = list(range(1, 9))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    y1 = ['P7_'+x for x in combined_list]
    #young 1 bottom half
    letters = list(string.ascii_uppercase[8:16])
    numbers = list(range(1, 8))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    y1+= ['P7_'+x for x in combined_list]

    ## young 1 plate 2
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(10,19))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    y1+=['P10_'+x for x in combined_list]

    # old 1
    letters = list(string.ascii_uppercase[:8])
    numbers = list(range(9,17))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    o1 = ['P7_'+x for x in combined_list]
    #young 1 bottom half
    letters = list(string.ascii_uppercase[8:16])
    numbers = list(range(8,17))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    o1+= ['P7_'+x for x in combined_list]

    # young 2
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(1, 9))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    y2 = ['P9_'+x for x in combined_list]

    # young 3
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(1, 9))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    y3 = ['P11_'+x for x in combined_list]

    # young 4
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(1, 9))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    y4 = ['P13_'+x for x in combined_list]

    # young 5
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(9,17))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    y5 = ['P14_'+x for x in combined_list]


    ################## old 2
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(1, 10))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    o2 = ['P10_'+x for x in combined_list]

    # old 3
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(9,18))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    o3 = ['P11_'+x for x in combined_list]

    # old 4
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(9,18))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    o4 = ['P13_'+x for x in combined_list]

    # old 5
    letters = list(string.ascii_uppercase[:16])
    numbers = list(range(1,9))
    combined_list = [f"{letter}{number}" for letter in letters for number in numbers]
    o5 = ['P14_'+x for x in combined_list]

    all_lists = [y1, y2, y3, y4, y5, o1, o2, o3, o4, o5]
    #all_lists

    fixed_lists = []
    for l in all_lists:
        fixed_lists.append([x+'_' for x in l])


    df_yo['group']=''
    df_yo['rep']=''
    df_yo['group_rep']=''

    ## label the group obs column as old or young and replicate as rep number
    for i in range(0,5):
        tmp_matching_indexes = [idx for idx in df_yo.index if any(search_str in idx for search_str in fixed_lists[i])]
        df_yo.loc[tmp_matching_indexes,'group'] = 'young'
        df_yo.loc[tmp_matching_indexes,'rep'] = i+1
        df_yo.loc[tmp_matching_indexes,'group_rep'] = 'y'+str(i+1)

    for i in range(5,10):
        tmp_matching_indexes = [idx for idx in df_yo.index if any(search_str in idx for search_str in fixed_lists[i])]
        df_yo.loc[tmp_matching_indexes,'group'] = 'old'
        df_yo.loc[tmp_matching_indexes,'rep'] = i+1
        df_yo.loc[tmp_matching_indexes,'group_rep']= 'o'+str(i+1)

    df_yo=df_yo.T
    return df_yo, fixed_lists


def label_hybrid_fibers(df, iso_pairs, label):
    ###label pure and mixed (2-type only) fiber types
    # TODO: fix warnings
    df_type = df.copy()
    df_type = df_type.reindex(df_type.index.values.tolist() + [label])

    isoforms = list({x for l in iso_pairs for x in l})

    # first check for hybrids that reach threshold
    # the iso_pairs list should be roughly in priority of order to check
    # as it will exit when a hit is found
    for s in df_type.columns:
        for j in iso_pairs:
            if (df_type[s].loc[j[0]] + df_type[s].loc[j[1]]) > 0.8:
                df_type[s].loc[label] = label[:-5] + str(j[0][3:]) + '/' + str(j[1][3:])
                break

                # if a single isoform is over threshold, set/replace it
    for s in df_type.columns:
        for i in isoforms:
            if df_type[s].loc[i] > 0.8:
                df_type[s].loc[label] = i
                break

    print(df_type.loc[label].value_counts())
    return df_type