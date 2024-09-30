import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.cm as cmx
import numpy as np

def plot_protein_counts_by_fiber(df, cutoff=0, cmap=None, alphas=0.5, label=None):
    '''plots protein counts above a given threshold  (cutoff) and prints summary statistics'''
    cutoff_counts = df[df > cutoff].count()

    plt.rcParams['figure.figsize'] = 5, 5
    plt.hist(cutoff_counts, alpha=0.5, label=label, color=[cmap[label] if cmap is not None else 'b'])

    plt.legend(fontsize=30)
    # TODO: determine quant method autmatically
    plt.title('Protein counts per fiber by iBAQ', fontsize=22)
    plt.grid(False)
    # TODO: determine quant method automatically
    plt.xlabel('# Proteins/Fiber', fontsize=20)
    plt.ylabel('# Fibers', fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)

    # numb samples (fibers) with less than 400 proteins by iBAQ
    print('# fibers w/ < 400 proteins by iBAQ: ', sum(cutoff_counts < 400), 'or ',
          np.round(100 * sum(cutoff_counts < 400) / len(cutoff_counts)), '%')

    # num samples (fibers) with 0 proteins by iBAQ
    print('# fibers w/ 0 proteins by iBAQ: ', sum(cutoff_counts == 0), 'or ',
          np.round(100 * sum(cutoff_counts == 0) / len(cutoff_counts)), '%')

    # num samples (fibers)total iBAQ
    print('# fibers total: ', len(cutoff_counts))

    print('mean # proteins/fiber by iBAQ: ', np.round(cutoff_counts.mean()), '+/-', np.round(cutoff_counts.std()))


def plot_pct_isoforms(df, sorter, isoforms=[], cmap=[], threshold=0.8):
    m = df.sort_values(by=sorter, axis=1)
    plt.rcParams['figure.figsize'] = 10, 5
    fig = plt.figure()
    # TODO: automatically define other attributes
    ax = plt.subplot(111)

    for isoform in isoforms:
        plt.plot(m.loc[isoform], label=isoform, linewidth=2, color=cmap[isoform])

    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])
    plt.axhline(y=threshold, color='black', linestyle='dashed')
    ax.legend(loc='center left', bbox_to_anchor=(0.75, 0.5), fontsize=20)
    plt.yticks(fontsize=20)
    plt.ylabel('Isoform Fraction', fontsize=22)
    plt.xlabel('Fibers', fontsize=22)
    ax.tick_params(labelbottom=False)
    ax.tick_params(axis='x', width=0, length=0)
    # ax.set_ylim(0,1)
    ax.set_facecolor('white')
    ax.grid(False)


def plot_isoform_correlations(df, protein_combos=[]):
    # TODO: automate.

    for combo in protein_combos:
        plt.scatter(df.loc[combo[0]], df.loc[combo[1]])
        plt.plot([0, 1], linestyle='dashed', color='r')
        plt.xlim([0, 1])
        plt.ylim([0, 1])
        plt.xlabel(combo[0])
        plt.ylabel(combo[1])
        prsn = np.round(df.loc[combo[0]].corr(df.loc[combo[1]]), 3)
        rho = np.round(df.loc[combo[1]].corr(df.loc[combo[0]], method='spearman'), 3)
        plt.text(0.1, 0.9, 'Pearson=' + str(prsn))
        plt.text(0.1, 0.8, 'Spearman=' + str(rho))
        plt.show()


def boxplots_by_age(data, protein):
    l_young = data.to_df()[(data.obs['group'] == 'young')][protein]
    l_old = data.to_df()[(data.obs['group'] == 'old')][protein]
    l_all = data.to_df()[protein]
    tmphitobs = cdata.obs.copy()

    merged_young = pd.merge(l_young.to_frame(), data.obs[['leiden']], left_index=True, right_index=True, how='inner')
    mm = merged_young.melt(id_vars=['leiden'])
    mm = mm[mm.leiden != 'none']
    mm['value'] = mm['value'].astype(float)
    sns.boxplot(data=mm, x='leiden', y='value')
    plt.title('young group')
    plt.show()

    merged_old = pd.merge(l_old.to_frame(), data.obs[['leiden']], left_index=True, right_index=True, how='inner')
    mm = merged_old.melt(id_vars=['leiden'])
    mm = mm[mm.leiden != 'none']
    mm['value'] = mm['value'].astype(float)
    sns.boxplot(data=mm, x='leiden', y='value')
    plt.title('old group')
    plt.show()

    merged_all = pd.merge(l_all.to_frame(), data.obs[['leiden']], left_index=True, right_index=True, how='inner')
    mm = merged_all.melt(id_vars=['leiden'])
    mm = mm[mm.leiden != 'none']
    mm['value'] = mm['value'].astype(float)
    sns.boxplot(data=mm, x='leiden', y='value')
    plt.title('all age groups')
    plt.show()


def boxplots_by_leiden(data, protein):
    l0 = data.to_df()[(data.obs['leiden'] == '0')][protein]
    l1 = data.to_df()[(data.obs['leiden'] == '1')][protein]
    l2 = data.to_df()[(data.obs['leiden'] == '2')][protein]
    l_all = data.to_df()[protein]

    tmphitobs = cdata.obs.copy()

    merged_0 = pd.merge(l0.to_frame(), data.obs[['group', 'group_rep']], left_index=True, right_index=True, how='inner')
    mm = merged_0.melt(id_vars=['group', 'group_rep'])
    mm = mm[mm.group != 'none']
    mm['value'] = mm['value'].astype(float)
    sns.boxplot(data=mm, x='group', y='value')
    plt.title('slow leiden')
    plt.show()

    merged_1 = pd.merge(l1.to_frame(), data.obs[['group', 'group_rep']], left_index=True, right_index=True, how='inner')
    mm = merged_1.melt(id_vars=['group', 'group_rep'])
    mm = mm[mm.group != 'none']
    mm['value'] = mm['value'].astype(float)
    sns.boxplot(data=mm, x='group', y='value')
    plt.title('hybrid leiden')
    plt.show()

    merged_2 = pd.merge(l2.to_frame(), data.obs[['group', 'group_rep']], left_index=True, right_index=True, how='inner')
    mm = merged_2.melt(id_vars=['group', 'group_rep'])
    mm = mm[mm.group != 'none']
    mm['value'] = mm['value'].astype(float)
    sns.boxplot(data=mm, x='group', y='value')
    plt.title('fast leiden')
    plt.show()

    merged_all = pd.merge(l_all.to_frame(), data.obs[['group', 'group_rep']], left_index=True, right_index=True,
                          how='inner')
    mm = merged_all.melt(id_vars=['group', 'group_rep'])
    mm = mm[mm.group != 'none']
    mm['value'] = mm['value'].astype(float)
    sns.boxplot(data=mm, x='group', y='value')
    plt.title('all leiden')
    plt.show()