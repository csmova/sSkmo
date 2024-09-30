import Bio.Seq
import Bio.SeqRecord
from Bio import SeqIO
import re
from re import findall as refindall
import pickle
from collections import Counter
from utils.misc import intersection
def seqs_from_fasta(fasta_file):
    '''
    Set up a list of seqIO objects to save our sequences for the next step
    fasta file is from human protein database, contains seq ID and seq.

    Parameters:
        fasta_file (.fasta)

    Returns:
        sequences (list)
    '''

    sequences = []

    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(seq_record)

    print('Found ' + str(len(sequences)) + ' sequences.')

    return sequences


def peps_from_fasta(sequences, min_len, max_len, num_missed):
    '''
    make master list all_seq = all possible peptides from fasta file
    that are 7-40 AA long including all cleavages + 1 missed cleavage

    Parameters:
        sequences (): object returned from from seqs_from_fasta
        min_len ():
        max_len ():
        num_missed ():

    Returns:
        all_seq
    '''

    all_seq = []
    for i in range(0, len(sequences)):

        # this part makes cleavages
        tmp_lst = []
        tmp = re.split(r'(?<=[R|K])', str(sequences[i].seq))

        # this part makes the missed cleavages
        if num_missed > 0:
            for j in range(num_missed):
                gap = j + 2
                for k in range(len(tmp) - (num_missed + 1)):
                    tmp_lst.append(''.join(tmp[k:k + gap]))

        # combine all the seqences of appropiate length
        tmp_lst.extend(tmp)
        for seq in tmp_lst:
            if min_len - 1 < len(seq) < max_len + 1:
                all_seq.append(seq)

    print('Generated ' + str(len(all_seq)) + ' possible peptides.')

    return all_seq


def get_protein_peptide_dict(sequences, min_len, max_len, num_missed, nonuniqpep_lst, save_file=r'my_len_dict.pkl'):
    '''returns dict of each protein and corresponding possible unique peptides
    from min_len to max_len AAs in length, with num_missed

    Parameters:
        sequences:

    Returns:
        len_dict (dict):
    '''

    len_dict = {}

    for i in range(0, len(sequences)):
        if i % 1000 == 0: print(i)  # counter to monitor progress
        tmp_lst = []
        tmp = re.split(r'(?<=[R|K])', str(sequences[i].seq))  # make fragments
        if num_missed > 0:  # make skipped fragments
            for j in range(num_missed):
                gap = j + 2
                for k in range(len(tmp) - (num_missed + 1)):
                    tmp_lst.append(''.join(tmp[k:k + gap]))
        tmp_lst.extend(tmp)
        sel_seq = []
        for seq in tmp_lst:
            if min_len - 1 < len(seq) < max_len + 1:
                if seq not in nonuniqpep_lst:
                    sel_seq.append(seq)

        len_dict[sequences[i].id.split('|')[1]] = len(sel_seq)

    # save len_dict
    with open(save_file, 'wb') as handle:
        pickle.dump(len_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

    return len_dict


def remove_pep_by_seq(df, seq_to_remove='DEELDQLK'):
    '''
    returns df with seq_to_removed dropped

    '''
    idx = 0
    for i in df["Stripped.Sequence"]:
        if seq_to_remove == i:
            print('Removing ', i, ' at index ', idx)
            df.drop([idx], inplace=True)
            df.reset_index(drop=True, inplace=True)
        elif seq_to_remove in i:
            print('Warning: found similar match but did not remove:', i)
        idx += 1
    return df


def remove_non_unique_peps(all_seq, df):
    '''
    counts and removes non-unique peptides
    '''
    # make list of peptides with value count > 1 (ie non-unique)
    c = Counter(all_seq)
    c

    nonuniqpep_lst = [k for k, v in c.items() if v > 1]
    print('Found ' + str(len(nonuniqpep_lst)) + ' possible non-unique peptides.')

    # list of non-unique peptides found
    nonunique = intersection(nonuniqpep_lst, list(df['Stripped.Sequence'].unique()))
    print('Found ' + str(len(nonunique)) + ' non-unique peptides in data.')

    indeces_drop = []

    # make a list of indices of peptides containing non-unique sequences
    for x in nonunique:
        indeces_drop.append(df[df['Stripped.Sequence'] == x].index)

    indeces_drop

    ind_drop_list = []

    # convert indeces_drop into list
    for x in range(0, len(indeces_drop)):
        # print(indeces_drop[x][0])
        ind_drop_list.append(indeces_drop[x][0])

    ind_drop_list

    # drop the proteins from the drop list (i.e., the non-unique ones)
    # TODO: removes additional proteins upon re-running?
    df.drop(df.index[ind_drop_list], inplace=True)

    print("Removed " + str(len(ind_drop_list)) + " non-unique peptides from data.")

    return df, nonuniqpep_lst


def remove_ptm_seqs(df, ptm='UniMod'):
    '''removes sequences with given PTM string identifier
    e.g., UniMod'''

    df_nomod = df[df["Precursor.Id"].str.contains(ptm) == False]
    print('Removed ', str(df.shape[0] - df_nomod.shape[0]), ' PTM containing sequences.')
    return df_nomod


aminoacid = {
    'I': 'C6H13NO2',
    'L': 'C6H13NO2',
    'K': 'C6H14N2O2',
    'M': 'C5H11NO2S',
    'F': 'C9H11NO2',
    'T': 'C4H9NO3',
    'W': 'C11H12N2O2',
    'V': 'C5H11NO2',
    'R': 'C6H14N4O2',
    'H': 'C6H9N3O2',
    'A': 'C3H7NO2',
    'N': 'C4H8N2O3',
    'D': 'C4H7NO4',
    'C': 'C3H7NO2S',
    'E': 'C5H9NO4',
    'Q': 'C5H10N2O3',
    'G': 'C2H5NO2',
    'P': 'C5H9NO2',
    'S': 'C3H7NO3',
    'Y': 'C9H11NO3'
}

monoisotopic = {
    'S': 31.972,
    'C': 12.0000,
    'H': 1.0078,
    'O': 15.9949,
    'N': 14.0031
}

#borrowed from stack overflow
#todo: double check function

from re import findall as refindall

def molecular_weight(molecule):
    return sum(
        monoisotopic[atom] * int(num or '1')
        for atom, num in refindall(r'([A-Z][a-z]*)(\d*)', molecule)
    )

def protein_mass(protein):
    return sum(molecular_weight(aminoacid[char]) for char in protein)