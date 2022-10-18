# Import necessary packages
import numpy as np
import pandas as pd
import scipy
import math
import os
import sys
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import StratifiedShuffleSplit

# Import local packages
sys.path.append('..')
sys.path.append('../../workflow')
from utils import get_diagnostic_df, get_extended_diagnostic_df
from config import *

# Start a seed
seed = 44

def get_cases(path):
    """
    Iterate through the provided path to collect all subject IDs, return
    data frame containing true positive infomration
    """
    
    # Read in all cases in selected output directory
    # path = '../output/2021-10-06_clinphen/summary'
    cases = []
    for entry in os.scandir(path):
        if entry.name.endswith('.lr.tsv'):
            cases.append(entry.name[:entry.name.find('.lr.tsv')])

    # Intialize case data frame
    case_df = pd.DataFrame(columns = ['copath'])
    for i, case in enumerate(cases):
        case_df.loc[i, 'copath'] = case

    # Merge in true positive values
    diag_df = get_diagnostic_df()
    case_df = case_df.merge(diag_df, how = 'left', left_on = 'copath', right_on = 'subjectId')
    
    # Calculate is true positive row
    def istp(row):
        return True if not pd.isna(row['gene']) else False

    case_df['isTp'] = case_df.apply(istp, axis = 1)
    demo_df = pd.read_excel(os.path.join(config['project_root'], config['demographic_data']))
    diag_df = pd.read_excel(os.path.join(config['project_root'], config['diagnostic_source']))
    demo_df = demo_df.merge(diag_df, right_on = 'MRN', left_on = 'MRN')[['CoPath.Number', 'MRN', 'RACE', 'ETHNICITY']].drop_duplicates().reset_index(drop = True)
    for col in ['RACE', 'ETHNICITY']:
        for nv in ['No Information', '<null>', 'Refuse to answer', pd.NA, np.NaN, 'nan']:
            demo_df.loc[demo_df[col] == nv, col] = 'Unknown'
    case_df = case_df.merge(demo_df, right_on = 'CoPath.Number', left_on = 'copath', how = 'left')
    for col in ['RACE', 'ETHNICITY']:
        case_df.loc[case_df[col].isna(), col] = 'Unknown'
    case_df = case_df.drop(columns = ['MRN', 'subjectId', 'CoPath.Number'])

    return case_df


def split_train_validation(case_df):
    """
    Set aside 30% of the data as a validation set
    """
    def get_cat(row):
        # return f"{row['isTp']}{row['RACE']}{row['ETHNICITY']}"
        return f"{row['isTp']}{row['RACE']}"

    # Parition training and validation data
    sss = StratifiedShuffleSplit(n_splits=1, test_size=0.3, random_state=seed)
    X = case_df[['copath', 'gene']]
    case_df['cat'] = case_df.apply(get_cat, axis = 1)
    Y = case_df[['cat']]
    for train_idx, val_idx in sss.split(X, Y):
        X_train, X_val = X.loc[train_idx], X.loc[val_idx]
        Y_train, Y_val = Y.loc[train_idx], Y.loc[val_idx]
    
    return X_train, X_val, Y_train, Y_val


def get_cohort_data(df: pd.DataFrame, res_dir: str):
    """
    Gather variant data from output/summary files ending in "lr.tsv" and combine
    into a single data frame
    INPUT:
        df -- One column data frame, header = ['copath'], indicating copath cases
        res_dir -- where to find case summary files ending in "lr.tsv"
    """

    def add_tp_column(row):
        # Return binary indicating true positive status of variant
        return 0 if pd.isna(row['gene']) else 1

    # Intialize list to collect case summary data frames
    case_dfs = []

    # Read in each case
    for entry in os.scandir(res_dir):
        case = entry.name[:entry.name.find('.lr.tsv')]
        if entry.name.endswith('lr.tsv') and case in df['copath'].to_list():
            summary_df = pd.read_csv(os.path.join(res_dir, entry.name), sep = '\t')
            summary_df = summary_df.merge(df.loc[df['copath'] == case], how = 'left', left_on = 'geneSymbol', right_on = 'gene')
            summary_df['isTp'] = summary_df.apply(add_tp_column, axis = 1)
            summary_df['subject'] = case
            case_dfs.append(summary_df)

    return pd.concat(case_dfs)


def scale_df(base_df: pd.DataFrame, pos: tuple):
    """"
    Return a scaled version of the base data frame (all scalars = 1)
    INPUT:
        base_df -- Base data frame indicating base LRs without scale
        pos -- tuple indicating scalars (geneLR, moiLR) 
    """
    scaled_base_df = base_df.copy()
    scaled_base_df['geneLR'] = scaled_base_df['geneLR'] * pos[0]
    scaled_base_df['moiLR'] = scaled_base_df['moiLR'] * pos[1]
    scaled_base_df['compositeLR'] = scaled_base_df['phenoLR'] + scaled_base_df['geneLR'] + scaled_base_df['moiLR']
    scaled_base_df = scaled_base_df.sort_values('compositeLR', ascending = False)
    return scaled_base_df


def get_roc_auc(current_df: pd.DataFrame):

    # Calculate the ROC area under the curve (AUC)
    df = current_df.copy()
    df = df.sort_values('compositeLR', ascending=False)
    df = df.reset_index(drop=True).reset_index().rename(columns = {'index': 'rank'})
    roc_auc = roc_auc_score(df['isTp'], df['compositeLR'])
    return roc_auc


def get_rank_improvement(current_df: pd.DataFrame, original_rankings: pd.DataFrame):
    """
    Calculate the p-value of a one-sided Wilcoxon signed rank test between a 
    rankings derived from unscaled LRs (original_rankings) and rankings derived
    from scaled LRs (current_df). Both input data frames should have both a 
    'subject' and 'compositeLR' column.
    """

    # Get ranks of true positive variants in the provided data frame
    tp_ranks = pd.DataFrame(index = current_df['subject'].unique())
    tp_ranks['rank'] = -1
    for idx, row in tp_ranks.iterrows():
        sub_df = current_df.loc[current_df['subject'] == idx]
        sub_df = sub_df.sort_values('compositeLR', ascending = False)
        sub_df = sub_df.reset_index(drop=True).reset_index().rename(columns = {'index': 'rank'})
        sub_df['rank'] += 1
        tp_rank = sub_df.loc[sub_df['isTp'] == 1].reset_index(drop=True)
        if len(tp_rank.index) != 0:
            tp_ranks.loc[idx, 'rank'] = tp_rank.loc[0, 'rank']
        else:
            tp_ranks.loc[idx, 'rank'] = pd.NA
    
    # Join new ranks with original ranks
    comp_df = original_rankings.merge(tp_ranks, how = 'left', left_index = True, right_index = True, suffixes = ("_original", "_new"))
    comp_df = comp_df[~comp_df['rank_new'].isna()]#.reset_index(drop=True)
    comp_df['rank_difference'] = comp_df['rank_original'] - comp_df['rank_new']

    # Get count of number one rankings
    number_ones = len(comp_df[comp_df['rank_new'] == 1].index)

    # Get average rank
    # avg_rank = comp_df['rank_new'].mean()

    # Return the p-value and comparison data frame
    return comp_df['rank_difference'].sum(), comp_df
    # return number_ones, comp_df



def get_wilcoxon(current_df: pd.DataFrame, original_rankings: pd.DataFrame):
    """
    Calculate the p-value of a one-sided Wilcoxon signed rank test between a 
    rankings derived from unscaled LRs (original_rankings) and rankings derived
    from scaled LRs (current_df). Both input data frames should have both a 
    'subject' and 'compositeLR' column.
    """

    # Get ranks of true positive variants in the provided data frame
    tp_ranks = pd.DataFrame(index = current_df['subject'].unique())
    tp_ranks['rank'] = -1
    for idx, row in tp_ranks.iterrows():
        sub_df = current_df.loc[current_df['subject'] == idx]
        sub_df = sub_df.sort_values('compositeLR', ascending = False)
        sub_df = sub_df.reset_index(drop=True).reset_index().rename(columns = {'index': 'rank'})
        sub_df['rank'] += 1
        tp_rank = sub_df.loc[sub_df['isTp'] == 1].reset_index(drop=True)
        if len(tp_rank.index) != 0:
            tp_ranks.loc[idx, 'rank'] = tp_rank.loc[0, 'rank']
        else:
            tp_ranks.loc[idx, 'rank'] = pd.NA
    
    # Join new ranks with original ranks
    comp_df = original_rankings.merge(tp_ranks, how = 'left', left_index = True, right_index = True, suffixes = ("_original", "_new"))
    comp_df = comp_df[~comp_df['rank_new'].isna()]#.reset_index(drop=True)

    # Calculate the one-sided Wilcoxon statistic
    data = scipy.stats.wilcoxon(comp_df['rank_original'], comp_df['rank_new'], alternative='greater')
    
    # Return the p-value and comparison data frame
    return data.pvalue, comp_df


def get_successors(current_pos, array, step_length):
    """
    In the greedy search, find neighbors along each dimension and step in amount
    indicated by step length (step_length). Start the step from the provided
    intial position (current_pos). Space dimension defined by (array).
    """

    # Intialize list to capture valid successor positions
    successor_pos = []
    successors = []

    for i in range(0, len(array.shape)):
        dim_positions = []

        # Stepping forward
        if current_pos[i] + step_length >= 1 and current_pos[i] + step_length < array.shape[i] + 1:
            dim_positions.append(current_pos[i] + step_length)

        # Stepping backward
        if current_pos[i] - step_length >= 1 and current_pos[i] - step_length < array.shape[i] + 1:
            dim_positions.append(current_pos[i] - step_length)
        
        # Add dimension positions
        dim_positions.sort()
        successor_pos.append(dim_positions)

    for i, dim in enumerate(successor_pos):
        for pos in dim:
            new_pos = list(current_pos)
            new_pos[i] = pos
            successors.append(tuple(new_pos))

    return successors


def build_integer_map(base_df, baseline_tp_ranks):

    path = '/Users/rsrxs003/Downloads/optimal_rank11.csv'
    if os.path.exists(path):
        return pd.read_csv(path)
    print('Building Grid')
    # Run the pipeline for each value
    # optimize_df = pd.DataFrame(columns = ['geneLR_scalar', 'moiLR_scalar', 'roc_auc', 'wilcoxon_pvalue'])
    optimize_df = pd.DataFrame(columns = ['geneLR_scalar', 'moiLR_scalar', 'roc_auc', 'rank_improvement'])

    for i in range(1, 11):
        for j in range(1, 11):
            rescaled_df = scale_df(base_df, (i, j))
            result = (get_roc_auc(rescaled_df), get_rank_improvement(rescaled_df, baseline_tp_ranks))
            optimize_df = optimize_df.append({
                'geneLR_scalar': i,
                'moiLR_scalar': j,
                'roc_auc': result[0],
                # 'wilcoxon_pvalue': result[1][0]
                'rank_improvement': result[1][0]
            }, ignore_index=True)
    
    # Add negative log term
    # optimize_df['wilcoxon_pvalue_negative_log'] = -1 * optimize_df['wilcoxon_pvalue'].apply(lambda x: math.log(x, 10))
    
    return optimize_df


def greedy_search(data: pd.DataFrame, res_dir: str, a: float):
    """
    Greedy search through two dimensional space in search of a value that
    satisfies both:
        1) Area under the receiver operating characteristic curve (AUC ROC)
        2) P-value of one sided wilcoxon signed rank test
    
    Define alpha (a) that functions to weight either AUC ROC (a = 1), the
    wilcoxon (a = 0), or somewhere in between (0 < a < 1).
    """

    # Read in data
    base_df = get_cohort_data(data, res_dir)
    
    # Adjust LR scaled values to be unscaled, creating a baseline data frame
    # base_df['geneLR'] = base_df['geneLR']/config['geneLR_scalar']
    # base_df['moiLR'] = base_df['moiLR']/config['moiLR_scalar']
    base_df['compositeLR'] = base_df['phenoLR'] + base_df['geneLR']

    # Obtain base rankings for wilcoxon comparison calculations
    baseline_tp_ranks = pd.DataFrame(index = base_df['subject'].unique())
    baseline_tp_ranks['rank'] = -1
    for idx, row in baseline_tp_ranks.iterrows():
        sub_df = base_df.loc[base_df['subject'] == idx]
        sub_df = sub_df.sort_values('compositeLR', ascending = False)
        sub_df = sub_df.reset_index(drop=True).reset_index().rename(columns = {'index': 'rank'})
        sub_df['rank'] += 1
        tp_rank = sub_df.loc[sub_df['isTp'] == True].reset_index(drop=True)
        if len(tp_rank.index) != 0:
            baseline_tp_ranks.loc[idx, 'rank'] = tp_rank.loc[0, 'rank']
        else:
            baseline_tp_ranks.loc[idx, 'rank'] = pd.NA

    # Define space
    dim_length, dim_count = 10, 2
    space = np.empty([dim_length] * dim_count)

    # Get integer values of space to determine AUC ROC and Wilcoxon ranges
    # print(f'Building Integer Map (a={a})')
    range_df = build_integer_map(base_df, baseline_tp_ranks)
    range_df.to_csv(f'~/Downloads/avg_rank_map{a}.csv')
    minimax = {
        'roc_min': range_df['roc_auc'].min(),
        'roc_max': range_df['roc_auc'].max(),
        # 'wilcoxon_min': range_df['wilcoxon_pvalue_negative_log'].min(),
        # 'wilcoxon_max': range_df['wilcoxon_pvalue_negative_log'].max()
        'rank_min': range_df['rank_improvement'].min(),
        'rank_max': range_df['rank_improvement'].max()
    }

    # Initalize variables
    pos = np.ones(dim_count)
    step_length = 1
    previous_result = 0.01
    max_steps = 15
    stop_percentage = 0.1

    # Intialize dictionary to store results
    results = {}

    # Keep going until we are within the bounds of the stop percentage
    for i in range(0, max_steps):

        # Get successors
        successors = get_successors(pos, space, step_length) + [pos]

        # Run the pipeline for each successor
        # successor_df = pd.DataFrame(columns = ['geneLR_scalar', 'moiLR_scalar', 'roc_auc', 'wilcoxon_pvalue'])
        successor_df = pd.DataFrame(columns = ['geneLR_scalar', 'moiLR_scalar', 'roc_auc', 'rank_improvement'])
        for successor in successors:
            rescaled_df = scale_df(base_df, successor)
            result = (get_roc_auc(rescaled_df), get_rank_improvement(rescaled_df, baseline_tp_ranks))

            # Add succesor to data frame
            successor_df = successor_df.append({
                'geneLR_scalar': successor[0],
                'moiLR_scalar': successor[1],
                'roc_auc': result[0],
                # 'wilcoxon_pvalue': result[1][0]
                'rank_improvement': result[1][0]
            }, ignore_index=True)

        """
        Normalize successors
        """
        def minimax_normalization(row, col, minimax):
            key = col[:col.find('_')] # Get minimax key
            return (row[col] - minimax[f'{key}_min']) / (minimax[f'{key}_max'] - minimax[f'{key}_min'])

        # Normal AUC ROC
        successor_df['roc_auc_norm'] = successor_df.apply(
            minimax_normalization,
            axis = 1,
            col = 'roc_auc',
            minimax = minimax
            )

        # Log transform wilcoxon p-value and normalize
        # successor_df['wilcoxon_pvalue_negative_log'] = -1 * successor_df['wilcoxon_pvalue'].apply(
        #                                         lambda x: math.log(x, 10))
        # successor_df['wilcoxon_norm'] = successor_df.apply(
        successor_df['rank_improvement_norm'] = successor_df.apply(
            minimax_normalization,
            axis = 1,
            # col = 'wilcoxon_pvalue_negative_log',
            col = 'rank_improvement',
            minimax = minimax
            )

        # Add a column that considers both AUC ROC and Wilcoxon P-value
        # successor_df['weighted_norm'] = a * successor_df['roc_auc_norm'] + (1-a) * successor_df['wilcoxon_norm']
        successor_df['weighted_norm'] = a * successor_df['roc_auc_norm'] + (1-a) * successor_df['rank_improvement_norm']

        # Get the best successor
        best_weighted_df = successor_df.loc[successor_df['weighted_norm'] == successor_df['weighted_norm'].max()].reset_index(drop=True)
        best_successor = {
            'geneLR_scalar': float(best_weighted_df.loc[0, 'geneLR_scalar']),
            'moiLR_scalar': float(best_weighted_df.loc[0, 'moiLR_scalar']),
            'roc_auc_norm': best_weighted_df.loc[0, 'roc_auc_norm'],
            # 'wilcoxon_pvalue_norm': best_weighted_df.loc[0, 'wilcoxon_norm'],
            'rank_improvement_norm': best_weighted_df.loc[0, 'rank_improvement_norm'],
            'weighted_norm': best_weighted_df.loc[0, 'weighted_norm'],
            }

        new_pos = (best_successor['geneLR_scalar'], best_successor['moiLR_scalar'])
        new_result = best_successor['weighted_norm']
    
        # Check to see if the best successor has already been assessed
        # print(f'{round((new_result / previous_result - 1) * 100, 2)}% increase in accuracy\n')
        if new_pos in results.keys():
            if new_result > previous_result:
                pos = new_pos
                previous_result = new_result
            step_length = step_length/2        # Half the step distance

        # Break if the result is within given percentage
        elif (new_result / previous_result - 1) * 100 < stop_percentage:
            if new_result > previous_result:
                pos = new_pos
                previous_result = new_result
                results.update({pos: previous_result}) # Lay a breadcrumb
            break

        else:
            if new_result > previous_result:
                pos = new_pos
                previous_result = new_result
                results.update({pos: previous_result}) # Lay a breadcrumb
            
    
    return pos, results


def optimize_scalars(path: str):
    """
    Take provided path, extract cases, optimize three points in scalar space
    corresponding to 1) (a = 0) minimum average rank, 2) (a = 1) maximize ROC
    AUC, 3) (a = 0.5) an evenly weighted average between the two.
    """

    # Get cases
    case_df = get_cases(path)

    # Split data into a training and validation set
    X_train, X_val, Y_train, Y_val = split_train_validation(case_df)
    # return Y_train, case_df

    # Partition the training data (5 iterrations)
    sss = StratifiedShuffleSplit(n_splits=5, test_size=0.2, random_state=seed)
    optimal_positions = []
    # for a in [0, 0.5, 1]:
    for a in [0.5]:
        print(f'ALPHA: {a}')
        current_fold = 1
        positions = []
        for train_idx, test_idx in sss.split(X_train, Y_train):
            df = X_train.iloc[train_idx].reset_index(drop=True)
            pos, res  = greedy_search(df, path, a)
            print(f'Pipeline complete for fold #{current_fold}')
            current_fold += 1
            positions.append(pos)
            print(pos, res)
            print()

        optimal_positions.append((pd.Series([x[0] for x in positions]).mean(), pd.Series([x[1] for x in positions]).mean()))
    return optimal_positions