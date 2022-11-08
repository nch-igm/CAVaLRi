import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import roc_curve, roc_auc_score


def get_roc_data(data):

    # Generate ROC curve metrics
    true_positives = data['result'].to_list()
    scores = data['postTestProbability'].to_list()
    fpr, tpr, thresholds = roc_curve(true_positives, scores, pos_label=1)
    roc_data = pd.DataFrame(list(zip(fpr, tpr, thresholds)), columns=['fpr', 'tpr', 'thresholds'])
    roc_auc = roc_auc_score(true_positives, scores)

    return roc_data, roc_auc


def build_roc(cohort, comp):

    config = cohort.config

    # Create figure
    fig, ax = plt.subplots(figsize = (12,10), dpi = 300)

    # Big Title        
    fig.suptitle("Pathogenic Variant Classfier ROC", fontsize="30", y=.93)

    # Plot midline
    midline = pd.DataFrame([[0,0], [1,1]])
    midline.columns = ['True Positive Rate', 'False Positive Rate']
    sns.lineplot(ax=ax, data=midline, x='False Positive Rate', y='True Positive Rate', color = 'grey', linewidth = 3.0, linestyle='--')

    # Get ROC input data
    roc_data, roc_auc = get_roc_data(cohort.cohort_variants)

    # Plot ROC Data
    roc = ax.plot(roc_data['fpr'].to_list(), roc_data['tpr'].to_list(), 
                    linewidth = 3.0, markersize=5, marker='d', color = 'blue', 
                    label='{} (AUC = {})'.format(config['plot_label'], round(roc_auc, 4)))
    
    # Get comparator data
    if comp:
        comparator_variants = pd.read_csv(os.path.join(config['comparator_results'], 'summary', 'cohort_variants.tsv'), sep='\t')
        comp_roc_data, comp_roc_auc = get_roc_data(comparator_variants)
        roc_c = ax.plot(comp_roc_data['fpr'].to_list(), comp_roc_data['tpr'].to_list(),
                        linewidth = 3.0, markersize=5, marker='d', color = 'black', 
                        label='{} (AUC = {})'.format(config['comparator_plot_label'], round(comp_roc_auc, 4)))

    # Show legend
    ax.legend(loc="lower right", fontsize="20")

    # Adjust font sizes
    ax.set_xlabel('False Positive Rate',fontsize = 25)
    ax.set_ylabel('True Positive Rate',fontsize = 25)
    ax.tick_params(labelsize = 20)
    

    return fig, roc_data


def build_topn(ax, **kwargs):

    # Set line width
    kwargs['linewidth'] = 5.0

    # for d in [data, filtered_data]:
    d = kwargs['data'].copy()

    # Rename rank column
    if 'geneRank' in d.columns.to_list():
        d = d.rename(columns={'geneRank': 'rank'})
    elif 'hitRank' in d.columns.to_list():
        d = d.rename(columns={'hitRank': 'rank'})
    else:
        d = d.rename(columns={'diseaseRank': 'rank'})

    # Initialize plot data frame
    p_data = pd.DataFrame(columns = ['top_n', 'percent_of_subjects'])

    # Iterate from 1-10 and populate plotting data frame
    d = d.loc[d['isTp'] == 1]
    worst_rank = d[(d['rank'] != '') & (d['rank'].notna())]['rank'].max()
    d.loc[(d['rank'] == '') | (d['rank'].isna()),'rank'] = worst_rank
    total_subjects = len(d.index)
    for i in range(1,51):
        n_rank_subjects = len(d.loc[d['rank'] <= i].index)
        p_data = pd.concat([
            p_data,
            pd.DataFrame({
                'top_n': i,
                'percent_of_subjects': n_rank_subjects/total_subjects*100
                }, index = [0])
        ])
    p_data.loc[p_data['top_n'] == 50, 'percent_of_subjects'] = 100
    kwargs['data'] = p_data
    kwargs['label'] = f"{kwargs['label']} (Average Rank: {round(d['rank'].mean(), 2)})"
    # Plot data
    sns.lineplot(ax=ax, x='top_n', y='percent_of_subjects', **kwargs)

    # Format x labels
    ax.set_xlim(0, 51)
    ax.set_ylim(-2, 102)
    ax.set_xticks([1, 10, 20, 30, 40, 50])
    ax.set_xticklabels(['1', '10', '20', '30', '40', '>50'])

    # Set labels
    ax.set_xlabel('Diagnostic Variant Rank', size = 45)
    ax.set_ylabel('Percent of Subjects', size = 45, x = 0.5)

    return p_data
    

def construct_topn(cohort, comp):

    config = cohort.config

    # Create figure
    topn_fig, topn_ax = plt.subplots(figsize = (12,10), dpi = 300)

    # Big Title        
    topn_fig.suptitle("Diagnostic Variants Found in Top N", fontsize="30", y=.93)

    # Plot the TOP-N curves
    kwargs = {
        'data': cohort.cohort_summary,
        'color': 'blue',
        'label': config['plot_label']
    }
    topn_data = build_topn(topn_ax, **kwargs)

    if comp:

        # Read in comp data
        comparator_data = pd.read_csv(os.path.join(cohort.root_path, config['comparator_results'], 'summary', 'cohort_summary.tsv'), sep = '\t')
    
        kwargs = {
            'data': comparator_data,
            'color': 'black',
            'label': config['comparator_plot_label']
        }
        topn_comp_data = build_topn(topn_ax, **kwargs)

        topn_data = topn_data.merge(topn_comp_data, on = 'top_n')

    # Style the plot
    # style(topn_fig, topn_ax)

    # Add legend
    legend = topn_ax.legend(loc="lower right", fontsize = 20)
    
    # Adjust font sizes
    topn_ax.set_xlabel('Diagnostic Variant Rank',fontsize = 25)
    topn_ax.set_ylabel('Percent of Subjects',fontsize = 25)
    topn_ax.tick_params(labelsize = 20)

    # Save the figure and topn data points
    return topn_fig, topn_data
    


def plot_figures(cohort):

    config = cohort.config

    try:
        x = config['comparator_results']
        comp = True
    except:
        comp = False

    ## BUILD ROC ##
    roc_fig, roc_data = build_roc(cohort, comp)

    ## BUILD TOP-N ##
    topn_fig, topn_data = construct_topn(cohort, comp)

    return roc_fig, roc_data, topn_fig, topn_data

