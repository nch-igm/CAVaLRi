import os
import sys
import pandas as pd
import numpy
import json
import yaml
import argparse


def build_cohort_summary(cohort):

    # Get the paths to all of the case results
    rows = []
    for case in cohort.cases:

        # Intialize row to return
        newrow = {
            'subjectId': case.case_id,
            'geneSymbol': '',
            'hitCompositeLR': '',
            'hitCompositeLR_log10': '',
            'hitPostTest': '',
            'hitRank': '',
            'isHit':'' ,
            'isTp': case.case_data['isTp'],
            'maxCompositeLR_log10': '',
            'maxPostTest': ''
        }

        # If the case is a true positive, find the hit data
        if newrow['isTp'] == 1:
            if len(case.case_data['tpHits']) == 0:
                newrow['isHit'] = 0
            else:
                newrow['hitRank'] = min(case.case_data['tpHits'])
                newrow['isHit'] = 1

            for d in case.case_data['diseases']:
                if d['geneRank'] == newrow['hitRank']:
                    newrow['geneSymbol'] = d['gene_data']['gene']
                    newrow['hitCompositeLR_log10'] = d['compositeLR']
                    newrow['hitCompositeLR'] = 10**d['compositeLR']
                    newrow['hitPostTest'] = d['postTestProbability']
                    break
        
        else:
            newrow['isHit'] = 0

        # Get top disease
        for d in case.case_data['diseases']:
            if d['diseaseRank'] == 1:
                newrow['maxCompositeLR_log10'] = d['compositeLR']
                newrow['maxPostTest'] = d['postTestProbability']
                break

        # Add newrow to result data frame
        rows.append(pd.DataFrame(newrow, index = [0]))

    return pd.concat(rows)
