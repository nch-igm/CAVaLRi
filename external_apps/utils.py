import os
import vcf
import pandas as pd
import copy

def system_exec(cmd: str, show_cmd: bool = False, warn: bool = False):
    if show_cmd is True:
        print(cmd)
    status = os.system(cmd)
    if status != 0:
        if warn is True:
            print(f"WARNING: Command '{cmd}' did not complete successfully: ({status}).")
        else:
            print(f"ERROR: Command '{cmd}' did not complete successfully ({status}).")
            raise RuntimeError

    return None


def excel_to_vcf(excel_path, vcf_template, output_folder): 
    """
    Use a template vcf to convert a excel to a vcf file
    """
    vcf_reader = vcf.Reader(open('/Users/rsrxs003/projects/LIRICAL/input/EXOME_RERUNS/M20-5833.vcf', 'r'))
    for record in vcf_reader:
        rec = record
        break

    vcf_reader = vcf.Reader(open(vcf_template, 'r'))

    df = pd.read_excel(excel_path, sheet_name='Filtered')

    # Same for all rows
    rec.ID = None
    rec.QUAL = None
    
    # Write to a new vcf
    vcf_writer = vcf.Writer(open(f'{output_folder}/spliceai.vcf', 'w'), vcf_reader)
    for idx, row in df.iterrows():

        # Make a copy of the rec and reassign attributes
        print(row)
        rec_current = copy.copy(rec)
        print(rec)
        rec_current.CHROM = f'chr{row["Chr"]}'
        rec_current.POS = row['Start']
        rec_current.REF = row['Ref']
        rec_current.ALT = row['Alt']
        print(rec_current)

        # Send to a line in the vcf
        vcf_writer.write_record(rec_current)

    return f'{output_folder}/spliceai.vcf'