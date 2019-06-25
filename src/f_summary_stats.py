# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, copy, pickle
import numpy as np
from collections import defaultdict
sys.path.append('/home/unix/maxwshen/')
from mylib import util, compbio
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'e_newgenotype_Cas9/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
lib_design = None


##
# f
##
def is_precise_1nt_ins(ins_nt, target_row):
  seq = target_row['Sequence']
  cutsite_idx = target_row['Cutsite index']
  return bool(ins_nt == seq[cutsite_idx - 1])

def get_stats(edit_dfs, target_row):
  simple_count_cats = [
    'wildtype',
    'hdr',
  ]
  all_cats = simple_count_cats + [
    'del_MH',
    'del_MHless',
    'ins_1nt_precise',
    'ins_1nt_imprecise',
    'ins_2ntplus',
  ]
  dd = {cat: 0 for cat in all_cats}
  for idx, row in edit_dfs.iterrows():
    if row['Category'] in simple_count_cats:
      dd[row['Category']] += row['Count']
    if row['Category'] == 'del':
      if row['Microhomology-Based'] == 'yes':
        dd['del_MH'] += row['Count']
      else:
        dd['del_MHless'] += row['Count']
    if row['Category'] == 'ins':
      if int(row['Length']) == 1:
        if is_precise_1nt_ins(row['Inserted Bases'], target_row):
          dd['ins_1nt_precise'] += row['Count']
        else:
          dd['ins_1nt_imprecise'] += row['Count']
      else:
        dd['ins_2ntplus'] += row['Count']
  return dd

##
# Primary
##

def summary_stats(exp_design_subset, lib_design, lib_nm):
  '''
    Edited outcomes only
  '''

  dd = defaultdict(list)
  timer = util.Timer(total = len(exp_design_subset))
  for idx, row in exp_design_subset.iterrows():
    condition = row['Name']
    target_nm = row['Target']
    target_row = target_design[target_design['Target'] == target_nm].iloc[0]

    df = pd.read_csv(inp_dir + f'{condition}_genotypes_0.csv', index_col = 0)
    peptide_nms = set(df['Peptide name'])

    for peptide_nm in peptide_nms:
      '''
        dataframe subsetting is slow. For larger libraries, parallelize or consider storing e_ output as dict of dfs.
      '''
      edit_dfs = df[df['Peptide name'] == peptide_nm]

      stats = get_stats(edit_dfs, target_row)
      for stat in stats:
        dd[stat].append(stats[stat])
      dd['Peptide name'].append(peptide_nm)
      dd['Condition'].append(condition)

    timer.update()

  df = pd.DataFrame(dd)
  df.to_csv(out_dir + f'{lib_nm}.csv')
  return


@util.time_dec
def main():

  libs = set(exp_design['Library'])

  for lib_nm in libs:
    print(lib_nm)
    exp_design_subset = exp_design[exp_design['Library'] == lib_nm]
    lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv')

    summary_stats(exp_design_subset, lib_design, lib_nm)

  return


if __name__ == '__main__':
  main()