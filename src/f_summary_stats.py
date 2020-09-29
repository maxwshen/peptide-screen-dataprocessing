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
      # Ignore long deletions
      if row['Length'] > 60:
        continue
      if row['Microhomology-Based'] == 'yes':
        dd['del_MH'] += row['Count']
      else:
        dd['del_MHless'] += row['Count']
    if row['Category'] == 'ins':
      # Ignore long insertions
      if row['Length'] > 20:
        continue
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

def summary_stats(condition, row):
  '''
    Edited outcomes only
  '''

  dd = defaultdict(list)
  condition = row['Name']
  target_nm = row['Target']
  target_row = target_design[target_design['Target'] == target_nm].iloc[0]

  n_members = len(lib_design)
  jump = 150
  timer = util.Timer(total = int(n_members / jump))
  for jdx in range(0, n_members, jump):
  # for jdx in range(0, 500, 500):
    start_jdx = jdx
    end_jdx = jdx + jump

    inp_pkl_fn = inp_dir + f'{condition}_genotypes_{start_jdx}.pkl'
    with open(inp_pkl_fn, 'rb') as f:
      ds = pickle.load(f)

    peptide_nms = set(ds.keys())
    for peptide_nm in peptide_nms:
      edit_dfs = ds[peptide_nm]

      stats = get_stats(edit_dfs, target_row)
      for stat in stats:
        dd[stat].append(stats[stat])
      dd['Peptide name'].append(peptide_nm)
      dd['Condition'].append(condition)

    timer.update()

  df = pd.DataFrame(dd)

  edited_cols = [col for col in df.columns if col not in ['wildtype', 'Peptide name', 'Condition']]
  df['Edited total'] = df.apply(
      lambda row: sum(row[edited_cols]), 
      axis = 'columns'
  )
  df['Edited fraction'] = df['Edited total'] / (df['wildtype'] + df['Edited total'])

  df['Total reads'] = df.apply(
    lambda row: sum(row[edited_cols + ['wildtype']]),
    axis = 'columns'
  )

  new_lib_design = lib_design.rename(columns = {'Name': 'Peptide name'})
  mdf = new_lib_design.merge(df, on = 'Peptide name', how = 'outer')

  mdf = mdf.drop_duplicates()
  mdf.to_csv(out_dir + f'{condition}.csv')
  return

##
# qsub
##
def gen_qsubs():
  # Generate qsub shell scripts and commands for easy parallelization
  print('Generating qsub scripts...')
  qsubs_dir = _config.QSUBS_DIR + NAME + '/'
  util.ensure_dir_exists(qsubs_dir)
  qsub_commands = []

  num_scripts = 0
  for condition in exp_design['Name']:
    exp_row = exp_design[exp_design['Name'] == condition].iloc[0]
    lib_nm = exp_row['Library']

    command = f'python {NAME}.py {condition}'
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + f'q_{script_id}_{condition}.sh'
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -j y -P regevlab -V -l h_rt=4:00:00,h_vmem=2G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))

  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return

@util.time_dec
def main(argv):

  cond = argv[0]

  exp_row = exp_design[exp_design['Name'] == cond].iloc[0]

  global lib_design
  lib_nm = exp_row['Dual library']
  lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv', index_col = 0)


  summary_stats(cond, exp_row)

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()