# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess, math, pickle, imp
sys.path.append('/home/unix/maxwshen/')
import fnmatch
import numpy as np
from collections import defaultdict
from mylib import util
from mylib import compbio
import pandas as pd

# Default params
inp_dir = _config.OUT_PLACE + 'c6_polish/'
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
lib_design = None

crispr_cutsite = None

##
# Data manipulation
##
def parse_header(header):
  # count = int(header.replace('>', '').replace('_', ''))
  hw = header.replace('>', '').split('_')
  count = int(hw[0])
  category = hw[1]
  return count, category

def get_indel_length(seq):
  # Assume exactly one indel in sequence
  central_indel_len = trim_start_end_dashes(seq).count('-')
  if central_indel_len != 0:
    return central_indel_len 
  else:
    # 3' end gap deletion
    return seq[::-1].index(seq[::-1].replace('-', '')[:1])

def trim_start_end_dashes(seq):
  alphabet = set(list(seq))
  if '-' in alphabet:
    alphabet.remove('-')
  alphabet = list(alphabet)
  start = min([seq.index(s) for s in alphabet])
  end = min([seq[::-1].index(s) for s in alphabet])
  if end > 0:
    return seq[start:-end]
  else:
    return seq[start:]

def is_main_category(category):
  '''
    Match: delX, insX
    Do not match: delnotcrispr, insnotcrispr
    delnotatcut, insnotatcut
  '''
  if category[:3] in ['ins', 'del']:
    if 'notcrispr' not in category and 'notatcut' not in category:
      return True
  return False

##
# Pandas dataframe manipulation
##
def compress(df):
  # Compress a dataframe along all columns except 'Count', summing it.
  cols = list(df.columns)
  cols.remove('Count')
  if len(cols) == 0:
    new_df = pd.DataFrame({'Count': [df['Count'].sum()]})
  else:
    grouped = df.groupby(cols)
    g = grouped.aggregate(sum)
    new_df = g.reset_index()
  return new_df

##
# Alignment properties
##
def calc_deletion_start_position(read, genome, del_len):
  # assumes exactly one deletion in alignment
  idx = genome.index(genome.replace('-', '')[:1]) + crispr_cutsite - del_len
  read_start = read.index(read.replace('-', '')[:1])
  for jdx in range(read_start, len(read)):
    if read[jdx] == '-':
      break
  del_start = jdx - idx
  return del_start, jdx

def calc_insertion_start_position(read, genome, ins_len):
  # assumes exactly one insertion in alignment
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  counter = 0
  for idx in range(genome_start_idx, len(genome)):
    if counter == crispr_cutsite:
      break
    if genome[idx] in 'ACGTN':
      counter += 1
  for jdx in range(genome_start_idx, len(genome)):
    if genome[jdx] == '-':
      break
  ins_start = jdx
  ins_end = ins_start + ins_len
  if abs(idx - ins_start) < abs(idx - ins_end):
    ins_position = ins_start - idx
  else:
    ins_position = ins_end - idx
  return ins_position, jdx

def check_mismatches(read, genome, del_start, del_len):
  # Check that 3 bp on both sides of deletion match perfectly
  match_5side = match(read[del_start - 3 : del_start], genome[del_start - 3 : del_start])
  del_end = del_start + del_len
  match_3side = match(read[del_end : del_end + 3], genome[del_end : del_end + 3])
  if bool(match_5side and match_3side):
    # is ok
    return 'no'
  else:
    return 'yes'

def match(s1, s2):
  for c1, c2 in zip(s1, s2):
    if c1 != c2:
      return False
  return True

def has_mh(read, genome, del_len, gt_pos):
  # assumes single deletion that arises from mmej mechanism
  if gt_pos < 0 or gt_pos > del_len:
    return 'na'
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  cut_loci = genome_start_idx + crispr_cutsite
  read_base_5side = read[cut_loci - del_len + gt_pos - 1]
  genome_base_3side = genome[cut_loci + gt_pos - 1]
  try:
    read_base_3side = read[cut_loci + gt_pos]
  except:
    read_base_3side = '*'
  genome_base_5side = genome[cut_loci -del_len + gt_pos]
  if gt_pos not in [0, del_len]:
    if bool(read_base_5side == genome_base_3side) or bool(read_base_3side == genome_base_5side):
      return 'yes'
  if gt_pos == 0:
    if bool(read_base_3side == genome_base_5side):
      return 'yes'
  if gt_pos == del_len:
    if bool(read_base_5side == genome_base_3side):
      return 'yes'
  return 'no'

def standardized_del_gtpos(genome, del_len, del_start):
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  cutsite = genome_start_idx + crispr_cutsite
  if del_start < 0 or del_start > del_len:
    return del_start
  left = genome[cutsite - del_len : cutsite]
  right = genome[cutsite : cutsite + del_len]
  start_idx = max(len(right) - len(left), 0)
  mhs = []
  mh = [start_idx]
  for idx in range(min(len(right), len(left))):
    if left[idx] == right[start_idx + idx]:
      mh.append(start_idx + idx + 1)
    else:
      mhs.append(mh)
      mh = [start_idx + idx +1]
  mhs.append(mh)
  for mh in mhs:
    if del_start in mh:
      return max(mh)

##
# Store 
##
def add_del(inp_package, del_dd):
  '''
    Output: Store info in ins_dd
  '''
  read = inp_package['Read']
  genome = inp_package['Genome']

  del_len = get_indel_length(read)
  del_start, ds_pos = calc_deletion_start_position(read, genome, del_len)

  del_dd['Count'].append(inp_package['Count'])
  del_dd['Category'].append('del')
  del_dd['Condition'].append(inp_package['Condition'])
  del_dd['Peptide name'].append(inp_package['Peptide name'])

  del_dd['Microhomology-Based'].append(has_mh(read, genome, del_len, del_start))
  del_dd['Indel with Mismatches'].append(check_mismatches(read, genome, ds_pos, del_len))
  del_dd['Genotype Position'].append(del_start)
  del_dd['Length'].append(del_len)

  return

def add_ins(inp_package, ins_dd):
  '''
    Output: Store info in ins_dd
  '''
  read = inp_package['Read']
  genome = inp_package['Genome']

  ins_len = get_indel_length(genome)
  ins_start, is_pos = calc_insertion_start_position(read, genome, ins_len)

  ins_dd['Count'].append(inp_package['Count'])
  ins_dd['Category'].append('ins')
  ins_dd['Condition'].append(inp_package['Condition'])
  ins_dd['Peptide name'].append(inp_package['Peptide name'])

  ins_dd['Inserted Bases'].append(read[is_pos : is_pos + ins_len])
  ins_dd['Indel with Mismatches'].append(check_mismatches(read, genome, is_pos, ins_len))
  ins_dd['Genotype Position'].append(ins_start)
  ins_dd['Length'].append(ins_len)

  return

##
# Primary 
##
def genotype_data(inp_dir, out_dir, nm, start, end):
  print(nm)
  start, end = int(start), int(end)
  master_df = pd.DataFrame()

  # Parse condition-specific settings
  exp_row = exp_design[exp_design['Name'] == nm].iloc[0]
  lib_nm = exp_row['Dual library']
  target_nm = exp_row['Target']

  # Library design
  global lib_design
  lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv')
  peptide_nms = list(lib_design['Name'])
  peptide_nms = peptide_nms[start : end]
  print(peptide_nms)

  # Target 
  target_row = target_design[target_design['Target'] == target_nm].iloc[0]
  target = target_row['Sequence']
  global crispr_cutsite
  crispr_cutsite = int(target_row['Cutsite index'])

  master_df = pd.DataFrame()
  master_d = dict()

  timer = util.Timer(total = len(peptide_nms))
  for peptide_nm in peptide_nms:
    peptide_fn = inp_dir + f'{peptide_nm}.txt'
    if not os.path.isfile(peptide_fn): continue

    count_dd = defaultdict(lambda: 0)
    ins_dd = defaultdict(list)
    del_dd = defaultdict(list)

    with open(peptide_fn) as f:
      for i, line in enumerate(f):
        if i % 4 == 0:
          header = line.strip()
          try:
            count, category = parse_header(header)
          except:
            count = 0
            category = 'other'
        if i % 4 == 1:
          read = line.strip()
        if i % 4 == 2:
          genome = line.strip()
        if i % 4 == 3:
          if not is_main_category(category):
            # Count only
            count_dd[category] += count
            continue

          inp_package = {
            'Count': count,
            'Category': category,
            'Read': read,
            'Genome': genome,
            'Condition': nm,
            'Peptide name': peptide_nm,
          }

          # CRISPR-induced insertion or deletion
          if 'ins' in category:
            add_ins(inp_package, ins_dd)
          elif 'del' in category:
            add_del(inp_package, del_dd)

    count_df = pd.DataFrame({
      'Category': list(count_dd.keys()),
      'Count': list(count_dd.values()),
    })
    count_df['Condition'] = nm
    count_df['Peptide name'] = peptide_nm

    ins_df = pd.DataFrame(ins_dd)
    del_df = pd.DataFrame(del_dd)

    fdf = count_df.append(ins_df, ignore_index = True, sort = False)
    fdf = fdf.append(del_df, ignore_index = True, sort = False)
    master_d[peptide_nm] = fdf

    for df in [del_df, ins_df, count_df]:
      master_df = master_df.append(df, ignore_index = True, sort = False)

    timer.update()

  master_df.to_csv(out_dir + '%s_genotypes_%s.csv' % (nm, start))

  with open(out_dir + f'{nm}_genotypes_{start}.pkl', 'wb') as f:
    pickle.dump(master_d, f)

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
    lib_nm = exp_row['Dual library']
    lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv')

    jump = 150
    for jdx in range(0, len(lib_design), jump):
      start_jdx = jdx
      end_jdx = jdx + jump
      command = 'python %s.py %s %s %s' % (NAME, condition, start_jdx, end_jdx)
      script_id = NAME.split('_')[0]

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, condition, start_jdx)
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
  print(NAME)

  nm = argv[0]
  start = int(argv[1])
  end = int(argv[2])

  # Function calls
  genotype_data(
    '%s%s/' % (inp_dir, nm), 
    out_dir, 
    nm, start, end
  )
  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()