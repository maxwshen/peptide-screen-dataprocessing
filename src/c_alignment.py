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
inp_dir = _config.OUT_PLACE + 'a_split/'
NAME = util.get_fn(__file__)
out_place = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_place)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
# lib_design = pd.read_csv(_config.DATA_DIR + f'lib_16_hdr_peptides_design.csv')

prefix_len = 8
suffix_len = 10

lib_design = None
prefixes = None
peptide_nms = None
prefix_to_peptide = None
suffixes = None
suffix_to_peptide = None

##
# Support
##
def find_peptide2_nm(read1):
  '''
  '''
  constant_region = 'GATCCTCCGCTAGA'
  if constant_region not in read1:
    return None, 'p2 missing constant'

  peptide_seq_start_idx = read1.index(constant_region) + len(constant_region)
  obs_suffix = read1[peptide_seq_start_idx : peptide_seq_start_idx + suffix_len]

  if obs_suffix in suffix_to_peptide:
    return suffix_to_peptide[obs_suffix], 'success'
  else:
    return None, 'p2 unmatched'


def find_peptide1_nm(read2):
  '''
  '''
  constant_region = 'CATGGCTAGC'
  if constant_region not in read2:
    return None, 'p1 missing_constant'

  peptide_seq_start_idx = read2.index(constant_region) + len(constant_region)
  obs_prefix = read2[peptide_seq_start_idx : peptide_seq_start_idx + prefix_len]

  if obs_prefix in prefix_to_peptide:
    return prefix_to_peptide[obs_prefix], 'success'
  else:
    return None, 'p1 unmatched'


##
# Alignments
##
def alignment(read, target):
  seq_align_tool = '/ahg/regevdata/projects/CRISPR-libraries/tools/seq-align/bin/needleman_wunsch'
  align = subprocess.check_output(seq_align_tool + ' --match 1 --mismatch -1 --gapopen -5 --gapextend -0 --freestartgap ' + read + ' ' + target, shell = True)
  align = align.decode('utf-8')
  if align[-2:] == '\n\n':
    align = align[:-2]
  return align


##
# IO
##
def store_alignment(alignment_buffer, peptide_nm, align_header, align, read_q):
  # Place quality scores in same configuration as read alignment
  read_align = align.split()[0]
  rq = ''
  curr_idx = 0
  for ch in read_align:
    if ch == '-':
      rq += ch
    else:
      rq += read_q[curr_idx]
      curr_idx += 1

  align_string = f'{align_header}\n{align}\n{rq}\n'
  if align_string.count('\n') != 4:
    import code; code.interact(local=dict(globals(), **locals()))
  alignment_buffer[peptide_nm].append(align_string)
  return


def init_alignment_buffer():
  alignment_buffer = defaultdict(list)
  return alignment_buffer


def flush_alignments(alignment_buffer, out_dir):
  print(f'Flushing... \n{datetime.datetime.now()}')
  for peptide_nm in alignment_buffer:
    with open(out_dir + f'{peptide_nm}.txt', 'a') as f:
      for align in alignment_buffer[peptide_nm]:
        f.write(align)
  print(f'Done flushing.\n{datetime.datetime.now()}')
  return


def prepare_outfns(out_dir, peptide_nms):
  for p1 in list(peptide_nms):
    for p2 in list(peptide_nms):
      out_fn = out_dir + f'{p1}-{p2}.txt'
      util.exists_empty_fn(out_fn)
  return


##
# Primary
##
def matchmaker(nm, split):
  print(split)
  stdout_fn = _config.SRC_DIR + f'nh_c_{nm}_{split}.out'
  util.exists_empty_fn(stdout_fn)
  out_dir = f'{out_place}{nm}/{split}/'
  util.ensure_dir_exists(out_dir)

  # Parse condition-specific settings
  exp_row = exp_design[exp_design['Name'] == nm].iloc[0]
  parent_fn = exp_row['Parent file']
  lib_nm = exp_row['Library']
  target_nm = exp_row['Target']

  # Library design
  global lib_design
  lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv')

  global prefixes
  global peptide_nms
  global prefix_to_peptide
  global suffixes
  global suffix_to_peptide
  prefixes = [s[:prefix_len] for s in lib_design['Sequence']]
  peptide_nms = list(lib_design['Name'])
  prefix_to_peptide = {prefix: nm for prefix, nm in zip(prefixes, peptide_nms)}
  suffixes = [compbio.reverse_complement(s[-suffix_len:]) for s in lib_design['Sequence']]
  suffix_to_peptide = {suffix: nm for suffix, nm in zip(suffixes, peptide_nms)}

  # Target 
  target_row = target_design[target_design['Target'] == target_nm].iloc[0]
  target = target_row['Sequence']
  target_strand = target_row['gRNA orientation']

  zf_split = str(split).zfill(3)
  read1_fn = inp_dir + f'{parent_fn}_R1_{zf_split}.fq'
  read2_fn = inp_dir + f'{parent_fn}_R2_{zf_split}.fq'

  count_stats = defaultdict(lambda: 0)
  count_stats['Success'] = 0

  alignment_buffer = init_alignment_buffer()
  prepare_outfns(out_dir, peptide_nms)

  tot_lines = util.line_count(read1_fn)
  timer = util.Timer(total = tot_lines)
  with open(read1_fn) as f1, open(read2_fn) as f2:
    for i, (line1, line2) in enumerate(zip(f1, f2)):
      if i % 4 == 0:
        h1 = line1.strip()
        h2 = line2.strip()
      if i % 4 == 1:
        read1 = line1.strip()
        read2 = line2.strip()
      if i % 4 == 3:
        q1, q2 = line1.strip(), line2.strip()
        count_stats['Read count'] += 1

        qs = [ord(s)-33 for s in q1 + q2]
        if np.mean(qs) < 25:
          count_stats['1a. Quality fail'] += 1
          continue

        res, msg = find_peptide1_nm(read2)
        if res is None:
          count_stats[f'2{msg}'] += 1
          continue
        p1_nm = res

        res, msg = find_peptide2_nm(read1)
        if res is None:
          count_stats[f'2{msg}'] += 1
          continue
        p2_nm = res

        peptide_nm = f'{p1_nm}-{p2_nm}'

        read1 = read1[6:]
        q1 = q1[6:]
        if target_strand == '-':
          read1 = compbio.reverse_complement(read1)
          q1 = q1[::-1]

        # Run alignment and store in buffer
        align_header = f'>1'
        align = alignment(read1, target)
        store_alignment(alignment_buffer, peptide_nm, align_header, align, q1)
        count_stats['Success'] += 1

      # flush_interval = 2000
      flush_interval = 200
      if i % int(tot_lines / flush_interval) == 1 and i > 1:
        # Flush alignment buffer
        flush_alignments(alignment_buffer, out_dir)
        alignment_buffer = init_alignment_buffer()

        # Stats for the curious
        with open(stdout_fn, 'a') as outf:
          outf.write(f'Time: {datetime.datetime.now()}\n')
          outf.write(f'Progress: {i / int(tot_lines / 100)}\n')
          outf.write(f'Line: {i}\n')
          for key in sorted(list(count_stats.keys())):
            outf.write(f'{key}, {count_stats[key]}\n')
        # break

      timer.update()
  
  # Final flush
  flush_alignments(alignment_buffer, out_dir)

  stats_df = pd.DataFrame(count_stats, index = [0])
  sorted_cols = sorted([s for s in stats_df.columns])
  stats_df = stats_df[sorted_cols]
  stats_df.to_csv(out_dir + f'stats_{nm}_{split}.csv')

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
  
  for nm in exp_design['Name']:
    for idx in range(0, _config.num_splits):

      command = 'python %s.py %s %s' % (NAME, nm, idx)
      script_id = NAME.split('_')[0]

      is_done = os.path.isfile(out_place + f'{nm}/{idx}/stats_{nm}_{idx}.csv')
      if is_done:
        continue

      # Write shell scripts
      sh_fn = qsubs_dir + 'q_%s_%s_%s.sh' % (script_id, nm, idx)
      with open(sh_fn, 'w') as f:
        f.write('#!/bin/bash\n%s\n' % (command))
      num_scripts += 1

      # Write qsub commands
      qsub_commands.append('qsub -j y -P regevlab -V -l h_rt=8:00:00 -wd %s %s &' % (_config.SRC_DIR, sh_fn))

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
  split = int(argv[1])

  # Function calls
  matchmaker(nm, split) 
  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()