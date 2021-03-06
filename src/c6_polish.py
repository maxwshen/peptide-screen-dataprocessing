# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess
sys.path.append('/home/unix/maxwshen/')
import numpy as np
from collections import defaultdict
import pandas as pd
from mylib import util
import pickle

# Default params
inp_place = _config.OUT_PLACE + 'c_alignment/'
NAME = util.get_fn(__file__)
out_place = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_place)

exp_design = pd.read_csv(_config.DATA_DIR + 'exp_design.csv')
target_design = pd.read_csv(_config.DATA_DIR + 'lib_targets_design.csv')
lib_design = None

expected_cutsite = None
hdr_outcome = None

##
# Left or right
##
def left_or_right_event(seq, cutsite):
  # Given a sequence with exactly one - event and possible end gaps, determine if it's on the left or right of the cutsite.
  left_test = len(seq[:cutsite + 1].replace('-', ' ').split())
  right_test = len(seq[cutsite - 1:].replace('-', ' ').split())
  if left_test == 2 and right_test == 1:
    return 'left'
  if left_test == 1 and right_test == 2:
    return 'right'

###
# IO
###
def save_alignments(data, out_dir, peptide_nm):
  out_fn = out_dir + f'{peptide_nm}.txt'
  sorted_repair_types = sorted(list(data.keys()))
  with open(out_fn, 'w') as f:
    for repair_type in sorted_repair_types:
      package = sort_combine_alignments(data[repair_type])
      aligns, qs = package

      sorted_keys = sorted(aligns, key = aligns.get, reverse = True)
      qsm = summarize_quality(qs)
      for alignment in sorted_keys:
        total_count = aligns[alignment]
        qf = qsm[alignment]
        f.write(f'>{total_count}_{repair_type}\n{alignment}\n{qf}\n')
  return

def sort_combine_alignments(aligns):
  store = defaultdict(lambda: 0)
  qs = defaultdict(list)
  for idx, line in enumerate(aligns):
    if idx % 4 == 0:
      header = line
      count_section = header.split('_')[0]
      try:
        count = int(count_section.replace('>', ''))
      except:
        print('point a')
        import code; code.interact(local=dict(globals(), **locals()))
    if idx % 4 == 1:
      read = line
    if idx % 4 == 2:
      genome = line
    if idx % 4 == 3:
      q = line

      alignment = read + '\n' + genome

      store[alignment] += count
      qs[alignment].append(q)

  return store, qs

def summarize_quality(qs):
  new_qs = dict()
  for align in qs:
    read = align.split()[1]
    convert = lambda x: np.array([ord(s)-33 for s in x])
    qm = np.array([convert(s) for s in qs[align]]).T
    summary = [int(round(np.mean(qm[s]))) for s in range(len(qm))]
    convertb = lambda x: ''.join([chr(s+33) for s in x])
    new_qs[align] = convertb(summary)
  return new_qs

##
# Basic sequence operations 
##
def match_rate(s1, s2):
  assert len(s1) == len(s2)
  return sum([bool(c1==c2) for c1,c2 in zip(s1, s2)]) / len(s1)

def get_cutsite_idx(read, genome):
  # Given cutsite position in genome, find idx such that 
  # read[:idx], read[idx:], &
  # genome[:idx], genome[idx:] 
  # are the CRISPR cut products...
  #   even with starting gaps and arbitrary indels
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  counter = 0
  for jdx in range(genome_start_idx, len(genome)):
    if counter == expected_cutsite:
      break
    if genome[jdx] in 'ACGTN':
      counter += 1
  cut_idx = jdx
  return cut_idx

def count_indels(s1, s2):
  num_dels, num_ins = 0, 0
  ts1, ts2 = trim_start_end_dashes(s1), trim_start_end_dashes(s2)
  num_dels = len(ts1.replace('-', ' ').split()) - 1
  num_ins = len(ts2.replace('-', ' ').split()) - 1
  return num_dels, num_ins

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

def findinf(seq, query):
  # if query not in seq, return large number. 
  # If query in seq, return its index.
  if query not in seq:
    return np.inf
  else:
    return seq.index(query)


##
# Measure alignment properties
##
def score_single_deletion(read, genome):
  # Find where cut should be, with possible starting gaps, but exactly one deletion and no insertions in alignment.
  cutsite = get_cutsite_idx(read, genome)
  # seq[:cutsite], seq[cutsite:] gives cut products

  if '-' not in read[cutsite - 2 : cutsite + 2]:
    # Indel in 10bp window requires looking at 12bp, just like how indel in 0bp window requires looking at 2bp
    return ('delnotcrispr', 0)
  if read[cutsite - 1] != '-' and read[cutsite] != '-':
    return ('delnotatcut', 1)
  else:
    del_len = trim_start_end_dashes(read).count('-')
    return ('del%s' % (del_len), 2)

def score_single_insertion(read, genome):
  # Find where cut should be, with possible starting gaps, but exactly one insertion and no deletions in alignment.
  cutsite = get_cutsite_idx(read, genome)
  # seq[:cutsite], seq[cutsite:] gives cut products

  if '-' not in genome[cutsite - 2 : cutsite + 2]:
    # Indel in 10bp window requires looking at 12bp, just like how indel in 0bp window requires looking at 2bp
    return ('insnotcrispr', 0)
  if genome[cutsite - 1] != '-' and genome[cutsite] != '-':
    return ('insnotatcut', 1)
  else:
    ins_len = trim_start_end_dashes(genome).count('-')
    return ('ins%s' % (ins_len), 2)

def detect_wildtype(read, genome):
  num_dels, num_ins = count_indels(read, genome)
  if num_dels == 0 and num_ins == 0:
    if '-' not in trim_start_end_dashes(read) and '-' not in trim_start_end_dashes(genome):
      return True
  return False

def detect_hdr(read, genome):
  curr_idx = 0
  hdr_diffs = 0
  read_hdr_match = 0
  for obs_nt, ref_nt in zip(read, genome):
    if ref_nt != '-':
      hdr_nt = hdr_outcome[curr_idx]
      if hdr_nt != ref_nt:
        hdr_diffs += 1
        if hdr_nt == obs_nt:
          read_hdr_match += 1
      curr_idx += 1
  return bool(read_hdr_match == hdr_diffs)

def check_matches_threshold(read, genome):
  s1, s2 = '', ''
  for idx in range(len(genome) - 1, -1, -1):
    if genome[idx] in 'ACGTN':
      if read[idx] != '-':
        s1 += read[idx]
        s2 += genome[idx]
  if len(s1) == 0:
    return False
  match_score = match_rate(s1, s2)
  if len(s1) <= 5 or match_score <= 0.80:
    return False
  return True

def detect_pcr_dimer(read, genome):
  return bool(len(read.replace('-', '')) <= 50)

def detect_pcr_recombination_anywhere(read, genome):
  # Alignment signature of PCR recombination:
  #   1. long insertion, chance overlap, long deletion
  #   2. long deletion, chance overlap, long insertion

  ir, ig = read + '*', genome + '*'
  ir, ig, element_nm, package = pcr_chomp_process_elements(ir, ig)

  match_threshold = 0.80
  elements = [element_nm]
  packages = [package]

  while element_nm != 'done':
    ir, ig, element_nm, package = pcr_chomp_process_elements(ir, ig)
    elements.append(element_nm)
    packages.append(package)

  # Check for signature with random matches
  for idx in range(2, len(elements)):
    check_recombination_trio = False
    if elements[idx-2] == 'ins' and elements[idx-1] == 'match' and elements[idx] == 'del':
      check_recombination_trio = True
      ins_len, del_len = packages[idx-2], packages[idx]
      
    if elements[idx-2] == 'del' and elements[idx-1] == 'match' and elements[idx] == 'ins':
      check_recombination_trio = True
      ins_len, del_len = packages[idx], packages[idx-2]
    
    if check_recombination_trio:
      s1, s2 = packages[idx-1]
      match_score = match_rate(s1, s2)
      match_len = len(s1)
      both_long = bool(ins_len >= 10 and del_len >= 10)
      one_verylong = bool(ins_len >= 30 or del_len >= 30)
      if both_long or one_verylong:
        if match_score <= match_threshold:
          return True
        if match_len <= 5:
          return True
        if len(set(s2)) == 1:
          return True
      if bool(ins_len >= 30 and del_len >= 30):
        if match_score <= 0.95:
          return True
        if match_len <= 10:
          return True
        if len(set(s2)) == 1:
          return True

  # Check for signature with 0 matches
  for idx in range(1, len(elements)):
    check_recombination_duo = False
    if elements[idx-1] == 'ins' and elements[idx] == 'del':
      check_recombination_duo = True
      ins_len, del_len = packages[idx-1], packages[idx]

    if elements[idx-1] == 'del' and elements[idx] == 'ins':
      check_recombination_duo = True
      ins_len, del_len = packages[idx], packages[idx-1]

    if check_recombination_duo:
      return True

  return False

def pcr_chomp_process_elements(read, genome):
  if read[0] == '*' or genome[0] == '*':
    return read, genome, 'done', None

  if read[0] == '-' and genome[0] != '-':
    element_nm = 'del'
    length = read.index(read.replace('-', '')[:1])
    package = length
  if read[0] != '-' and genome[0] == '-':
    element_nm = 'ins'
    length = genome.index(genome.replace('-', '')[:1])
    package = length
  if read[0] != '-' and genome[0] != '-':
    element_nm = 'match'
    read_idx = min(findinf(read, '-'), findinf(read, '*'))
    genome_idx = min(findinf(genome, '-'), findinf(genome, '*'))
    length = min(read_idx, genome_idx)
    package = (read[:length], genome[:length])

  new_read = read[length:]
  new_genome = genome[length:]
  return new_read, new_genome, element_nm, package

def detect_combination_indel(read, genome):
  # Counts indels in 10 bp window centered at cutsite
  cutsite = get_cutsite_idx(read, genome)
  num_dels, num_ins = count_indels(read, genome)
  total_indels = num_dels + num_ins

  left_idx = 0
  counter = 0
  for idx in range(cutsite -1, -1, -1):
    if genome[idx] != '-' and read[idx] != '-':
      counter += 1
    if counter == 5:
      break
  left_idx = idx

  right_idx = 0
  counter = 0
  for idx in range(cutsite, len(genome),):
    if counter == 5:
      break
    if genome[idx] != '-' and read[idx] != '-':
      counter += 1
  right_idx = idx

  # Regions between and surrounding indels should have match score at least 80% if they are real combination indels, otherwise they're probably recombination products -- label as other.  
  ir, ig = read + '*', genome + '*'
  ir, ig, element_nm, package = pcr_chomp_process_elements(ir, ig)
  elements = [element_nm]
  packages = [package]
  while element_nm != 'done':
    ir, ig, element_nm, package = pcr_chomp_process_elements(ir, ig)
    elements.append(element_nm)
    packages.append(package)

  for element, package in zip(elements, packages):
    if element == 'match':
      s1, s2 = package
      match_score = match_rate(s1, s2)
      if match_score <= 0.80:
        return 'other'

  # If match score test passes, 
  #   if 0 indels in 10bp window: label as forgiven wt (not here)
  #   if 1 indel in 10bp window: label as forgiven indel
  #   if > 1 indel in 10bp window: label as combination indel
  readwindow = '*' + read[left_idx : right_idx] + '*'
  genomewindow = '*' + genome[left_idx : right_idx] + '*'
  num_dels, num_ins = count_indels(readwindow, genomewindow)
  
  num_indels_in_window = num_dels + num_ins
  if num_indels_in_window == 0:
    return 'combinationindelnotcrispr'
  
  if num_indels_in_window == 1:
    return 'forgivenindel'
  
  if num_indels_in_window >= 2:
    if num_indels_in_window == total_indels:
      return 'combinationindel'
    else:
      return 'forgivencombinationindel'
  return 'other'

##
# Categorize
##
def categorize_alignment(read, genome):
  # Alignment endgaps are assumed to be meaningless
  
  # check if homopolymer
  if len(set(read.replace('-', ''))) == 1:
    return 'homopolymer'

  if 'N' in read:
    return 'hasN'

  '''
    Detect PCR recombination within context (should be rare).
    Only anticipated for gRNA-target library assays. For editing at a fixed target, ignore this class.
  '''
  # if detect_pcr_recombination_anywhere(read, genome):
    # return 'pcr_recombination'

  # Check if context matches too poorly
  if not check_matches_threshold(read, genome):
    return 'poormatches'

  if detect_pcr_dimer(read, genome):
    return 'pcrdimer'

  # Categorize HDR
  if detect_hdr(read, genome):
    return 'hdr'

  # Categorize wildtype
  if detect_wildtype(read, genome):
    return 'wildtype'

  num_dels, num_ins = count_indels(read, genome)

  # Categorize single deletion
  if num_dels == 1 and num_ins == 0:
    category, score = score_single_deletion(read, genome)
    return category

  # Categorize single insertion
  elif num_dels == 0 and num_ins == 1:
    category, score = score_single_insertion(read, genome)
    return category

  # Categorize combination indel
  elif num_dels + num_ins >= 2:
    return detect_combination_indel(read, genome)
    # can be "forgiven_indel", 
    # "combination_indel",
    # "combination_indel_notcrispr",
    # "other"
  return 'other'

###
# Explore equal-scoring alignments: deletion shifting
###
def shift_single_deletion(read, genome, old_category):
  # Find where cut should be, with possible starting gaps, but exactly one deletion and no insertions in alignment.
  cutsite = get_cutsite_idx(read, genome)
  # seq[:cutsite], seq[cutsite:] gives cut products

  lr = left_or_right_event(read, cutsite)
  if lr == 'left':
    return rightshift_single_deletion(read, genome, old_category, cutsite)
  if lr == 'right':
    return leftshift_single_deletion(read, genome, old_category, cutsite)

def leftshift_single_deletion(read, genome, old_category, cutsite):
  # Deletion is on 3' side of cutsite. Try shifting deletion left
  del_start = cutsite + read[cutsite:].index('-')
  del_len = trim_start_end_dashes(read).count('-')
  del_end = del_start + del_len

  best_read = read
  best_category = old_category
  best_score = -1
  for jdx in range(del_start - 1, -1, -1):
    shift_len = del_start - jdx
    read_shift = read[jdx:del_start]
    genome_side = genome[del_end - shift_len : del_end] 
    if read_shift == genome_side:
      new_read = read[:jdx] + '-'*del_len + read_shift + read[del_end:]
      new_category, score = score_single_deletion(new_read, genome)
      if score >= best_score:
        best_score = score
        best_read = new_read
        best_category = new_category
    else:
      break
  if del_len in [1, 2] and 'not' in best_category:
    best_category += '_%s' % (del_len)
  return best_read, genome, best_category

def rightshift_single_deletion(read, genome, old_category, cutsite):
  # Deletion is on 5' side of cutsite. Try shifting deletion right
  read_start_idx = read.index(read.replace('-', '')[:1])
  del_start = read_start_idx + read[read_start_idx:].index('-')
  del_len = trim_start_end_dashes(read).count('-')
  del_end = del_start + del_len

  best_read = read
  best_category = old_category
  best_score = -1
  for jdx in range(del_end + 1, len(read)):
    shift_len = jdx - del_end
    read_shift = read[del_end:jdx]
    genome_side = genome[del_start : del_start + shift_len]
    if read_shift == genome_side:
      new_read = read[:del_start] + read_shift + '-'*del_len + read[del_end + shift_len:]
      new_category, score = score_single_deletion(new_read, genome)
      if score >= best_score:
        best_score = score
        best_read = new_read
        best_category = new_category
    else:
      break
  if del_len in [1, 2] and 'not' in best_category:
    best_category += '_%s' % (del_len)
  return best_read, genome, best_category

##
# Explore equal-scoring alignments: insertion shifting
def shift_single_insertion(read, genome, old_category):
  # Find where cut should be, with possible starting gaps, but exactly one deletion and no insertions in alignment.
  cutsite = get_cutsite_idx(read, genome)
  # seq[:cutsite], seq[cutsite:] gives cut products

  lr = left_or_right_event(genome, cutsite)
  if lr == 'left':
    return rightshift_single_insertion(read, genome, old_category, cutsite)
  if lr == 'right':
    return leftshift_single_insertion(read, genome, old_category, cutsite)

def leftshift_single_insertion(read, genome, old_category, cutsite):
  # Insertion is on 3' side of cutsite. Try shifting insertion left
  ins_start = cutsite + genome[cutsite:].index('-')
  ins_len = trim_start_end_dashes(genome).count('-')
  ins_end = ins_start + ins_len

  best_genome = genome
  best_category = old_category
  best_score = -1
  for jdx in range(ins_start - 1, -1, -1):
    shift_len = ins_start - jdx
    genome_shift = genome[jdx:ins_start]
    read_side = read[ins_end - shift_len : ins_end] 
    if genome_shift == read_side:
      new_genome = genome[:ins_start - shift_len] + '-'*ins_len + genome_shift + genome[ins_end:]
      new_category, score = score_single_insertion(read, new_genome)
      if score >= best_score:
        best_genome = new_genome
        best_score = score
        best_category = new_category
    else:
      break
  if ins_len == 1 and 'not' in best_category:
    best_category += '_%s' % (ins_len)
  return read, best_genome, best_category  

def rightshift_single_insertion(read, genome, old_category, cutsite):
  # Insertion is on 5' side of cutsite. Try shifting insertion left
  genome_start_idx = genome.index(genome.replace('-', '')[:1])
  ins_start = genome_start_idx + genome[genome_start_idx:].index('-')
  ins_len = trim_start_end_dashes(genome).count('-')
  ins_end = ins_start + ins_len

  best_genome = genome
  best_category = old_category
  best_score = -1
  for jdx in range(ins_end + 1, len(genome)):
    shift_len = jdx - ins_end
    genome_shift = genome[ins_end:jdx]
    read_side = read[ins_start : ins_start + shift_len]
    if genome_shift == read_side:
      new_genome = genome[:ins_start] + genome_shift + '-'*ins_len +  genome[ins_end + shift_len:]
      new_category, score = score_single_insertion(read, new_genome)
      if score >= best_score:
        best_genome = new_genome
        best_score = score
        best_category = new_category
    else:
      break
  if ins_len == 1 and 'not' in best_category:
    best_category += '_%s' % (ins_len)
  return read, best_genome, best_category

###
# Main function
###
def remaster_aligns(inp_fn, data):
  # timer = util.Timer(total = util.line_count(inp_fn))
  with open(inp_fn) as f:
    for i, line in enumerate(f):
      if i % 4 == 0:
        header = line.strip()
      if i % 4 == 1:
        read = line.strip()
      if i % 4 == 2:
        genome = line.strip()

      if i % 4 == 3:
        q = line.strip()

        # Main -- Find indel category, assuming end gaps are meaningless
        try:
          category = categorize_alignment(read, genome)
        except:
          print('point c')
          print(inp_fn)
          print(f'line {i}')
          import code; code.interact(local=dict(globals(), **locals()))

        if category in ['delnotatcut', 'delnotcrispr']:
          read, genome, category = shift_single_deletion(read, genome, category)
        if category in ['insnotatcut', 'insnotcrispr']:
          read, genome, category = shift_single_insertion(read, genome, category)

        alignment = [header, read, genome, q]
        # if header == '-----':
        # if '\x00\x00\x00\x00\x00\x00' in header:
          # print('point b')
          # import code; code.interact(local=dict(globals(), **locals()))

        if check_alignment(alignment):
          data[category] += alignment

      # timer.update()
  return


def check_alignment(alignment):
  [header, read, genome, q] = alignment
  if header[0] != '>':
    return False
  if len(alignment) != 4:
    return False
  if len(read) != len(genome):
    return False
  if len(read) != len(q):
    return False

  return True


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
    # if condition[-1] == 'B':
      # continue
    # 190625: library is test12
    # if larger, generate qsubs for subsets of library for parallelization

    lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv')

    jump = 10
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
  nm = argv[0]
  start = int(argv[1])
  end = int(argv[2])

  print(NAME)
  print(nm)

  # Parse condition-specific settings
  exp_row = exp_design[exp_design['Name'] == nm].iloc[0]
  lib_nm = exp_row['Dual library']
  target_nm = exp_row['Target']

  # Library design
  global lib_design
  lib_design = pd.read_csv(_config.DATA_DIR + f'lib_{lib_nm}_design.csv')
  peptide_nms = list(lib_design['Name'])
  peptide_nms = peptide_nms[start : end]

  # Target 
  target_row = target_design[target_design['Target'] == target_nm].iloc[0]
  target = target_row['Sequence']
  global expected_cutsite
  expected_cutsite = int(target_row['Cutsite index'])
  global hdr_outcome
  hdr_outcome = target_row['Sequence - HDR outcome']

  out_dir = out_place + nm +'/'
  util.ensure_dir_exists(out_dir)

  inp_dir = inp_place + nm + '/'

  timer = util.Timer(total = len(peptide_nms))
  for peptide_nm in peptide_nms:
    data = defaultdict(list)
    for split in os.listdir(inp_dir):
      if split == 'aligns':
        continue
      inp_fn = inp_dir + '%s/%s.txt' % (split, peptide_nm)
      remaster_aligns(inp_fn, data)
    print(peptide_nm)
    save_alignments(data, out_dir, peptide_nm)
    timer.update()

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(sys.argv[1:])
  else:
    gen_qsubs()