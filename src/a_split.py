# 
from __future__ import division
import _config
import sys, os, fnmatch, datetime, subprocess
sys.path.append('/home/unix/maxwshen/')
import numpy as np
from collections import defaultdict
from mylib import util
import pandas as pd

# Default params
inp_dir = _config.DATA_DIR
NAME = util.get_fn(__file__)
out_dir = _config.OUT_PLACE + NAME + '/'
util.ensure_dir_exists(out_dir)

exp_design = pd.read_csv(inp_dir + f'exp_design.csv')

##
# Functions
##
def split(inp_fn, out_nm):
  inp_fn_numlines = util.line_count(inp_fn)

  num_splits = 60
  split_size = int(inp_fn_numlines / num_splits)
  if num_splits * split_size < inp_fn_numlines:
    split_size += 1
  while split_size % 4 != 0:
    split_size += 1
  # print 'Using split size %s' % (split_size)

  split_num = 0
  timer = util.Timer(total = num_splits)
  for idx in range(1, inp_fn_numlines, split_size):
    start = idx
    end = start + split_size  
    out_fn = out_dir + out_nm + '_%s.fq' % (split_num)
    
    skip = False
    if os.path.isfile(out_fn):
      size_mb = os.path.getsize(out_fn) / 1e6
      if size_mb > 0:
        skip = True

    if not skip:
      command = 'tail -n +%s %s | head -n %s > %s' % (start, inp_fn, end - start, out_fn)
      subprocess.check_output(command, shell = True)

    split_num += 1
    # print(command)
    timer.update()

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
  fns = [fn for fn in os.listdir(inp_dir) if '.fq' in fn]

  for fn in fns:
    command = 'python %s.py %s' % (NAME, fn)
    script_id = NAME.split('_')[0]

    # Write shell scripts
    sh_fn = qsubs_dir + 'q_%s_%s.sh' % (script_id, fn)
    with open(sh_fn, 'w') as f:
      f.write('#!/bin/bash\n%s\n' % (command))
    num_scripts += 1

    # Write qsub commands
    qsub_commands.append('qsub -V -l h_rt=2:00:00,h_vmem=1G -wd %s %s &' % (_config.SRC_DIR, sh_fn))

  # Save commands
  commands_fn = qsubs_dir + '_commands.sh'
  with open(commands_fn, 'w') as f:
    f.write('\n'.join(qsub_commands))
  subprocess.check_output('chmod +x %s' % (commands_fn), shell = True)

  print('Wrote %s shell scripts to %s' % (num_scripts, qsubs_dir))
  return

##
# Main
##
@util.time_dec
def main(nm = ''):
  print(NAME)
  
  split(inp_dir + nm, nm.replace('.fq', ''))

  # # Function calls
  # for fn in os.listdir(inp_dir):
  #   if '.fastq' in fn:
  #     split(inp_dir + fn, fn.replace('.fastq', ''))

  return


if __name__ == '__main__':
  if len(sys.argv) > 1:
    main(nm = sys.argv[1])
  else:
    gen_qsubs()