#!/usr/bin/env python3


import os
import sys
import subprocess
import argparse
from datetime import datetime, timezone
import time
import json


def wipe_tmp():
  """
  Description:
    Wipe all paths from the ephemeral directory.
  """

  subprocess.run(
    f"rm -f ./tmp/interim*",
    shell=True
  )

def get_delta_t_h(t: float) -> str:
  """
  Description:
    Converts time from seconds (and a fraction of seconds), into human readable format with appropriate unit.

  Parameters:
    `t`. Time.

  Returns:
    `t_h`. Human readable time.
  """

  if t is None:
    return None

  h = int(t//60**2)
  m = int((t-h*60**2)//60)
  s = int(t-m*60-h*60**2)
  fr = round(1000*(t-s-m*60-h*60**2))

  if 24*60**2 < t:
    t_h = f"{h}h"
  elif 60**2 <= t and t < 24*60**2:
    t_h = f"{h:02d}h:{m:02d}m"
  elif 60 <= t and t < 60**2:
    t_h = f"{m:02d}m:{s:02d}s"
  elif 1 <= t and t < 60:
    t_h = f"{s:02d}s"
  elif 0 <= t and t < 1:
    t_h = f"{fr:03d}ms"

  return t_h

def subprocess_finalization(id: str):
  t_per_subprocess_f = time.time()
  delta_t_per_subprocess = t_per_subprocess_f - t_per_subprocess_i
  delta_t_h_per_subprocess = get_delta_t_h(t=delta_t_per_subprocess)
  session_report["delta_t_per_subprocess"].append(delta_t_per_subprocess)
  session_report["delta_t_h_per_subprocess"].append(delta_t_h_per_subprocess)
  print("[Task Completed]: Subprocess", id, "running time:", delta_t_h_per_subprocess)

# ===== Initialize =====

t_i = time.time()
datetime_now = datetime.now(timezone.utc).strftime("d%Y%m%dt%H%M%S")

# Directories.
dp_root_abs = os.path.abspath(os.path.dirname(__file__))
dp_interim = os.path.join(dp_root_abs, "data/interim")
dp_cluster = os.path.join(dp_root_abs, "data/cluster")
dp_tmp = os.path.join(dp_root_abs, "tmp")

acceptable_algorithms = ("mcl", "hipmcl")

parser = argparse.ArgumentParser(
  description="\
bioscflow: Perform community detection on a sequence similarity network.\n\
The vertex weight is set as the bit score estimated from an all-vs-all local sequence \n\
alignment (BLASTp), applied on some set of user provided database of FASTA sequences.\
  ",
  formatter_class=argparse.RawTextHelpFormatter  # Preserve newline characters in stdout of `--help`.
)

parser.add_argument(
  "--db",
  required=True,
  type=str,
  help="\
Filepath of FASTA database (canonical directory path in :/data/raw).\n\
Warning: Each record header must start with a unique identifier."
)
parser.add_argument(
  "-I",
  required=True,
  type=float,
  help="Inflation parameter of MCL."
)
parser.add_argument(
  "--algorithm",
  choices=["mcl", "hipmcl"],
  required=True,
  type=str,
  help="Clustering algorithm to use.\n|-> mcl: Reference MCL implementation\n|-> hipmcl: HipMCL (single-node mode)"
)
parser.add_argument(
  "--nopreproc",
  action="store_const",
  const=1,
  default=0,
  help=f"\
Disable preprocessing; apply community detection on an already prepared\n\
adjacency matrix.\n\
  "
)
parser.add_argument(
  "--debug",
  action="store_const",
  const=1,
  default=0
)

ARGS = parser.parse_args()

fp_interim0 = os.path.join(dp_tmp, "interim0")
debug = ARGS.debug

# Config.
evalue_ceil = 5e-3  # Acceptable E-value interval [0, EVALUE_CEIL].
inflation_prm = ARGS.I
n_threads_palladium = 96
algorithm = ARGS.algorithm

# Data.
fp_vertices = ARGS.db
basename_edges = ".".join(os.path.basename(fp_vertices).split(".")[:-1])
fp_casted_vertices = os.path.join(dp_interim, ".".join([basename_edges, "casta"]))

# Graph.
fp_edges_list = os.path.join(dp_interim, f"{basename_edges}.adj")
fp_edges_mci = os.path.join(dp_interim, f"{basename_edges}.mci")
fp_edges_tab = os.path.join(dp_interim, f"{basename_edges}.tab")
fp_cluster_list = os.path.join(dp_cluster, f"{basename_edges}_{algorithm}_I{inflation_prm}.cl")

if not(debug):
  fp_direct_session_report = os.path.join(dp_cluster, f"session_report_{basename_edges}_{algorithm}_I{inflation_prm}_{datetime_now}.json")

fp_session_report = os.path.join(dp_cluster, "session_report_last.json")

# Log.
session_report = {
  "algorithm": algorithm,
  "inflation": inflation_prm,
  "dataset_name": basename_edges,
  "dataset_size": None,
  "n_significant_edges": None,
  "n_clusters": None,
  "datetime": datetime_now,
  "delta_t": None,
  "delta_t_h": None,
  "subprocess_id": [],
  "delta_t_per_subprocess": [],
  "delta_t_h_per_subprocess": [],
  "exit_code_per_subprocess": [],
}

wipe_tmp()

# ===== Pre-process =====

if not os.path.exists(fp_vertices):
  raise FileNotFoundError(f'E: "{fp_vertices}" does not exist.')

session_report["dataset_size"] = int(
  subprocess.run(
    f'grep "^>" {fp_vertices} | wc -l',
    capture_output=True,
    text=True,
    shell=True
  ).stdout.strip()
)

# Cast.
if (os.path.exists(fp_casted_vertices) and os.path.getsize(fp_casted_vertices) > 0) or ARGS.nopreproc:
  print("[Warning]: Skipping cast.")
else:
  t_per_subprocess_i = time.time()
  session_report["exit_code_per_subprocess"].append(
    subprocess.run(
      f"libexec/cast {fp_vertices} > {fp_casted_vertices}",
      shell=True
    ).returncode
  )
  session_report["subprocess_id"].append("CAST")
  subprocess_finalization(id=session_report["subprocess_id"][-1])

if (os.path.exists(fp_edges_list) and os.path.getsize(fp_edges_list) > 0) or ARGS.nopreproc:
  print("[Warning]: Skipping sequence alignment.")
else:

  if not os.path.exists(fp_casted_vertices):
    raise FileNotFoundError(f'E: "{fp_casted_vertices}" does not exist.')

  # Sequence alignment.
  t_per_subprocess_i = time.time()
  session_report["exit_code_per_subprocess"].append(
    subprocess.run(
      f"libexec/diamond blastp --query {fp_casted_vertices} --db {fp_vertices} --out {fp_interim0} --outfmt 6 qseqid sseqid pident bitscore evalue",
      shell=True
    ).returncode
  )
  session_report["subprocess_id"].append("BLAST")
  subprocess_finalization(id=session_report["subprocess_id"][-1])

  # Eliminate low significance edges.
  # - Warning: The follow-up graph edges, of these vertex-pairs, will have 0 weight.
  t_per_subprocess_i = time.time()
  session_report["exit_code_per_subprocess"].append(
    subprocess.run(
      f"awk '$(NF)<={evalue_ceil} {{print $1,$2,$(NF-1)}}' {fp_interim0} > {fp_edges_list}",
      shell=True
    ).returncode
  )
  session_report["subprocess_id"].append("FILTER")
  subprocess_finalization(id=session_report["subprocess_id"][-1])

if not os.path.exists(fp_edges_list):
  raise FileNotFoundError(f'E: "{fp_edges_list}" does not exist.')

session_report["n_significant_edges"] = int(
  subprocess.run(
    f"wc -l {fp_edges_list} | awk 'NR==1 {{print $1}}'",
    capture_output=True,
    text=True,
    shell=True
  ).stdout.strip()
)

# ===== Detect Communities =====

# Community detection.
t_per_subprocess_i = time.time()
if algorithm == "mcl":
  if (os.path.exists(fp_edges_mci) and os.path.getsize(fp_edges_mci) > 0) or ARGS.nopreproc:
    print("[Warning]: Skipping MCI conversion.")
  else:
    subprocess.run(
      f"libexec/mcxload --stream-mirror -abc {fp_edges_list} -o {fp_edges_mci} -write-tab {fp_edges_tab}",
      shell=True
    )
  session_report["exit_code_per_subprocess"].append(
    subprocess.run(
      f"libexec/mcl {fp_edges_list} --abc -te {n_threads_palladium} -I {inflation_prm} -o {fp_cluster_list} -scheme 4",
      shell=True
    ).returncode
  )
elif algorithm == "hipmcl": 
  session_report["exit_code_per_subprocess"].append(
    subprocess.run(
      f"libexec/hipmcl -M {fp_edges_list} -I {inflation_prm} -per-process-mem 16 -o {fp_cluster_list}",
      shell=True
    ).returncode
  )
session_report["subprocess_id"].append("COMDET")
subprocess_finalization(id=session_report["subprocess_id"][-1])

if not os.path.exists(fp_cluster_list):
  raise FileNotFoundError(f'E: "{fp_cluster_list}" does not exist.')

# ===== Finalize =====

session_report["n_clusters"] = int(
  subprocess.run(
    f"wc -l {fp_cluster_list} | awk '{{print $1}}'",
    capture_output=True,
    text=True,
    shell=True
  ).stdout.strip()
)
t_f = time.time()
session_report["delta_t"] = t_f - t_i
session_report["delta_t_h"] = get_delta_t_h(t=session_report["delta_t"])
del t_i, t_f

with open(file=fp_session_report, mode="w") as f:
  json.dump(obj=session_report, fp=f, indent=2)

with open(file=fp_direct_session_report, mode="w") as f:
  json.dump(obj=session_report, fp=f, indent=2)

wipe_tmp()
sys.exit(any(session_report["exit_code_per_subprocess"]))
