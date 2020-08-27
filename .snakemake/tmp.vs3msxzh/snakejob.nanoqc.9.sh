#!/bin/sh
# properties = {"type": "single", "rule": "nanoqc", "local": false, "input": ["Data/X_nanopore_runs1234.fastq"], "output": ["Analyses/nanoqc_report/X_nanopore_runs1234", "Analyses/nanoqc_logs/X_nanopore_runs1234_NanoQC.log"], "wildcards": {"sample": "X_nanopore_runs1234"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 9, "cluster": {"partition": "batch", "nodes": 1, "ntasks": 1, "ncpus": 8, "job-name": "Variant", "output": "slurm-%j-%x.out", "error": "slurm-%j-%x.err"}}
cd /home/user6/Mini-project-group-06 && \
/export/apps/snakemake/5.7.0/bin/python3.6 \
-m snakemake Analyses/nanoqc_logs/X_nanopore_runs1234_NanoQC.log --snakefile /home/user6/Mini-project-group-06/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /home/user6/Mini-project-group-06/.snakemake/tmp.vs3msxzh Data/X_nanopore_runs1234.fastq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules nanoqc --nocolor --notemp --no-hooks --nolock \
--mode 2  --use-conda  && touch "/home/user6/Mini-project-group-06/.snakemake/tmp.vs3msxzh/9.jobfinished" || (touch "/home/user6/Mini-project-group-06/.snakemake/tmp.vs3msxzh/9.jobfailed"; exit 1)

