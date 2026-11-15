#!/bin/bash
#SBATCH --job-name=roks-frozen
#SBATCH --partition=cpu-skylake
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%A_%a.out

# Diagnostic only: CP2K intentionally stops after its first SCF energy.
# Success here means that energy was evaluated on the immutable seed density.
set -euo pipefail
campaign_dir="$HOME/agent-calcdir/cp2k-delta-scf/roks-st-gap-2026-11-13"
cd "$campaign_dir"
mapfile -t trial_ids < trials.txt
trial_id="${trial_ids[$SLURM_ARRAY_TASK_ID]}"
cd "trials/$trial_id"
test ! -e probe.out

saved_srun=$(command -v srun)
saved_tmpdir="${TMPDIR:-/tmp}"
if type module &>/dev/null; then module purge; fi
set +u
source "$HOME/opt/cp2k/main/cp2k_env"
set -u
mpi_launcher=$(command -v mpirun)
cp2k_binary=$(command -v cp2k.psmp)
export PATH="$(dirname "$mpi_launcher"):$(dirname "$saved_srun"):$PATH"
export TMPDIR="$saved_tmpdir"
export PRTE_MCA_prte_tmpdir_base="$TMPDIR"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
export OMP_PLACES=cores OMP_PROC_BIND=true
export OMPI_MCA_coll='^hcoll' OMPI_MCA_pml=ob1 OMPI_MCA_btl=self,sm
export SKALA_MODEL="${SKALA_MODEL:-${GAUXC_SKALA_MODEL:-}}"
export GAUXC_SKALA_MODEL="$SKALA_MODEL"
test -r "$SKALA_MODEL"
export PYTHONPATH="$campaign_dir/code"
project_python="$HOME/work/cp2k-deltascf/.venv/bin/python"
test -x "$project_python"

{
  date -u +%FT%TZ
  hostname
  printf 'job=%s array=%s index=%s trial=%s\n' "$SLURM_JOB_ID" "$SLURM_ARRAY_JOB_ID" "$SLURM_ARRAY_TASK_ID" "$trial_id"
  printf 'TMPDIR=%s\nCP2K=%s\nMPI=%s\nMODEL=%s\n' "$TMPDIR" "$cp2k_binary" "$mpi_launcher" "$SKALA_MODEL"
  sha256sum "$cp2k_binary" "$SKALA_MODEL" "$HOME/opt/cp2k/main/cp2k_env" seed.wfn probe.inp
  sha256sum "$HOME/opt/cp2k/main/share/cp2k/data/EMSL_BASIS_SETS" "$HOME/opt/cp2k/main/share/cp2k/data/POTENTIAL"
  "$cp2k_binary" --version
  "$mpi_launcher" --version
} > runtime.txt 2>&1
set +e
"$mpi_launcher" -np "$SLURM_NTASKS" --bind-to core "$cp2k_binary" -i probe.inp -o probe.out > launcher.log 2>&1
cp2k_status=$?
set -e
printf '%s\n' "$cp2k_status" > cp2k_exit_code.txt
"$project_python" -c 'from roks_st_gap_common import validate; validate(".")'  > validation.log 2>&1
gzip -n -c density.Log > density.Log.gz.tmp
gzip -t density.Log.gz.tmp
mv density.Log.gz.tmp density.Log.gz
rm density.Log
