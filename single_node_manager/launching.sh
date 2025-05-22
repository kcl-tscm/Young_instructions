#!/bin/bash -l
#$ -l h_rt=48:00:00
#$ -l mem=4.5G
#$ -pe mpi 40
#$ -N mpi_pool
#$ -A KCL_Cedric
#$ -P Gold
#$ -cwd                                                                                                                                                                                   

module purge
module load beta-modules
module load gcc-libs/10.2.0
module load compilers/gnu/10.2.0
module load gerun

echo "==================================================================="
echo "Starting per-core job submission (Manual Batching) on node: $(hostname)"
echo "PBS Job ID: $PBS_JOBID" # Use PBS_JOBID for PBS
echo "Current working directory: $(pwd)"
echo "PBS temporary directory: $PBS_O_WORKDIR" # PBS equivalent for original working directory
echo "Job started on: $(date)"
echo "==================================================================="

# --- Determine the number of physical cores using lscpu ---
CORES_PER_SOCKET=$(lscpu | grep "Core(s) per socket:" | awk '{print $NF}')
NUM_SOCKETS=$(lscpu | grep "Socket(s):" | awk '{print $NF}')

# Calculate total physical cores
NPROC=$((CORES_PER_SOCKET * NUM_SOCKETS))

echo "Number of parallel processes (NPROC) set to: $NPROC"

# --- Rest of the job execution logic (manual batching loop remains the same) ---

MAX_PARAM_FILE_INDEX=450 # <--- IMPORTANT: SET THIS TO YOUR HIGHEST par_ FILE INDEX
CURRENT_TASK_ID=1        # Start processing from par_1
RUNNING_JOBS=0           # Counter for currently active background jobs
declare -A PIDS          # Associative array to store PIDs and their corresponding TASK_ID

echo "Preparing to process $MAX_PARAM_FILE_INDEX parameter files."
echo "Launching up to $NPROC independent 'CoS' processes concurrently..."

while [ "$CURRENT_TASK_ID" -le "$MAX_PARAM_FILE_INDEX" ]; do
    if [ "$RUNNING_JOBS" -lt "$NPROC" ]; then
        PARAM_FILE="par_$CURRENT_TASK_ID"
        if [ ! -f "$PARAM_FILE" ]; then
            echo "Warning: Parameter file '$PARAM_FILE' not found. Skipping."
            CURRENT_TASK_ID=$((CURRENT_TASK_ID + 1))
            continue
        fi
        echo "Launching CoS with $PARAM_FILE (Task ID: $CURRENT_TASK_ID)..."
        ./CoS "$PARAM_FILE" &
        PID=$!
        PIDS[$PID]=$CURRENT_TASK_ID
        RUNNING_JOBS=$((RUNNING_JOBS + 1))
        CURRENT_TASK_ID=$((CURRENT_TASK_ID + 1))
    else
        echo "Maximum concurrent jobs ($NPROC) reached. Waiting for a slot to free up..."
        wait -n
        RUNNING_JOBS=$((RUNNING_JOBS - 1))
        echo "A job completed. $RUNNING_JOBS slots occupied, $NPROC total."
    fi
done

echo "All $MAX_PARAM_FILE_INDEX jobs have been launched. Waiting for any remaining jobs to complete..."
wait

echo "All CoS processes completed."
echo "==================================================================="
echo "Job finished on: $(date)"
echo "==================================================================="
