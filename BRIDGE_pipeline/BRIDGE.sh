#!/bin/bash
conda activate bridge_env
# Start a new tmux session and run commands
tmux new-session -d -s mysession 'bash -c "
    rm -f pipeline.png;
    /home/ubuntu/datapond/Dependencies/conda/bin/conda init bash;
    /home/ubuntu/datapond/Dependencies/conda/bin/conda activate bridge_env;
    #sudo mount /dev/xvdh /workdir_nextflow;
    sudo mount /dev/nvme1n1 /workdir_nextflow;

    
    # chmod -R 777 /workdir_nextflow;  # Uncomment if you want to change permissions
    nextflow run bridge_local.nf -params-file '"$1"' -resume -work-dir /workdir_nextflow/work;
"'

# Attach to the tmux session
tmux attach-session -t mysession