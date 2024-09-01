#!/bin/bash

train_extxyz="./train.extxyz"  # Path to the training extxyz file
valid_extxyz="./valid.extxyz"  # Path to the validation extxyz file

train_yaml="./pre_train.yaml"  # Path to the training script

n_graph_build_cores=4  # Number of CPU cores for data preprocessing

n_nodes=1  # Number of nodes for multi-GPU training
n_gpus=1  # Number of GPUs per nodes for multi-GPU training

echo "This is an example script for using SevenNet, from data preprocessing to training."
echo "It is recommended to adjust n_gpus in train_sevennet.sh and batch_size in pre_train.yaml for real training."
echo "Refer to convert_MPTrj2xyz.py to use MPtrj data. Note that *all* data points in MPtrj are used for SevenNet-0 (July112024)"


# check whether command 'sevenn' exist or not. If not exist, raise error and quit
if ! command -v sevenn >/dev/null 2>&1; then
    echo "SevenNet is not installed. Please check 'https://github.com/MDIL-SNU/SevenNet' for installation."
    exit 1
fi

for required in $train_extxyz $valid_extxyz; do
    if [ ! -f $required ]; then
        echo "No such file ${required}, training requires *.extxyz files as dataset."
        exit 1
    fi
done

cutoff=$(grep 'cutoff:' $train_yaml | awk '{print $2}')

# Build training graph data if it doesn't exist
if [ ! -f train.sevenn_data ]; then
    sevenn_graph_build -f ase -n $n_graph_build_cores -o train.sevenn_data $train_extxyz $cutoff
    mv graph_build_log train_graph_build.log
else
    echo "train.sevenn_data already exists, skipping graph build."
fi

# Build validation graph data if it doesn't exist
if [ ! -f valid.sevenn_data ]; then
    sevenn_graph_build -f ase -n $n_graph_build_cores -o valid.sevenn_data $valid_extxyz $cutoff
    mv graph_build_log valid_graph_build.log
else
    echo "valid.sevenn_data already exists, skipping graph build."
fi


if [ 1 -eq $n_gpus ] && [ 1 -eq $n_nodes ]; then
    sevenn pre_train.yaml -s
else
    # multi GPU training
    torchrun --standalone --nnodes=$n_nodes --nproc_per_node $n_gpus --no_python sevenn pre_train.yaml -d -s
fi
