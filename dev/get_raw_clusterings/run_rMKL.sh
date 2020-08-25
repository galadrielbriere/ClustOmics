#!/bin/bash

mkl_path=/usr/run_MKL_DR/application
matlab_path=/usr/local/MATLAB/MATLAB_Runtime/v90
num_dim=5
num_neigh=9
result_dir=./dev/get_raw_clusterings/rMKL_similarity_matrices
base_out_path=./data

cd $result_dir

for cancer_dir in $result_dir/* ; do
        cd $cancer_dir
        cancer=$(basename "$cancer_dir")
        for data in $cancer_dir/*; do

                if [ -d $data ]
                then
                        datatype=$(basename "$data")
                        ids_file="ids_"$datatype".txt"
                        ids=$cancer_dir/$ids_file
                        out_path=$base_out_path"/"$cancer
                        echo $data $ids $out_path
                        sh $mkl_path/run_run_MKL_DR.sh $matlab_path $data $ids $out_path $num_neigh $num_dim

                fi

        done

done
