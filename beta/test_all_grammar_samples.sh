#!/bin/bash

for folder in /Users/danielbergman/grammar_samples/user_projects/*; do
    if [ -d "$folder" ]; then
        folder_name=$(basename $folder)
        if [ "$folder_name" == "epi_caf_invasion" ]; then
            config_files="config/PhysiCell_settings_*.xml"
        elif [ "$folder_name" == "neuro_dev" ]; then
            config_files="config/PhysiCell_settings_*.xml"
        elif [ "$folder_name" == "tam_egf" ]; then
            config_files="config/PhysiCell_settings_*.xml"
        else
            config_files="config/PhysiCell_settings.xml"
        fi
        make load PROJ=$folder_name
        make clean
        make -j 20
        for config_file in $config_files; do
            python ./beta/test_run_sample.py project $config_file 60
        done
        make data-cleanup
        make reset
    fi
done