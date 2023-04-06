#!/bin/bash

script=../../src/sftp_multi/sftp_multi.py

python3 ${script} --filelist filename.list \
                --username USER \
                --host example.com \
                --remote_dir /remote/data/share/ \
                --local_dir /local/example/data \
                --n_split 2
