#!/bin/bash
files=("file6.gz" "file7.gz" "file8.gz" "file9.gz")
for idx in ${!files[@]}
do
	if [[ $idx == 0 ]]; then
		echo "cd /remote/data/share/"
		echo "lcd /local/example/data"
	fi
	echo "get ${files[$idx]}"
done | sftp USER@example.com
