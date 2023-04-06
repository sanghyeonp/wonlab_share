#!/bin/bash
files=("file1.gz" "file2.gz" "file3.gz" "file4.gz" "file5.gz")
for idx in ${!files[@]}
do
	if [[ $idx == 0 ]]; then
		echo "cd /remote/data/share/"
		echo "lcd /local/example/data"
	fi
	echo "get ${files[$idx]}"
done | sftp USER@example.com
