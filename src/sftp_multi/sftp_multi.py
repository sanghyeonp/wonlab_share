import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util import *
from packages import *



def parse_args():
    parser = argparse.ArgumentParser(description=":: SFTP multiple file downloader script generator ::")

    parser.add_argument('--filelist', required=True, type=str,
                    help='Path to file having list of file names. One file name per line.')
    parser.add_argument('--username', required=True, type=str,
                        help='SFTP user name.')
    parser.add_argument('--host', required=True, type=str,
                        help='SFTP host URL.')
    parser.add_argument('--remote_dir', required=True, type=str,
                        help='Remote server directory path where file exists.')
    parser.add_argument('--local_dir', required=True, type=str,
                        help='Local server directory path where the file will be copied.')
    parser.add_argument('--n_split', required=False, type=int, default=1,
                        help='Number splits in filelist')
    
    args = parser.parse_args()
    return args


def generate_cmd(filenames, username, host, remote_dir, local_dir):
    filenames = ' '.join(['"{}"'.format(f) for f in filenames])

    cmd = ['#!/bin/bash', 
            'files=({})'.format(filenames), 
            'for idx in ${!files[@]}', 
            'do', 
            '\tif [[ $idx == 0 ]]; then', 
            '\t\techo "{}"'.format(remote_dir), 
            '\t\techo "{}"'.format(local_dir),
            '\tfi',
            '\techo "get ${files[$idx]}"',
            'done | sftp {}@{}'.format(username, host)
            ]
    
    return [v + "\n" for v in cmd]


def main(filelist, username, host, remote_dir, local_dir, n_split):
    filenames = read_filelist(filelist)
    filename_chunks = split_into_chunks(filenames, n_split)

    for idx, filename_chunk in enumerate(filename_chunks):
        cmd = generate_cmd(filename_chunk, username, host, remote_dir, local_dir)
        with open('run_{}.sh'.format(idx + 1), 'w') as f:
            f.writelines(cmd)
        
        run_bash("chmod u+x run_{}.sh".format(idx + 1))


if __name__ == "__main__":
    args = parse_args()

    main(filelist=args.filelist, 
            username=args.username, 
            host=args.host, 
            remote_dir=args.remote_dir, 
            local_dir=args.local_dir, 
            n_split=args.n_split
            )
