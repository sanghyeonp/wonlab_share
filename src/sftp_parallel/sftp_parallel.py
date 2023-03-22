import sys
import os
current = os.path.dirname(os.path.realpath(__file__))
parent = os.path.dirname(current)
sys.path.append(parent)

from util_sftp import *
from packages import *


def parse_args():
    parser = argparse.ArgumentParser(description=":: SFTP parallelization ::")

    parser.add_argument('--username', required=True, type=str,
                        help='SFTP user name.')
    parser.add_argument('--host', required=True, type=str,
                        help='SFTP host URL.')
    parser.add_argument('--password', required=True, type=str,
                        help='SFTP password.')
    parser.add_argument('--remote_dir', required=True, type=str,
                        help='Remote server directory path where file exists.')
    parser.add_argument('--local_dir', required=True, type=str,
                        help='Local server directory path where the file will be copied.')
    parser.add_argument('--filelist', required=True, type=str,
                        help='Path to file that contains the name of files to be copied from remote server.')
    parser.add_argument('--cores', required=False, type=int, default=1,
                        help='Number of cores to use.')
    
    args = parser.parse_args()
    return args


def main(username, host, password, remote_dir, local_dir, filelist, cores):
    file_list = read_filelist(f=filelist)
    remote_file_path_list = [os.path.join(remote_dir, f) for f in file_list]

    log_list = download_parallel(remote_file_path_list=remote_file_path_list, 
                                    password=password, 
                                    username=username, 
                                    hostname=host, 
                                    local_dir_path=local_dir, 
                                    fnc=download_, 
                                    cores=cores
                                    )
    
    with open("download.log", "w") as f:
        f.writelines(log_list)


if __name__ == "__main__":
    args = parse_args()

    main(username=args.username, 
        host=args.host, 
        password=args.password, 
        remote_dir=args.remote_dir, 
        local_dir=args.local_dir, 
        filelist=args.filelist, 
        cores=args.cores
        )