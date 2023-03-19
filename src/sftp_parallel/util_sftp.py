from src.packages import *
from src.util import run_bash

def download_(log_list, idx, password, username, hostname, remote_file_path, local_dir_path):
    cmd = 'sshpass -p "{}" scp {}@{}:{} {}'.format(password, username, hostname, remote_file_path, local_dir_path)
    start = datetime.now()
    log_ = "Download index: {:,}\n\t{}\n\tStart: {}\n".format(idx + 1, cmd, start); print(log_)
    run_bash(cmd)
    end = datetime.now()
    logs_ += "\tEnd: {}\n\tElapsed (hh:mm:ss.ms): {}\n".format(end, end - start)
    log_list.append(log_)


def download_parallel(remote_file_path_list, password, username, hostname, local_dir_path, fnc, cores=1):
    m = Manager()
    log_list = m.list()

    pool = Pool(processes=cores)
    inputs = [(log_list, idx, password, username, hostname, remote_file_path, local_dir_path) for idx, remote_file_path in enumerate(remote_file_path_list)]
    pool.starmap(fnc, inputs)
    pool.close()
    pool.join()

    return list(log_list)


def read_filelist(f):
    with open(f, 'r') as f1:
        rows = f1.readlines()
        file_list = [v.strip() for v in rows]
    return file_list
