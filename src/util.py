from src.packages import *

def run_bash(bash_cmd):
    """
    Run bash command.
    Return a list containing standard output, line by line.
    """
    popen = subprocess.Popen(bash_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, _ = popen.communicate()
    return str(stdout, 'utf-8').strip().split('\n')


def logger(log_list, log, verbose):
    log_list.append(log)
    if verbose:
        print(log)
    return log_list


def save_log(log_list, out):
    with open(out, 'w') as f:
        f.writelines([l+'\n' for l in log_list])


def map_delim(delim):
    mapping = {'tab':'\t', 
                'comma':',',
                'whitespace':' '
                }
    return mapping[delim]

