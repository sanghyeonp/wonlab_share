from packages import *

def FileExists(file):
    return True if os.path.exists(file) else False


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
                'whitespace':' ',
                'formatted':'formatted'
                }
    return mapping[delim]


def read_formatted_file(file):
    with open(file, 'r') as f:
        rows = [row.strip() for row in f.readlines()]
        data = []
        for row in rows:
            row = row.split(sep=" ")
            row = [ele for ele in row if ele]
            data.append(row)
    
    df = pd.DataFrame(data[1:], columns=data[0])
    return df


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(io.StringIO(''.join(lines)),
                        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
                            'QUAL': str, 'FILTER': str, 'INFO': str},
                        sep='\t'
                    ).rename(columns={'#CHROM': 'CHROM'})