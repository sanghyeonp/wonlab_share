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


def logger(log_list, log, verbose=False):
    if isinstance(log, str):
        log_list.append(log)
    else:
        log_list += log
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


def split_into_chunks(list1, n_chunk):
    # Split list into n chunks
    k, m = divmod(len(list1), n_chunk)
    return list(list1[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n_chunk))


def read_filelist(filelist):
    with open(filelist, 'r') as f:
        files = [v.strip() for v in f.readlines()]
    return files


def isin_list(str1, list1):
    for l in list1:
        if l in str1:
            return True
    return False
