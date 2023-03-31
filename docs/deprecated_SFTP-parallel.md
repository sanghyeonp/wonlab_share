
```diff
+ 이 script로 가능함. 하지만 SFTP 연결을 5개 이상이 되면 자주 느려지면서 멈춤.
```

## **설명** -> [README version 2](#sftp-file-transfer-in-parallel-version-2)

---

Public 데이터를 다루다보면 여러 대용량 파일들을 외부기관에서 받거나 또는 보내는 일이 생깁니다.  
그럴때, SFTP를 활용해서 데이터를 주고 받는 일이 많습니다.  
파일이 한두개라면 shell에 command를 입력해서 바로 다운로드 받을 수 있습니다.  
하지만 받고자 하는 파일이 여러개 일 때는 하나하나 다운받으려면 계속 확인하고 신경써야합니다.  
또, 파일을 하나하나씩 받기 때문에 시간이 오래 걸립니다.

그래서 여러 파일들을 한번에 다운로드 받을 수 있으면서, parallel로 받을 수 있게 script를 짜보았습니다.

---

# SFTP file transfer : Parallel

### **To-do**

### **Version control**
*Version 1:*
- Remote directory에 있는 다운로드 받을 파일 이름들을 가지는 파일 리스트를 input으로 넣어서 parallel하게 다운로드.

## How to use?
1. `git clone`을 활용해 관련 script를 다운로드 하기.
```
git clone https://github.com/sanghyeonp/wonlab_share.git
```

2. `src/sftp_parallel/sftp_parallel.py`를 이용하여 다운로드 진행하기.
```
python3 src/sftp_parallel/sftp_parallel.py --username <User name> \
                        --host <Host URL> \
                        --password <Password> \
                        --remote_dir <Remote server directory path> \
                        --local_dir <Local server directory path> \
                        --filelist <File path containing list of file names> \
                        --cores <Number of cores>
```

<ins>**아래는 현재 구현되어있는 arguments에 대한 설명입니다.**</ins>

- 아래 arguments들은 필수 입니다.  
`--username` :  SFTP username.  
`--host` :  SFTP host URL.  
`--password` :  SFTP password.  
`--remote_dir` :  Remote server에 있는 다운로드 하고자 하는 파일들이 있는 directory path.  
`--local_dir` :  Local server에 여기로 다운로드 하고싶은 directory path.  
`--filelist` : Remote server에서 다운로드 하고자 하는 파일들의 파일 이름들을 한줄에 하나씩 적혀진 파일 path.  

- 아래 arguments들은 **optional arguments** 입니다.  
`--cores` : Number of files to download simultaneously. Default = 1.


---
# SFTP command 설명
참고: https://docs.oracle.com/cd/E26502_01/html/E29001/remotehowtoaccess-14.html

Example:  
- username: USER
- password: PWD
- host: example.com
- remote directory: /remote/data/share/
- local directory: /local/example/data
- Remote directory에 있는 다운로드 할 파일 이름: FILE1.txt, FILE2.txt, FILE3.txt

## - SFTP에 접속해서 file transfer 하기 (Interactively)
참고: https://unix.stackexchange.com/questions/167266/sftp-command-to-get-download-tar-gz-file

### Step 1. SFTP로 remote에 연결하기
`sftp <USER NAME>@<HOST>`
```
$ sftp USER@example.com
```

### Step 2. Passsword 입력

### Step 3. Remote에서 다운로드 할 파일이 있는 directory로 들어가기
`cd <REMOTE DIRECTORY>`
```
> cd /remote/data/share/
```

### Step 4. 다운로드 할 파일이 저장 될 local directory 설정하기
`lcd <LOCAL DIRECTORY>`
```
> lcd /local/example/data/
```

### Step 5. 파일 다운로드 하기
`get <FILE NAME>`
```
> get FILE1.txt
```


## - Single-line command: Single file transfer 하기
참고: https://stackoverflow.com/questions/50096/how-to-pass-password-to-scp

`sshpass -p "<PASSWORD>" scp <USER NAME>@<HOST>:<REMOTE FILE PATH> <LOCAL DIRECTORY PATH>`

```
$ sshpass -p "PWD" scp USER@example.com:/remote/data/share/FILE1.txt /local/example/data
```

## - Single-line command: Single directory transfer 하기

[SFTP에 접속해서 file transfer 하기 (Interactively)](#--sftp에-접속해서-file-transfer-하기-interactively)에 Step 4. 까지는 동일.  
Step 5에서 파일 하나가 아닌, directory 전체를 다운로드.

`get -r <REMOTE DIRECTORY>`
```
> get -r /remote/data/share/
```

**Note:** Directory 전체를 다운로드 받아서, 여러 파일들을 다 다운로드 받을 수 있음. 하지만 Directory 전체가 한번에 다운로드 진행되지 않고, directory내에 파일 하나하나가 다운로드 진행됨. 그래서 다운로드 실행하고 신경 안 써도 되는 건 맞지만, 모든 파일을 다 다운로드 하는데 걸리는 시간이 파일 하나하나 다운로드 받는 시간과 동일.

---

# SFTP file transfer in parallel version 2

## Version 1이  안된 이유.
`sshpass` 및 `scp`를 이용하여 다운로드 받을 파일들을 one-liner command로 만들어서 python의 subprocess로 parallel하게 진행해보았지만, 실패했습니다.  
이유는 이렇게 접속을 하면 불안정해서 접속이 계속 끊어지게 되고, 다운로드가 중간에 계속 멈추게 됩니다.  

그래서 interactive하게 SFTP 접속을 해놓고 ([방법 참고 - sftp에 접속해서 file transfer 하기 interactively](./README_sftp.md#--sftp에-접속해서-file-transfer-하기-interactively)), 여러 파일을 다운로드 받아야겠다고 생각하고 진행해보았습니다.

---

## Version 2 방안

`screen` 만들어서 SFTP 접속 후, path 설정 후, `mget` command로 특정 여러 파일 다운로드하기.  

[sftp에 접속해서 file transfer 하기 interactively](./README_sftp.md#--sftp에-접속해서-file-transfer-하기-interactively)의 Step 4 까지는 동일.  

그리고, `ls`를 사용해서 다운로드 하고자 하는 파일 이름들을 terminal에서 얻고,

```
> ls
```

파일 이름들을 pattern을 토대로 여러 set들로 나누고,

*예를 들어, X001.gz, X002.gz, X003.gz, X101.gz, X102.gz, X20.gz, X203.gz 이런식으로 파일이 있으면,*  
*Set 1: X001.gz, X002.gz, X003.gz*  
*Set 2: X101.gz, X102.gz*  
*Set 3: X20.gz, X203.gz*


Step 5에서 하나의 파일 set에 대해 아래와 같이 진행.

`mget <FILENAME specified with pattern>`
```
# 예를 들어, set 1을 다운로드 하려면,
> mget X0*.gz

# 예를 들어, set 2을 다운로드 하려면,
> mget X1*.gz
```

그리고 나머지 파일 set들도 새로운 screen 만들어서 다시 Step 1 부터 5까지 진행.


### -> **단점**
- 귀찮음. screen을 여러개 켜야함.

### -> **참고**
- screen을 여러개 켜서 한다고, 하나당 core 하나 먹고 그러지는 않는 것 같음. (File transfer여서 resource를 사용하지 않나?)
- 하지만 screen 4개 이상 열어서, 같은 host에 4개 보다 많이 연결하면 전체적으로 다운로드 받는 속도가 느려짐.
- Wonlab에 `/data2/MGI/`를 이런식으로 받음.
    - 총 1,542개의 GWAS summary statistics를 받아야 했는데,
    - `/data2/MGI/00_download/01_generate_filelist.py`를 보면 위 예제처럼 이름이 되어있어서, X0*.output.gz 부터 X9*.output.gz을 나눠서 다운로드 진행함.

---

```diff
+ 혹시 다른 좋은 방법을 아시게 되면, 공유 해주시면 감사하겠습니다
```

---

### -> **확인해보기**
- 근데 SFTP로 접속해서 해도, 중간중간 "stalled"라고 다운로드가 느려지면서 멈췄다가 다시 진행하는걸 확인 함 --> 그러면 version 1 방법으로 했을 때, 다운로드가 중간에 멈췄는게 이렇게 "stalled"된 상태였던게 아닐까? 그걸 나는 잠깐 보고 안되네?라고 생각했을 수 있지 않을까?  
![SFTP stalled](https://github.com/sanghyeonp/wonlab_share/raw/master/img/sftp_stalled.png)
    - 되긴 됨. 하지만 4개 이상 연결하지 말 것.

