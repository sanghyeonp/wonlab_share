# SFTP file transfer in parallel

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
git clone https://github.com/sanghyeonp/sftp_parallel.git
```

2. `sftp_parallel.py`를 이용하여 다운로드 진행하기.
```
python3 sftp_parallel.py --username <User name> \
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