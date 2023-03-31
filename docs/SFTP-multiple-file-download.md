

Public 데이터를 다루다보면 여러 대용량 파일들을 외부기관에서 받거나 또는 보내는 일이 생깁니다.  
그럴때, SFTP를 활용해서 데이터를 주고 받는 일이 많습니다.  
파일이 한두개라면 SFTP에 command를 입력해서 바로 다운로드 받을 수 있습니다.  
하지만 받고자 하는 파일이 여러개 일 때는 하나하나 다운받으려면 계속 확인하고 신경써야합니다.  
또, 파일을 하나하나씩 받기 때문에 시간이 오래 걸립니다.

그래서 여러 파일들을 한번의 command 수행으로 다운로드 받을 수 있게 script를 짜보았습니다.

여러 파일을 한번의 commend 수행으로 다운로드 받을 수 있는 방법들은 아래와 같이 많습니다.  

방법 1: 다운로드 받을 파일들이 directory에 포함되어 있다면  

`get -R`  [참고](#sftp-접속-후-single-directory-transfer-하기-interactively)

방법 2: 다운로드 받을 파일들 이름의 일정한 pattern이 있다면   
`mget `  [참고](#sftp-접속-후-multiple-file-transfer-하기-interactively)

하지만 위 방법 1과 방법 2로 하다보면 중간에 SFTP 연결이 끊켜서 다시 이어서 진행해야 하는데, `get -R`로는 이미 다운로드 받은 파일들도 다시 다운로드 받게 되고, `get`만 활용하면, 한번에 하나의 파일만 specify 할 수 있고, 그리고 `mget`은 파일 이름들 패턴을 줘야하는 단점들이 있습니다.  

그래서 차라리 파일 이름 list를 만들고, list를 이용해서 다운로드 할 수 있다면, 전체 directory를 다운로드 받는 다던가, 파일 이름들에 특정 패턴이 있어야 한다던가, 아니면 중간에 다운로드가 끊켜도 list만 수정해서 다시 여러 파일을 commend 실행 한번으로 받을 수가 있겠습니다.

방법 3: 파일 이름에 일정한 패턴이 없고, 여러 파일을 다운로드 받고 싶다면  
`use script`

---

# SFTP multiple file transfer script generator

**원리**: `sftp_multi.py`를 이용해서 여러 파일을 한번의 bash script execution으로 다운로드 할 수 있게 bash script를 만들어 줌.  
bash script가 만들어지면, bash script를 실행 후, password를 입력하면 자동으로 해당 bash script에 list된 파일들이 순서대로 다운로드 진행됨.  




### **To-do**

### **Version control**
*Version 1:*


## How to use?
1. `git clone`을 활용해 관련 script를 다운로드 하기.
```
git clone https://github.com/sanghyeonp/wonlab_share.git
```

2. `src/sftp_multi/sftp_multi.py`를 이용하여 bash script 만들기.
```
python3 src/sftp_multi/sftp_multi.py --filelist <File path containing list of file names> \
                        --username <User name> \
                        --host <Host URL> \
                        --remote_dir <Remote server directory path> \
                        --local_dir <Local server directory path> \
                        --n_split <Number of runnable file to generate>
```

<ins>**아래는 현재 구현되어있는 arguments에 대한 설명입니다.**</ins>

- 아래 arguments들은 필수 입니다.  
`--filelist` : Remote server에서 다운로드 하고자 하는 파일들의 파일 이름들을 한줄에 하나씩 적혀진 파일 path.  
`--username` :  SFTP username.  
`--host` :  SFTP host URL.  
`--remote_dir` :  Remote server에 있는 다운로드 하고자 하는 파일들이 있는 directory path.  
`--local_dir` :  Local server에 여기로 다운로드 하고싶은 directory path.  


- 아래 arguments들은 **optional arguments** 입니다.  
`--n_split` : Number of runnable files to generate. Default = 1.  

3. `screen`을 켠 후.

4. 2.에서 만들어진 bash script를 실행하고, prompt에 password 치라고 나오면, password 입력하면 자동으로 다운로드 진행됨.

위와 같은 방법으로 하면, 여러 파일을 한번의 bash script 실행으로 다운로드 자동으로 됨.  
하지만 여러파일을 **동시에** 다운로드가 되는 것이 아님. 한번에 한 파일씩.  
여러파일을 동시에 다운로드를 하려면, `--n_split`로 다운로드 하고자 하는 파일들을 여러 set로 나누면, 각 set에 포함된 파일들을 다운로드 할 수 있는 여러 bash script가 만들어짐.  
그러면 `screen`을 여러개 켜서, 한 `screen` 당 하나의 bash script를 실행 시켜놓으면, 여러 파일을 동시에 다운로드 할 수 있음.  

---
# Example
`/data1/sanghyeon/wonlab_contribute/combined/examples/sftp_multi`

---

# - SFTP 접속 후, single file transfer 하기 (Interactively)
참고: https://unix.stackexchange.com/questions/167266/sftp-command-to-get-download-tar-gz-file

## Step 1. SFTP로 remote에 연결하기
`sftp <USER NAME>@<HOST>`
```
$ sftp USER@example.com
```

## Step 2. Passsword 입력

## Step 3. Remote에서 다운로드 할 파일이 있는 directory로 들어가기
`cd <REMOTE DIRECTORY>`
```
> cd /remote/data/share/
```

## Step 4. 다운로드 할 파일이 저장 될 local directory 설정하기
`lcd <LOCAL DIRECTORY>`
```
> lcd /local/example/data/
```

## Step 5. 파일 다운로드 하기
`get <FILE NAME>`
```
> get FILE1.txt
```

# - SFTP 접속 후, single directory transfer 하기 (Interactively)

[SFTP에 접속해서 file transfer 하기 (Interactively)](#sftp-접속-후-single-file-transfer-하기-interactively)에 Step 4. 까지는 동일.  
Step 5에서 파일 하나가 아닌, directory 전체를 다운로드.

`get -r <REMOTE DIRECTORY>`
```
> get -r /remote/data/share/
```

Note: *Directory 전체를 다운로드 받아서, 여러 파일들을 다 다운로드 받을 수 있음. 하지만 Directory 전체가 한번에 다운로드 진행되지 않고, directory내에 파일 하나하나가 다운로드 진행됨. 그래서 다운로드 실행하고 신경 안 써도 되는 건 맞지만, 모든 파일을 다 다운로드 하는데 걸리는 시간이 파일 하나하나 다운로드 받는 시간과 동일.*


# - SFTP 접속 후, multiple file transfer 하기 (Interactively)

[SFTP에 접속해서 file transfer 하기 (Interactively)](#sftp-접속-후-single-file-transfer-하기-interactively)에 Step 4. 까지는 동일.  
Step 5에서 파일 이름들의 pattern을 이용해서 여러 파일 다운로드 받기

*예를 들어, X001.gz, X002.gz, X003.gz, X101.gz, X102.gz, X20.gz, X203.gz 이런식으로 파일이 있으면,*  

`mget <FILENAME specified with pattern>`
```
> mget X*.gz
```

---

# SCP command 설명
# - Single-line command: Single file transfer 하기 (Non-interactively)
참고: https://stackoverflow.com/questions/50096/how-to-pass-password-to-scp

`sshpass -p "<PASSWORD>" scp <USER NAME>@<HOST>:<REMOTE FILE PATH> <LOCAL DIRECTORY PATH>`

```
$ sshpass -p "PWD" scp USER@example.com:/remote/data/share/FILE1.txt /local/example/data
```
