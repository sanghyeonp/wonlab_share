# SFTP file transfer in parallel version 2

## Version 1이  안된 이유.
`sshpass` 및 `scp`를 이용하여 다운로드 받을 파일들을 one-liner command로 만들어서 python의 subprocess로 parallel하게 진행해보았지만, 실패했습니다.  
이유는 이렇게 접속을 하면 불안정해서 접속이 계속 끊어지게 되고, 다운로드가 중간에 계속 멈추게 됩니다.  

그래서 interactive하게 SFTP 접속을 해놓고 ([방법 참고 - sftp에 접속해서 file transfer 하기 interactively](./README.md#--sftp에-접속해서-file-transfer-하기-interactively)), 여러 파일을 다운로드 받아야겠다고 생각하고 진행해보았습니다.

---

## Version 2 방안

`screen` 만들어서 SFTP 접속 후, path 설정 후, `mget` command로 특정 여러 파일 다운로드하기.  

[sftp에 접속해서 file transfer 하기 interactively](./README.md#--sftp에-접속해서-file-transfer-하기-interactively)의 Step 4 까지는 동일.  

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
![SFTP stalled](./imgs/sftp_stalled.png)
    - 되긴 됨. 하지만 4개 이상 연결하지 말 것.

