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

파일 이름들을 여러 set들로 나누고,

Step 5에서 하나의 파일 set에 대해 아래와 같이 진행.

`mget <FILE NAME 1> <FILE NAME 2> <FILE NAME 3> <FILE NAME 4> ...`
```
> mget FILE1.txt FILE2.txt FILE3.txt FILE4.txt
```

그리고 나머지 파일 set들도 새로운 screen 만들어서 다시 Step 1 부터 5까지 진행.


### -> **단점**
- 귀찮음. screen을 여러개 켜야함.

### -> **참고**
- screen을 여러개 켜서 한다고, 하나당 core 하나 먹고 그러지는 않는 것 같음. (File transfer여서 resource를 사용하지 않나?)
- Wonlab에 `/data2/MGI/`를 이런식으로 받음.
    - 총 1,542개의 GWAS summary statistics를 받아야 했는데,
    - `/data2/MGI/00_download/01_generate_filelist.py`를 사용해서 기존에 받은 파일 (`/data2/MGI/00_download/exclude.list`)을 제외한 1,532개의 파일 이름들을 100개를 하나의 set로 구성해서 쪼갬.
    - screen 하나 당 파일 set 하나를 `mget`으로 다운로드 진행.