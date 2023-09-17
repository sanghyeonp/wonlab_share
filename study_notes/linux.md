
# 서버에 새로운 group 및 user setting하고 directory assign 하기

### 서버에 만들어져 있는 group list 확인하는 법.
```
$ less /etc/group
```

### 새로운 group 만드는 법.
`groupadd`
```
$ sudo groupadd <group name>
```

### 새로운 group에 user 추가하는 법.
`usermod`
```
$ sudo usermod -a -G <group name> <user name>
```

### 특정 directory의 소유그룹을 변경.  
`chgrp` : 파일이나 directory의 소유그룹을 변경하는 명령어.
```
$ sudo chgrp -R <group name> <directory path>
```

### Directory 내의 파일의 read, write, execute 설정.  
`chmod` : 파일에 대한 read, write, execute 권한 설정하는 명령어.

8진수 형식으로 파일 모드 지정 방법
- 000 : 모든 사용자가 읽고 쓰고 실행할 수 없도록 지정.
- 777 : 모든 사용자가 읽고 쓰고 실행할 수 있는 권한 지정.
- 755 : 소유자는 모든 권한, 그룹 및 그 외 사용자는 읽기와 실행만 가능.
- 440 : 소유자 및 그룹은 읽기 가능, 그외 사용자는 권한 없음.

참고: https://recipes4dev.tistory.com/175
```
$ sudo chmod -R <8진수 모드> <directory path>
```

### 권한에 대해서.
기본적으로 10개의 letter로 표현함.  
0123456789  
- 0 : file인지 directory 인지. File이면 - / directory이면 d
- 123 : user (u)에 대한 권한
    - 1 : read
    - 2 : write
    - 3 : execute
- 456 : group (g)에 대한 권한
    - 4 : read
    - 5 : write
    - 6 : execute
- 789 : other (o)에 대한 권한
    - 7 : read
    - 8 : write
    - 9 : execute

예제.  
- -rw-rw-r-- : 파일이고, user는 read + write, group은 read + write, 그리고 other는 read
- drwxr-xr-x : directory이고, user는 read + write + execute, group은 read + execute, 그리고 other는 read + execute.

권한 변경은 `chmod` 활용.


### 소유자 변경.
`chown` : 소유자 변경 명령어
```
$ sudo chown -R <user name:user group> <directory path>
```

# Storage check

## <Directory path> 의 storage를 사용량으로 sorting해서.
```
# du -d 1 -h <Directory path> | sort -rh > storage_20230917.txt
```
