# Multiple comparisons correction

안녕하세요? 박상현입니다.

분석을 하다보면 multiple comparisons correction (MCC)을 하게 될 경우가 많이 생깁니다.

공유드린 script를 활용하면, MCC 중 FDR과 Bonferroni correction에 대한 결과를 돌려줍니다.

---

### **To-do**


### **Version control**
*Version 1:*
- FDR, Bonferroni correction 할 수 있음.

---

## How to use?
1. `git clone`을 활용해 관련 script를 다운로드 하기.
```
git clone https://github.com/sanghyeonp/MCC.git
```

2. `mcc.py`를 이용하여 multiple comparisons correction 진행하기.
```
python3 mcc.py --file <file path> \
                --p_col <P-value column 이름>
```
&nbsp;

<ins>**아래는 현재 구현되어있는 arguments에 대한 설명입니다.**</ins>

- 아래 arguments들은 필수 입니다.  
`--file` :  Path to the input file.  
`--p_col` :  Name of the P-value column.  

- 아래 arguments들은 input file과 관련된 **optional arguments** 입니다.  
`--delim` : Input file의 delimiter. Choices = ['tab', 'comma', 'whitespace']. Default = 'tab'.  
`--skip_rows` : 처음 몇개 줄을 skip할 건지.  

- 아래 arguments들은 output file과 관련된 **optional arguments** 입니다.  
`--outf` : Name of the output file. Default = mcc.InputFileName.  
`--outd` : Directory path of the output file. Default = Current working directory.  
`--quoting` : (flag) Specify to quote the data in the output file. Default = False.  

- 아래 arguments들은 **optional arguments** 입니다.  
`--log` : (flag) Save the log. Default = False.  
`--verbose` : (flag) Verbose output. Default = False.

---

## Example

`/example` 에 있는 예제를 참고하세요.

---
## Contributors

---
## Reference