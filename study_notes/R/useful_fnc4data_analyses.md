
# skimr::skim
```
library(skimr)
df <- data.frame(A = c(1:10), B = c(rep("A", 5), rep("B", 5)))
skim(df)

```

# dput
```
df <- data.frame(A = c(1:10), B = c(rep("A", 5), rep("B", 5)))
dput(names(df))
# Output
"A", "B"
```
