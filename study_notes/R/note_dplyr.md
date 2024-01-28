

# Retain all duplicates
```
df2 <- df1[duplicated(df1$rsid) | duplicated(df1$rsid, fromLast = TRUE), ]
```

# Leave only unique
```
df2 <- df1 %>%
    distinct(SNP, .keep_all = TRUE)
```

# Resolve conflict
```
# install.packages("conflicted")
library(conflicted)
# set conflict preference
conflict_prefer("select", "dplyr")
```
