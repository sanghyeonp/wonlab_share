# Simple parallel computation
```
library(foreach)
library(doParallel)

n_parallel <- 2

x_list <- 1:5
y_list <- 10:15

func_addition <- function(x, y){
    z <- x + y
    return (z)
}


cl <- makeCluster(n_parallel)
registerDoParallel(cl)

foreach(idx = 1:length(x_list)
        ) %dopar% {
            
            res <- func_addition(x = x_list[idx],
                                y = y_list[idx])
            print(res)
        }

stopCluster(cl)
```

# .packages, .export, .verbose
```
library(dplyr)
library(foreach)
library(doParallel)

a1 <- 50

df1 <- data.frame(x = c(1:10),
                  y = c(11:20))

func <- function(dat, idx, a1){
    row <- dat[idx, ]
    row <- row %>% mutate(z = (x + y) * a1)
    return (row)
}

cl <- makeCluster(n_parallel)
registerDoParallel(cl)

foreach(idx = 1:nrow(df1),
        # Packages to be called across subprocesses
        .packages = c("dplyr"),
        # Name of variables to be used across subprocesses  
        .export = c("a1"),
        # Useful for troubleshooting
        .verbose = TRUE) %dopar% {
    
    res <- func(dat = df1,
                idx = idx,
                a1 = a1)
    print(res)
}
stopCluster(cl)
```

# .combine

```
library(dplyr)
library(foreach)
library(doParallel)

a1 <- 50

df1 <- data.frame(x = c(1:10),
                  y = c(11:20))

func <- function(dat, idx, a1){
    row <- dat[idx, ]
    row <- row %>% mutate(z = (x + y) * a1)
    return (row)
}

cl <- makeCluster(n_parallel)
registerDoParallel(cl)

df.res <- foreach(idx = 1:nrow(df1),
        # Packages to be called across subprocesses
        .packages = c("dplyr"),
        # Name of variables to be used across subprocesses  
        .export = c("a1"),
        # To combine the output of function as dataframe
        .combine = "bind_rows",
        # Useful for troubleshooting
        .verbose = TRUE) %dopar% {
    
    func(dat = df1,
                idx = idx,
                a1 = a1)
}
stopCluster(cl)
```
