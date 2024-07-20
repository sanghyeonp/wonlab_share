```
elapsed <- function(start, end) {
    # https://stackoverflow.com/questions/32100133/print-the-time-a-script-has-been-running-in-r
    dsec <- as.numeric(difftime(end, start, unit = "secs"))
    hours <- floor(dsec / 3600)
    minutes <- floor((dsec - 3600 * hours) / 60)
    seconds <- dsec - 3600*hours - 60*minutes
    paste0(
        sapply(c(hours, minutes, seconds), function(x) {
            formatC(x, width = 2, format = "d", flag = "0")
        }), collapse = ":")
}

start <- Sys.time()

# CODE

end <- Sys.time()

print(elapsed(start, end))
```
