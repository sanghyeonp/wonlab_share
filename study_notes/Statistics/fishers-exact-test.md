

# Fisher's exact test

```
##########
# Fisher's exact test
# https://statsandr.com/blog/fisher-s-exact-test-in-r-independence-test-for-a-small-sample/
##########

dat <- data.frame(
  "smoke_no" = c(7, 0),
  "smoke_yes" = c(2, 5),
  row.names = c("Athlete", "Non-athlete"),
  stringsAsFactors = FALSE
)
colnames(dat) <- c("Non-smoker", "Smoker")

dat

mosaicplot(dat,
           main = "Mosaic plot",
           color = TRUE
)


chisq.test(dat)$expected


fisher.test(dat, alternative = "greater")

fisher.test(dat, alternative = "two.sided")

# one-sided test p-value
p <- (factorial(9)*factorial(5)*factorial(7)*factorial(7)) / (factorial(7)*factorial(2)*factorial(0)*factorial(5)*factorial(14))



```
