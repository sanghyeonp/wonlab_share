
## Odds ratio to Beta
Beta = ln(OR)
```
# In R
beta <- log(OR)
```
References:
- https://www.biostars.org/p/276869/
- https://www.biostars.org/p/336405/

## Beta to Odds ratio
```
or <- exp(beta)
```

## 95% confidence interval of Beta
```
# In R
beta_lo_ci <- beta - 1.96 * se
beta_up_ci <- beta + 1.96 * se
```

## Z-score from Beta
```
# In R
Z <- beta/se
```
