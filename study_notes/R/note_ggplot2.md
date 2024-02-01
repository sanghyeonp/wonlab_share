# Visualize specified colors
```
library(scales)
show_col(c("#9C89B8", "#F0A6CA", "#EFC3E6", "#F0E6EF"))
```
# Center the plot by origin
```
library(ggpmisc)

p1 +
  scale_x_continuous(limits = symmetric_limits) +
  scale_y_continuous(limits = symmetric_limits)
```
