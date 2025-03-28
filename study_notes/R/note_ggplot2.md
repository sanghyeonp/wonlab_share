# How expand in scale_?_continuous() works?
Reference: https://stackoverflow.com/questions/44170871/how-does-ggplot-scale-continuous-expand-argument-work

```
SOME GGPLOT OBJECT +
scale_x_continuous(limits = c(1, 7), 
                   expand = expansion(mult = c(0, 0.5), 
                                         add = c(2, 0))
# left most position will be 1 - (7-1) * 0.0  -2 = -1, 
# right most position will be 7 + (7-1) * 0.5 = 10
```

# Axis title with markdown
```
library(ggtext)
SOME GGPLOT OBJECT +
  labs(y = "-log<sub>10</sub>(*P*)") +
  theme(axis.title.y = ggtext::element_markdown(size=12, face="bold", family="Helvetica", color="black")

### 그냥 PDF로 저장해도 됨..?
# As PDF
ggsave(..., device=pdf)

### When saving the plot using ggsave(), save it with Cario-based (단점: Helvetica font가 없음 cario는)
# As PNG
ggsave(..., device=png, type="cairo")
# As PDF
ggsave(..., device=cairo_pdf)
```

Subscript: "\<sub>X\</sub>"\
Superscript: "\<sup>X\</sup>"\
Italic: "\*X*"\
Bold: "\**X**"\
Break line: "\<br>"\
Different text style (e.g., font, color): "\<span style = 'color:red;font-size:8pt'>X\</span>"\
En-dash: "--"


# Visualize specified colors
```
library(scales)
show_col(c("#9C89B8", "#F0A6CA", "#EFC3E6", "#F0E6EF"))
```
# Color-blind friendly color palette
Source: https://github.com/JLSteenwyk/ggpubfigs
```
library(ggpubfigs)
library(scales)
scales::show_col(ggpubfigs::friendly_pal("muted_nine"))
```

# Center the plot by origin
```
library(ggpmisc)

p1 +
  scale_x_continuous(limits = symmetric_limits) +
  scale_y_continuous(limits = symmetric_limits)
```

# Heatmap with circles

<img src="heatmap_with_circles.png" width=50% height=50%>

```
library(ggplot2)

df <- readRDS("table.heatmap_with_circles.rds")

trait1.list <- unique(df$trait1)
trait2.list <- unique(df$trait2)
trait.list <- unique(c(trait1.list, trait2.list))

df <- rbind(df, df %>% rename(trait1=trait2, trait2=trait1))

df <- rbind(df, data.frame(trait1=trait.list, trait2=trait.list, rg=1, se=NA, Z=NA, P=NA))

df$trait1 <- factor(df$trait1, levels=sort(unique(df$trait1)))
df$trait2 <- factor(df$trait2, levels=rev(sort(unique(df$trait2))))

max_circle_size <- 30 # 저장하는 figure 사이즈에 따라 유동적으로 결정됨. Figure 파일 저장 후, 확인하면서 조절하기.
min_rg <- min(abs(df$rg))
ggplot(df, aes(trait1, trait2)) +
    geom_tile(color = "black", fill = NA) +
    geom_point(aes(size = abs(rg), fill = rg), shape = 21, color = "black") +
    geom_text(aes(label = round(rg, 3)), color = "black", size = 4) +
    scale_size_continuous(range = c(max_circle_size*min_rg, max_circle_size), guide = "none") + # Adjust the size range
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, 
                         breaks=c(0, 0.5, 1), limits=c(0, 1)) +
    theme_light() +
    labs(x = "", y = "", fill = "rg") +
    theme(axis.text.x = element_text(size=10, angle=45, hjust=1, color="black", family = "Helvetica"),
          axis.text.y = element_text(size=10, color="black", family = "Helvetica"),
          axis.ticks = element_blank(),
          legend.title = element_text(color="black", family = "Helvetica"),
          legend.text = element_text(color="black", family = "Helvetica"),
          panel.grid.major = element_blank(),
          panel.border = element_blank()) +
    coord_fixed(ratio = 1) # Keep the aspect ratio fixed
```
