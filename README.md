# Keep Me Around (`kma`): Intron Retention Detection

`kma` is an R package that performs intron retention estimation and detection
using biological replicates and resampling.

# Installation

To install, first ensure you have the required packages:

```{r}
required_packages <- c("devtools", "data.table", "reshape2", "dplyr")
install.packages(required_packages)
library("devtools")
install_github("pachterlab/kma")
```

Assuming all goes well, load `kma`:

```{r}
library("kma")
```

# Tutorial

After it has been installed, please see the vignette in R:

```{r}
vignette("kma", "kma")
```
