# Keep Me Around (`kma`): Intron Retention Detection

`kma` is an R package that performs intron retention estimation and detection
using biological replicates and resampling. Updated code can always be found at
https://github.com/pachterlab/kma

# Installation

To install, first ensure you have the required packages:

```{r}
required_packages <- c("devtools", "data.table", "reshape2", "dplyr")
install.packages(required_packages)
```

You can then install the package using `devtools`:

```{r}
devtools::install_github("pachterlab/kma")
```

Assuming all goes well, load `kma`:

```{r}
library("kma")
```

# Tutorial

After it has been installed, please see the vignette in R:

```{r}
vignette("kma")
```

# Bugs and feature requests

Please file these on [Github](https://github.com/pachterlab/kma/issues).

# Future work

- Additional exploratory analysis plotting tools
- Provide differential intron usage analysis between experimental conditions
    - We currently have some ideas on how to do this and will likely be
      implementing it soon
- Provide time series analysis

# Authors

Software was developed by Harold Pimentel. Methods were developed with Lior
Pachter and John Conboy.

# Related open source tools

Below you will find a list of related tools and how they differ from `kma`.

### DEXSeq

[DEXSeq](http://bioconductor.org/packages/release/bioc/html/DEXSeq.html) is interested in differential usage across genic regions. As a result,
it does not determine whether an intron is being "used" (relative to transript
expression), simply that it is being "differentially used."

### MISO

[MISO](http://genes.mit.edu/burgelab/miso/) can calculate the intronic percent spliced in (PSI), though it currently
requires a modified annotation [from their website](https://miso.readthedocs.org/en/fastmiso/annotation.html). `kma` can currently work with any
annotation, as the annotation will be processed during the pre-processing step.
Also, MISO does not currently provide [built-in suppport for
replicates](http://miso.readthedocs.org/en/fastmiso/#answer13).
