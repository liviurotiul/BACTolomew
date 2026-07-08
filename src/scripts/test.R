# 1) Make sure ggplot2 is recent
library(ggplot2)
packageVersion("ggplot2")
# should show 4.0.1 (or at least >= 4.0.0)

# 2) Install remotes if missing
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# 3) Install the latest ggtree from GitHub (NOT BiocManager)
remotes::install_github("YuLab-SMU/ggtree")

# 4) Check versions
packageVersion("ggtree")

# 5) Try loading it
library(ggtree)
sessionInfo()
