# mixopt

Empirical comparisons of algorithms for solving the "mixture
distribution" optimization problem. See the code and results
[here](https://pcarbo.github.io/mixopt).

## How to build the webpages

Run the following commands in R from the analysis directory:

```
library(rmarkdown)
render("index.Rmd",output_dir = "../docs")
```
