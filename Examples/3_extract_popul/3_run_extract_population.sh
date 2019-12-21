
## measure the population of each cluster, definited by the
## x/y-axis ranges

../../3_rama_extract_popul.py                          \
  -in ../2_md_heatmap/fgf21-wt.A208.rama.out.bz2       \
  -list fgf21-wt.A208.rama.cluster.list >              \
        fgf21-wt.A208.rama.cluster.txt


## fgf21-wt.A208.rama.cluster.list
#  list of x/y-axis ranges to count population
#
## fgf21-wt.A208.rama.cluster.txt
#  result of the populations. Best view with 'cat'
#
## fgf21-wt.A208.rama_histo.png
#  original A208 distribution heatmap
#
## fgf21-wt.A208.rama_histo.popul.png
#  modified (in PowerPoint) heatmap to indicate the population
#  in each cluster, defined by the boxes
#
#
#  19.12.20
