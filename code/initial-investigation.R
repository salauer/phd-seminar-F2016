## Flu data investigation by Stephen Lauer and Sara Nunez
## See: https://predict.phiresearchlab.org/legacy/flu/targets.html and https://github.com/cdcepi/flu-tools-2015-16

library(cdcfluview)
library(ggplot2)
library(readr)

usflu <- get_flu_data("national", "ilinet", years=1997:2010)
write.csv(usflu, "data/usflu.csv")
usflu <- read_csv("data/usflu.csv", na="X")

qplot(data=usflu, x=YEAR+(WEEK-1)/53,y=X.UNWEIGHTED.ILI, geom="line")
