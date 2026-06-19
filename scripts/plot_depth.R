#!/usr/bin/env Rscript
depth_filename <- commandArgs(TRUE)[1]

d <- read.delim(depth_filename, header = FALSE)[[3]]
npos <- length(d) / 2
depth <- d[rep(c(TRUE, FALSE), each = npos)] + d[rep(c(FALSE, TRUE), each = npos)]

png_filename <- sub("\\.tsv$", ".png", depth_filename)

png(png_filename, width = 900, height = 700)
plot(depth, ylab = "Number of reads")
invisible(dev.off())
