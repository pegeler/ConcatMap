#!/usr/bin/env Rscript
depth_filename <- commandArgs(TRUE)[1]

df <- read.csv(depth_filename, header = FALSE, sep = "\t")
npos <- nrow(df) / 2
depth <- df[rep(c(TRUE, FALSE), each = npos), 3] + df[rep(c(FALSE, TRUE), each = npos), 3]

png_filename <- sub("\\.tsv$", ".png", depth_filename)

png(png_filename, width = 900, height = 700)
plot(depth, ylab = "Number of reads")
invisible(dev.off())
