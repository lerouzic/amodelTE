#!/usr/bin/env Rscript

source("../src/amodel.R")

pdf("figA.pdf", width=8, height=4)

layout(t(1:2))
plot(mod1.nosel(0.2, 0.03)["n",], xlab="Generations", ylab="Copy number", type="l")
plot(mod1.nosel(0.2, 0.03)["p",], xlab="Generations", ylab="Frequency of regulatory clusters", type="l")

dev.off()
