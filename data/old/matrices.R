#!/usr/bin/env Rscript

# This is a strange way to store a dataset, but it does the job.
# sourcing this file provides six 2x2 matrices (3 matrices for each population).
# G is the G matrix
# P is the P matrix
# Gv is the estimation variance of each element of G

G.tovar <- matrix(c(1.04, 0.57, 0.57, 1.40)/100, ncol=2)
Gv.tovar <- matrix(c(0.0376,0.0219, 0.0219, 0.0337)/10000, ncol=2)
P.tovar <- matrix(c(2.45, 1.18, 1.18, 2.8)/100, ncol=2)
G.tulum <- matrix(c(0.72, 0.48, 0.48, 0.83)/100, ncol=2)
Gv.tulum <- matrix(c(0.0427, 0.0190, 0.0190, 0.0183)/10000, ncol=2)
P.tulum <- matrix(c(4.49,1.47,1.47,2.52)/100, ncol=2)
