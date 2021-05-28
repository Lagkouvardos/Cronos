# install package rstudioapi and set file path as current directory
if (!require("rstudioapi")) { install.packages("rstudioapi", dependencies = T) }
library("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# install missing packages, load packages, sources
source('environment.r')
environment.start()

# load otu table and tree
data.folder <- "data/otus/infants/"
workload <- otus.load(data.folder, "OTUs-Table.tab", "OTUs-NJTree.tre")

# load metadata
metadata <- as.matrix(read_excel('data/meta.xlsx'))
metadata[, 'ID'] <- str_pad(str_trim(metadata[, 'ID']), 3, pad = "0")

# extract individuals, timepoints
samples <- row.names(workload$otus)
individuals <- unlist(lapply(samples, function(i){substr(i, start = 2, stop = 4)}))
timepoints <- unlist(lapply(samples, function(i){substr(i, start = 5, stop = 7)}))

# build unique and ordered timepoint list
unique.timepoints <- as.vector(sort(unique(timepoints)))

# cluster data per timepoint
clusterings <- time.series.generate.clusterings(workload$otus, workload$otus.tree)

# generate time series from computed clusters of every timepoint 
time.series <- time.series.generate(clusterings, samples)
#time.series <- cbind( metadata[, 'Group'], time.series)
#colnames(time.series)[1] <- '00'

write.table(time.series, file = paste(data.folder, 'time-series.txt', sep=''))

# generate transition matrix for selected time points
timepoint.chain <- c("01", "03")
timepoint.chain <- c("01", "03", "05", "07", "09", "12", "24", "MM")
timepoint.chain <- c("01", "03", "05", "07", "09", "12", "24")
transition.matrix <- markov.transition.matrix(time.series, timepoint.chain)
plots.network(transition.matrix, timepoint.chain)

timepoint.chain <- c("01","05")
transition.matrix <- markov.transition.matrix(time.series, timepoint.chain)
plots.network(transition.matrix, timepoint.chain)

timepoint.chain <- c("01","07")
transition.matrix <- markov.transition.matrix(time.series, timepoint.chain)
plots.network(transition.matrix, timepoint.chain)

timepoint.chain <- c("01","09")
transition.matrix <- markov.transition.matrix(time.series, timepoint.chain)
plots.network(transition.matrix, timepoint.chain)

timepoint.chain <- c("03","05")
transition.matrix <- markov.transition.matrix(time.series, timepoint.chain)
plots.network(transition.matrix, timepoint.chain)

timepoint.chain <- c("03","07")
transition.matrix <- markov.transition.matrix(time.series, timepoint.chain)
plots.network(transition.matrix, timepoint.chain)

g <- graph_from_adjacency_matrix(transition.matrix, weighted = T, mode = 'directed')
bn <- as.bn(g)
plot(bn)
#heatmap(transition.matrix, Rowv = NA, Colv = NA, scale = 'none', symm = T)

# generate svm classifiers for every timepoint
svm.classifiers <- svm.generate.classifiers(workload$otus, clusterings) 
# rf.classifiers <- rf.generate.classifiers(workload$otus, clusterings) 

