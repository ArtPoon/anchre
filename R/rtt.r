#!/usr/bin/env Rscript
args <- commandArgs(TRUE)

if (length(args) != 1) {
	stop('Usage: rtt.r <input NWK>')
}

input.nwk <- args[1]
#if (!file.exists(input.nwk)) {
#	stop('input file does not exist')
#}

# load tree object from file
library(ape)
tree <- read.tree(text=input.nwk)


# extract tip dates from labels
tips <- tree$tip.label
tip.dates <- sapply(tips, function(s) {
	items <- strsplit(s, split='_')[[1]]
	return (as.integer(items[length(items)]))
})


#source('chronos.R')

rtt <- function (t, tip.dates, ncpu = 1, objective = "correlation", 
    			 opt.tol = .Machine$double.eps^0.25)  {
    if (objective == "correlation") 
        objective <- function(x, y) cor.test(y, x)$estimate
    else if (objective == "rsquared") 
        objective <- function(x, y) summary(lm(y ~ x))$r.squared
    else if (objective == "rms") 
        objective <- function(x, y) -summary(lm(y ~ x))$sigma^2
    else stop("objective must be one of \"correlation\", \"rsquared\", or \"rms\"")

	# unroot the tree
    ut <- unroot(t)
    
    # pad zero branch lengths
    #padded.brlens <- sapply(t$edge.length, function(x) ifelse(x>0, x, 1e-5))
    #compute.brlen(ut, padded.brlens)
    
    # compute pairwise distances between tips and all nodes
    dist <- dist.nodes(ut)[ ,1:length(tip.dates)]
    #dist <- dist.nodes(ut)[, 1:(ut$Nnode + 2)]
    
    # Save the tip labels, they might get reordered during the re-root
	# which would muss up th ordering of tip.dates
	saved.tips <- t$tip.label
	ut$tip.label <- unlist(lapply(1:length(ut$tip.label), toString)) # Give them 1..n, as strings

    f <- function (x, parent, child) {
        edge.dist <- x * dist[parent, ] + (1 - x) * dist[child,]
        objective(tip.dates, edge.dist)
    }

    obj.edge <- if (ncpu > 1) 
        unlist(parallel::mclapply(1:nrow(ut$edge), function (e) {
            opt.fun <- function (x) f(x, ut$edge[e,1], ut$edge[e,2])
            optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
        }, mc.cores=ncpu))
    else apply(ut$edge, 1, function (e) {
        opt.fun <- function (x) f(x, e[1], e[2])
        optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$objective
    })

    best.edge <- which.max(obj.edge)

    best.edge.parent <- ut$edge[best.edge, 1]
    best.edge.child <- ut$edge[best.edge, 2]
    best.edge.length <- ut$edge.length[best.edge]

    opt.fun <- function (x) f(x, best.edge.parent, best.edge.child)
    best.pos <- optimize(opt.fun, c(0, 1), maximum = TRUE, tol = opt.tol)$maximum

    new.root <- list(edge = matrix(c(2L, 1L), 1, 2), tip.label = "new.root", 
        edge.length = 1, Nnode = 1L, root.edge = 1)
    class(new.root) <- "phylo"
    ut <- bind.tree(ut, new.root, where = best.edge.child, position = best.pos * 
        best.edge.length)
    ut <- collapse.singles(ut)
    ut <- root(ut, "new.root")
    rt <- drop.tip(ut, "new.root")
    
    # Reorder the tip dates, and restore the labels
  	permutation <- as.integer(rt$tip.label)
  	tip.dates <- tip.dates[permutation]
	rt$tip.label <- saved.tips[permutation]
	
	# Get inferredtime at root
    distances <- node.depth.edgelength(rt)  # distance from root to node
    tip.dists <- distances[1:length(tip.dates)]
    model <- lm(tip.dists ~ tip.dates)
    a <- model$coefficients[1]
    b <- model$coefficients[2]
    root.time <- as.double(-a/b)  # x-intercept
    
    ## use chronos() to rescale node heights
    ## FIXME: this is unstable
    # tip.times <- data.frame(index=1:length(tip.dates), name=rt$tip.label, time=tip.dates, row.names=1:length(tip.dates))
    # max.time <- max(tip.times$time)
    # min.tip.time <- min(tip.times$time)
    
    # calib <- makeChronosCalib(rt)
    # calib$age.min <- min.tip.time - root.time
    # calib$age.max <- max.time - root.time
    # calib <- rbind(calib, data.frame(node=tip.times$index,
    	# age.min=min.tip.time - tip.times$time,
    	# age.max=max.time - tip.times$time,
    	# soft.bounds=FALSE))
    # dated.tree <- RLchronos(rt, lambda=1, 
    	# model='discrete', 
    	# calibration=calib, 
    	# control=chronos.control(nb.rate.cat=1), quiet=TRUE)
    
    # package results as list object to return
    res <- {}
    res$tree <- rt
    res$origin <- root.time
    #res$dated.tree <- dated.tree
    
    return(res)
}

res <- rtt(tree, tip.dates)

write.tree(res$tree)  # to stdout
#write.tree(res$dated.tree)
