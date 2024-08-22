
## usage
# Rscript <this script> <input matrix (gzipped or not)> <out pdf> <out newick> <clustering method>

# available clustering method: average|ward.D|ward.D2|single|complete|mcquitty|median|centroid

## library
library(gplots)
library(phylogram)

## args
args <- commandArgs(trailingOnly = TRUE)

fin  = args[1] ## input identity matrix
fout = args[2] ## pdf heatmap output
fdnd = args[3] ## dendrogram newick output
meth = args[4] ## clustering method --> [average|ward.D|ward.D2|single|complete|mcquitty|median|centroid]

## [optional args]
fphy = args[5] ## phylogenetic tree file (newick or nexus) if given
useC = args[6] ## generate chronogram if given

## parse identity matrix (separater: ';')
get_matrix = function(f) {
	### can read txtfile, gzipped file, bz2 file, ...
	## summary(con)$class == 'gzfile'
	con = file(f)
	m = as.matrix(read.delim(con, sep=";", row.names = 1, as.is = T))

	return(m)
}

m = get_matrix(fin)

## check if upper triangle distance or not
if (m[1,2] != m[2,1] & m[2,1] == 0) { ### upper triangle
	m = m + t(m)
}

### min max check
if (min(m) < 0) {
	print("minimum value of matrix is less than 0. The range of value should be 0 to 1. aborting...") ; quit()
}
if (max(m) > 1) {
	print("minimum value of matrix is more than 1. The range of value should be 0 to 1. aborting...") ; quit()
}

## distance --> similarity
m = 1 - m

## ladderize dendrogram by dendextend
if (is.na(fphy)) { ## calculate dendrogram (when phylogenetic tree is not given)
	## additional library
	library(dendextend)

	m.dist = dist(m)
	m.dend = as.dendrogram(hclust(m.dist, method=meth))
	m.dend = rev(dendextend::ladderize(m.dend))

} else { ## parse phylogenetic tree and use as dendrogram
	library(ape)
	library(phangorn)

	lab = tolower(fphy)
	if (length(grep("\\.nex$", lab)) > 1 || length(grep("\\.nxs$", lab)) || length(grep("\\.nexus", lab))) { ## parse extension
		## nexus: .nex, .nxs, .nexus
		tr = ape::read.nexus(fphy)
	} else {
	 	## Newick or New Hampshire
		tr = ape::read.tree(fphy)
	}

	## midpoint rooting and ladderize
	tr = phangorn::midpoint(tr)
	tr = ape::ladderize(tr, right=F)

	## use chronogram, might take a time.
	if (!is.na(useC) && useC == "chronogram") {
		tr = ape::chronos(tr)
	}

	## generate dendrogram using phylogram package
	m.dend = phylogram::as.dendrogram.phylo(tr)


	## [!!!] very complicated steps to arrange the order of heatmap... this is necesary.
	## (ref: https://www.polarmicrobes.org/merging-a-phylogenetic-tree-with-a-heatmap-in-r/)
	c.order    = order.dendrogram(m.dend)
	c.name     = labels(m.dend)
	c.position = data.frame(c.name, c.order)
	c.position = c.position[order(c.position$c.order),]
	new.order  = match(c.position$c.name, row.names(m))
	m          = m[new.order, new.order]
}

## write dendrogram as a newick file
phylogram::write.dendrogram(m.dend, fdnd)

## [TODO] multiple color scale?

## viridis
# n  = 20
# h1 = 300
# c1 = 40
# l1 = 15
# h2 = 75
# c2 = 95
# l2 = 90
# hmcol = hcl(seq(h1, h2, length.out=n), seq(c1, c2, length.out=n), seq(l1, l2, length.out=n))

## rainbow
n = 10
d = 240 / (2 * n - 1) * 2
hmcol = c(hsv(seq(240/360,     (120+d)/360, length.out=n), seq(90/100, 60/100, length.out=n), seq(60/100, 90/100, length.out=n)),
          hsv(seq((120-d)/360,       0/360, length.out=n), seq(60/100, 90/100, length.out=n), seq(90/100, 60/100, length.out=n)))


## pdf size control
# size = ceiling(sqrt(nrow(m))) + 3
size = log10(nrow(m)) ** 2 * 15 + 3
factor = .1+1/log10(nrow(m))

## make output
pdf(file=fout, width=size, height=size, useDingbats=F)
heatmap.2(m, dendrogram="both", col=hmcol, trace="none",
					breaks=seq(0, 1, 0.05),     ## from 0 to 1, inverval: 0.05
					Rowv=m.dend, Colv=m.dend,   ## use dendrogram order
					# cexRow=.5, cexCol=.5,     ## label font size
					cexRow=factor, cexCol=factor,     ## label font size
					offsetRow=0, offsetCol=0,   ## label offset
					adjCol=c(NA, .5),           ## adjust position of column labels
					keysize=factor,
					tracecol="cyan",
					densadj=.5,
					key.title="",
					key.xlab="simka similarity",
					key.ylab="Count",
					)
dev.off()
