### use three modified function below ----------------------------------
# chromosome visualization by Genomic Region
vis("chr17",7661779,7687538, all, platform='EPIC', genome = "hg38",
show.probeNames = TRUE, draw = T)
# chromosome visualization by Gene name
visgene("CDKN1B", all, platform = "EPIC", genome = 'hg38')+
  WLegendV("betas", TopRightOf("betas"))
# chromosome visualization by Probe ID 
visProbes("cg03487391", all, platform='EPIC', genome = "hg38") +
  WLegendV("betas", TopRightOf("betas"))

# function code block ----------------------------------------------------
ass <- function(
    betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband,
    heat.height = NULL, mapLine.height = 0.2,
    show.probeNames = TRUE, show.samples.n = NULL,
    show.sampleNames = TRUE, sample.name.fontsize = 10,
    dmin = 0, dmax = 1) {
  
  if (is.null(show.samples.n)) { show.samples.n <- ncol(betas); }
  if (is.null(heat.height) && length(txns) > 0) {
    heat.height <- 10 / length(txns); }
  w <- WGrob(plt.txns, name = 'txn')
  w <- w + WGrob(plt.mapLines, Beneath(pad=0, height=mapLine.height))
  w <- w + WHeatmap(
    t(betas), Beneath(height = heat.height),
    name = 'betas',
    cmp = CMPar(dmin=dmin, dmax=dmax,colorspace.name = "diverge_hsv", 
                stop.points=c("blue","yellow")),
    xticklabels = show.probeNames,
    xticklabel.rotat = 45,
    yticklabels = show.sampleNames,
    yticklabel.fontsize = sample.name.fontsize,
    yticklabels.n = show.samples.n,
    xticklabels.n = length(probes))
  w <- w + WGrob(plt.cytoband, TopOf('txn', height=0.15))
  w
}

vis <- function(chrm, beg, end, betas, platform = NULL,
                genome = NULL, draw = TRUE, cluster.samples = FALSE,
                na.rm = FALSE, nprobes.max = 1000, txn.types = "protein_coding",
                txn.font.size = 6, ...) {
  
  if (is.null(dim(betas))) { betas <- as.matrix(betas) }
  platform <- sesameData_check_platform(platform, rownames(betas))
  genome <- sesameData_check_genome(genome, platform)
  
  reg <- GRanges(chrm, IRanges::IRanges(beg, end))
  
  genomeInfo <- sesameData_getGenomeInfo(genome)
  txns <- subsetByOverlaps(genomeInfo$txns, reg)
  
  probes <- sesameData_getManifestGRanges(platform, genome = genome)
  probes <- subsetByOverlaps(probes, reg)
  probes <- probes[names(probes) %in% rownames(betas)]
  if (na.rm) { probes <- probes[apply(
    betas[names(probes), ], 1, function(x) !all(is.na(x)))] }
  
  if (length(probes) == 0) {
    stop("No probe overlap region ", sprintf(
      '%s:%d-%d', chrm, beg, end)) }
  if (length(probes) > nprobes.max) {
    stop(sprintf('Too many probes (%d). Shrink region?', length(probes))) }
  
  ## plot transcripts
  plt.txns <- plotTranscripts(txns, reg, beg, end,
                              txn.types = txn.types, txn.font.size = txn.font.size)
  plt.mapLines <- plotMapLines(probes, beg, end)
  plt.cytoband <- plotCytoBand(chrm, beg, end, genomeInfo)
  
  ## clustering
  betas <- betas[names(probes),,drop=FALSE]
  if (cluster.samples) {
    betas <- column.cluster(betas[names(probes),,drop=FALSE])$mat }
  
  if (draw) { ass(betas, txns, probes,
                  plt.txns, plt.mapLines, plt.cytoband, ...)
  } else { return(betas); }
}

visProbes <- function(
    probeNames, betas,
    platform = NULL, genome = NULL,
    upstream = 1000, dwstream = 1000, ...) {
  
  if (is.null(dim(betas))) { betas <- as.matrix(betas); }
  platform <- sesameData_check_platform(platform, rownames(betas))
  genome <- sesameData_check_genome(genome, platform)
  
  probes <- sesameData_getManifestGRanges(platform, genome)
  probeNames <- probeNames[probeNames %in% names(probes)]
  
  if (length(probeNames)==0)
    stop('Probes specified are not well mapped.')
  
  target.probes <- probes[probeNames]
  
  regBeg <- min(GenomicRanges::start(target.probes)) - upstream
  regEnd <- max(GenomicRanges::end(target.probes)) + dwstream
  
  vis(
    as.character(GenomicRanges::seqnames(
      target.probes[1])), regBeg, regEnd,
    betas, platform = platform, genome = genome, ...)
}

# visualizeGene
visgene <- function(gene_name, betas,
                    platform = NULL, genome = NULL,
                    upstream = 2000, dwstream = 2000, ...) {
  
  if (is.null(dim(betas))) { betas <- as.matrix(betas); }
  platform <- sesameData_check_platform(platform, rownames(betas))
  genome <- sesameData_check_genome(genome, platform)
  
  txns <- sesameData_getGenomeInfo(genome)$txns
  target.txns <- txns[GenomicRanges::mcols(txns)$gene_name == gene_name]
  stopifnot(length(target.txns) > 0)
  target.strand <- as.character(GenomicRanges::strand(target.txns[[1]][1]))
  if (target.strand == '+') {
    pad.start <- upstream
    pad.end <- dwstream
  } else {
    pad.start <- dwstream
    pad.end <- upstream
  }
  
  merged.exons <- GenomicRanges::reduce(unlist(target.txns))
  vis(
    as.character(GenomicRanges::seqnames(merged.exons[1])),
    min(GenomicRanges::start(merged.exons)) - pad.start,
    max(GenomicRanges::end(merged.exons)) + pad.end,
    betas, platform = platform, genome = genome, ...)
}

exonToCDS <- function(exons, cdsStart, cdsEnd) {
    if (is.na(cdsStart) || is.na(cdsEnd) || cdsEnd <= cdsStart) {
        return(NULL); } 
    cds <- exons[(
        (GenomicRanges::start(exons) < cdsEnd) &
            (GenomicRanges::end(exons) > cdsStart))]
    GenomicRanges::start(cds) <- pmax(
        GenomicRanges::start(cds), cdsStart)
    GenomicRanges::end(cds) <- pmin(
        GenomicRanges::end(cds), cdsEnd)
    cds
}

plotTranscript1 <- function(txn, reg, i, beg, end,
    isoformHeight, padHeight, txn.font.size) {

    txn_name <- names(txn)[1]; exons <- txn[[1]]
    meta <- as.data.frame(GenomicRanges::mcols(txn))
    plt.width <- end - beg
    txn.beg <- max(beg, min(GenomicRanges::start(exons))-2000)
    txn.end <- min(end, max(GenomicRanges::end(exons))+2000)
    exons <- subsetByOverlaps(exons, reg)
    txn.strand <- as.character(GenomicRanges::strand(exons[1]))
    lined <- (c(txn.beg, txn.end)-beg) / plt.width # direction is in arrow ends
    
    y.bot <- (i-1) * isoformHeight + padHeight
    y.bot.exon <- y.bot + padHeight
    y.hei <- isoformHeight - 2 * padHeight

    ## transcript name
    g <- gList(grid.text(sprintf('%s (%s)', meta$gene_name, txn_name),
        x=mean(lined), y=y.bot + y.hei + padHeight * 0.5,
        just=c('center','bottom'),
        gp = gpar(fontsize = txn.font.size), draw=FALSE))

    ## plot transcript line
    g <- gList(g, gList(grid.lines(
        x=lined, y=y.bot+y.hei/2, arrow=arrow(length=unit(0.06, "inches"),
            ends=ifelse(txn.strand == "+", "last", "first")), draw=FALSE)))

    g <- gList(g, gList(grid.lines(x=c(0,1), y=y.bot+y.hei/2,
        gp=gpar(lty='dotted'), draw=FALSE)))

    ## plot exons
    g <- gList(g, gList(
        grid.rect((GenomicRanges::start(exons)-beg)/plt.width,
            y.bot + y.hei/2 - y.hei/3, GenomicRanges::width(exons)/plt.width,
            y.hei/3*2, gp=gpar(fill='grey10', lwd=0),
            just=c('left','bottom'), draw=FALSE)))

    ## plot cds
    cds <- exonToCDS(exons, as.integer(meta$cdsStart), as.integer(meta$cdsEnd))
    if (length(cds) > 0) {
        g <- gList(g, gList(
            grid.rect((GenomicRanges::start(cds)-beg)/plt.width,
                y.bot + y.hei/2 - y.hei/6, GenomicRanges::width(cds)/plt.width,
                y.hei/6*2, gp=gpar(fill='red', lwd=0),
                just=c('left','bottom'), draw=FALSE)))
    }
    g
}

## helper function to plot transcript
plotTranscripts <- function(
    txns, reg, beg, end,
    txn.types = c("protein_coding"), txn.font.size = 6) {

    if (!is.null(txn.types)) {
        txns <- txns[
            GenomicRanges::mcols(txns)$transcript_type %in% txn.types] }
    
    if (length(txns) == 0) {
        return(gList(
            grid.rect(0,0.1,1,0.8, just = c('left','bottom'), draw=FALSE),
            grid.text('No transcript found', x=0.5, y=0.5, draw=FALSE)))
    }
    
    isoformHeight <- 1/length(txns)
    padHeight <- isoformHeight*0.2

    do.call(gList, lapply(seq_along(txns), function(i) {
        plotTranscript1(txns[i], reg, i, beg, end,
            isoformHeight, padHeight, txn.font.size)
    }))
}

plotMapLines <- function(probes, beg, end) {
    nprobes <- length(probes)
    x00 <- ((GenomicRanges::start(probes) - beg) / (end - beg))
    y0 <- rep(0.5, length.out=length(probes))
    x1 <- ((seq_len(nprobes) - 0.5)/nprobes)
    y1 <- rep(0, length.out=nprobes)
    x0 <- c(x00, x00)
    x1 <- c(x1, x00)
    y0 <- c(y0, rep(0.5, length.out=length(probes)))
    y1 <- c(y1, rep(1, length.out=length(probes)))
    grid.segments(x0, y0, x1, y1, draw=FALSE)
}

plotCytoBand <- function(
    chrom, beg, end, genomeInfo) {

    cytoBand <- genomeInfo$cytoBand

    ## set cytoband color
    requireNamespace("pals")
    cytoBand2col <- setNames(
        pals::ocean.gray(10)[seq(9,3)],
        c('stalk', 'gneg', 'gpos25', 'gpos50', 'gpos75', 'gpos100'))
    cytoBand2col['acen'] <- 'red'
    cytoBand2col['gvar'] <- cytoBand2col['gpos75']

    ## chromosome range
    cytoBand.target <- cytoBand[cytoBand$chrom == chrom,]
    chromEnd <- max(cytoBand.target$chromEnd)
    chromBeg <- min(cytoBand.target$chromStart)
    chromWid <- chromEnd - chromBeg
    bandColor <- cytoBand2col[as.character(cytoBand.target$gieStain)]

    pltx0 <- (c(beg, end)-chromBeg)/chromWid
    gList(
        grid.text( # coordinate name
            sprintf("%s:%d-%d", chrom, beg, end), 0, 0.9,
            just = c('left','bottom'), draw = FALSE),
        ## cytoband box
        grid.rect(0, 0.35, 1, 0.35, just = c("left", "bottom"),
            gp = gpar(col = "black", lwd=2, lty="solid"), draw = FALSE),
        grid.rect( # cytoband
            vapply(cytoBand.target$chromStart,
                function(x) (x-chromBeg)/chromWid, 1),
            0.35,
            (cytoBand.target$chromEnd - cytoBand.target$chromStart)/chromWid,
            0.35, gp = gpar(fill = bandColor, col = bandColor),
            just = c('left','bottom'), draw = FALSE),
        grid.segments( # sentinel bar
            x0 = pltx0, y0 = 0.1, x1 = pltx0, y1 = 0.9,
            gp = gpar(col = "red"), draw = FALSE))
}

#' assemble plots
#'
#' @param betas beta value
#' @param txns transcripts GRanges
#' @param probes probe GRanges
#' @param plt.txns transcripts plot objects
#' @param plt.mapLines map line plot objects
#' @param plt.cytoband cytoband plot objects
#' @param heat.height heatmap height (auto inferred based on rows)
#' @param mapLine.height height of the map lines
#' @param show.probeNames whether to show probe names
#' @param show.samples.n number of samples to show (default: all)
#' @param show.sampleNames whether to show sample names
#' @param sample.name.fontsize sample name font size
#' @param dmin data min
#' @param dmax data max
#' @return a grid object
assemble_plots <- function(
    betas, txns, probes, plt.txns, plt.mapLines, plt.cytoband,
    heat.height = NULL, mapLine.height = 0.2,
    show.probeNames = TRUE, show.samples.n = NULL,
    show.sampleNames = TRUE, sample.name.fontsize = 10,
    dmin = 0, dmax = 1) {
    
    if (is.null(show.samples.n)) { show.samples.n <- ncol(betas); }
    if (is.null(heat.height) && length(txns) > 0) {
        heat.height <- 10 / length(txns); }
    w <- WGrob(plt.txns, name = 'txn')
    w <- w + WGrob(plt.mapLines, Beneath(pad=0, height=mapLine.height))
    w <- w + WHeatmap(
        t(betas), Beneath(height = heat.height),
        name = 'betas',
        cmp = CMPar(dmin=dmin, dmax=dmax, label2color = col_f),
        xticklabels = show.probeNames,
        xticklabel.rotat = 45,
        yticklabels = show.sampleNames,
        yticklabel.fontsize = sample.name.fontsize,
        yticklabels.n = show.samples.n,
        xticklabels.n = length(probes))
    w <- w + WGrob(plt.cytoband, TopOf('txn', height=0.15))
    w
}

MapToContinuousColors <- function(data, cmp=CMPar(), given.cm=NULL) {
  
  if (!is.null(given.cm)) {
    given.cm$colors <- apply(
      given.cm$mapper(given.cm$scaler(data)), 1,
      function(x) do.call(rgb, c(as.list(x), maxColorValue=255)))
    return(given.cm)
  }
  
  if (is.null(cmp$stop.points)) {
    if (is.null(cmp$cmap) &&
        is.null(cmp$brewer.name) &&
        is.null(cmp$colorspace.name)) {
      cmp$cmap <- 'jet'
    }
    if (!is.null(cmp$cmap)) {
      ## get(cmp$cmap)
      ## data(cmp$cmap)
      cmp$stop.points <- get(paste0(cmp$cmap,'.stops'))
    } else if (!is.null(cmp$brewer.name)) {
      ## use display.brewer.all for the brewer colors
      ## note that brewer.n cannot be >8 typically
      if (cmp$brewer.n < 3)
        cmp$brewer.n <- 3
      cmp$stop.points <- brewer.pal(cmp$brewer.n, cmp$brewer.name)
      if (cmp$rev)
        cmp$stop.points <- rev(cmp$stop.points)
    } else {
      ## colorspace.name can be
      ## diverge_hcl, diverge_hsv, terrain_hcl, heat_hcl, sequential_hcl and rainbow_hcl
      ## colorspace.n can be very large
      cmp$stop.points <- get(cmp$colorspace.name)(cmp$colorspace.n)
    }
  }
  
  ## cap data
  if (!is.null(cmp$dmax))
    data[data>=cmp$dmax] <- cmp$dmax
  if (!is.null(cmp$dmin))
    data[data<=cmp$dmin] <- cmp$dmin
  
  .dmax <- max(cmp$dmax, data, na.rm=TRUE)
  .dmin <- min(cmp$dmin, data, na.rm=TRUE)
  if (.dmax==.dmin) # when range==0
    .dmax <- .dmax+1
  data <- (data - .dmin) / (.dmax-.dmin)
  
  cm <- ColorMap(
    dmin = .dmin, dmax = .dmax,
    scaler = function(x) {(x-.dmin)/(.dmax-.dmin)},
    mapper = colorRamp(cmp$stop.points, alpha=TRUE))
  cm$colors = apply(cm$mapper(data), 1, function(x) {
    if (any(is.na(x)))
      x <- col2rgb(cmp$na.color, alpha=TRUE)
    do.call(rgb, c(as.list(x), maxColorValue=255))
  })
  cm
}



