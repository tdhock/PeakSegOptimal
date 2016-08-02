\name{PeakSegPDPAchrom}
\alias{PeakSegPDPAchrom}
\title{PeakSegPDPAchrom}
\description{Find the optimal change-points using the Poisson loss and the
PeakSeg constraint. This function is a user-friendly interface to
the PeakSegPDPA function.}
\usage{PeakSegPDPAchrom(count.df, max.peaks = NULL)}
\arguments{
  \item{count.df}{data.frame with columns count, chromStart, chromEnd.}
  \item{max.peaks}{integer > 0: maximum number of peaks.}
}

\value{List of data.frames: segments can be used for plotting the
segmentation model, loss describes model loss and feasibility,
modelSelection.feasible describes the set of all linear penalty
(lambda) values which can be used to select the feasible models,
modelSelection.decreasing selects from all models that decrease
the Poisson loss relative to simpler models (same as PeakSegFPOP).}

\author{Toby Dylan Hocking}




\examples{

## samples for which pdpa recovers a more likely model, but it is
## not feasible for the PeakSeg problem (some segment means are
## equal).
sample.id <- "McGill0322"
sample.id <- "McGill0079"
sample.id <- "McGill0106"
n.peaks <- 3
data("H3K4me3_XJ_immune_chunk1")
H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
one.sample <- by.sample[[sample.id]]
pdpa.fit <- PeakSegPDPAchrom(one.sample, 9L)
pdpa.segs <- subset(pdpa.fit$segments, n.peaks == peaks)
both.segs.list <- list(pdpa=data.frame(pdpa.segs, algorithm="PDPA"))
if(require(PeakSegDP)){
  dp.fit <- PeakSegDP(one.sample, 9L)
  dp.segs <- subset(dp.fit$segments, n.peaks == peaks)
  both.segs.list$dp <- data.frame(dp.segs, algorithm="cDPA")
}
both.segs <- do.call(rbind, both.segs.list)
both.breaks <- subset(both.segs, 1 < first)
library(ggplot2)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  facet_grid(algorithm ~ ., scales="free")+
  geom_step(aes(chromStart/1e3, coverage),
            data=one.sample, color="grey")+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               color="green",
               data=both.segs)+
  geom_vline(aes(xintercept=chromStart/1e3),
             color="green",
             linetype="dashed",
             data=both.breaks)

## samples for which pdpa recovers some feasible models that the
## heuristic dp does not.
sample.id.vec <- c(
  "McGill0091", "McGill0107", "McGill0095",
  "McGill0059", "McGill0029", "McGill0010")
sample.id <- sample.id.vec[1]
one.sample <- by.sample[[sample.id]]
pdpa.fit <- PeakSegPDPAchrom(one.sample, 9L)
gg.loss <- ggplot()+
  scale_size_manual(values=c(cDPA=2, PDPA=3))+
  scale_fill_manual(values=c(cDPA="white", PDPA="black"))+
  geom_point(aes(peaks, PoissonLoss,
                 size=algorithm, fill=algorithm),
             shape=21,
             data=data.frame(pdpa.fit$loss, algorithm="PDPA"))
if(require(PeakSegDP)){
  dp.fit <- PeakSegDP(one.sample, 9L)
  gg.loss <- gg.loss+
    geom_point(aes(peaks, error,
                   size=algorithm, fill=algorithm),
               shape=21,
               data=data.frame(dp.fit$error, algorithm="cDPA"))
}
gg.loss

}
