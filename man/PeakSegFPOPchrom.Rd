\name{PeakSegFPOPchrom}
\alias{PeakSegFPOPchrom}
\title{PeakSegFPOPchrom}
\description{Find the optimal change-points using the Poisson loss and the
PeakSeg constraint. This function is a user-friendly interface to
the PeakSegFPOP function.}
\usage{PeakSegFPOPchrom(count.df, penalty = NULL)}
\arguments{
  \item{count.df}{data.frame with columns count, chromStart, chromEnd.}
  \item{penalty}{integer > 0: maximum number of peaks.}
}

\value{List of data.frames: segments can be used for plotting the
segmentation model, loss summarizes the penalized PoissonLoss and
feasibilty of the computed model.}

\author{Toby Dylan Hocking}




\examples{

library(coseg)
data("H3K4me3_XJ_immune_chunk1")
sample.id <- "McGill0106"
H3K4me3_XJ_immune_chunk1$count <- H3K4me3_XJ_immune_chunk1$coverage
by.sample <-
  split(H3K4me3_XJ_immune_chunk1, H3K4me3_XJ_immune_chunk1$sample.id)
one.sample <- by.sample[[sample.id]]
penalty.constant <- 2000
fpop.fit <- PeakSegFPOPchrom(one.sample, penalty.constant)
fpop.breaks <- subset(fpop.fit$segments, 1 < first)

library(ggplot2)
ggplot()+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "lines"))+
  geom_step(aes(chromStart/1e3, coverage),
            data=one.sample, color="grey")+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               color="green",
               data=fpop.fit$segments)+
  geom_vline(aes(xintercept=chromStart/1e3),
             color="green",
             linetype="dashed",
             data=fpop.breaks)

max.peaks <- as.integer(fpop.fit$segments$peaks[1]+1)
pdpa.fit <- PeakSegPDPAchrom(one.sample, max.peaks)
models <- pdpa.fit$modelSelection
models$PoissonLoss <- rev(pdpa.fit$loss$PoissonLoss)
models$algorithm <- "PDPA"
fpop.fit$loss$algorithm <- "FPOP"
ggplot()+
  geom_abline(aes(slope=segments-1, intercept=PoissonLoss, color=peaks),
              data=pdpa.fit$loss)+
  geom_point(aes(penalty.constant, penalized.loss, fill=algorithm),
             shape=21,
             data=fpop.fit$loss)+
  geom_point(aes(min.lambda, min.lambda*model.complexity + PoissonLoss,
                 fill=algorithm),
             shape=21,
             data=models)

}
