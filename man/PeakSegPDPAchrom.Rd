\name{PeakSegPDPAchrom}
\alias{PeakSegPDPAchrom}
\title{PeakSegPDPAchrom}
\description{Compute the PeakSeg constrained, Poisson loss, Segment Neighborhood
model using a constrained version of the Pruned Dynamic
Programming Algorithm.}
\usage{PeakSegPDPAchrom(count.df, max.peaks = NULL)}
\arguments{
  \item{count.df}{data.frame with columns count, chromStart, chromEnd.}
  \item{max.peaks}{integer > 0: maximum number of peaks.}
}

\value{List of data.frames. }

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
dp.fit <- PeakSegDP(one.sample, 9L)
dp.segs <- subset(dp.fit$segments, n.peaks == peaks)
both.segs <- rbind(
  data.frame(dp.segs, algorithm="cDPA"),
  data.frame(pdpa.segs, algorithm="PDPA"))
both.breaks <- subset(both.segs, 1 < first)

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
dp.fit <- PeakSegDP(one.sample, 9L)

ggplot()+
  scale_size_manual(values=c(cDPA=3, PDPA=1))+
  geom_point(aes(peaks, error,
                 size=algorithm, color=algorithm),
             data=data.frame(dp.fit$error, algorithm="cDPA"))+
  geom_point(aes(peaks, PoissonLoss,
                 size=algorithm, color=algorithm),
             data=data.frame(pdpa.fit$loss, algorithm="PDPA"))

}
