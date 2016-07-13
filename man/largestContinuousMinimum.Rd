\name{largestContinuousMinimum}
\alias{largestContinuousMinimum}
\title{largestContinuousMinimum}
\description{Find the run of minimum cost with the largest size.}
\usage{largestContinuousMinimum(cost, size)}
\arguments{
  \item{cost}{
}
  \item{size}{
}
}



\author{Toby Dylan Hocking}




\examples{
data(chr11ChIPseq)
one <- subset(chr11ChIPseq$coverage, sample.id=="McGill0322")
fit <- PeakSegDP(one, 5L)
## First compute the optimal number of peaks function.
exact.df <- with(fit$error, exactModelSelection(error, segments, peaks))
## Then compute the PeakError of these models with respect to the
## annotated regions.
regions <- subset(chr11ChIPseq$regions, sample.id=="McGill0322")
peak.list <- split(fit$segments, fit$segments$peaks)
require(PeakError)
all.error <- NULL
error.regions <- NULL
for(peaks.chr in names(peak.list)){
  peaks <- subset(peak.list[[peaks.chr]], status=="peak")
  error <- PeakErrorChrom(peaks, regions)
  error.regions <- rbind(error.regions, data.frame(peaks=peaks.chr, error))
  all.error <- rbind(all.error, {
    data.frame(peaks=as.integer(peaks.chr),
               errors=with(error, sum(fp+fn)))
  })
}
## plot the annotation error of the 6 models.
ann.colors <- c(noPeaks = "#f6f4bf", peakStart = "#ffafaf", 
                peakEnd = "#ff4c4c", peaks = "#a445ee")
library(ggplot2)
ggplot()+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    fill=annotation),
                data=error.regions, alpha=1/2)+
  geom_step(aes(chromStart/1e3, count), data=one, color="grey40")+
  geom_tallrect(aes(xmin=chromStart/1e3, xmax=chromEnd/1e3,
                    linetype=status),
                data=error.regions, fill=NA, color="black", size=2)+
  scale_fill_manual(values=ann.colors)+
  scale_linetype_manual(values=c("false negative"=3,
                        "false positive"=1,
                        "correct"=0))+
  geom_segment(aes(chromStart/1e3, mean,
                   xend=chromEnd/1e3, yend=mean),
               data=fit$segments, color="green")+
  geom_segment(aes(chromStart/1e3, 0,
                   xend=chromEnd/1e3, yend=0),
               data=subset(fit$segments, status=="peak"),
               size=3, color="deepskyblue")+
  theme_bw()+
  theme(panel.margin=grid::unit(0, "cm"))+
  facet_grid(peaks ~ ., scales="free", labeller=function(var, val){
    s <- ifelse(val==1, "", "s")
    paste0(val, " peak", s)
  })  
rownames(all.error) <- all.error$peaks
exact.df$errors <-
  all.error[as.character(exact.df$model.complexity), "errors"]
indices <- with(exact.df, {
  largestContinuousMinimum(errors, max.log.lambda-min.log.lambda)
})
## The target interval (min.log.lambda, max.log.lambda) is the
## largest continuous interval such that the error is minimal.
target.interval <- with(indices, {
  data.frame(min.log.lambda=exact.df$min.log.lambda[start],
             max.log.lambda=exact.df$max.log.lambda[end])
})
print(target.interval)
}
