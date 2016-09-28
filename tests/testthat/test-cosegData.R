context("cosegData")

data.name.vec <- c(
  "H3K36me3_AM_immune_McGill0079_chr3_60000_66170270"
)
for(data.name in data.name.vec){
  data(list=data.name, package="cosegData")
  data.list <- get(data.name)
  fit <- with(data.list, PeakSegFPOPchrom(coverage, as.numeric(penalty)))
}
