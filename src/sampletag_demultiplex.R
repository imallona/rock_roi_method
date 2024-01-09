#!/usr/bin/env R
##
## Assigns sampletags to CBs in a BD-compliant way. Untested using real data.
##
## 09th Jan 2024
## Izaskun Mallona

## Attempts to follow:
## https://bd-rhapsody-bioinfo-docs.genomics.bd.com/steps/steps_sample_tag.html
## comments quoted come from the URL above

## also simulates
## 30 cells,
##   cb1-10 are (truly) sampletag1, cb11-20 are sampletag2, cb20-30 are pure noise; all of them have noise

set.seed(1)
fd <- t(rbind(data.frame(st1 = rbinom(10, 10, 0.6),
                         st2 = rbinom(10, 10, 0.1),
                         st3 = rbinom(10, 10, 0.05)),
              data.frame(st1 = rbinom(10, 10, 0.1),
                         st2 = rbinom(10, 10, 0.6),
                         st3 = rbinom(10, 10, 0.05))))
fd <- cbind(fd, t(data.frame(st1 = rbinom(10, 10, 0.4),
                             st2 = rbinom(10, 10, 0.3),
                             st3 = rbinom(10, 10, 0.4))))

colnames(fd) <- paste0('cb', 1:30)

fd

## "A high quality singlet is a putative cell where more than 75% of Sample Tag reads are from a single tag"

demux <- data.frame(cb = colnames(fd))

## "When a singlet is identified, the counts for all the other tags are considered Sample Tag noise."

foo <- t(apply(fd, 2, function(x) {
    tmp <- proportions(x)[proportions(x) >= 0.75]
    sts <- names(tmp)
    counts <- as.numeric(x[sts][1])
    if (length(sts) == 0) {
        sts <- NA
        counts <- 0
    }
    return(c(sts, counts))
}))


demux$st <- foo[,1]
demux$highqual_counts <- as.numeric(foo[,2])
demux$total_counts <- colSums(fd)
demux$noise <- colSums(fd) - demux$highqual_counts
demux$status <- ifelse(!is.na(demux$st), 'highqual', 'undetermined')

## "The minimum Sample Tag read count for a putative cell to be positively identified with a
##  Sample Tag is defined
##  as the lowest read count of a high quality singlet for that Sample Tag"
thres <- list()
for (st in rownames(fd)) {
    highquals <- demux[demux$status == 'highqual' & demux$st == st, 'cb']
    if (length(highquals) > 0)
        thres[[st]] <- min(fd[st, highquals])
    else
        thres[[st]] <- NA
}

for (st in names(thres)) {
    demux[,st] <- ifelse(fd[st, ] >= thres[[st]], yes = fd[st,], no = NA)
}

## "
## The percentage of noise contribution of each Sample Tag for all cells is calculated by dividing
## the total per tag noise by the total overall noise.

## In addition, the total amount of noise versus the total Sample Tag count per putative cell is recorded
## so that a trend line can be established to estimate the total per-cell noise
## given an observed number of total Sample Tag count for a cell.
## The figure below shows an overall noise profile where each dot is a cell.
## A trend line (in red) is fitted and used to establish the expected amount of noise given
## a total Sample Tag count. Cells that are off the trend line are likely multiplets.
## The level of antigen expression across cells can vary, contributing to variation in
## Sample Tag count per cell. Generally, cells with higher total Sample Tag
## counts have higher noise Sample Tag counts.
## "

(tags_noise <- data.frame(aggregate(demux$noise, list(st = demux$st), function(x) sum(x)/sum(demux$noise)),
                          row.names = 1))

plot(demux$highqual_counts, demux$noise)
(m <- lm(demux$noise ~ demux$total_counts))

## let's hope the slope is positive!
stopifnot(m$coefficients[2] > 0)

## predoct
demux$predicted_noise <- predict(object = m, data = demux$total_counts)

for (st in names(thres)) {
    demux[paste0('weighted_noise_', st)] <- demux[, 'predicted_noise'] * tags_noise [st, 'x']

    demux[paste0('corrected_', st)] <- demux[, st] - demux[, paste0('weighted_noise_', st)]
    demux[paste0('corrected_', st)][is.na(demux[paste0('corrected_', st)])] <- 0
}


demux$called_st <- gsub('corrected_', '',
                           apply(demux[, paste0('corrected_', names(thres))],
                                 1,
                                 function(x) paste0(names(demux[, paste0('corrected_', names(thres))])[x > min(x)], collapse = ',')))


demux$status[demux$status != 'highqual' & nchar(demux$called_st) > 0] <- 'called'

table(demux$status)

## clean a bit
demux <- demux[,c('cb', 'st', 'status', 'called_st')]

demux
