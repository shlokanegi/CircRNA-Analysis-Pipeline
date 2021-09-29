#---- CLP Analysis ----

# Set the working directory
setwd("~/data/analyses")

# Load packages
require(devtools)
install_github('dieterich-lab/CircTest')
library(CircTest)

# Prepare count tables
### CIRCULAR READS
circ.bks.circxpr <- combined.df[, c(1, 2, 6)]
colnames(circ.bks.circxpr)[1] <- "sample_id"
colnames(circ.bks.circxpr)[2] <- "circ_id"
colnames(circ.bks.circxpr)[3] <- "circ.reads"
circ.circ.xpr <- data.table(circ.bks.circxpr) # Formed table to perform reshaping

Circ <- data.frame(dcast(circ.circ.xpr, formula=circ_id ~ sample_id, 
                           value.var="circ.reads", fill=NA), row.names="circ_id")
Circ$circ_id <- rownames(Circ)
Circ <- na.omit(Circ)


### LINEAR READS
input <- list()
for(i in 1:length(Samples)){
  input[[i]] <- paste0("../", Samples[i],"_bks_linear_counts.tab")
}
names(input) = Samples

# Combining counts.tab file for all samples
bks_linear_counts.tab.gz <- rbindlist(lapply(input, fread, data.table=T,
            showProgress=F), use.names=T, idcol="sample_id")

# Just replicating the V4 column and naming it as circ_id
bks_linear_counts.tab.gz[, `:=`(circ_id = sub('.*"([^"]+)".*', "\\cr1", V4))]

# Summing up the left and right flanked reads of the back-spliced junction
ciri_bks_linexp <- bks_linear_counts.tab.gz[, .(lin.reads=sum(V7)), 
                                            by = .(sample_id, circ_id)]
circ.lin.xpr <- data.table(ciri_bks_linexp)

Linear <- data.frame(dcast(circ.lin.xpr, formula=circ_id ~ sample_id, 
        value.var="lin.reads", fill=NA), row.names="circ_id")
Linear$circ_id <- rownames(Linear)
Linear <- na.omit(Linear) #269 circRNAs

Linear <- Linear[Circ$circ_id, ]

# Run Circ.test
test <- Circ.test(Circ, Linear, group=c(rep(1,3), rep(2,3)), 
                  circle_description = 7)   # We have to provide numbers 
                                              #in groupings (not strings)
View(test$summary_table)

# PLOTS
for (i in rownames(test$summary_table)) {
  Circ.lineplot(Circ, Linear, plotrow=i, groupindicator1=c(rep('B-Cell',3),rep('Monocyte',3)),
                circle_description = 7)
  n = readline(prompt = ">>>")
  if(n == "n") break
}

for (i in rownames(test$summary_table)) { 
  Circ.ratioplot(Circ, Linear, plotrow=i, groupindicator1=c(rep('B-Cell',3),rep('Monocyte',3)), 
                 lab_legend='Condition', circle_description = 7)
  n = readline(prompt = ">>>")
  if(n == "n") break
}
