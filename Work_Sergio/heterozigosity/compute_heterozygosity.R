schwein<-read.table("Pt_schwein_noNA.frq",h=T)
troglo = read.table("Pan_troglodytes_noNA.fr", h=T)
verus<-read.table("Pt_verus_noNA.frq",h=T)

het<-function(x){2*x*(1-x)}

verus <- verus[verus[,"MAF"]>0,]
schwein <- schwein[schwein[,"MAF"]>0,]
troglo <- troglo[troglo[,"MAF"]>0,]

verus <- cbind(verus,position= as.numeric(gsub("22:",'',verus[,"SNP"])))
verus <- cbind(verus, pi=het(verus$MAF)
*(length(verus$MAF)/(verus[length(verus[,"position"]),"position"] - verus
[1,"position"])))

schwein <- cbind(schwein,position= as.numeric(gsub("22:",'',schwein[,"SNP"])))
schwein <- cbind(schwein, pi=het(schwein$MAF) *(length(schwein$MAF)/(schwein
[length(schwein[,"position"]),"position"] - schwein [1,"position"])))

troglo <- cbind(troglo,position= as.numeric(gsub("22:",'',troglo[,"SNP"])))
troglo <- cbind(troglo, pi=het(troglo$MAF) *(length(troglo$MAF)/(troglo
[length(troglo[,"position"]),"position"] - troglo [1,"position"])))


####################Second Part, go through windows and calculate heterozygosity ################################
## Function for generating sliding windows
slidingwindowplot <- function(mainv, xlabv, ylabv, ylimv, window.size,
step.size,input_x_data,input_y_data)
{
if (window.size > step.size)
step.positions <- seq(window.size/2 + 1, length(input_x_data)- window.size/2,
by=step.size)
else
step.positions <- seq(step.size/2 + 1, length(input_x_data)- step.size,
by=step.size)
n <- length(step.positions)
means_x <- numeric(n)
means_y <- numeric(n)
for (i in 1:n) {
chunk_x <- input_x_data[(step.positions[i]-
window.size/2):(step.positions[i]+window.size-1)]
means_x[i] <- mean(chunk_x,na.rem=TRUE)
chunk_y <- input_y_data[(step.positions[i]-
window.size/2):(step.positions[i]+window.size-1)]
means_y[i] <- mean(chunk_y,na.rem=TRUE)
}
plot(means_x,means_y,type="b",main=mainv,xlab=xlabv,ylab=ylabv,ylim=ylimv,cex=0.25,
pch=20, cex.main=0.75)
vec <- c(0.025,0.5,0.975)
zz <- means_y[!is.na(means_y)]
abline(h=quantile(zz,0.025,na.rem=TRUE),col="blue")
abline(h=quantile(zz,0.925,na.rem=TRUE),col="blue")
abline(h=mean(input_y_data))
}


pdf ("nucleotide_diversity_in_3_subspecies")
par(mfrow=c(2,2))
windowsize<- 3000
steps<- 100
# Pan troglodytes verus
mainvv = paste("verus pi = ",format(mean(verus$pi,na.rem=TRUE), digits=3), "SNPs =",
length(verus$pi), "Win: ", windowsize, "Step: ", steps)
slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
ylab=expression(paste("pi")),ylimv=c(0.00,0.0016), window.size=windowsize/4,
step.size=steps, input_x_data=verus$position/1000000,input_y_data=verus$pi)

# Pan troglodytes schweinfurthii
mainvv = paste("schwein pi = ",format(mean(schwein$pi,na.rem=TRUE), digits=3),"SNPs
=", length(schwein$pi),"Win: ", windowsize, "Step: ", steps )
slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
ylab=expression(paste("pi")),ylimv=c(0.000,0.0016), window.size=windowsize,
step.size=steps, input_x_data=schwein$position/1000000,input_y_data=schwein$pi)

# Pan troglodytes troglodytes
mainvv = paste("troglo pi = ",format(mean(troglo$pi,na.rem=TRUE), digits=3),"SNPs
=", length(troglo$pi),"Win: ", windowsize, "Step: ", steps )
slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
ylab=expression(paste("pi")),ylimv=c(0.00,0.0016), window.size=windowsize,
step.size=steps, input_x_data=troglo$position/1000000,input_y_data=troglo$pi)
