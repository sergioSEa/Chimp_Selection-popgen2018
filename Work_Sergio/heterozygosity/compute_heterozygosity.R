setwd("/home/mjq483/exam/Chimp_Selection-popgen2018/Work_Sergio/heterozygosity/")
#Load Populations Frequency files
schwein<-read.table("Pt_schwein_noNA.frq",h=T)
troglo = read.table("Pt_troglo_noNA.frq", h=T)
verus<-read.table("Pt_verus_noNA.frq",h=T)
#Compute heterozygosity as 2pq, assuming H-W proportions
het<-function(x){2*x*(1-x)}
#Take the SNPs which are heterozygotes 
verus <- verus[verus[,"MAF"]>0,]
schwein <- schwein[schwein[,"MAF"]>0,]
troglo <- troglo[troglo[,"MAF"]>0,]
#Creates a list per each subspecie, the last column is the heterozygosity
#Format --> CHR    SNP A1 A2    MAF NCHROBS position           pi
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
# ls() should show the three lists and a function , typeof() shows which type it is 

#Optional: Plot total mean heterozygosity
# Making a barplot with the nucleotide diversity
pdf ("nucleotide_global_diversity")
par(mfrow=c(1,1))
val = c(mean(schwein$pi), mean(troglo$pi),mean(verus$pi))
barplot(val,ylim=c(0.000,0.00015), ylab="pi", xlab="Population",
names.arg=c("schwein","trogl","verus"))
dev.off()


#The heterozygosity observed is smaller than in the exercise. This can be due to the fact 
#Mutations are affecting exons, which normally have negative selection and changes are not preserved


#Plot heterozygosity per Chromosome 
ver = 0
tro = 0
sch = 0
for (chr in seq(1:22)){
    v  = verus[verus[,"CHR"]==chr,]
    ver = c(ver,mean(v$pi))  
    t  = troglo[troglo[,"CHR"]==chr,]
    tro = c(tro,mean(t$pi))  
    s  = schwein[schwein[,"CHR"]==chr,]
    sch = c(sch,mean(s$pi))  
    }
pdf ("chrom_diversity")
ver = ver[2:length(ver)]
plot(seq(1:22), ver, xaxt="n", xlab="Number of Chromosome", ylab="Expected Heterozygosity",
   ylim=c(0.00001, 0.00008), type = "b")
#lines(seq(1:22), ver, xaxt="n", )
axis(1, at=seq(1, 22, by= 1), las=2)
title(main="Chromosomal Heterozygosity")

tro = tro[2:length(tro)]
lines(seq(1:22), tro, xaxt="n", col= "red", type= "b")

sch = sch[2:length(sch)]
lines(seq(1:22), sch, xaxt="n", col= "blue", type = "b")


legend(15, 0.00008 , legend=c("Verus", "Troglodites", "Schwein"),pch=15,col= c("black","red", "blue"))
dev.off()




####################Second Part, go through windows and calculate heterozygosity ################################ Ï'll use this to fine-tune the signals observed previously
## Function for generating sliding windows
slidingwindowplot <- function(mainv, xlabv, ylabv, ylimv, window.size,
step.size,input_x_data,input_y_data)
{
if (window.size > step.size){

step.positions <- seq(window.size/2 + 1, length(input_x_data)- window.size/2, by=step.size)
}else{
step.positions <- seq(step.size/2 + 1, length(input_x_data)- step.size, by=step.size) }
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

#pdf ("nucleotide_diversity_in_3_subspecies")
par(mfrow=c(2,2))
windowsize = 400
steps = 100

# Pan troglodytes verus
#mainvv = paste("verus pi = ",format(mean(verus$pi,na.rem=TRUE), digits=3), "SNPs =",
#length(verus$pi), "Win: ", windowsize, "Step: ", steps)
#slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
#ylab=expression(paste("pi")),ylimv=c(0.00,0.0016), window.size=windowsize/4,
#step.size=steps, input_x_data=verus$position/1000000,input_y_data=verus$pi)

# Pan troglodytes schweinfurthii
#mainvv = paste("schwein pi = ",format(mean(schwein$pi,na.rem=TRUE), digits=3),"SNPs
#=", length(schwein$pi),"Win: ", windowsize, "Step: ", steps )
#slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
#ylab=expression(paste("pi")),ylimv=c(0.000,0.0016), window.size=windowsize,
#step.size=steps, input_x_data=schwein$position/1000000,input_y_data=schwein$pi)

# Pan troglodytes troglodytes
#mainvv = paste("troglo pi = ",format(mean(troglo$pi,na.rem=TRUE), digits=3),"SNPs
#=", length(troglo$pi),"Win: ", windowsize, "Step: ", steps )
#slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")),
#ylab=expression(paste("pi")),ylimv=c(0.00,0.0016), window.size=windowsize,
#step.size=steps, input_x_data=troglo$position/1000000,input_y_data=troglo$pi)


#I am doing it all with  the Python scripts
#Interest : Troglodytes Chr 18
#pdf("nucleotide_diversity_tro_18")
#troglo = troglo[troglo[,"CHR"]==18,]

#mainvv = paste("troglo pi = ",format(mean(troglo$pi,na.rem=TRUE), digits=3),"SNPs=", length(troglo$pi),"Win: ", windowsize, "Step: ", steps )

#print(troglo)
#slidingwindowplot(mainv=mainvv, xlab=expression(paste("Position (x ", 10^6,")")), ylab=expression(paste("pi")),ylimv=c(0.000045,0.000055), window.size=windowsize,step.size=steps, input_x_data=troglo$position/1000000,input_y_data=troglo$pi)

#troglo = troglo[troglo[,"CHR"]==2,]
options(max.print=1000000)
print(troglo)

