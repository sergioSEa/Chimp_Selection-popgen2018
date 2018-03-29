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
#The pi value should be specific by chrom 
verus <- cbind(verus,position = as.numeric(verus[,"SNP"]))
schwein <- cbind(schwein,position= as.numeric(schwein[,"SNP"]))
troglo <- cbind(troglo,position= as.numeric(troglo[,"SNP"]))
ver = 0
tro= 0
sch = 0
for (chr in seq(1:22)){
    v  = verus[verus[,"CHR"]==chr,]
    v_l = (v$position[length(v$position)] - v$position[1])
    ver = c(ver,het(v$MAF)/v_l)

    t  = troglo[troglo[,"CHR"]==chr,]
    t_l = (t$position[length(t$position)] - t$position[1])
    tro = c(tro,het(t$MAF)/t_l) 

    s  = schwein[schwein[,"CHR"]==chr,]
    s_l = (s$position[length(troglo$position)] - troglo$position[1])
    sch = c(sch,het(s$MAF)/v_l) 
    }

verus <- cbind(verus, pi=ver[2:length(ver)] )

schwein <- cbind(schwein, pi=sch[2:length(sch)])

troglo <- cbind(troglo, pi = tro[2:length(tro)])

#Optional: Plot total mean heterozygosity
# Making a barplot with the nucleotide diversity
options(scipen=7)
pdf ("nucleotide_global_diversity.pdf")
par(mfrow=c(1,1))
val = c(sum(schwein$pi), sum(troglo$pi),sum(verus$pi))
barplot(val,ylim=c(0.000,0.00015), ylab="pi", xlab="Population",
names.arg=c("schwein","trogl","verus"),xaxs="i",yaxs="i")
dev.off()


#The heterozygosity observed is smaller than in the exercise. This can be due to the fact 
#Mutations are affecting exons, which normally have negative selection and changes are not preserved


#Plot heterozygosity per Chromosome 
ver = 0
tro = 0
sch = 0
for (chr in seq(1:22)){
    v  = verus[verus[,"CHR"]==chr,]
    ver = c(ver,sum(v$pi))  
    t  = troglo[troglo[,"CHR"]==chr,]
    tro = c(tro,sum(t$pi))  
    s  = schwein[schwein[,"CHR"]==chr,]
    sch = c(sch,sum(s$pi))  
    }
pdf ("chrom_diversity.pdf")
ver = ver[2:length(ver)]
plot(seq(1:22), ver, xaxt="n", xlab="Number of Chromosome", ylab="Expected Heterozygosity",
   ylim=c(0.00000, 0.000015), type = "b",xaxs="i",yaxs="i")
#lines(seq(1:22), ver, xaxt="n", )
axis(1, at=seq(1, 22, by= 1), las=2)
title(main="Chromosomal Heterozygosity")

tro = tro[2:length(tro)]
lines(seq(1:22), tro, xaxt="n", col= "red", type= "b")

sch = sch[2:length(sch)]
lines(seq(1:22), sch, xaxt="n", col= "blue", type = "b")


legend(15, 0.000015 , legend=c("Verus", "Troglodites", "Schwein"),pch=15,col= c("black","red", "blue"))
dev.off()




####################Second Part, go through windows and calculate heterozygosity ################################ 
## Function for generating sliding windows


#troglo = troglo[troglo[,"CHR"]==3,]
#options(max.print=1000000)
#print(troglo)
pbs =read.table("PBS_interesting.txt", header =T)


window_size = 12000 #3000000 
step = 7000 #1500000
for (chr in seq(1:22)){
    t  = troglo[troglo[,"CHR"]==chr,]
    p = pbs[pbs[,"CHR"]==chr,]$POS/1000000
   
    t = cbind(t, h = het(t$MAF))
		
    end_chr = t$position[length(t$position)]
    init_chr = t$position[1]
    len_chr = end_chr - init_chr
    
    n_windows = round(len_chr/step)
    vector = 0
    vector_y = 0

       for (i in seq(1,n_windows)){
        if (i == 1){
            id = window_size/2
            b_window = init_chr
            e_window = init_chr + window_size
       }else{
            b_window = id - window_size/2
            e_window = id + window_size/2
	}
        n = 0
      
       #print(t$position[t$position> b_window & t$position <= e_window])
	in_window = t[t[,"position"] > b_window & t$position <= e_window,]
	#x = sum(in_window$h) / (in_window$position[length(in_window$position)] - in_window$position[1])
	x = sum(in_window$h) / (e_window -b_window)
	vector <- c(vector, x)
	if (length(x) == 0){
	vector  <- c(vector, 0)
	}
        vector_y <- c(vector_y,id/1000000)
        id =id +  step
}

pdf(paste(toString(chr),"diversity.pdf",sep = "-"))
plot(vector_y[2:length(vector_y)],vector[2:length(vector)], xlab="Position(MB)", ylab="Diversity", type = "l",xaxs="i",yaxs="i")
abline(v= p, col= "red")
title(main=paste(toString(chr),"Chromosomal Heterozygosity",sep = "-"))
dev.off()



}
