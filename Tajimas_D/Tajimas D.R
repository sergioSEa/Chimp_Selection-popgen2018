chimp_data = read.table('F:\\courses\\Population Genetics\\Project\\ChimpExome\\pedfile.ped', header = F, sep="\t")
schwein_data = chimp_data[1:11,7:46679]
trogl_data = chimp_data[12:23,7:46679]
verus_data = chimp_data[24:29,7:46679]    # outgroup
rm(chimp_data)

# define ancestral alleles
ances_allele = rep('', 46673)
for(i in 1:46673){
  count = list('A'=0, 'T'=0, 'C'=0, 'G'=0)
  for(j in 1:6){
    if(verus_data[j,i] == '0 0')
      next
    allele = unlist(strsplit(as.character(verus_data[j,i]), ' '))
    count[allele[1]] = as.numeric(count[allele[1]]) + 1
    count[allele[2]] = as.numeric(count[allele[2]]) + 1
  }
  max_allele = 'A'
  if(count$T > count[max_allele]){
    max_allele = 'T'
  }
  if(count$C > count[max_allele]){
    max_allele = 'C'
  }
  if(count$G > count[max_allele]){
    max_allele = 'G'
  }
  ances_allele[i] = max_allele
}

cal_derived_AF = function(genotype){
  site_num = dim(genotype)
  derived_AF = rep(0, site_num[2])
  for(i in 1:site_num[2]){
    total = 0
    derived = 0
    for(j in 1:site_num[1]){
      total = total + 2
      if(genotype[j,i] == '0 0'){
        derived = total = 1
        break
      }
      allele = unlist(strsplit(as.character(genotype[j,i]), ' '))
      if(allele[1] != ances_allele[i]){
        derived = derived + 1
      }
      if(allele[2] != ances_allele[i]){
        derived = derived + 1
      }
    }
    derived_AF[i] = derived/total
  }
  return(derived_AF)
}

cal_tajima_D = function(sfs, n){
  a1 = 0
  a2 = 0
  for(i in 1:(n-1)){
    a1 = a1 + 1/i
    a2 = a2 + 1/(i^2)
  }
  b1 = (n + 1)/3/(n - 1)
  b2 = 2 * (n^2 + n + 3)/9/n/(n - 1)
  c1 = b1 - 1/a1
  c2 = b2 - (n + 2)/a1/n + a2/(a1^2)
  e1 = c1/a1
  e2 = c2/(a1^2 + a2)
  S = sum(sfs)
  watterson = S/a1
  para_t = 0
  for(eta in 1:length(sfs)){
    i = as.numeric(names(sfs[eta])) * n
    para_t = para_t + i * (n - i) * sfs[eta]
  }
  tajima = 1/choose(n, 2) * para_t
  tajima_D = (tajima - watterson)/sqrt(e1*S + e2*S*(S-1))
  return(tajima_D)
}

# for schwein group
plot_tajimasD_whole_genome = function(pop_data, n, pop_name){
  pop_tajima_D = rep(0, 46)
  for(i in 1:46){
    pop_derived_AF = cal_derived_AF(pop_data[,((i-1)*1000+1):(i*1000)])
    pop_sfs = table(pop_derived_AF)
    # get rid of fixed sites
    if(names(pop_sfs[1]) == '0')
      pop_sfs = pop_sfs[-1]
    # get rid of SNPs with missing genotype
    if(names(pop_sfs[length(pop_sfs)]) == '1')
      pop_sfs = prop.table(pop_sfs[-length(pop_sfs)])
    #sch_sfs = c(0.39, 0.19, 0.13, 0.1, 0.08, 0.06, 0.06)
    #names(sch_sfs) = c('0.125', '0.25', '0.375', '0.5', '0.625', '0.75', '0.875')
    pop_tajima_D[i] = cal_tajima_D(pop_sfs, n)
  }
  plot(1:46, pop_tajima_D, ylim=c(-2, 2), xlab='genome', ylab="Tajima's D", main=paste(pop_name, 'population whole genome'))
  abline(h=0, col='red')
  fit = lm(pop_tajima_D ~ 1)
  abline(fit, col='blue')
}

plot_tajimasD_whole_genome(schwein_data, 22, 'Schwein')
plot_tajimasD_whole_genome(trogl_data, 24, 'Trogl')


# chromosome-wise
plot_tajimasD_among_chr = function(pop_data, n, pop_name){
  index = c(5211,3188,2544,2038,2226,2575,1973,1618,2031,2134,2874,2387,1015,
            1465,1444,1859,2611,894,2916,1310,657,982,721)
  pop_tajima_D = rep(0, 23)
  loc = 1
  for(i in 1:23){
    pop_derived_AF = cal_derived_AF(pop_data[,loc:(loc+index[i]-1)])
    loc = loc + index[i]
    pop_sfs = table(pop_derived_AF)
    # get rid of fixed sites
    pop_sfs = pop_sfs[-1]
    # get rid of SNPs with missing genotype
    pop_sfs = prop.table(pop_sfs[-length(pop_sfs)])
    pop_tajima_D[i] = cal_tajima_D(pop_sfs, n)
  }
  plot(1:23, pop_tajima_D, ylim=c(-2, 2), xlab='genome', ylab="Tajima's D", main=paste(pop_name,'population among chr'))
  abline(h=0, col='red')
  fit = lm(pop_tajima_D ~ 1)
  abline(fit, col='blue')
}

plot_tajimasD_among_chr(schwein_data, 22, 'Schwein')
plot_tajimasD_among_chr(trogl_data, 24, 'Trogl')


# target loci
pbs_data = read.table('F:\\courses\\Population Genetics\\Project\\ChimpExome\\PBSrelevant.txt', header = F, sep=" ")
pos_data = read.table('F:\\courses\\Population Genetics\\Project\\ChimpExome\\pruneddata.bim', header = F, sep="\t")
pbs_pos = matrix(0, 50, 2)
for(i in 1:50){
  pos = which(pos_data == pbs_data[i,3])
  pbs_pos[i,2] = pos[1] - 46673
  pbs_pos[i,1] = pos_data[pbs_pos[i,2], 1]
}

# for each chromosome
plot_tajimasD_within_chr = function(pop_data, n, pop_name, dots){
  index = c(5211,3188,2544,2038,2226,2575,1973,1618,2031,2134,2874,2387,1015,
            1465,1444,1859,2611,894,2916,1310,657,982,721)
  pbs_tajima = c()
  loc = 1
  for(i in 1:23){
    chr_data = pop_data[,loc:(loc+index[i]-1)]
    loc = loc + index[i]
    chr_step = floor(index[i]/dots)
    pop_tajima_D = rep(0, dots)
    pbs_points = c()
    for(j in 1:dots){
      pop_derived_AF = cal_derived_AF(chr_data[,((j-1)*chr_step+1):(j*chr_step)])
      pop_sfs = table(pop_derived_AF)
      # get rid of fixed sites
      if(names(pop_sfs[1]) == '0')
        pop_sfs = pop_sfs[-1]
      # get rid of SNPs with missing genotype
      if(names(pop_sfs[length(pop_sfs)]) == '1')
        pop_sfs = prop.table(pop_sfs[-length(pop_sfs)])
      pop_tajima_D[j] = cal_tajima_D(pop_sfs, n)
      # label out PBS SNPs
      pbs_chr = which(pbs_pos == i)
      for(k in pbs_chr){
        pos = pbs_pos[k, 2]
        if(i > 1){
          pos = pos - sum(index[1:i-1])
        }
        if(pos >= (j-1)*chr_step+1 && pos <= j*chr_step){
          pbs_tajima = c(pbs_tajima, pop_tajima_D[j])
          pbs_points = c(pbs_points, j)
        }
      }
    }
    pbs_points = unique(pbs_points)
    # plotting
    # title = paste(pop_name,' population within chr',i, sep='')
    # png(filename=paste("F:\\courses\\Population Genetics\\Project\\Tajima's_D_within_chr\\", title, '.png', sep=''))
    # plot(1:dots, pop_tajima_D, ylim=c(-2, 2), xlab='chromosome', ylab="Tajima's D", main=title)
    # abline(h=0, col='green')
    # fit = lm(pop_tajima_D ~ 1)
    # abline(fit, col='blue')
    # points(pbs_points, pop_tajima_D[pbs_points], col='red')
    # dev.off()
  }
  new_pbs_data = cbind(pbs_data, round(pbs_tajima, digits=7))
  write.table(new_pbs_data,"F:\\courses\\Population Genetics\\Project\\ChimpExome\\PBS_Tajima's_D.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

plot_tajimasD_within_chr(schwein_data, 22, 'Schwein', 15)
plot_tajimasD_within_chr(trogl_data, 24, 'Trogl', 15)

