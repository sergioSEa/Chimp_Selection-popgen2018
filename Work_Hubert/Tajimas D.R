chimp_data = read.table('F:\\courses\\Population Genetics\\Project\\ChimpExome\\pedfile.ped', header = F, sep="\t")
schwein_data = chimp_data[1:11,7:46679]
trogl_data = chimp_data[12:23,7:46679]
verus_data = chimp_data[24:29,7:46679]    # outgroup

# define ancestral alleles
ances_allele = rep('', 46673)
for(i in 1:100){
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
      if(genotype[j,i] == '0 0')
        next
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
sch_tajima_D = rep(0, 46)
for(i in 1:46){
  sch_derived_AF = cal_derived_AF(schwein_data[,((i-1)*1000+1):(i*1000)])
  sch_sfs = prop.table(table(sch_derived_AF))
  #sch_sfs = c(0.39, 0.19, 0.13, 0.1, 0.08, 0.06, 0.06)
  #names(sch_sfs) = c('0.125', '0.25', '0.375', '0.5', '0.625', '0.75', '0.875')
  sch_tajima_D[i] = cal_tajima_D(sch_sfs, 22)
}
plot(1:46, sch_tajima_D, ylim=c(-2, 2), xlab='genome', ylab="Tajima's D", main='Schwein population')
abline(h=0, col='red')
fit1 = lm(sch_tajima_D ~ 1)
abline(fit1, col='blue')

# for trogl group
tro_tajima_D = rep(0, 46)
for(i in 1:46){
  tro_derived_AF = cal_derived_AF(trogl_data[,((i-1)*1000+1):(i*1000)])
  tro_sfs = prop.table(table(tro_derived_AF))
  #sch_sfs = c(0.39, 0.19, 0.13, 0.1, 0.08, 0.06, 0.06)
  #names(sch_sfs) = c('0.125', '0.25', '0.375', '0.5', '0.625', '0.75', '0.875')
  tro_tajima_D[i] = cal_tajima_D(tro_sfs, 22)
}
plot(1:46, tro_tajima_D, ylim=c(-2, 2), xlab='genome', ylab="Tajima's D", main='Trogl population')
abline(h=0, col='red')
fit2 = lm(tro_tajima_D ~ 1)
abline(fit2, col='blue')

# 5211 1
# 3188 2
# 2544 3
# 2038 4
# 2226 5
# 2575 6
# 1973 7
# 1618 8
# 2031 9
# 2134 10
# 2874 11
# 2387 12
# 1015 13
# 1465 14
# 1444 15
# 1859 16
# 2611 17
# 894 18
# 2916 19
# 1310 20
# 657 21
# 982 22
# 721 23
# chromosome-wise
index = c(5211,3188,2544,2038,2226,2575,1973,1618,2031,2134,2874,2387,1015,
          1465,1444,1859,2611,894,2916,1310,657,982,721)
tro_tajima_D = rep(0, 23)
loc = 1
for(i in 1:23){
  tro_derived_AF = cal_derived_AF(trogl_data[,loc:(loc+index[i]-1)])
  loc = loc + index[i]
  tro_sfs = prop.table(table(tro_derived_AF))
  #sch_sfs = c(0.39, 0.19, 0.13, 0.1, 0.08, 0.06, 0.06)
  #names(sch_sfs) = c('0.125', '0.25', '0.375', '0.5', '0.625', '0.75', '0.875')
  tro_tajima_D[i] = cal_tajima_D(tro_sfs, 22)
}
plot(1:23, tro_tajima_D, ylim=c(-2, 2), xlab='genome', ylab="Tajima's D", main='Trogl population')
abline(h=0, col='red')
fit1 = lm(tro_tajima_D ~ 1)
abline(fit1, col='blue')






