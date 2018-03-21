import sys
import matplotlib.pyplot as plt
import numpy as np 

data =np.genfromtxt(sys.argv[1], dtype =None, names=('Num', 'CHR', 'SNP', 'A1', 'A2', 'MAF', "NCHROBS","position","pi"), skip_header=1)

end_chr = data["position"][-1] +100
init_chr = data["position"][0] - 100
len_chr = end_chr - init_chr
window_size = 4000000 
step = 2500000 
n_windows = len_chr/step

x = []
y= []
def average(l):
    if len(l) != 0:
        total = 0
        for i in l:
            total += i
        return total/len(l)
    else:
        return 0
        
dont = False 
if dont != True:        
    for i in range(n_windows):
        if i == 0:
            id = window_size/2
            b_window = init_chr
            e_window = init_chr + window_size
        else:
            b_window = id - window_size/2
            e_window = id + window_size/2
        mean = []
        n = 0
        for pos in data["position"]:
            if e_window>pos>= b_window:
                mean.append(data["pi"][n])
                n+=1
            elif pos > e_window:
                break 
        y.append(average(mean))
        x.append(id/1000000)
        
        id += step
    fig = plt.figure(figsize=(8.0, 5.0))
    plt.plot(x,y) 
    plt.xlabel("Chromosomal Position (Mb)")
    plt.ylabel("Expected Heterozygosity")
    plt.savefig('Troglo_18_het.png')
    plt.show()
    plt.close()

dont = True  
if dont != True:    
    #No windows 
    y =  [] 
    x = []
    stick = []
    cummulative = 0
    previous = 1 
    for i in range(len(data["position"])): 
        chr = data["CHR"][i]
        position = data["position"][i] 
        if chr != previous:
            fig = plt.figure(figsize=(8.0, 5.0))
            plt.scatter(x,y, s =1) 
            plt.ylim([0,0.00015])
            plt.xlabel("Chromosomal Position (Mb)")
            plt.ylabel("Expected Heterozygosity")
            plt.savefig('het'+str(previous)+'.png')
            plt.close()   
            y =  [] 
            x = []
            
        else:
            x.append(position/1000000)
            y.append(data["pi"][i])
        previous = chr
    fig = plt.figure(figsize=(8.0, 5.0))
    plt.scatter(x,y, s=1) 
    plt.ylim([0,0.00015])

    plt.xlabel("Chromosomal Position (Mb)")
    plt.ylabel("Expected Heterozygosity")
    plt.savefig('het'+str(previous)+'.png')

    plt.close()    