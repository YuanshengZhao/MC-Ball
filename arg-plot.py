import matplotlib.pyplot as plt
import sys

if(sys.argv[1]=='x'):
    sz=3
    nplt=(int)(sys.argv[2])
else: 
    sz=1

dat=[]
for z in range(len(sys.argv)-sz):
    dat.append(float(sys.argv[z+sz]))

if(sz==3):
    if(len(dat)%(nplt+1) != 0):
        print('X Y mismatch!')
    else:
        kz=int(len(dat)/(nplt+1))
        ip=1
        while(ip<=nplt):
            plt.plot(dat[:kz],dat[kz*ip:kz*(ip+1)],marker="*")
            ip+=1
else:
    plt.plot(dat)
plt.show()