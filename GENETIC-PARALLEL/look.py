import os

icount=0

while icount<128:
    os.system("cat ./"+str(icount)+"/config >> sample.xyz")
    icount+=1
