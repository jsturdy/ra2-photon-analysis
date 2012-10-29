import sys,os
for i in range(110):
    cmd = "python genLooperFull.py -m -s gjets -f %d -t directPhotonsID -i ID -c 100 > & gjets.id.100.%d.out &"%(i,i)
    print cmd
