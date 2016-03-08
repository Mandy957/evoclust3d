import os
import sys

need = "/work/j7adams/Summer2015/august-rerun/NeedToRerun"
DIR = "/work/j7adams/Summer2015/selectome-Summer2015/results/"

todo = []

os.system("ssh j7adams@saw.sharcnet.ca ls -1 /scratch/j7adams/results/*/PValues.txt > doneonsaw")

doneOnSaw = [line.replace("/scratch/j7adams/results/","").replace("/PValues.txt\n","") for line in open("doneonsaw","r").readlines()]


for line in open(need,"r").readlines():
    f = line.replace("\n","")
    
    potentialpath = DIR+f+"/PValues.txt"
    if os.path.exists(potentialpath) or f in set(doneOnSaw):
        pass
    else:
        todo.append(f)




for f in todo[int(sys.argv[1]):int(sys.argv[2])]:
    
    suffix = "python adaptation3d.py --data-dir /home/j7adams/Adaptation3D/ --aln /scratch/j7adams/notimportant --tree /scratch/j7adams/notimportant --out-dir /work/j7adams/Summer2015/selectome-Summer2015/results/"+f
    
    command = suffix
    
    os.system(command)

