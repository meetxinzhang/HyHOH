import sys
import numpy as np
import pandas as pd
from itertools import groupby

def sort_txt(in_file,out_file):
    in_lines = [line for line in open(in_file,"r",encoding="utf-8")]
 #   a = in_lines.split(" ")
 #   b = list(filter(lambda x : x, a))
  #  out_lines = sorted(in_lines,key = lambda line:(int(line.split(" ")[1])),reverse = True)
    out_lines = sorted(in_lines,key = lambda line:(float(list(filter(lambda x : x ,line.split(" ")))[1].replace('\n',''))),reverse = False)
  
 #   print(float(list(filter(lambda x : x ,out_lines[0].split(" ")))[1].replace('\n','')))
    listY = []
    for perline in out_lines:
        a = round(float(list(filter(lambda x : x ,perline.split(" ")))[1].replace('\n','')),2)
        listY.append(a)

  #  print(listY)

    with open("listY.txt","w",encoding="utf-8") as fw:
        fw.writelines(str(listY))

    for key,group in groupby(listY,key = lambda x : x / 0.01 ):
        print('{}：{}'.format(key/100,len(list(group))))
     #   print(k,len(list(g)))

   
  
    with open(out_file,"w",encoding="utf-8") as fw:
        fw.writelines(out_lines)


if __name__ == "__main__":
    in_file = "rmsd.xvg"
    out_file = "reOrderRMSD.xvg"
    sort_txt(in_file,out_file)

#gmx trjconv -f md_0_1.xtc  -o mostFre.xtc -drop rmsd.xvg -dropunder 0.14 -dropover 0.15
