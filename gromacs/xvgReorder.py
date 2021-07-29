import sys
import numpy as np
import pandas as pd
from itertools import groupby

def sort_txt(in_file,out_file):
    pre_in_lines = [line for line in open(in_file,"r",encoding="utf-8")]

    in_lines = []
    
    #不要那些以 # 和 @ 开头的行
    for i in range(0,len(pre_in_lines)):
        curline = pre_in_lines[i]
        if (not curline.strip().startswith("#")) and (not curline.strip().startswith("@")):
           in_lines.append(curline)

    #排序
    out_lines = sorted(in_lines,key = lambda line:(float(list(filter(lambda x : x ,line.split(" ")))[1].replace('\n',''))),reverse = False)
  
    with open(out_file,"w",encoding="utf-8") as fw:
        fw.writelines(out_lines)


if __name__ == "__main__":
    in_file = sys.argv[1]
    prefix = in_file.split(".") 
    out_file =  prefix[0] + "-reorder.xvg"
    sort_txt(in_file,out_file)


#命令行
# python3  xvgReorder.py  /home/wurp/workspace/antibody/6YZ5/rmsd_vs_average.xvg
#在输入文件的路径下输出文件  rmsd_vs_average-reorder.xvg
