import ase
import os
from ase.io import read,write
import argparse
import numpy as np

parser=argparse.ArgumentParser()
parser.add_argument("--in_path",type=str,default="B8_minimal")
parser.add_argument("--n",type=int,default=3)
parser.add_argument("--cubic",action="store_true")
parser.add_argument("--out_path",type=str,default="B3_sub_minimal")
args=parser.parse_args()

assert(args.in_path is not None)
assert(os.path.exists(args.in_path))
assert(args.out_path is not None)
if not os.path.exists(args.out_path):
    os.mkdir(args.out_path)
if args.cubic and not os.path.exists(args.out_path+"_cubic"):
    os.mkdir(args.out_path+"_cubic")

files=os.listdir(args.in_path)

EDGE_INDEX=[[0,0,1,1,2,2],[1,2,0,2,0,1]]

if __name__=="__main__":
    for file in files:
        if file.endswith(".res"):
            atoms=read(os.path.join(args.in_path,file))
            positions=np.array(atoms.get_positions())
            count=0
            with open(os.path.join(args.in_path,file),"r") as f:
                lines=f.readlines()
            while count<args.n:
                indices=np.random.randint(0,len(positions),size=3)
                pos=positions[indices]
                distances=np.linalg.norm(pos[EDGE_INDEX[0]]-pos[EDGE_INDEX[1]],axis=1)
                if np.all(distances>1.0) and np.all(distances<3.75):
                    indices=indices+4
                    line_indices=[0,1,2,3]+indices.tolist()+[12]
                    with open(os.path.join(args.out_path,file.replace(".res","_"+str(count)+".res")),"w") as f:
                        for i in line_indices:
                            f.write(lines[i])
                        count+=1
    if args.cubic:
        files=os.listdir(args.out_path)
        for file in files:
            if file.endswith(".res"):
                atoms=read(os.path.join(args.out_path,file))
                cell_line="CELL 1.0 20 20 20 90.00 90.00 90.00\r"
                positions=np.array(atoms.get_positions())
                positions=positions+10
                frac_pos=positions/20
                with open(os.path.join(args.out_path,file),"r") as f:
                    lines=f.readlines()
                with open(os.path.join(args.out_path+"_cubic",file),"w") as f:
                    for i in range(len(lines)):
                        if i == 1:
                            f.write(cell_line)
                        elif i == 4:
                            f.write("B 1 {} {} {} 1.0\r".format(frac_pos[0][0],frac_pos[0][1],frac_pos[0][2]))
                        elif i == 5:
                            f.write("B 1 {} {} {} 1.0\r".format(frac_pos[1][0],frac_pos[1][1],frac_pos[1][2]))
                        elif i == 6:
                            f.write("B 1 {} {} {} 1.0\r".format(frac_pos[2][0],frac_pos[2][1],frac_pos[2][2]))
                        else:
                            f.write(lines[i])
                        
                
                
                
            