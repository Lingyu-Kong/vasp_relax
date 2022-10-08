import ase
from ase.io import read, write,Trajectory
from ase.io.trajectory import TrajectoryWriter
from ase.calculators.vasp import Vasp
import argparse
import os
import time
import wandb

parser=argparse.ArgumentParser()
parser.add_argument("--in_path",type=str,default="B3_sub_minimal")
parser.add_argument("--interleave",type=int,nargs="+",default=[0,2])  ## [0,1422]
parser.add_argument("--n_core",type=int,default=6)
parser.add_argument("--kspacing",type=float,default=0.2)
parser.add_argument("--gamma",action="store_true")
parser.add_argument("--wandb",action="store_true")
args=parser.parse_args()

ASE_VASP_COMMAND="mpirun -np "+str(args.n_core)+" vasp_std"

assert(args.in_path is not None)
assert(os.path.exists(args.in_path))
assert(os.path.exists(os.getcwd()+"/VASP_PP"))  ## VASP_PP must be under the working directory
assert(len(args.interleave)==2)
assert(args.interleave[0]<args.interleave[1])

if args.wandb:
    wandb.login(key="37f3de06380e350727df28b49712f8b7fe5b14aa")
    wandb.init(project="vasp run",name=args.in_path+str(args.interleave),config=args)

## set VASP_PP_PATH
os.environ["VASP_PP_PATH"]=os.path.join(os.getcwd(),"VASP_PP")

outfile=args.in_path+str(args.interleave)+".traj"

files=os.listdir(args.in_path)
files.sort()
files=files[args.interleave[0]:args.interleave[1]]

calc = Vasp(xc='PBE',
                    kspacing=args.kspacing,
                    gamma=args.gamma,
                    restart=None,            
                    command=ASE_VASP_COMMAND,
                    directory="./vasp_run",)

if __name__=="__main__":
    traj_writer=TrajectoryWriter(outfile,mode="a",properties=["energy","forces","stress"])
    for i,file in enumerate(files):
        atoms=read(os.path.join(args.in_path,file))
        atoms.set_calculator(calc)
        atoms.get_potential_energy()
        traj_writer.write(atoms)
        print("write {}/{}".format(i+1,len(files)))
        os.system("rm -rf ./vasp_run/*")
    if args.wandb:
        wandb.save(outfile)
    