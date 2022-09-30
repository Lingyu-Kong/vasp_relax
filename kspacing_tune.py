import ase
from ase.io import read, write,Trajectory
from ase.calculators.vasp import Vasp
from ase.constraints import ExpCellFilter
from ase.optimize.precon import PreconLBFGS, Exp
import argparse
import os
import time
import wandb

parser=argparse.ArgumentParser()
parser.add_argument("--path",type=str,default=None)  ## file path
parser.add_argument("--name",type=str,default=None)  ## system name
parser.add_argument("--n_core",type=int,default=6)
parser.add_argument("--kspacing_interleave",type=float,nargs="+",default=[0.02,0.5])
parser.add_argument("--step_size",type=float,default=0.1)
parser.add_argument("--encut",type=float,default=400)
parser.add_argument("--ediff",type=float,default=1e-4)
parser.add_argument("--ismear",type=int,default=0)
parser.add_argument("--sigma",type=float,default=0.02)
parser.add_argument("--gamma",action="store_true")
parser.add_argument("--nelm",type=int,default=200)
parser.add_argument("--wandb",action="store_true")
parser.add_argument("--print_log",action="store_true")
args=parser.parse_args()

if args.wandb:
    wandb.login(key="37f3de06380e350727df28b49712f8b7fe5b14aa")
    wandb.init(project="vasp run",name="kspacing_tune "+args.name,config=args)

ASE_VASP_COMMAND="mpirun -np "+str(args.n_core)+" vasp_std"

## set VASP_PP_PATH
os.environ["VASP_PP_PATH"]=os.path.join(os.getcwd(),"VASP_PP")

atoms=read(args.path)

if __name__=="__main__":
    kspacing=args.kspacing_interleave[1]
    while kspacing>=args.kspacing_interleave[0]:
        calc = Vasp(xc='PBE',
            encut=args.encut,
            ediff=args.ediff,
            ismear=args.ismear,
            sigma=args.sigma,
            kspacing=kspacing,
            gamma=args.gamma,
            nelm=args.nelm,
            restart=None,            
            command=ASE_VASP_COMMAND,
            directory="./vasp_run",
            txt="-" if args.print_log else None)
        atoms.set_calculator(calc)
        energy=atoms.get_potential_energy()
        wandb.log({"kspacing":kspacing,"energy":energy})
        print("kspacing: ",kspacing," energy: ",energy)
        kspacing-=args.step_size
        if kspacing<=(0.1+1e-2):
            args.step_size=0.02
        
        