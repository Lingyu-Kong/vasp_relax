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
parser.add_argument("--path",type=str,default=None)  ## path and path/data must exist
parser.add_argument("--n_core",type=int,default=8)
parser.add_argument("--relax_steps",type=int,default=1000)
parser.add_argument("--shake_steps",type=int,default=0)     
parser.add_argument("--sample_freq",type=int,default=1)
parser.add_argument("--interleave",nargs="+",type=int,default=[0,1500])  ## interleave [low,high)
parser.add_argument("--fmax",type=float,default=0.001)
parser.add_argument("--smax",type=float,default=0.02)
parser.add_argument("--encut",type=float,default=400)
parser.add_argument("--ediff",type=float,default=1e-4)
parser.add_argument("--ismear",type=int,default=0)
parser.add_argument("--sigma",type=float,default=0.02)
parser.add_argument("--amplitude",type=float,default=0.2)
parser.add_argument("--kspacing",type=float,default=0.02)
parser.add_argument("--gamma",action="store_true")
parser.add_argument("--nelm",type=int,default=200)
parser.add_argument("--wandb",action="store_true")
parser.add_argument("--print_log",action="store_true")
args=parser.parse_args()
ASE_VASP_COMMAND="mpirun -np "+str(args.n_core)+" vasp_std"

assert(args.path is not None)
assert(os.path.exists(args.path))
assert(os.path.exists(os.path.join(args.path,"data")))
assert(args.relax_steps>0)
assert(args.shake_steps>=0)
assert(os.path.exists(os.getcwd()+"/VASP_PP"))  ## VASP_PP must be under the working directory

if args.wandb:
    wandb.login(key="37f3de06380e350727df28b49712f8b7fe5b14aa")
    wandb.init(project="vasp run",name="relax "+args.path+" ["+str(args.interleave[0])+","+str(args.interleave[1])+"]",config=args)

## set VASP_PP_PATH
os.environ["VASP_PP_PATH"]=os.path.join(os.getcwd(),"VASP_PP")

## relax trajectories are stored in path/relax
if not os.path.exists(os.path.join(args.path,"relax")):
    os.mkdir(os.path.join(args.path,"relax"))
if not os.path.exists(os.path.join(args.path,"vasp_run")):
    os.mkdir(os.path.join(args.path,"vasp_run"))
    
input_path=os.path.join(args.path,"data")
relax_path=os.path.join(args.path,"relax")
shake_path=os.path.join(args.path,"shake")
input_files=os.listdir(input_path)
input_files.sort()
print(len(input_files))
input_files=input_files[args.interleave[0]:args.interleave[1]]


calc = Vasp(xc='PBE',
            encut=args.encut,
            ediff=args.ediff,
            ismear=args.ismear,
            sigma=args.sigma,
            kspacing=args.kspacing,
            gamma=args.gamma,
            nelm=args.nelm,
            restart=None,            
            command=ASE_VASP_COMMAND,
            directory=args.path+"/vasp_run",
            txt="-" if args.print_log else None)

if __name__=="__main__":
    os.system("rm -rf "+args.path+"/vasp_run/*")
    os.system("rm -rf "+args.path+"/relax/*")
    for file in input_files:
        if file.endswith(".res"):
            try:
                start_time=time.time()
                ## relax begins:
                atoms=read(os.path.join(input_path,file))
                atoms.set_calculator(calc)
                ecf=ExpCellFilter(atoms)
                traj = Trajectory(relax_path+"/traj_"+file.replace(".res",".traj"), 'w', atoms, properties=["energy","forces","stress"])
                optimizer = PreconLBFGS(ecf, precon=Exp(3), use_armijo=True, master=True)
                optimizer.attach(traj.write, interval=args.sample_freq)
                optimizer.run(fmax=args.fmax, smax=args.smax, steps=args.relax_steps)
                traj.close()
                end_time=time.time()
                ## relax ends
                print("relax finished for {} in {} seconds".format(file,end_time-start_time))
                if args.wandb:
                    wandb.log({"relaxed energy":atoms.get_potential_energy()})
            except:
                print("relax failed for {}".format(file))
            ## clean up the vasp_run directory
            os.system("rm -rf "+args.path+"/vasp_run/*")
    if args.wandb:
        wandb.save(args.path+"/relax/*")
        


