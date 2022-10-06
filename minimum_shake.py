import argparse
import os
import time
import wandb
from ase.io import read, write
from ase.io.trajectory import TrajectoryWriter
from ase.calculators.vasp import Vasp
import zipfile38 as zipfile

parser=argparse.ArgumentParser()
parser.add_argument("--path",type=str,default=None)
parser.add_argument("--n_core",type=int,default=6)
parser.add_argument("--amplitude",type=float,default=0.02)
parser.add_argument("--shake_steps",type=int,default=20)  
parser.add_argument("--kspacing",type=float,default=0.2)
parser.add_argument("--interleave",nargs="+",type=int,default=[0,2])
parser.add_argument("--gamma",action="store_true")
parser.add_argument("--wandb",action="store_true")
args=parser.parse_args()

ASE_VASP_COMMAND="mpirun -np "+str(args.n_core)+" vasp_std"

## set VASP_PP_PATH
os.environ["VASP_PP_PATH"]=os.path.join(os.getcwd(),"VASP_PP")

if args.wandb:
    wandb.login(key="37f3de06380e350727df28b49712f8b7fe5b14aa")
    wandb.init(project="vasp run",name="shake "+args.path+str(args.interleave),config=args)

assert(args.path is not None)

if not os.path.exists(args.path+"_shake"):
    os.mkdir(args.path+"_shake")
else:
    os.system("rm -rf "+args.path+"_shake/*")
outpath=args.path+"_shake"

calc = Vasp(xc='PBE',
            kspacing=args.kspacing,
            gamma=args.gamma,
            restart=None,            
            command=ASE_VASP_COMMAND,
            directory=args.path+"/vasp_run",)

files=os.listdir(args.path)
files.sort()
files=files[args.interleave[0]:args.interleave[1]]


def zipDir(dirpath, outFullName):
    """
    压缩指定文件夹
    :param dirpath: 目标文件夹路径
    :param outFullName: 压缩文件保存路径+xxxx.zip
    :return: 无
    """
    zip = zipfile.ZipFile(outFullName, "w", zipfile.ZIP_DEFLATED)
    for path, dirnames, filenames in os.walk(dirpath):
        # 去掉目标跟路径，只对目标文件夹下边的文件及文件夹进行压缩
        fpath = path.replace(dirpath, '')
 
        for filename in filenames:
            zip.write(os.path.join(path, filename), os.path.join(fpath, filename))

if __name__=="__main__":
    for file in files:
        traj_writer=TrajectoryWriter(os.path.join(outpath,file.replace(".res",".traj")),mode="a",properties=["energy","forces","stress"])
        atoms=read(os.path.join(args.path,file))
        atoms.set_calculator(calc)
        traj_writer.write(atoms)
        for i in range(args.shake_steps):
            atoms.rattle(stdev=args.amplitude)
            atoms.set_calculator(calc)
            traj_writer.write(atoms)
            print("step: ",i," energy: ",atoms.get_potential_energy())
        os.system("rm -r "+os.path.join(args.path,"vasp_run"))
        traj_writer.close()
    zipDir(outpath,"shake_"+args.path+".zip")
    if args.wandb:
        wandb.save("shake_"+args.path+".zip")
