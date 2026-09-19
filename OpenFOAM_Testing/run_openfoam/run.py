import os
import shutil
import subprocess
from pathlib import Path

case_dir = Path(__file__).parent.resolve()
os.chdir(case_dir)

subprocess.run("topoSet > log.init 2>&1", shell=True, check=True)
subprocess.run("setFields >> log.init 2>&1", shell=True, check=True)
subprocess.run("decomposePar >> log.init 2>&1", shell=True, check=True)
subprocess.run("mpirun -np 6 interFoam -parallel > log.interFoam 2>&1", shell=True, check=True) 