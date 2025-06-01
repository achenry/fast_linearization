# import weis, ROSCO, pyFAST
import os

weis_dir = "../toolboxes/WEIS" # "#os.path.abspath(weis.__file__)
rosco_dir = "../toolboxes/WEIS/ROSCO/ROSCO" #os.path.abspath(ROSCO.__file__)
pyfast_dir = "../toolboxes/python-toolbox/" #os.path.abspath(pyFAST.__file__)
fast_install_dir = '/Users/aoifework/Documents/OpenFAST-src/openfast_usflowt/install'

# start conda
os.system("ssh scompile")
os.system("source /curc/sw/anaconda3/latest")

# compile ROSCO
if not os.path.exists(os.path.join(rosco_dir, 'install')):
    cwd = os.getcwd()
    os.chdir(rosco_dir)
    if os.path.exists('./build'):
        os.system("rm -rf build")
    os.mkdir('build')
    os.chdir('build')
    os.system(f"cmake .. && make -j4 install")
    os.chdir(cwd)

os.system(f"conda create -y --name weis-env python=3.8")
os.system(f"conda activate weis-env")
os.system("conda install -y cmake cython geopy git jsonschema make matplotlib-base numpy numpydoc openmdao openpyxl "
          "pandas pip pytest pyyaml ruamel_yaml scipy setuptools shapely six sympy swig xlrd")
os.system("conda install -y petsc4py mpi4py compilers")
os.system("pip install simpy marmot-agents jsonmerge")

# setup pyFAST
os.chdir(pyfast_dir)
os.system(f"python {os.path.join(os.path.dirname(pyfast_dir), 'setup.py')} develop")

# setup weis
os.chdir(weis_dir)
os.system(f"python {os.path.join(os.path.dirname(weis_dir), 'setup.py')} develop")

# copy model files and code
os.system("scp -r '../models/OpenFAST-Dev3.0_Models_Dist_v3' aohe7145@login.rc.colorado.edu:/projects/aohe7145/OpenFAST-models")
os.system("scp -r ./* aohe7145@login.rc.colorado.edu:/projects/aohe7145/linearization")
