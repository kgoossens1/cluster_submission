#!/usr/bin/env python3

"""
This is a tool with the purpose of submitting jobs to the UAntwerp clusters.
Currently supporting job submission of Gaussian, Gromacs and Amber jobs to
the tier-2 cluster Hopper, Leibniz and Vaughan (both GPU and CPU nodes) and
the tier-1 cluster Breniac. 
"""

import os
import sys
import subprocess
import shlex
import time
import yaml
from argparse import ArgumentParser, HelpFormatter

__author__ = "Kenneth Goossens"
__version__ = "1.1.1"
__email__ = "goossens_kenny@hotmail.com"

"""
Core information about the different clusters and modules is gathered in 
yml dictionaries. Any updates to modules or clusters can easily be modified
there and should be consistent through the tool. 
"""

#Smart parser for newlines
class SmartFormatter(HelpFormatter):
    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        return HelpFormatter._split_lines(self, text, width)

#Classes to store required parameters and paths for amber, gaussian and gromacs.
class Amber_variables:
    def __init__(self,
                 dirname,
                 outpath,
                 runtype,
                 parampath,
                 simtype,
                 inputfiles,
                 toppath,
                 restart,
                 plumed):
        self.path = dirname
        self.runtype = runtype       
        self.outputpath = outpath
        self.simtype= simtype
        self.inpath = parampath
        self.toppath = toppath
        self.restartflag = restart
        if self.runtype == "cpu":
            self.run_exec = f"mpirun -np {tasks} pmemd.MPI"
        else:
            self.run_exec = "pmemd.cuda"

        if plumed:
            self.plumed = plumed
            if not os.path.exists(f"{self.path}/{self.plumed}"):
                print(f"WARNING: file '{self.path}/{plumed}' not found. \nTerminating...")
                sys.exit(0)
#currently no explicit plumed functionality implemented (is there an Amber version with plumed installed?)
        for item in inputfiles:
            if item.split(".")[-1] in ["prmtop", "prm"]:
                if not os.path.exists(f"{self.toppath}/{item}"):
                    print(f"WARNING: file '{self.toppath}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.topology = item
            if item.split(".")[-1] in ["rst7", "crd", "inpcrd", "ncrst"]:
                if not os.path.exists(f"{self.path}/{item}"):
                    print(f"WARNING: file '{self.toppath}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.inputfile = item
        return

    def min(self):
        if not os.path.exists(f"{self.inpath}/min.in"):
            print(f"WARNING: file '{self.inpath}/min.in' not found. \nTerminating...")
            sys.exit(0)
        min_command = f"""
mkdir em; cd em
{self.run_exec} -O -i {self.inpath}/min.in -p {self.toppath}/{self.topology} -c {self.path}/{self.inputfile} -ref {self.path}/{self.inputfile} -o {self.outputpath}/min.log -r {self.outputpath}/min.rst7
cd ..
        """
        return min_command

    def eq(self):
        for eq in ["nvt", "npt"]:
            if not os.path.exists(f"{self.inpath}/{eq}.in"):
                print(f"WARNING: file '{self.inpath}/{eq}.in' not found. \nTerminating...")
                sys.exit(0)
        if "min" in self.simtype:
            self.path = "../em"
            self.inputfile = "em.rst7"
        eq_command = f"""
mkdir nvt; cd nvt
{self.run_exec} -O -i {self.inpath}/nvt.in -p {self.toppath}/{self.topology} -c {self.path}/{self.inputfile} -ref {self.path}/{self.inputfile} -o {self.outputpath}/nvt.log -r {self.outputpath}/nvt.rst7
cd ..\n
mkdir npt; cd npt
{self.run_exec} -O -i {self.inpath}/npt.in -p {self.toppath}/{self.topology} -c ../nvt/nvt.rst7 -ref ../nvt/nvt.rst7 -o {self.outputpath}/npt.log -r {self.outputpath}/npt.rst7
cd ..
        """
        return eq_command

    def prod(self):
        if not os.path.exists(f"{self.inpath}/md.in"):
            print(f"WARNING: file '{self.inpath}/md.in' not found. \nTerminating...")
            sys.exit(0)
        if "min" in self.simtype:
            self.path = "../em"
            self.inputfile = "em.rst"
        if "eq" in self.simtype:
            self.path = "../npt"
            self.inputfile = "npt.gro"
            self.checkpoint = "npt.cpt"
        prod_command = f"""
mkdir md; cd md
{self.run_exec} -O -i {self.inpath}/md.in -p {self.toppath}/{self.topology} -c {self.path}/{self.inputfile} -ref {self.path}/{self.inputfile} -o {self.outputpath}/md.log -r {self.outputpath}/md.rst7 -x md.nc 
cd ..
            """
        return prod_command

    def sed_commands(self):
        commands = []
        if self.restartflag:
            restart_command = shlex.split(f"sed -i 's/-O/-A/g' {homedir}/{args.jobname}.{extension}")
            commands.append(restart_command)
        for command in commands:
            process = subprocess.Popen(command)
        return

class Gromacs_variables:
    def __init__(self,
                 dirname,
                 outpath,
                 runtype,
                 simtype,
                 inputfiles,
                 toppath,
                 parampath,
                 restart,
                 plumed):
        self.path = dirname
        self.runtype = runtype
        self.outputpath = outpath
        self.simtype = simtype    
        self.mdppath = parampath
        self.toppath = toppath
        self.restartflag = restart
        if plumed:
            self.plumed = plumed
            if not os.path.exists(f"{self.path}/{self.plumed}"):
                print(f"WARNING: file '{self.path}/{plumed}' not found. \nTerminating...")
                sys.exit(0)
        for item in inputfiles:
            if item.split(".")[-1] == "top":
                if not os.path.exists(f"{self.toppath}/{item}"):
                    print(f"WARNING: file '{self.toppath}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.topology = item
            elif item.split(".")[-1] in ["gro", "pdb"]:
                if not os.path.exists(f"{self.path}/{item}"):
                    print(f"WARNING: file '{self.path}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.inputfile = item
            elif item.split(".")[-1] == "tpr":
                if not os.path.exists(f"{self.path}/{item}"):
                    print(f"WARNING: file '{self.path}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.tprfile = item.rsplit(".", 1)[0]
            elif item.split(".")[-1] == "ndx":
                if not os.path.exists(f"{self.toppath}/{item}"):
                    print(f"WARNING: file '{self.toppath}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.index = item
            elif item.split(".")[-1] == "cpt":
                if not os.path.exists(f"{self.path}/{item}"):
                    print(f"WARNING: file '{self.toppath}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.checkpoint = item
        if self.runtype == "cpu":
            self.run_exec = f"mpirun -np {tasks} gmx_mpi mdrun"
        else:
            self.run_exec = f"mpirun -np {tasks} gmx_mpi mdrun -ntomp {omp}"

        return
    def min(self):
        if not os.path.exists(f"{self.mdppath}/em.mdp"):
            print(f"WARNING: file '{self.mdppath}/em.mdp' not found. \nTerminating...")
            sys.exit(0)
        min_command = f"""
mkdir em; cd em
gmx grompp -f {self.mdppath}/em.mdp -c {self.path}/{self.inputfile} -p {self.toppath}/{self.topology} -o em.tpr
{self.run_exec} -deffnm em -pin on
cd ..
        """
        return min_command

    def eq(self):
        for eq in ["nvt", "npt"]:
            if not os.path.exists(f"{self.mdppath}/{eq}.mdp"):
                print(f"WARNING: file '{self.mdppath}/{eq}.mdp' not found. \nTerminating...")
                sys.exit(0)
        if "min" in self.simtype:
            self.path = "../em"
            self.inputfile = "em.gro"
        eq_command = f"""
mkdir nvt; cd nvt
gmx grompp -f {self.mdppath}/nvt.mdp -c {self.path}/{self.inputfile} -p {self.toppath}/{self.topology} -r {self.path}/{self.inputfile} -o nvt.tpr -maxwarn 2
{self.run_exec} -deffnm nvt -pin on
cd ..\n
mkdir npt; cd npt
gmx grompp -f {self.mdppath}/npt.mdp -c ../nvt/nvt.gro -p {self.toppath}/{self.topology} -r ../nvt/nvt.gro -o npt.tpr -maxwarn 2
{self.run_exec} -deffnm npt -pin on
cd ..\n
        """
        return eq_command

    def prod(self):
        if not os.path.exists(f"{self.mdppath}/md.mdp"):
            print(f"WARNING: file '{self.mdppath}/md.mdp' not found. \nTerminating...")
            sys.exit(0)
        if "min" in self.simtype:
            self.path = "../em"
            self.inputfile = "em.gro"
        if "eq" in self.simtype:
            self.path = "../npt"
            self.inputfile = "npt.gro"
            self.checkpoint = "npt.cpt"
        prod_command = f"""
mkdir md; cd md
gmx grompp -f {self.mdppath}/md.mdp -c {self.path}/{self.inputfile} -p {self.toppath}/{self.topology} -t {self.path}/{self.checkpoint} -r {self.path}/{self.inputfile} -o md.tpr -maxwarn 2
{self.run_exec} -deffnm md -pin on
cd ..\n
        """
        return prod_command

    def tpr(self):
        tpr_command = f"{self.run_exec} -s {self.path}/{self.tprfile}.tpr -deffnm {self.tprfile} -pin on"
        return tpr_command

    def restart(self):
        if not hasattr(self, "tprfile"):
            print("WARNING: tpr file is required along with chk file for restarting a simulation. \nTerminating...")
            sys.exit(0)
        if not hasattr(self, "checkpoint"):
            print("WARNING: chk file is required along with tpr file for restarting a simulation. \nTerminating...")
            sys.exit(0)
        restart_command = f"{self.run_exec} -deffnm {self.path}/{self.tprfile} -cpi {self.checkpoint} -append yes -pin on"
        return restart_command

    def sed_commands(self):
        commands = []
        if hasattr(self, "index"):
            index_command = shlex.split(f"sed -i 's|-r|-n\ {self.toppath}/{self.index}\ -r|g' {homedir}/{args.jobname}.{extension}")
            commands.append(index_command)
        if hasattr(self, "plumed"):
            plumed_command = shlex.split(f"sed -i 's|-t|-plumed\ {self.path}/{self.plumed}\ -t|g' {homedir}/{args.jobname}.{extension}")
            commands.append(plumed_command)
        for command in commands:
            process = subprocess.Popen(command)
        return

class Gaussian_variables:
    def __init__(self,
                 queue,
                 nodes,
                 path,
                 outpath,
                 runtype,
                 inputfiles,
                 restart):
        self.path = wd
        self.outputpath = outpath
        self.runtype = runtype  #for gpu compatibility. Not included yet
        self.procs = tasks
        self.mem = int(float(resources[queue]["memory"])*int(resources[queue]["cores"]))
        self.restartflag = restart
        self.nodes = nodes
        for item in inputfiles:
            if item.split(".")[-1] == "com":
                if not os.path.exists(f"{self.path}/{item}"):
                    print(f"WARNING: file '{self.path}/{item}' not found. \nTerminating...")
                    sys.exit(0)
                self.filename = item
                self.filepath = os.path.join(self.path, self.filename)
            elif item.split(".")[-1] == "chk":
                self.checkpoint = item
        return

    def gaussian_commands(self):
        commands = []
        change_proc = f"sed -i 's/[cC][pP][uU]=.*/CPU=0-{int((self.procs/self.nodes)-1)}/g ; s/[nN][pP][rR][oO][cC]=.*/CPU=0-{in((self.procs/self.nodes)-1)}/g' {self.path}/${{job_name}}.com"
        if hasattr(self, "checkpoint"):
            change_chk_path = f"sed -i '/[cC][hH][kK]=.*/d' {self.path}/${{job_name}}.com; sed -i '/CPU=0/a %Chk={self.outputpath}\/{self.checkpoint}' {self.path}/${{job_name}}.com"
        else:
            change_chk_path = ""
        if self.restartflag and hasattr(self, "checkpoint"):
            set_restart = f"sed -i 's/[oO][pP][tT]=(/opt=(restart,/' {self.path}/${{job_name}}.com"
#only restarts geometry optimization
        else:
            set_restart = ""
        if self.runtype == "gpu":
            if self.nodes == 1:
                set_gpus = f"sed -i '/[gG][pP][uU][cC][pP][uU]/d' {self.path}/${{job_name}}.com; sed -i '/%[cC][pP][uU]=/a %GPUCPU=0=0' {self.path}/${{job_name}}.com"
            elif self.nodes == 2:
                set_gpus = f"sed -i '/[gG][pP][uU][cC][pP][uU]/d' {self.path}/${{job_name}}.com; sed -i '/%[cC][pP][uU]=/a %GPUCPU=0,1=0,16' {self.path}/${{job_name}}.com"
        else:
            set_gpus = ""
        change_mem = f"sed -i '/[mM][eE][mM]=.*[Bb]/d' {self.path}/${{job_name}}.com; sed -i '/CPU=0/a %mem={self.mem}GB' {self.path}/${{job_name}}.com"
        commands = f"""
job_name='{self.filename.rsplit(".",1)[0]}'\n
export GAUSS_SCRDIR=$VSC_SCRATCH_NODE
{set_restart}
{change_proc}
{change_chk_path}
{set_gpus}
{change_mem}
g16 {self.path}/${{job_name}}.com\n
        """
        return commands

def queue_type(queue, nodes, openMP, memory, walltime, jobname, project_name, mail, reportpath):
    """Check queue type and create appropriate header for submission script."""
    if queue == "leibniz_slurm": #leibniz partition, 148 nodes, 28 cores per node, 128 GB RAM, max. job time of 3 days, Slurm
        header = f"""#!/bin/bash\n
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task={omp}
#SBATCH --time={walltime}
#SBATCH --job-name={jobname}
#SBATCH --output={reportpath}/{jobname}.out
#SBATCH --error={reportpath}/{jobname}.err
        """
        if mail == True:
            header += "\n#SBATCH --mail-type=BEGIN,FAIL,END\n"
    elif queue == "leibniz_pbs": #leibniz partition, completely transfered to Slurm.
        header = f"""#!/bin/bash\n
#PBS -L tasks={tasks}:lprocs={omp}:memory={memory}gb
#PBS -W x=template:singleswitch
#PBS -l walltime={walltime}
#PBS -N {jobname}
#PBS -o {reportpath}/{jobname}.out
#PBS -e {reportpath}/{jobname}.err
#PBS -q leibniz
        """
        if mail == True:
            header += "\n#PBS -m abe\n"
    elif queue == "vaughan": #vaughan partition, 152 nodes, 64 cores per node, 256 GB RAM, max. job time of 3 days, Slurm 
        header = f"""#!/bin/bash\n
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task={omp}
#SBATCH --time={walltime}
#SBATCH --job-name={jobname}
#SBATCH --output={reportpath}/{jobname}.out
#SBATCH --error={reportpath}/{jobname}.err
        """
        if mail == True:
            header += "\n#SBATCH --mail-type=BEGIN,FAIL,END\n"
    elif queue == "hopper": #leibniz partition, 23 nodes, 20 cores per node, 256 GB RAM, max. job time of 7 days, Slurm
        header = f"""#!/bin/bash\n
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task={omp}
#SBATCH --time={walltime}
#SBATCH --job-name={jobname}
#SBATCH --output={reportpath}/{jobname}.out
#SBATCH --error={reportpath}/{jobname}.err
#SBATCH --partition ivybridge
        """
        if mail == True:
            header += "\n#SBATCH --mail-type=BEGIN,FAIL,END\n"
    elif queue == "breniac": #Tier-1 cluster, 409 nodes, 28 cores per node, RAM varies from 128 to 256 GB, max job time of 3 days, Torque
        header = f"""#!bin/bash\n
#PBS -l nodes={nodes}:ppn={tasks_per_node}:singleisland
#PBS -l pmem=4gb
#PBS -l walltime={walltime}
#PBS -N {jobname}
#PBS -o {reportpath}/{jobname}.out
#PBS -e {reportpath}/{jobname}.err
#PBS -A {project_name}
        """
        if mail == True:
            header += "\n#PBS -m abe\n"
#GPU moved to slurm on leibniz, also gpu's on vaughan
    elif queue in ["ampere_gpu", "pascal_gpu"]: 
        header = f"""#!/bin/bash\n
#SBATCH --nodes={nodes}
#SBATCH --ntasks-per-node={tasks_per_node}
#SBATCH --cpus-per-task={omp}
#SBATCH --time={walltime}
#SBATCH --job-name={jobname}
#SBATCH --output={reportpath}/{jobname}.out
#SBATCH --error={reportpath}/{jobname}.err
#SBATCH --partition {queue}
#SBATCH --gres=gpu:{nodes}
        """
        if mail == True:
            header += "\n#SBATCH --mail-type=BEGIN,FAIL,END\n"
    return header

def load_modules(queue, software, runtype, plumed):
    """Load required modules for requested task, cluster and resources"""
    module_load = "\n"
    if queue == "hopper":
        module_load += f"module load {modules['hopper_modules']['general']}\n"
        if plumed:
            module_load += f"module load {modules['hopper_modules'][software]['plumed']}\n"
        else:
            module_load += f"module load {modules['hopper_modules'][software][runtype]}\n"
    elif queue == "breniac":
        module_load += f"module load {modules['breniac_modules'][software][runtype]}\n"
    else:
        module_load += f"module load {modules['modules']['general']}\n"
        if plumed:
            module_load += f"module load {modules['modules'][software]['plumed']}\n"
        else:
            module_load += f"module load {modules['modules'][software][runtype]}\n"
# specific modules to load for gpu usage
#        if runtype == "gpu":
#            module_load += f"module load {modules['modules'][runtype]}\n"
        module_load += "\n"
    return module_load

def get_commands(software, queue, nodes, path, outpath, runtype,
        simtype, inputfiles, toppath, parampath, restart, plumed):
    """Write the execution commands for requested task and resources."""
    commands = f"\nmkdir -p {outpath}; cd {outpath}\n"
    if software == "gaussian":
        gaussian = Gaussian_variables(queue, nodes, path, outpath, runtype, inputfiles, restart)
        commands += gaussian.gaussian_commands()
    elif software == "gromacs":
        gromacs = Gromacs_variables(path, outpath, runtype, simtype, inputfiles, toppath, parampath, restart, plumed)
        if gromacs.restartflag:
            commands += gromacs.restart()
            return commands
        if hasattr(gromacs, "tprfile"):
            commands += gromacs.tpr()
            return commands
        for sim in simtype:
            commands += getattr(gromacs, sim)()
    elif software == "amber":
        amber = Amber_variables(path, outpath, runtype, parampath, simtype, inputfiles, toppath, restart, plumed)
        for sim in simtype:
            commands += getattr(amber, sim)()
    return commands

def write_jobscript(path, jobname, queue, software, runtype, nodes):
    """Combine strings into submission script."""
    print(f"writing script for submitting {software} job using {runtype}s on {nodes} {queue} node(s)... ('{jobname}.{extension}')")
    script_content = header + modules + commands
    with open(f"{homedir}/{jobname}.{extension}", "w") as file:
        file.write(script_content)

def mod_jobscript(software, path, outpath, runtype, simtype, inputfiles, toppath, parampath, restart, plumed): 
    """Modify job submission script with sed commands if required."""
    if software == "gromacs":
        gromacs = Gromacs_variables(path, outpath, runtype, simtype, inputfiles, toppath, parampath, restart, plumed)
        gromacs.sed_commands()
    if software == "amber":
        amber = Amber_variables(path, outpath, runtype, parampath, simtype, inputfiles, toppath, restart, plumed)
        amber.sed_commands()
    return

def submit_job(queue, jobname, nosubmit, keep):
    """Submit job to queue if submit flag is on. Keep or delete after submission."""
    if nosubmit:
        return
    if queue in ["leibniz_pbs", "breniac"]:
        job = f"{jobname}.pbs"
        submit = shlex.split(f"qsub {job}")
    else:
        job = f"{jobname}.sh"
        if queue in ["pascal_gpu", "ampere_gpu"]:
            submit = shlex.split(f"sbatch -p {queue} {homedir}/{job}")
        else:
            submit = shlex.split(f"sbatch {homedir}/{job}")
    process = subprocess.Popen(submit)
    print_partition_info(queue)
    if keep:
        sys.exit()
    time.sleep(1)
    os.remove(f"{homedir}/{job}")
    return 

def print_partition_info(queue):
    if queue in ["pascal_gpu", "leibniz_slurm", "leibniz_pbs"]:
        queue = "leibniz"
    elif queue == "ampere_gpu":
        queue = "vaughan"        
    print("job submitted.")
    print('Queue state:\nPARTITION\tAVAIL\t     NODES\t   STATE')
    if queue in ["breniac", "leibniz_pbs"]:
        command1 = "showq"
        command2 = shlex.split("grep ' active'")
        p1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
        subprocess.Popen(command2, stdin=p1.stdout)
    else:    
        with open(os.path.join(script_dir,"partitions.yml"), "r") as file:
            partitions = yaml.safe_load(file)
            for partition in partitions[queue]:
                command1 = shlex.split(f"sinfo -p {partition}")
                command2 = shlex.split(f"grep {partition}")
                command3 = shlex.split("awk '{print $1, $2, $4, $5}'")
                p1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
                p2 = subprocess.Popen(command2, stdin=p1.stdout, stdout=subprocess.PIPE)
                p3 = subprocess.Popen(command3, stdin=p2.stdout, stdout=subprocess.PIPE)
                output = p3.stdout.read().decode("UTF-8")
                print(output.replace(" ", "\t".expandtabs(12)))
    return

if __name__ == "__main__":
    basedir = os.environ["VSC_SCRATCH"]
    homedir = os.environ["VSC_HOME"]
    script_dir = os.path.abspath(os.path.dirname(__file__))
    with open((os.path.join(script_dir, "modules.yml")), "r") as file:
        modules = yaml.safe_load(file)
    with open((os.path.join(script_dir, "resources.yml")), "r") as file:
        resources = yaml.safe_load(file)
    parser = ArgumentParser(formatter_class=SmartFormatter)
    requiredNamed = parser.add_argument_group('required arguments')                        
    requiredNamed.add_argument("-q", "--queue",
                        dest="queue",
                        type=str,
                        help=("The name of the queue to which to submit the job. "
                        "Options are 'leibniz_pbs', 'leibniz_slurm', 'vaughan', "
                        "'hopper', 'breniac', 'pascal_gpu', 'ampere_gpu', 'arcturus_gpu'."),
                        choices=["leibniz_pbs",
                                 "leibniz_slurm",
                                 "vaughan",
                                 "hopper",
                                 "breniac",
                                 "pascal_gpu",
                                 "ampere_gpu"
                                 "arcturus_gpu"],
                        metavar="QUEUE")
    requiredNamed.add_argument("-p", "--path",
                        dest="dirname",
                        type=str,
                        default=basedir,
                        help=("The path to the directory with input files." 
                        " Looks for the directory in the scratch folder." 
                        " If not provided, scratch folder is used."),
                        metavar="PATH")
    parser.add_argument("-po", "--path_out",
                        dest="outpath",
                        type=str,
                        help=("The path to where output is written. Default"
                        " is the same as input path. "
                        "Looks for the directory in the scratch folder."),
                        metavar="OUTPATH")
    parser.add_argument("-pp", "--path_parameters",
                        dest="parampath",
                        type=str,
                        default=os.path.join(homedir, "parameter_files"),
                        help=("The path to where the simulation parameters are. Default"
                        " is the 'parameter_files' folder within the folder of this script. Option 'input' " 
                        "uses same directory as input files. Only required for gromacs and amber jobs."),
                        metavar="PARAMPATH")
    parser.add_argument("-pt", "--path_topology",
                        dest="toppath",
                        type=str,
                        help=("The path to where the topology file (and optionally the index file) is. "
                        "Default is the folder with input files. Only required for "
                        "gromacs and amber jobs."),
                        metavar="TOPPATH")
    parser.add_argument("-pr", "--path_reports",
                        dest="reportpath",
                        type=str,
                        default=os.path.join(homedir, "reports"),
                        help=("The path to where the stderr and stout are written. Default"
                        " is the 'reports' folder in the home directory."),
                        metavar="REPORTPATH")
    requiredNamed.add_argument("-i", "--input",
                        dest="input",
                        type=str,
                        nargs="+", 
                        help=("R|Names of the input files.\nRequired for gaussian: "
                        ".com file\nOptional for gaussian: .chk file for restart or "
                        "checkpoint writing.\nRequired for gromacs: .gro/.pdb file "
                        "and .top file or a .tpr file.\nOptional for gromacs: .ndx "
                        "file, .cpt file (if restart or if only production run)\n"
                        "Required for amber: .rst7/.crd/.ncrst/.inpcrd, .prm/.prmtop"),
                        metavar="INPUT")
    parser.add_argument("-o", "--output",
                        dest="jobname",
                        type=str,
                        default=None,
                        help=("The job name. Will also be used as the name of "
                        "the submission script."),
                        metavar="OUTPUT")
#    parser.add_argument("-r", "--runtype",
#                        dest="runtype",
#                        type=str,
#                        default="cpu",
#                        help="Run type (cpu or gpu).",
#                        choices=["cpu", "gpu"],
#                        metavar="RUNTYPE")
    requiredNamed.add_argument("-x", "--software",
                        dest="software",
                        type=str,
                        default="gromacs",
                        help="Software to run (gromacs, amber or gaussian).",
                        choices=["gromacs", "amber", "gaussian"],
                        metavar="SOFTWARE")
    parser.add_argument("-n", "--nodes",
                        dest="nodes",
                        type=int,
                        default=1,
                        help="The number of nodes to use for the job. Default is 1.",
                        metavar="NODES")
    parser.add_argument("-m", "--memory",
                        dest="memory",
                        type=str,
                        default="default",
                        help=("Amount of memory used for the job in GB (default uses "
                        "most of the available memory)."),
                        metavar="MEMORY")
    parser.add_argument("-wt", "--walltime",
                        dest="walltime",
                        type=str,
                        default="24:00:00",
                        help="Wall time in format HH:MM:SS. Default is 24h. 'max' "
                        "gives the maximum wall time for the cluster.",
                        metavar="WALLTIME")
    parser.add_argument("-omp", "--openMP",
                        dest="openMP",
                        type=int,
                        default=1,
                        help=("Defines the number of openMP threads to be used per "
                        "process. Default is 1 (no openMP parallelisation)."),
                        metavar="OPENMP")
    parser.add_argument("-s", "--simtype",
                        dest="simtype",
                        type=str,
                        default=["min", "eq", "prod"],
                        nargs="+",
                        help=("Provide for amber or gromacs simulation. options are "
                        "'min', 'eq', 'prod' (multiple are possible). Default is all."),
                        metavar="SIMTYPE")
    parser.add_argument("-A", "--project_name",
                        dest="project_name",
                        type=str,
                        default=resources["breniac"]["project_name"],
                        help="Project name for tier 1 submission.",
                        metavar="PROJECT_NAME")

    parser.add_argument("-restart",
                        dest="restart",
                        default=False,
                        action='store_true',
                        help=("If true, a process will be restarted from a previous "
                        "job. For Gaussian, a checkpoint file should be present."))
    parser.add_argument("-plumed",
                        dest="plumed",
                        default=False,
                        help=("If true, plumed is used in the job. The file name of "
                        "the settings file should be 'plumed.dat'."))
    parser.add_argument("-mail",
                        dest="mail",
                        default=False,
                        action='store_true',
                        help="If true, sends an email when job starts, ends or fails.")
    parser.add_argument("-keep",
                        dest="keep",
                        default=False,
                        action='store_true',
                        help="If true, don't remove job script after submission.")
    parser.add_argument("-nosubmit",
                        dest="nosubmit",
                        default=False,
                        action='store_true',
                        help="If true, don't submit the job to the queue.")

    args = parser.parse_args()
#Check for default settings
    simtype = list(args.simtype)
    tasks = int(resources[args.queue]["cores"]/int(args.openMP)*int(args.nodes))
    tasks_per_node = int(resources[args.queue]["cores"]/int(args.openMP))
    omp = args.openMP
    if args.queue in ["leibniz_pbs", "breniac"]:
        extension = "pbs"
    else:
        extension = "sh"
    if args.memory == "default":
        memory = resources[args.queue]["memory"]
    if args.walltime == "max":
        args.walltime = resources[args.queue]["time"]
#set paths
    if args.dirname is not basedir:
        wd = os.path.join(basedir, args.dirname)
    else:
        wd = basedir
    if args.outpath is not None:
        args.outpath = os.path.join(basedir, args.outpath)
    else:
        args.outpath = wd
    if args.parampath == "input":
        args.parampath = wd
    elif args.parampath == os.path.join(homedir, "parameter_files"):
        pass 
    else:
        args.parampath = os.path.join(basedir, args.parampath)
    if args.toppath is not None:
        args.toppath = os.path.join(basedir, args.toppath)
    else:
        args.toppath = wd

    if args.queue in ["pascal_gpu", "ampere_gpu"]:
        args.runtype = "gpu"
    else:
        args.runtype = "cpu"
    if args.jobname is None:
        if args.software == "gaussian":
            for item in args.input:
                if item.split(".")[-1] == "com":
                    args.jobname = item.rsplit(".", 1)[0]
        elif args.software == "gromacs":
            args.jobname = "gromacs_job"
        elif args.software == "amber":
            args.jobname = "amber_job"

#Run functions
    if args.queue in ["pascal_gpu", "ampere_gpu"]:
        if args.software == "gromacs":
            print("WARNING: Gromacs support for GPU calculations is currently not available. \nTerminating...")
            sys.exit(0)
        if args.software == "amber":
            pass
        if args.software == "gaussian":
            print("WARNING: GPU calculations are not possible with curren NVIDIA Tesla A100/P100 GPUs (Vaughan) because they do not have enough memory. \nContinuing...")
    if args.software == "amber":
        if args.plumed:
            print("WARNING: Amber support for plumed calculations is currently not available. \nTerminating...")
            sys.exit(0)
        if args.openMP != 1:
            print("WARNING: Amber does not support the use of hybrid MPI/openMP. \nContinuing...")

    header = queue_type(args.queue, args.nodes, args.openMP, args.memory, args.walltime, args.jobname, args.project_name, args.mail, args.reportpath)
    modules = load_modules(args.queue, args.software, args.runtype, args.plumed)
    try:
        commands = get_commands(args.software, args.queue, args.nodes, wd, args.outpath, args.runtype, simtype, args.input, args.toppath, args.parampath, args.restart, args.plumed)
    except AttributeError:
        print("WARNING: Standard Gromacs job require at least a .gro and a .top file as input. If "
        "only doing a production run, a .cpt file is also required. If doing another type of "
        "calculation, not all required input files were provided. \nTerminating...")
        sys.exit(0)
    write_jobscript(wd, args.jobname, args.queue, args.software, args.runtype, args.nodes)
    mod_jobscript(args.software, wd, args.outpath, args.runtype, simtype, args.input, args.toppath, args.parampath, args.restart, args.plumed)
    submit_job(args.queue, args.jobname, args.nosubmit, args.keep)
