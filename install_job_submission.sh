#!/bin/bash

echo -e "Installing job submission script...\n"
mkdir -p $VSC_HOME/scripts
mv ./submit_job.py $VSC_HOME/scripts/
mv ./autocomplete_submit_job.sh $VSC_HOME/scripts/
source $VSC_HOME/scripts/autocomplete_submit_job.sh

shopt -s nocasematch
read -p "Do you want the alias 'submit_job' to be created automatically?: " create_alias

if [[ $create_alias == y || $create_alias == "yes" ]] ; then
    echo "Creating alias..."
    echo "alias submit_job='python $VSC_HOME/scripts/submit_job.py'" >> $VSC_HOME/.bashrc
    echo "source '$VSC_HOME/scripts/autocomplete_submit_job.sh'" >> $VSC_HOME/.bashrc
else
    echo "Not creating an alias."
fi
sleep 2

read -p $'About to create directory "$VSC_HOME/reports" where job reports are written by default.\x0aThe directory can also be chosen using command line option "-pr" during execution.\x0aDo you wish to use a different default directory?(yes/no)\x0a' alt_reports
if [[ $alt_reports == y || $alt_reports == "yes" ]] ; then
    read -p "Please input the full path to the directory for writing job reports: " reports_dir
    sed -i "s|os.path.join(homedir, 'reports')|'${reports_dir}'|" $VSC_HOME/scripts/submit_job.py
    echo -e "Report path changed to '$reports_dir'.\n"
else
    mkdir -p $VSC_HOME/reports
fi
sleep 2

echo -e "Creating directory '$VSC_HOME/parameter_files' where simulation input files are looked for by default for Gromacs and Amber simulations.\nInput files should have names 'min', 'nvt', 'npt' and 'md' with *in or *mdp extensions. Files can be found in another directory using command line option '-pp' during execution..."
mkdir -p $VSC_HOME/parameter_files
sleep 2

read -p "Do you have an active project on the Breniac cluster? If yes, answer with the project name: " prj_name
    if [[ $prj_name == n || $prj_name == "no" ]] ; then
        echo "No project added. If you have an active project in the future, you can provide the name in the future through the command line option '-A' during execution."
    else
        sed -i "s/prj_name/${prj_name}/g" $VSC_HOME/scripts/submit_job.py
        echo -e "The project name has been implemented. This project name will automatically be used when submitting jobs to the Breniac cluster.\nYou can also provide a different name in the future through the command line option '-A' during execution."
    fi
sleep 2

echo "Installation complete."
sleep 2

read -p "Do you want to remove this installing script? " remove_script
if [[ $remove_script == y || $remove_script == "yes" ]] ; then
    echo "Removing install script..."
    chmod 777 .git/objects/*/*
    rm -r *
fi
