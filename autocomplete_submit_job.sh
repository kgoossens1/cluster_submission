_submit_compl() {
    local cur prev opts1 opts2
    COMREPLY=()
    cur="${COMP_WORDS[COMP_CWORD]}"
    prev="${COMP_WORDS[COMP_CWORD-1]}"
    opts1="-q -p -i -x -h -po -pp -pt -pr -o -n -A -m -wt -omp -s -restart -plumed -mail -keep -nosubmit"
    opts2="vaughan leibniz_pbs leibniz_slurm hopper breniac pascal_gpu ampere_gpu"
    opts3="gromacs amber gaussian"
    opts4="min eq prod"
    opts5="-q -p -i -x -h -po -pp -pt -pr -o -n -A -m -wt -omp -s -restart -plumed -mail -keep -nosubmit min eq prod"
    if [[ ${COMP_CWORD} -eq 1 ]]; then
        COMPREPLY=( $(compgen -W "${opts1}" -- ${cur}))
    elif [[ ${prev} == -q ]] ; then
        COMPREPLY=( $(compgen -W "${opts2}" -- ${cur}))
    elif [[ ${prev} == -x ]] ; then
	COMPREPLY=( $(compgen -W "${opts3}" -- ${cur}))
    elif [[ ${prev} == -s ]] ; then
	COMPREPLY=( $(compgen -W "${opts4}" -- ${cur}))
    elif [[ ${prev} == "min" || ${prev} == "eq" ]] ; then
	COMPREPLY=( $(compgen -W "${opts5}" -- ${cur}))
    elif [[ ${prev} == -p* ]] ; then
        COMPREPLY=($(compgen -o plusdirs -f -X '!*/' -- ${cur}))
    elif [[ ${prev} == -o ]] ; then
        COMPREPLY=($(compgen -o plusdirs -f -X '!*/' -- ${cur}))
    else
	COMPREPLY=( $(compgen -W "${opts1}" -- ${cur}))
    fi
}

complete -o filenames -F _submit_compl submit_job

