PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
export PATH="/usr/ds/bin:$PATH"
export PATH="$PATH:/usr/ds/anaconda3/bin"
export PATH="$PATH:/bin:/usr/bin"

# increase length of qstat
export SGE_LONG_QNAMES=40

# parse git branch
parse_git_branch() {
     git branch 2> /dev/null | sed -e '/^[^*]/d' -e 's/* \(.*\)/ (\1)/'
}

# set a fancy prompt (non-color, unless we know we "want" color)
case "$TERM" in
    xterm-color|*-256color) color_prompt=yes;;
esac

# uncomment for a colored prompt, if the terminal has the capability; turned
# off by default to not distract the user: the focus in a terminal window
# should be on the output of commands, not on the prompt
#force_color_prompt=yes

if [ -n "$force_color_prompt" ]; then
    if [ -x /usr/bin/tput ] && tput setaf 1 >&/dev/null; then
        # We have color support; assume it's compliant with Ecma-48
        # (ISO/IEC-6429). (Lack of such support is extremely rare, and such
        # a case would tend to support setf rather than setaf.)
        color_prompt=yes
    else
        color_prompt=
    fi
fi

if [ "$color_prompt" = yes ]; then
    PS1="${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\[\033[31m\]\$(parse_git_branch)\[\033[00m\]\$ "
else
    PS1="${debian_chroot:+($debian_chroot)}\u@\h:\w\$ "
fi

# if [ "$color_prompt" = yes ]; then
#     PS1='${debian_chroot:+($debian_chroot)}\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ '
# else
#     PS1='${debian_chroot:+($debian_chroot)}\u@\h:\w\$ '
# fi
unset color_prompt force_color_prompt







file="/core/uge/LFPN/common/settings.sh"
if test -f "$file"; then
    source "$file"
fi


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/lfleddermann/mySoftware/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/lfleddermann/mySoftware/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/lfleddermann/mySoftware/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/lfleddermann/mySoftware/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<



scpToMpi(){
    # This is necessary if ~ was expaned.
    local remote_file="${2/home/lfleddermann}"
    command scp -r -o 'ProxyJump lfleddermann@10.218.100.51' "$1" lfleddermann@sirona01:"$remote_file"
}


scpFromMpi(){
    # replace /local/home with /home/mpi-username/. This is necessary if ~ was expaned.
    local remote_file="${1/home/lfleddermann}"
    command scp -r -o 'ProxyJump lfleddermann@10.218.100.51' lfleddermann@sirona01:"$remote_file" "$2"
}

rsyncToCaosDB(){
    # replace /local/home with /home/caosdb-username/. This is necessary if ~ was expaned.
    local remote_file="${2/home/caosdb-lfleddermann}"
    command rsync -av --progress -rltpPuzh --append "$1" caosdb-lfleddermann@134.76.17.187:"$remote_file"
}

scpToGwdg(){
    # First is From, Second is To
    local remote_file="${2}"
    echo $remote_file 
    echo ${1}
    echo ${2}
    echo scp -rp ${1} lfledde@transfer-scc.gwdg.de:$remote_file
    command scp -rp ${1} lfledde@transfer-scc.gwdg.de:$remote_file
}



#keep directory when using ssh to sirona:
# from local sshcd sironaXX will change computer to
# sirona and open sirona in the directory you were
# in on the local computer
sshcd() { ssh -A -t "$@" "cd '${PWD}'; bash -l"; }
# aliases
# change to sirona and stay in the same directory
alias camulus001='sshcd camulus001'
alias gwdglogin='ssh login-mdc.hpc.gwdg.de -l lfledde -i ~/.ssh/gwdg-id -J lfledde@login.gwdg.de'
alias mpilogin='ssh lfleddermann@10.218.100.51'


# helper
alias open='xdg-open'
alias DiskUsage="df -h | sed -n '1p;/$USER/p'"

# alias to qstat:
alias bashrc='vim ~/.bashrc'
alias Qfree='qstat  -g c' # Commented out: only for grannus queues -q "grannus.q,mvapich2-grannus.q" -g c'
alias Qall="qstat -u '*'"
alias QjobsGrannusparallelQ='qstat -pe mvapich2-teutates.q -u "*"'
alias Qothers='qstat -u "*" | grep -v $USER' # Commented out only for grannus queues: -q "grannus.q,mvapich2-grannus.q" 
alias Qlistallqs="qstat -g c"
alias QnumMyjobs="qstat | grep ' r ' | wc -l"

# Path handeling
alias ".."="cd ..; "
alias cdHybrid="cd /home/lfleddermann/NumericalExperiments/Gitlab/HybridReservoir"
alias cdbmp="cd /data.bmp/lfleddermann/"
alias cdlorenz="cd /data.bmp/lfleddermann/DataAnalysis/2022_Lorenz63Reservoir/"
alias cdhind="cd /data.bmp/lfleddermann/DataAnalysis/2022_HindmarshRoseReservoir/"
alias cdcovid="cd /data.bmp/lfleddermann/DataAnalysis/2023_Paper_SIRS/"
alias cdres="cd /data.bmp/lfleddermann/DataAnalysis/2023_Paper_ReservoirDifferentTimescales"
alias cdresdr="cd /data.bmp/lfleddermann/DataAnalysis/2023_Paper_Reservoir_DimensionReduction/DRRC"
alias cdMyCode="cd /data.bmp/lfleddermann/DataAnalysis/MyCode/"

# git 
alias gitlog="git log --graph --oneline --all"
