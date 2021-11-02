export A='/Users/tylermasthay/Documents/Documents/UT/Engquist/wasserstein/Wasserstein/src/1D/forward/scalar'

export MAD_COMPILE='/Users/tylermasthay/Documents/Documents/UT/Year4/Fomel/CWasserstein/compile'

export SCHOOL='/Users/tylermasthay/Documents/Documents/UT'

export FOMEL='/Users/tylermasthay/Documents/Documents/UT/Year4/Fomel/CWasserstein'

export MAD='/Users/tylermasthay/Documents/Documents/UT/Year4/Fomel/madagascar'

export TMP='/var/tmp/Tasmanian-Devils'
export PROJ="$FOMEL/project/src"
export PROJTMP="/var/tmp/Tasmanian-Devils/project/src"

export C_INCLUDE_PATH=$MAD/include
export LIBRARY_PATH=$MAD/lib

source /Users/tylermasthay/Documents/Documents/UT/Year4/Fomel/madagascar/share/madagascar/etc/env.sh
source $MAD_COMPILE/compile_aliases

PATH=$PATH:$(head -1 .conda/environments.txt)/bin


# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/Users/tylermasthay/opt/anaconda3/bin/conda' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/Users/tylermasthay/opt/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/Users/tylermasthay/opt/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/Users/tylermasthay/opt/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

alias python='python3'
fpath=(~/.zfunc $path)
autoload -U $fpath[1]/*(.:t)

#git aliases
alias ga='git add'
alias gau='git add -u'
alias gst='git status'
alias gsl='git status -uno'
alias gp='git push'
alias gc='git commit'
alias gcm='git commit -m'
alias ge='gau; gcm $1; gp'
alias gb='git branch'

#bash alias
alias bash=bash-5.1.8
