setenv PATH " ~/nc/bin ~/bin $PATH"
set   path=( ~/nc/bin ~/bin $path)
alias h history
set history=30
set shell=/bin/csh
alias t textmod
alias lc ls -FC
alias l ls -l
alias hn hostname
alias rs 'eval `resize`'
alias psg 'ps -ef | grep \!* | grep -v grep'
if ($?COLUMNS == 0) exit
stty intr ^c
stty erase ^h

