#!/bin/bash

GITUSERNAME="fvidal"
GITEMAIL="fabrice.vidal@obspm.fr"

cd /home/$USER

FILE="/home/$USER/git-completion.bash"

if [ ! -f $FILE ]; then
    wget https://raw.githubusercontent.com/git/git/master/contrib/completion/git-completion.bash
    echo 'source ~/git-completion.bash' >> ~/.bashrc
fi

git config --global user.name $GITUSERNAME
git config --global user.email $GITEMAIL

git config --global credential.helper "cache --timeout=3600"

git config --global alias.co checkout
git config --global alias.br branch
git config --global alias.ci commit
git config --global alias.st status

git config --global alias.unstage 'reset HEAD --'
git config --global alias.last 'log -1 HEAD'

git config --global color.ui true

git config --global branch.autosetuprebase always

git config --global alias.laststat '!f() { git diff --stat $1@{1}..$1 --dirstat=cumulative,files; }; f'
