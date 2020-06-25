#!/bin/bash

function check_var() {
    if [[ $2 != "" ]]
    then
        echo $1 defined to $2
    else
        echo $1 not defined !!!
    fi
}

GITUSERNAME="Sevin Arnaud"
GITEMAIL="arnaud.sevin@obspm.fr"

check_var GITUSERNAME $GITUSERNAME
check_var GITEMAIL $GITEMAIL

read  -n 1 -p "Are you really sure to be $GITUSERNAME <$GITEMAIL> ? [n]/y" ANS
if [ "A$ANS" != "Ay" ]
then
    echo "Installation cancelled"
    exit 0
fi
echo

cd /home/$USER

FILE="/home/$USER/git-completion.bash"

if [ ! -f $FILE ]; then
    echo "downloading git-completion.bash and source it in ~/.bashrc"
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

git config --global alias.assume-unchanged 'update-index --assume-unchanged'
git config --global alias.no-assume-unchanged 'update-index --no-assume-unchanged'

echo "configuration done"
