#!/bin/bash

source /opt/exp_soft/llr/root/v6.06.00-el6-gcc48/etc/init.sh

cd ${pwd_}

# check what we're doing
echo '###############################################'
echo '#  __        __   _                           #'
echo '#  \ \      / /__| | ___ ___  _ __ ___   ___  #'
echo '#   \ \ /\ / / _ \ |/ __/ _ \| '\''_ ` _ \ / _ \ #'
echo '#    \ V  V /  __/ | (_| (_) | | | | | |  __/ #'
echo '#     \_/\_/ \___|_|\___\___/|_| |_| |_|\___| #'
echo '###############################################'
                                                   
echo ""
                                        
echo "Hello, my name is" `uname -n`
echo "The date is" `date`
echo "My uptime is" `uptime`
echo "I am now in" $PWD
echo ""
echo "The environment variables are the following:"
env

echo
echo '###############################################'
echo '#      _       _           _             _    #'
echo '#     | | ___ | |__    ___| |_ __ _ _ __| |_  #'
echo '#  _  | |/ _ \| '\''_ \  / __| __/ _` | '\''__| __| #'
echo '# | |_| | (_) | |_) | \__ \ || (_| | |  | |_  #'
echo '#  \___/ \___/|_.__/  |___/\__\__,_|_|   \__| #'
echo '###############################################'
                                                      
echo root -l -b -q runLimit_RaaNS_Workspace.C+'("'${it}'", "RFrac2Svs1S_PbPbvsPP", "simPdf_syst", "workspace", "dOS_DATA", "'${ACTag}'", '${CL}', '${calculatorType}', '${testStatType}', '${useCLs}')'
time root -l -b -q runLimit_RaaNS_Workspace.C+'("'${it}'", "RFrac2Svs1S_PbPbvsPP", "simPdf_syst", "workspace", "dOS_DATA", "'${ACTag}'", '${CL}', '${calculatorType}', '${testStatType}', '${useCLs}')'

echo "ok, I'm done"
