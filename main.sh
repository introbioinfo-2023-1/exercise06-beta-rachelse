#!/bin/bash
# File: main.sh                                                                #
# Description:
#     This script contains exercise01 in "Introduction to Bioinformatics" lecture
#     Do not edit this script
################################################################################

echo "--------------------------------------------------"
# Execute commands
for command in ./command/*.sh
do
    bash $command
echo "--------------------------------------------------"
done

################################################################################