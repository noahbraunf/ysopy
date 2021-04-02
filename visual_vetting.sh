#!/bin/bash

# INSTRUCTIONS + POINTERS
#sudo apt-get install fim # if fim is not already installed
#chmod u+x PATH_TO_THIS FILE 
# Execute that ^^^ prior to running this for the first time
# Run this file by pasting it's path into terminal session
# (when you're in the directory with the plots, that is)

# ACTUAL CODE
echo Enter report file path: 
read report_file
exec > $report_file
for i in $(ls); do
    fim $i
    read varname
    echo $i,$varname
    done
sed -i '1s/^/filename,class\n/' $report_file 
# ^^^ Adds column labels so you can read file into pandas df later
