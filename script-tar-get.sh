#!/bin/bash
FC=/home/ibrik/fc/
ssh -l ibrik lazaro ls -Cr -1 $FC
ssh -l ibrik lazaro /home/ibrik/fc/./script-tar.sh
sftp ibrik@148.206.45.51:/home/ibrik/fc/latest.tar.bz2
ssh -l ibrik lazaro 'rm -rf /home/ibrik/fc/latest.tar.bz2'
