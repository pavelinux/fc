#!/bin/bash
args=("$@")
/opt/intel/bin/ifort -o ${args[1]} ${args[0]}
chmod 755 ${args[1]}
./${args[1]}
