#!/bin/bash
args=("$@")
ifort -o ${args[1]} ${args[0]}
chmod 755 ${args[1]}
./${args[1]}
