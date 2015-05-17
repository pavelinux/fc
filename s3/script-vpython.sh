#!/bin/bash
args=("$@")
MESSAGE='from visual import *'
nl -ba -s '=sphere(' ${args[0]} > C-py.xyz
sed 's/     //g' C-py.xyz > ${args[0]}
sed 's/^/A/g' ${args[0]} > C-py.xyz
sed '/A1/ i from visual import *' C-py.xyz > ${args[0]}
sed 's/)/),radius=0.5)/g' ${args[0]} > C-py.xyz
mv C-py.xyz Coordenadas.py
rm -rf ${args[0]}
#scp Coordenadas.py pavel@atenea:.
