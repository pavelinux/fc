#!/bin/bash
args=("$@")
/opt/intel/bin/ifort -o ${args[1]} ${args[0]}
chmod 755 ${args[1]}
./${args[1]}
#JMOL='/home/pavel/Descargas/jmol-14.2.13_2015.03.23/./jmol.sh'
#$JMOL POSICIONES.XYZ
