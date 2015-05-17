#!/bin/bash
ifort -o coordenadas md1.f90
chmod 755 coordenadas
./coordenadas
