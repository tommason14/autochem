#!/bin/sh
command -v gfortran > /dev/null && fort=gfortran || fort=ifort
[ -f ~/.zshrc ] && rc=~/.zshrc || rc=~/.bashrc
echo "Compiling autochem/core/thermo.f"
$fort autochem/core/thermo.f -o autochem/core/thermo.exe
echo "Setting PYTHONPATH in $rc"
echo "export PYTHONPATH=\"$(pwd):\$PYTHONPATH\"" >> $rc
