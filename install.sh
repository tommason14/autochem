#!/bin/sh
echo "Compiling autochem/core/thermo.f"
gfortran autochem/core/thermo.f -o autochem/core/thermo.exe
[ -f ~/.zshrc ] && rc=~/.zshrc || rc=~/.bashrc
echo "Setting PYTHONPATH in $rc"
echo "export PYTHONPATH=\"$(pwd):\$PYTHONPATH\"" >> $rc
