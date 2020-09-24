#!/bin/sh
command -v gfortran > /dev/null && fort=gfortran || fort=ifort
[ -f ~/.zshrc ] && rc=~/.zshrc || rc=~/.bashrc
echo "Compiling autochem/core/thermo.f"
$fort autochem/core/thermo.f -o autochem/core/thermo.exe
echo "Adding to PATH and PYTHONPATH in $rc"
echo "export PATH=\"$(pwd)/bin:\$PATH\"" >> $rc
echo "export PYTHONPATH=\"$(pwd):\$PYTHONPATH\"" >> $rc
