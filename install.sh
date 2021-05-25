#!/usr/bin/env bash
python3 -m pip install numpy pandas dfply --user
command -v gfortran > /dev/null && fort=gfortran || fort=ifort
[[ $SHELL =~ zsh ]] && rc=~/.zshrc || rc=~/.bashrc
echo "Compiling autochem/core/thermo.f"
$fort autochem/core/thermo.f -o autochem/core/thermo.exe
echo "Adding to PATH and PYTHONPATH in $rc"
echo "export PATH=\"$(pwd)/bin:\$PATH\"" >> $rc
echo "export PYTHONPATH=\"$(pwd):\$PYTHONPATH\"" >> $rc
