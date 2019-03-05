#!/usr/bin/env sh

if [[ $# != 2 && $1 != '-h' && $1 != '--help' ]] ; then
  echo 'ERROR: Incorrect arguments passed'
  echo
  echo 'Must run `qcp -t 8 | tee bonds.txt` first'
  echo 'Then run this script with the two atoms to look for'
  echo 'i.e. bond_lengths.sh C N'
elif [[ $1 == '-h' || $1 == '--help' ]] ; then 
  echo
  echo 'Must run `qcp -t 8 | tee bonds.txt` first'
  echo 'Then run this script with the two atoms to look for'
  echo 'i.e. bond_lengths.sh C N'
else
  cat bonds.txt | grep -E '^[[:lower:]]|^\s*'$1'[0-9]*\s*'$2'[0-9]*|\s*'$2'[0-9]*\s*'$1'[0-9]*'
fi
