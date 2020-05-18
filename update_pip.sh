#!/bin/sh
rm -r autochem.egg-info build dist
printf "Reminder to increment version number in setup.py. Press any key to continue"
read
vim setup.py
python3 setup.py bdist_wheel sdist
printf "Upload to pip? [y/n] "
read 
[ $REPLY = 'y' ] || [ $REPLY = 'Y' ] && twine upload dist/*
