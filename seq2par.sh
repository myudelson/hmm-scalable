#!/bin/bash

# PAR to SEQ

# find . -name '*.cpp' -type f -print | xargs grep -e '^\(.*\)//PAR$' | wc -l

# //PAR comment, //SEQ uncomment
find . -name '*.cpp' -type f -print | xargs perl -i -p -e  's/^(.*)\/\/SEQ$/\/\/\1\/\/SEQ/g'
find . -name '*.cpp' -type f -print | xargs perl -i -p -e  's/^\/\/(.*)\/\/PAR$/\1\/\/PAR/g'

find . -name '*.h' -type f -print | xargs perl -i -p -e  's/^(.*)\/\/SEQ$/\/\/\1\/\/SEQ/g'
find . -name '*.h' -type f -print | xargs perl -i -p -e  's/^\/\/(.*)\/\/PAR$/\1\/\/PAR/g'