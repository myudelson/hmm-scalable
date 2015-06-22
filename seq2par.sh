#!/bin/bash

# PAR to SEQ
# convert sequential (SEQ) code to parallel (PAR)

# find . -name '*.cpp' -type f -print | xargs grep -e '^\(.*\)//PAR$' | wc -l


# implementation
find . -name '*.cpp' -type f -print | xargs perl -i -p -e  's/^(\s*[^\/\/])(.*)\/\/SEQ$/\/\/\1\2\/\/SEQ/g'
find . -name '*.cpp' -type f -print | xargs perl -i -p -e  's/^(\s*)\/\/(.*)\/\/PAR$/\1\2\/\/PAR/g'

# headers
find . -name '*.h' -type f -print | xargs perl -i -p -e  's/^(\s*[^\/\/])(.*)\/\/SEQ$/\/\/\1\2\/\/SEQ/g'
find . -name '*.h' -type f -print | xargs perl -i -p -e  's/^(\s*)\/\/(.*)\/\/PAR$/\1\2\/\/PAR/g'