#!/bin/bash

# convert parallel code (PAR) to sequential (SEQ)

# //PAR comment, //SEQ uncomment

# implementation
find . -name '*.cpp' -type f -print | xargs perl -i -p -e  's/^(\s*[^\/\/])(.*)\/\/PAR$/\/\/\1\2\/\/PAR/g'
find . -name '*.cpp' -type f -print | xargs perl -i -p -e  's/^(\s*)\/\/(.*)\/\/SEQ$/\1\2\/\/SEQ/g'

# headers
find . -name '*.h' -type f -print | xargs perl -i -p -e  's/^(\s*[^\/\/])(.*)\/\/PAR$/\/\/\1\2\/\/PAR/g'
find . -name '*.h' -type f -print | xargs perl -i -p -e  's/^(\s*)\/\/(.*)\/\/SEQ$/\1\2\/\/SEQ/g'