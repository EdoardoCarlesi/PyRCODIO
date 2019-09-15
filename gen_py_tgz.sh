#!/bin/bash

file_tgz='py_lgtools.tgz'
sources='*.py libcosmo libics *.md *.sh'

tar -cjf $file_tgz $sources
