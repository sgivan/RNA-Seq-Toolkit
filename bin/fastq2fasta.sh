#!/bin/bash
grep -A 1 "^@" $1 | grep -v "^--" | sed "s/^@/>/" > $1.fa
