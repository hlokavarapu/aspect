#!/bin/bash

fileName=$1

sed -i "s/subsection Density/subsection Analytical density/" $1
sed -i "s/subsection Pressure/subsection Analytical pressure/" $1
