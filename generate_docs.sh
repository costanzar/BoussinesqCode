#!/bin/bash

# Ensure the output folder exists
mkdir -p docs/api

# Run FORD with source directory and output directory explicitly set
ford -d scr -o docs/api ford.yml

