#!/bin/bash

python3 create_metadata.py

./wham -20 -8 50 1e-8 119.8 0 metadata.dat PMF.dat
