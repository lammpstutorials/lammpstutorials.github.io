#!/bin/bash

python3 create_metadata.py

./wham -25 25 50 1e-8 119.8 0 metadata.dat PMF.dat
