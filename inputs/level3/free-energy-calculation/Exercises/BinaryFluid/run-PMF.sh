#!/bin/bash

python3 create_metadata.py

./wham -25 25 50 1e-8 119.8 0 metadata-particle1.dat PMF-particle1.dat

./wham -25 25 50 1e-8 119.8 0 metadata-particle2.dat PMF-particle2.dat