#!/usr/bin/env python
# coding: utf-8

import numpy as np
import re, yaml

docs = ""
with open("log.lammps") as f:
    for line in f:
        pattern = (
            r"^(keywords:.*$|"
            r"data:$|"
            r"---$|"
            r"\.\.\.$|"
            r"  - \[.*\]$)"
        )
        m = re.search(pattern, line)
        if m:
            docs += m.group(0) + '\n'
thermo = list(
    yaml.load_all(docs,
                    Loader=yaml.CSafeLoader))

Etot_vs_step = []
for line in thermo[0]['data']:
    Etot_vs_step.append([line[0], line[4]])
for line in thermo[1]['data']:
    Etot_vs_step.append([line[0], line[4]])
np.savetxt("etot_vs_step.dat", Etot_vs_step)
