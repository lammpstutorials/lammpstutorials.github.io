#!/usr/bin/env python
# coding: utf-8

import sys, os, git
import numpy as np

from utilities import ReadRST, WriteTex, FixDocument

current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")
sys.path.append(git_path+"/docs/inputs/shared-pyplot-files/")
if os.path.exists(git_path+'/ebook/tutorials') is False:
    os.mkdir(git_path+'/ebook/tutorials')

tutorials = {"level0": ["lennard-jones-fluid"]}
             #"level1": ["breaking-a-carbon-nanotube"],
             #"level2": ["polymer-in-water", "nanosheared-electrolyte"],
             #"level3": ["water-adsorption-in-silica", "free-energy-calculation", "reactive-silicon-dioxide"]}

for level in tutorials.keys():
    if os.path.exists(git_path+'/ebook/tutorials/'+level) is False:
        os.mkdir(git_path+'/ebook/tutorials/'+level)
    for tutorial in tutorials[level]:
        rst_file_name = git_path+'/docs/sphinx/source/tutorials/'+level+'/'+tutorial+'.rst'
        tex_file_name = git_path+'/ebook/tutorials/'+level+'/'+tutorial+'.tex'
        RST = ReadRST(rst_file_name)
        RST.convert_file()
        TEX = WriteTex(tex_file_name, RST, git_path)
        TEX.convert_file()
        FIX = FixDocument(tex_file_name, RST, git_path)
        FIX.fix_document()

