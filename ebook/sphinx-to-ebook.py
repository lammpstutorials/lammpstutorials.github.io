#!/usr/bin/env python
# coding: utf-8

import sys, os, git

from ReadRST import ReadRST
from WriteTEX import WriteTex
from FixDocument import FixDocument

current_path = os.getcwd()
git_repo = git.Repo(current_path, search_parent_directories=True)
git_path = git_repo.git.rev_parse("--show-toplevel")
sys.path.append(git_path+"/docs/inputs/shared-pyplot-files/")
if os.path.exists(git_path+'/ebook/tutorials') is False:
    os.mkdir(git_path+'/ebook/tutorials')

tutorials = {"level1": ["lennard-jones-fluid",
                        "breaking-a-carbon-nanotube"],
             "level2": ["polymer-in-water",
                        "nanosheared-electrolyte"],
             "level3": ["water-adsorption-in-silica",
                        "free-energy-calculation",
                        "reactive-silicon-dioxide"],
             "vmd": ["vmd-tutorial"]}

for level in tutorials.keys():
    if os.path.exists(git_path+'/ebook/tutorials/'+level) is False:
        os.mkdir(git_path+'/ebook/tutorials/'+level)
    for tutorial in tutorials[level]:
        print(level, "tutorial", tutorial)
        print("-----------------------------------------")
        rst_file_name = git_path+'/docs/sphinx/source/tutorials/'+level+'/'+tutorial+'.rst'
        tex_file_name = git_path+'/ebook/tutorials/'+level+'/'+tutorial+'.tex'
        RST = ReadRST(rst_file_name)
        RST.convert_file()
        assert len(RST.label_positions) == 1, """Careful, more than one label"""
        TEX = WriteTex(tex_file_name, RST, git_path)
        TEX.convert_file()
        FIX = FixDocument(tex_file_name)
        FIX.fix_document()

#non_tutorials = {"solutions": ["solutions"],
#                 "before-you-start": ["before-you-start"]}

if os.path.exists(git_path+'/ebook/non-tutorials/') is False:
    os.mkdir(git_path+'/ebook/non-tutorials/')

print("before-you-start")
print("-----------------------------------------")

rst_file_name = git_path+'/docs/sphinx/source/non-tutorials/before-you-start.rst'
tex_file_name = git_path+'/ebook/non-tutorials/before-you-start.tex'
RST = ReadRST(rst_file_name)
RST.convert_file()
assert len(RST.label_positions) == 1, """Careful, more than one label"""
TEX = WriteTex(tex_file_name, RST, git_path, nonumber=True)
TEX.convert_file()
FIX = FixDocument(tex_file_name)
FIX.fix_document()

print("solutions")
print("-----------------------------------------")

rst_file_name = git_path+'/docs/sphinx/source/non-tutorials/solutions.rst'
tex_file_name = git_path+'/ebook/non-tutorials/solutions.tex'
RST = ReadRST(rst_file_name)
RST.convert_file()
assert len(RST.label_positions) == 1, """Careful, more than one label"""
TEX = WriteTex(tex_file_name, RST, git_path, nonumber=True)
TEX.convert_file()
FIX = FixDocument(tex_file_name)
FIX.fix_document()
