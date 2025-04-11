# LAMMPS tutorials

This is the repository of the [LAMMPS tutorials](https://lammpstutorials.github.io/)
webpage. All the LAMMPS input scripts and data files can be found in a separate
repository named [lammpstutorials-inputs](https://github.com/lammpstutorials/lammpstutorials-inputs).

The tutorials are compatible with the XXXX2024 release of LAMMPS.

## About LAMMPS tutorials

The LAMMPStutorials website is made of seven tutorials that are ordered by increasing difficulty.
[Lennard-Jones fluid](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/lennard-jones-fluid.html)
is meant for absolute LAMMPS and molecular dynamics beginners, and the complexity of the simulation is
progressively increased for [Pulling on a carbon nanotube](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/breaking-a-carbon-nanotube.html), 
[Polymer in water](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level2/polymer-in-water.html), 
[Nanosheared electrolyte](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level2/nanosheared-electrolyte.html),
and [Reactive silicon dioxide](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/reactive-silicon-dioxide.html). 
Finally, [Water adsorption in silica](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/water-adsorption-in-silica.html) and
[Free energy calculation](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/free-energy-calculation.html) use some more advanced simulation methods that are commonly used when studying soft matter systems, respectively grand canonical Monte Carlo simulations and a free energy method named umbrella sampling.

<p float="left">
    <a href="https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/lennard-jones-fluid.html">
        <img src="https://raw.githubusercontent.com/lammpstutorials/lammpstutorials.github.io/2Aug2023/docs/avatars/level1/lennard-jones-fluid/avatar-Lennard-Jones-LAMMPS.png" width="32%" />
    </a>
    <a href="https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/breaking-a-carbon-nanotube.html">
        <img src="https://raw.githubusercontent.com/lammpstutorials/lammpstutorials.github.io/2Aug2023/docs/avatars/level1/breaking-a-carbon-nanotube/CNT.png" width="32%" />
    </a>
    <a href="https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/water-adsorption-in-silica.html">
        <img src="https://raw.githubusercontent.com/lammpstutorials/lammpstutorials.github.io/2Aug2023/docs/avatars/level3/water-adsorption-in-silica/water-adsorption.png" width="32%" />
    </a>
</p>

## Access the files

You can access all the files by cloning this repository with its submodules:

```
git clone https://github.com/lammpstutorials/lammpstutorials.github.io.git --recurse-submodule
```

Alternatively, you can download the [inputs](https://github.com/lammpstutorials/lammpstutorials-inputs) only:

```
git clone https://github.com/lammpstutorials/lammpstutorials.github.io.git
```

The Matplotlib Pyplot functions for the figures are shared [here](https://github.com/simongravelle/pyplot-perso).

### Template ###

The template from the first page has been adapted from [HTML5 UP](https://html5up.net/).
The other pages use the [Sphinx](https://www.sphinx-doc.org/) generator with the
[furo style](https://github.com/pradyunsg/furo).

