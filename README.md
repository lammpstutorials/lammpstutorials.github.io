# LAMMPS tutorials

This is the repository of the [LAMMPS tutorials](https://lammpstutorials.github.io/)
webpage. All the LAMMPS input scripts and data files can be found
in a separate [repository](https://github.com/lammpstutorials/lammpstutorials-inputs).

## About LAMMPS tutorials

LAMMPS tutorials is made of seven tutorials that are ordered by increasing difficulty.
[Lennard-Jones fluid](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/lennard-jones-fluid.html)
is meant for absolute LAMMPS and molecular dynamics beginners, and the complexity of the simulation is
progressively increased for [Pulling on a carbon nanotube](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/breaking-a-carbon-nanotube.html), 
[Polymer in water](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level2/polymer-in-water.html), 
[Nanosheared electrolyte](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level2/nanosheared-electrolyte.html),
and [Reactive silicon dioxide](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/reactive-silicon-dioxide.html). 
Finally, [Water adsorption in silica](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/water-adsorption-in-silica.html),
[Free energy calculation](https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/free-energy-calculation.html), 
are using some more advanced simulation methods that are commonly used when studying soft matter systems, respectively grand canonical Monte Carlo simulations and a free energy method named umbrella sampling.

<p float="left">
    <a href="https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/lennard-jones-fluid.html">
        <img src="https://raw.githubusercontent.com/lammpstutorials/lammpstutorials.github.io/version2.0/docs/avatars/level1/lennard-jones-fluid/avatar-Lennard-Jones-LAMMPS.png" width="32%" />
    </a>
    <a href="https://lammpstutorials.github.io/sphinx/build/html/tutorials/level1/breaking-a-carbon-nanotube.html">
        <img src="https://raw.githubusercontent.com/lammpstutorials/lammpstutorials.github.io/version2.0/docs/avatars/level1/breaking-a-carbon-nanotube/CNT.png" width="32%" />
    </a>
    <a href="https://lammpstutorials.github.io/sphinx/build/html/tutorials/level3/water-adsorption-in-silica.html">
        <img src="https://raw.githubusercontent.com/lammpstutorials/lammpstutorials.github.io/version2.0/docs/avatars/level3/water-adsorption-in-silica/water-adsorption.png" width="32%" />
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
The other pages use the [Sphinx](https://www.sphinx-doc.org/) generator](https://www.sphinx-doc.org/) with the [furo style](https://github.com/pradyunsg/furo).

### About me & Contact ###

I am a computer physicist in soft matter and fluids at interfaces. You can 
find more information on my [personal webpage](https://simongravelle.github.io/).

See the [contact page](https://lammpstutorials.github.io/sphinx/build/html/non-tutorials/contact-me.html). 
You can report issues here on Github, or send me an [email](https://simongravelle.github.io/). Your feedback is always appreciated.

### License and Acknowledgments ###

All the LAMMPS inputs/data/parameter files and Python scripts are released under the 
GNU general public license v3.0: feel free to adapt and re-publish them.  


This project has received funding from the European
Union's Horizon 2020 research and innovation programme
under the Marie Skłodowska-Curie grant agreement No 101065060.

![Acknowledgments-logos](https://raw.githubusercontent.com/simongravelle/credits/1c44b5ae76a33c5bbbd33a54243365c6abdc24b2/cnrs-uga-liphy-msca.png)