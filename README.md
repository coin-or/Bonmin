# Bonmin 1.8

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

Projects such as this one are maintained by a small group of volunteers under
the auspices of the non-profit [COIN-OR Foundation](https://www.coin-or.org)
and we need your help! Please consider [sponsoring our
activities](https://github.com/sponsors/coin-or) or [volunteering](mailto:volunteer@coin-or.org) to help!

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Bonmin?sort=semver)](https://github.com/coin-or/Bonmin/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

Bonmin (*B*asic *O*pen-source *N*onlinear *M*ixed *IN*teger programming) is
an open-source code for solving general MINLP (Mixed Integer NonLinear
Programming) problems. It builds on top of [Cbc](https://github.com/coin-or/Cbc)
and [Ipopt](https://github.com/coin-or/Ipopt).
It can be used as a callable library, through its C++ or C interfaces, or as
a stand-alone executable that reads models written in [AMPL](https://ampl.com).

Bonmin features several algorithms, selected with the `algorithm` option:
 * B-BB, an NLP-based branch-and-bound algorithm (the default),
 * B-OA, an outer-approximation decomposition algorithm,
 * B-QG, an implementation of Quesada and Grossmann's branch-and-cut algorithm,
 * B-Hyb, a hybrid outer-approximation based branch-and-cut algorithm,
 * B-Ecp, an ECP cuts based branch-and-cut algorithm a la FilMINT,
 * B-iFP, an iterated feasibility pump for MINLP.

These algorithms are exact when the problem is convex. For nonconvex
problems, they are heuristics.

Bonmin can also be used through
 * AMPL (natively, using the `bonmin` executable),
 * GAMS (through the [GAMSlinks](https://github.com/coin-or/GAMSlinks) project),
 * [JuMP](https://github.com/jump-dev/AmplNLWriter.jl) (Julia, via AmplNLWriter),
 * [Pyomo](https://github.com/Pyomo/pyomo) (Python),
 * [NEOS](https://neos-server.org/neos/solvers/minco:Bonmin/AMPL.html).

This project was initiated in 2004 by IBM and Carnegie Mellon University as
part of a joint effort to study algorithms for MINLP. You may find additional
information at http://egon.cheme.cmu.edu/ibm/page.htm, in particular
 * A publicly available library of test instances of convex MINLPs,
 * Research papers on new algorithmic procedures for MINLP.

Bonmin was developed by Pierre Bonami and the other contributors listed in
[AUTHORS](Bonmin/AUTHORS). 

Bonmin is written in C++ and is released as open source under the [Eclipse Public License 2.0](http://www.opensource.org/licenses/EPL-2.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org).

The Bonmin development site is https://github.com/coin-or/Bonmin.

## CITE

Paper: https://doi.org/10.1016/j.disopt.2006.10.011

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/Bonmin/actions/workflows/windows-ci.yml/badge.svg?branch=stable/1.8)](https://github.com/coin-or/Bonmin/actions/workflows/windows-ci.yml?query=branch%3Astable/1.8)

[![Linux and MacOS Builds](https://github.com/coin-or/Bonmin/actions/workflows/linux-ci.yml/badge.svg?branch=stable/1.8)](https://github.com/coin-or/Bonmin/actions/workflows/linux-ci.yml?query=branch%3Astable/1.8)

## DOWNLOAD

What follows is a quick start guide for obtaining or building
Bonmin on common platforms. More detailed information is
available [here](https://coin-or.github.io/user_introduction.html).

### Docker image

There is a Docker image that provides Bonmin, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/repository/docker/coinor/coin-or-optimization-suite)

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/Bonmin/releases). Older binaries
are archived as part of Bonmin
[here](https://www.coin-or.org/download/binary/Bonmin).

 * *Linux* (see https://repology.org/project/coin-or-bonmin/versions for a complete listing): 
   * arch:
     ```
     $ sudo pacman -S  coin-or-bonmin
     ```
   * Debian/Ubuntu:
     ```
     $ sudo apt-get install  coinor-libbonmin-dev
     ```
   * Fedora/Redhat/CentOS:
     ```
     $ sudo dnf install  coin-or-Bonmin coin-or-Bonmin-devel
     ```
   * freebsd:
     ```
     $ sudo pkg install math/bonmin
     ```
 * *Windows*: The easiest way to get Bonmin on Windows is to download an archive as described above.

* *conda* (cross-platform):
     ```
     $ conda install coin-or-bonmin
     ```

Due to license incompatibilities, pre-compiled binaries lack some 
functionality. If binaries are not available for your platform for the latest 
version and you would like to request them to be built and posted, feel free 
to open a discussion on Github. 

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Bonmin from the
 [releases](https://github.com/coin-or/Bonmin/releases) page,
 * Cloning this repository from [Github](https://github.com/coin-or/Bonmin), or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

### Dependencies

Bonmin has a number of dependencies, which are detailed in
[config.yml](.coin-or/config.yml). Dependencies on other COIN-OR projects are
automatically downloaded when obtaining the source with `coinbrew`. For some
of the remaining third-party dependencies, automatic download scripts and
build wrappers are provided (and will also be automatically run for required
and recommended dependencies), while other libraries that are aeasy to obtain
must be installed using an appropriate package manager (or may come with your
OS by default). 

## BUILDING from source

These quick start instructions assume you are in a bash shell. 

### Using `coinbrew`

To download and build Bonmin from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Bonmin@1.8
./coinbrew build Bonmin
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/Bonmin
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

### Linear solvers for Ipopt

Bonmin solves its NLP subproblems with Ipopt, which needs a sparse linear
solver. When building with `coinbrew`, [MUMPS](https://github.com/coin-or-tools/ThirdParty-Mumps)
is built by default. The HSL codes (e.g., MA27) are often faster. To use them,
obtain the source from [HSL](https://www.hsl.rl.ac.uk/ipopt/) and follow the
instructions in [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL)
before building.

## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxydoc` 

in the build directory. If Bonmin was built via `coinbrew`, then the build
directory will be `./build/Bonmin/1.8` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/Bonmin/Doxygen).

## DOCUMENTATION

 * [User's Manual](https://coin-or.github.io/Bonmin) (LaTeX source in the [doc](Bonmin/doc) directory)
 * [List of options](https://coin-or.github.io/Bonmin/options_list.html)
 * Source code [examples](Bonmin/examples)


## Project Links

 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [COIN-OR general discussion forum](https://github.com/orgs/coin-or/discussions)
 * [Bonmin Discussion forum](https://github.com/coin-or/Bonmin/discussions)
 * [Report a bug](https://github.com/coin-or/Bonmin/issues/new)
 * [Doxygen generated documentation](http://coin-or.github.io/Bonmin/Doxygen)

