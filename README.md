
# Master Thesis: Adaptive timestepping for the Simulation of Electric Machines

This repository has been forked from the LehrFEM++ repository, all the files relevant to the Master thesis can be found in the folder [cylinder_test](projects/cylinder_test)

This repository contains:

1. The LehrFEM++ library
2. A program used to solve the Eddy current equation with iterative methods [BDF-1/2](projects/cylinder_test/solve_non-linear.cc)
3. Programs used to solve the Eddy current equation with ROW methods on a static mesh [ROW-static](projects/cylinder_test/solve_ROW_no_rotation_main.cc) and a rotating mesh [ROW-static](projects/cylinder_test/solve_ROW_complete.cc)
4. 

1. [The specification](spec.md) for how a standard README should look.
2. A link to [a linter](https://github.com/RichardLitt/standard-readme-preset) you can use to keep your README maintained ([work in progress](https://github.com/RichardLitt/standard-readme/issues/5)).
3. A link to [a generator](https://github.com/RichardLitt/generator-standard-readme) you can use to create standard READMEs.
4. [A badge](#badge) to point to this spec.
5. [Examples of standard READMEs](example-readmes/) - such as this file you are reading.

Standard Readme is designed for open source libraries. Although it’s [historically](#background) made for Node and npm projects, it also applies to libraries in other languages and package managers.


## Table of Contents

- [Background](#background)
- [Install](#install)
- [Usage](#usage)


## Background

Standard Readme started with the issue originally posed by [@maxogden](https://github.com/maxogden) over at [feross/standard](https://github.com/feross/standard) in [this issue](https://github.com/feross/standard/issues/141), about whether or not a tool to standardize readmes would be useful. A lot of that discussion ended up in [zcei's standard-readme](https://github.com/zcei/standard-readme/issues/1) repository. While working on maintaining the [IPFS](https://github.com/ipfs) repositories, I needed a way to standardize Readmes across that organization. This specification started as a result of that.

> Your documentation is complete when someone can use your module without ever
having to look at its code. This is very important. This makes it possible for
you to separate your module's documented interface from its internal
implementation (guts). This is good because it means that you are free to
change the module's internals as long as the interface remains the same.

> Remember: the documentation, not the code, defines what a module does.

~ [Ken Williams, Perl Hackers](http://mathforum.org/ken/perl_modules.html#document)

Writing READMEs is way too hard, and keeping them maintained is difficult. By offloading this process - making writing easier, making editing easier, making it clear whether or not an edit is up to spec or not - you can spend less time worrying about whether or not your initial documentation is good, and spend more time writing and using code.

By having a standard, users can spend less time searching for the information they want. They can also build tools to gather search terms from descriptions, to automatically run example code, to check licensing, and so on.

The goals for this repository are:

1. A well defined **specification**. This can be found in the [Spec document](spec.md). It is a constant work in progress; please open issues to discuss changes.
2. **An example README**. This Readme is fully standard-readme compliant, and there are more examples in the `example-readmes` folder.
3. A **linter** that can be used to look at errors in a given Readme. Please refer to the [tracking issue](https://github.com/RichardLitt/standard-readme/issues/5).
4. A **generator** that can be used to quickly scaffold out new READMEs. See [generator-standard-readme](https://github.com/RichardLitt/generator-standard-readme).
5. A **compliant badge** for users. See [the badge](#badge).

## Install

This project uses [node](http://nodejs.org) and [npm](https://npmjs.com). Go check them out if you don't have them locally installed.

```sh
$ npm install --global standard-readme-spec
```

## Usage

This is only a documentation package. You can print out [spec.md](spec.md) to your console:

```sh
$ standard-readme-spec
# Prints out the standard-readme spec
```


[![Build Status](https://github.com/craffael/lehrfempp/workflows/Continuous%20Integration/badge.svg?branch=master)](https://github.com/craffael/lehrfempp/actions)

# LehrFEM++
Simple C++ Finite Element Framework for research and eduction optimzed for clarity and
flexibility with some trade-off concerning performance. This libary is used for the course _Numerical Methods for Partial Differential Euqations_ taught by Prof. R. Hiptmair at ETH Zurich.

* LehrFEM++ follows the [Google C++ Style
Guide](https://google.github.io/styleguide/cppguide.html#Naming).
* Adhere to the LehrFEM [coding style
  guidelines](https://github.com/craffael/lehrfempp/wiki/Contribute).
* Whenever adding core functionality, thorough testing is mandatory, following the
  instructions in the [testing
  guidelines](https://github.com/craffael/lehrfempp/wiki/Contribute).
* [Doxygen Class Documentation](https://craffael.github.io/lehrfempp)
* [Getting Started Guide](https://craffael.github.io/lehrfempp/getting_started.html)

## Contributors
- Raffael Casagrande (core developer)
- Ralf Hiptmair (core developer)
- Tobias Rohner (`projects/ipdg_stokes`, hp-fem in `lf::fe`)
- Anian Ruoss (Second order Geometry, Mesh Generators)
- Philippe Peter (`projects/dpg`)
- Amélie Justine Loher (`projects/FisherKPP`)
- Gina Magnussen (TIKZ output)
- Julien Gacon (`lf::base::comm`)
- Simon Meierhans

