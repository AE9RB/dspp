# dspp - DSP library for C++ [![Build Status](https://travis-ci.org/AE9RB/dspp.png?branch=master)](https://travis-ci.org/AE9RB/dspp)

The dspp library is designed for high performance signal processing.
It requires a C++11 or newer compiler and makes good use of modern
techniques and the Standard Template Library.

This README is for anyone wanting to develop the dspp library itself.
If you just want to use the dspp library, begin with the [main documentation]
(http://AE9RB.github.io/dspp).

## Development Environment

You'll need a few tools to build the tests and documentation. These are
common tools available for all operating systems so you shouldn't have
too much trouble getting them installed.

 * [CMake](http://www.cmake.org)
 * [Doxygen](http://www.doxygen.org)
 * [Graphviz](http://www.graphviz.org)

Cmake can create a typical Makefile as well as project files for Xcode,
Visual Studio, and many others. Here's a quick start for Makefile users:

```
$ cmake .
$ make && ctest
```
