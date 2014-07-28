# dspp Documentation

The dspp library is designed for high performance signal processing.
It requires a C++11 or newer compiler and makes good use of modern
techniques and the Standard Template Library.

A good math library should have formulas describing how the functions work.
So you'll find formulas, examples, and descriptions for everything.
Math is displayed using [MathJax](http://www.mathjax.org).
You can right-click or ctrl-click on any formula to bring up a menu.
Give it a try:

\f[\begin{gather}
\nabla \cdot \mathbf{D} = \rho_\text{f} \\
\nabla \cdot \mathbf{B} = 0 \\
\nabla \times \mathbf{E} = -\frac{\partial \mathbf{B}} {\partial t} \\
\nabla \times \mathbf{H} = \mathbf{J}_\text{f} + \frac{\partial \mathbf{D}} {\partial t}
\end{gather}\f]

A testing framework is part of the dspp library and all commits are
verified with a continuous integration server. You can also use this
testing framework for your own applications, if you so desire. The testing
framework includes a microbenchmarking tool. This makes it difficult to
accidentally regress performance and easy to test different implementations.
You can be confident that performance is suitable for real-time applications
like Software Defined Radio.

## Installation

 1. Download or clone from GitHub. <https://github.com/AE9RB/dspp>
 2. Add the include directory to your compiler include path.
 3. #include <dspp.hpp>
 
The dspp library is a collection of headers like the Standard Template
Library. There's nothing to compile except your own application.
There are no dependencies to install.

## Contribute

View the [README](https://github.com/AE9RB/dspp) on GitHub to learn more.
Submit a pull request when you have something to contribute.
