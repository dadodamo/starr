# starr

R package for fitting spatio-temporal AR(1) models.

## Installation and setup

The `starr` package is compiling and executing C++ code that has several dependencies: 

- [Cmake](https://cmake.org): A powerful, comprehensive solution for managing the software build process. 
- [Google Protocol Buffers](https://github.com/protocolbuffers/protobuf): A language-neutral, platform-neutral extensible mechanism for serializing structured data.
- [Boost](https://www.boost.org): Provides free peer-reviewed portable C++ source libraries.
- [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page): A C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.

If not already installed, the package will download and install Boost and Eigen automatically during configuration using 
[Cmake](https://cmake.org). `Cmake` and `protobuf` will have to be installed by the user manually. See  [Cmake](https://cmake.org)
and [protobuf](https://github.com/protocolbuffers/protobuf) for more info. 

After `cmake` and `protobuf` are installed, open the terminal and `cd` to your desired package location and run:
`````
git clone https://github.com/dadodamo/starr.git
`````

Now the package can be installed via `R`. Simply open R and set your working directory to 
`````
setwd('[your chosen path]/starr/R/starpackage')
`````
and run
`````
devtools::install(quick = TRUE)
library(starr)
`````
Before using `starr::run_code()`, we have to compile the C++ code once in each R session. We can do that through R by simply calling the function 
`starr::compile_code()`, which will build the C++ executable and load the dependencies. 
