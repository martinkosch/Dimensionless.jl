[![Build Status](https://travis-ci.com/martinkosch/Dimensionless.jl.svg?branch=master)](https://travis-ci.com/martinkosch/Dimensionless.jl)
[![codecov](https://codecov.io/gh/martinkosch/Dimensionless.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/martinkosch/Dimensionless.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://martinkosch.github.io/Dimensionless.jl/dev)

# Dimensionless.jl
Dimensionless is a package built on top of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). It contains tools to conduct dimensional analysis and solve similitude problems.

## Installation
Install the package using Julia's package manager:
```julia
julia>]add Dimensionless
```

Currently, Unitful's functionality is not exported from Dimensionless. To use Dimensionless, Unitful should be imported as well:
```julia
julia> using Dimensionless, Unitful
```
