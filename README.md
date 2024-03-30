[![CI](https://github.com/martinkosch/Dimensionless.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/martinkosch/Dimensionless.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/github/martinkosch/Dimensionless.jl/graph/badge.svg?token=G1SAGH9JWN)](https://codecov.io/github/martinkosch/Dimensionless.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://martinkosch.github.io/Dimensionless.jl/dev)

# Dimensionless.jl
Dimensionless is a package built on top of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). It contains tools to conduct dimensional analysis and solve similitude problems. An [examplary use case](https://martinkosch.github.io/Dimensionless.jl/dev/example/) can be found in the documentation. 

## Installation
Install the package using Julia's package manager:
```julia
pkg> add Dimensionless
```

Currently, Unitful's functionality is not exported from Dimensionless. To use Dimensionless, Unitful should be imported as well:
```julia
julia> using Dimensionless, Unitful
```
