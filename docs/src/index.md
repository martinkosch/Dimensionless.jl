# Introduction
Dimensionless is a package built on top of [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). It contains tools to conduct dimensional analysis and solve similitude problems based on [the Buckingham Π theorem](https://en.wikipedia.org/wiki/Buckingham_%CF%80_theorem).

## Installation
Install the package using Julia's package manager:
```julia
pkg> add Dimensionless
```

Currently, Unitful's functionality is not exported from Dimensionless. To use Dimensionless, Unitful should be imported as well:
```julia
julia> using Dimensionless, Unitful
```
