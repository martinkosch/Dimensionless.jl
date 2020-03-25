# Introduction
This package provides tools for dimensional analysis and similitude problems. Dimensionless.jl uses units and dimensions from [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).

## Installation
Install the package using Julia's package manager:
```julia
pkg> add Dimensionless
```

Currently, Unitful's functionality is not exported from Dimensionless. To use Dimensionless, Unitful should be imported as well:
```julia
julia> using Dimensionless, Unitful
```
