# Introduction
This package provides tools for dimensional analysis and similitude problems. Dimensionless.jl uses units and dimensions from [Unitful.jl](https://github.com/PainterQubits/Unitful.jl). 

## Installation
Clone the package repository from GitHub using Julia's package manager:
```julia
pkg> add https://github.com/martinkosch/Dimensionless.jl
```

Currently, Unitful's functionality is not exported from Dimensionless. To use Dimensionless, Unitful needs to be imported as well:
```julia
julia> using Dimensionless, Unitful
```
