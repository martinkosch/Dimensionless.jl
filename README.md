[![CI](https://github.com/martinkosch/Dimensionless.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/martinkosch/Dimensionless.jl/actions/workflows/CI.yml)
[![codecov](https://codecov.io/github/martinkosch/Dimensionless.jl/graph/badge.svg?token=G1SAGH9JWN)](https://codecov.io/github/martinkosch/Dimensionless.jl)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://martinkosch.github.io/Dimensionless.jl/dev)

# Dimensionless.jl
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

## Usage
The main functionality of this package is to find and make use of bases made of problem-specific variables and their units. Such bases can then be used to transform additional variables that describe a physical problem to corresponding dimensionless values and back to a dimensional form.

For example, if a problem is, among other variables, characterized by a mass ``m = 15 kg``, a length ``L = 75 cm``, and the gravitational constant ``g = 9.81 m/s^2``, a corresponding basis for the problem can be calculated:
```julia
julia> using Dimensionless, Unitful

julia> basis = DimBasis("m"=>15u"kg", "L"=>75u"cm", "g"=>9.81u"m/s^2")
DimBasis{...}
```
This basis can now be used to coherently tranform any variable within the spanned space of this basis to a dimensionless value and back to a dimensional value again:  
```julia
julia> p_dimless = dimensionless(0.2u"bar", basis)
76.45259938837918

julia> p_dimful = dimensionful(p_dimless, u"Pa", basis)
0.2 bar
```
This procedure enables reducing the number of simulations or experiments as well as streamlining problem formulations. A full [examplary use case](https://martinkosch.github.io/Dimensionless.jl/dev/example/) can be found in the documentation. 
