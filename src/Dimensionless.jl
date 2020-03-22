module Dimensionless

using Unitful, LinearAlgebra

export DimBasis, dim_matrix
export number_of_dimensions, number_of_dimensionless, dimensionless, dimensionful, change_basis
export dimensionless_string, print_dimensionless

include("dim_basis.jl")
include("dim_basis_ops.jl")
include("utils.jl")

end # module
