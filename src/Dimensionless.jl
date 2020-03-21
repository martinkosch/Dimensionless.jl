module Dimensionless

using Unitful, LinearAlgebra

export dim_basis, dim_matrix
export count_dims, num_dimensionless, dimensionless, dimensionful, change_basis
export dimensionless_string

include("dim_bases.jl")
include("dim_basis_ops.jl")
include("utils.jl")

end # module
