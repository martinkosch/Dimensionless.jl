module Dimensionless

using Unitful, LinearAlgebra

export basis_dimensions, dimensionless, dimensionful

struct basis_dimensions
    dimensions
    dimensional_matrix
end

function dimensionless()

end

function dimensionful()

end

end # module
