using Dimensionless, Unitful
using Test

@testset "DimensionBasis tests" begin
    @test_throws ErrorException DimensionBasis(9.81u"m/s^2", 1u"mm", 1u"s")
    @test_throws ErrorException DimensionBasis(9.81u"m/s^2")

    basis = DimensionBasis(9.81u"m/s^2", 6371u"km", 1420788u"kg", 1u"mol/g")
    @test_throws ErrorException dimensional_matrix(basis.basis_dims, 1u"A")

    @test basis.dim_mat == Rational.([-2 0 0 0; 1 1 0 0; 0 0 0 1; 0 0 1 -1])
    for dim in [:Mass, :Length, :Time, :Amount]
        @test Unitful.Dimension{dim} in basis.dim_set
    end
    for dim in [:Temperature, :Current, :Luminosity]
        @test !(Unitful.Dimension{dim} in basis.dim_set)
    end

end
