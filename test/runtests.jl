using Dimensionless, Unitful
using Test

basis = BasisQuantities(9.81u"m/s^2", 6371u"km", 1420788u"kg", 1u"mol/g")

@testset "Basic tests" begin
    @test basis.dimensional_matrix == Rational.([-2 0 0 0; 1 1 0 0; 0 0 0 1; 0 0 1 -1])
end
