using Dimensionless, Unitful

applications_basis  =    dim_basis(5.5u"m", 9e5u"kg", 2000u"m/s")
model_basis         =    dim_basis(20u"cm", 48u"kg", 1500u"m/s")
change_basis(1600u"N*m", applications_basis, model_basis)

applications_basis  =    dim_basis("L"=>1u"m",      "ρ"=>1028u"kg/m^3",  "μ"=>1.88e-3u"Pa*s")
model_basis         =    dim_basis("L"=>0.025u"m",  "ρ"=>998u"kg/m^3",   "μ"=>1e-3u"Pa*s")
change_basis(5u"m/s", applications_basis, model_basis)
change_basis(u"N", model_basis, applications_basis)

applications_basis  =    dim_basis("L"=>dimension(u"m^0.5"),      "ρ"=>dimension(u"kg/m^3"),  "μ"=>dimension(u"Pa*s"), "i"=>dimension(u"A"))
# applications_basis  =    dim_basis("L"=>u"m",      "ρ"=>u"kg/m^3",  "μ"=>u"Pa*s", "i"=>u"A")
dimensionless_string("v"=>dimension(u"A*m/s"), applications_basis)
dimensionless_string("v"=>u"A*m/s", applications_basis)

problem_vars = ("P"=>dimension(u"W"), "n"=>dimension(u"s^-1"), "D"=>dimension(u"m"), "μ"=>dimension(u"Pa*s"), "ρ"=>dimension(u"kg/m^3"))
num_dimensionless(problem_vars...)
count_dims(problem_vars...)
dimensionless_string("P"=>dimension(u"W"), dim_basis(problem_vars[[2,3,5]]...))
dimensionless_string("P"=>dimension(u"A"), dim_basis("P"=>dimension(u"A")))
dimensionless_string("n"=>dimension(u"1/s"), dim_basis(problem_vars[[3,4,5]]...))

function empiric_barometric(h)
    a0 = 7.001985e-02
    a1 = -4.336216e-03
    a2 = -5.009831e-03
    a3 = 1.621827e-04
    a4 = -2.471283e-06
    a5 = 1.904383e-08
    a6 = -7.189421e-11
    a7 = 1.060067e-13
    return 10^(((((((a7*h/1000 + a6)*h/1000 + a5)*h/1000 + a4)*h/1000 + a3)*h/1000 + a2)*h/1000 + a1)*h/1000 + a0)
end
