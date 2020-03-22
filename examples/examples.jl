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
