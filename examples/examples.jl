using Dimensionless, Unitful

applications_basis  =    dimension_basis(5.5u"m", 9e5u"kg", 2000u"m/s")
model_basis         =    dimension_basis(20u"cm", 48u"kg", 1500u"m/s")
change_basis(1600u"N*m", applications_basis, model_basis)

applications_basis  =    dimension_basis(1u"m",      1028u"kg/m^3",  1.88e-3u"Pa*s")
model_basis         =    dimension_basis(0.025u"m",   998u"kg/m^3",   1e-3u"Pa*s")
change_basis(5u"m/s", applications_basis, model_basis)
old_to_new_multiplier(u"N", model_basis, applications_basis)
