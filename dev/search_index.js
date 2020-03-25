var documenterSearchIndex = {"docs":
[{"location":"library/#Library-1","page":"Library","title":"Library","text":"","category":"section"},{"location":"library/#Types-1","page":"Library","title":"Types","text":"","category":"section"},{"location":"library/#","page":"Library","title":"Library","text":"DimBasis","category":"page"},{"location":"library/#Dimensionless.DimBasis","page":"Library","title":"Dimensionless.DimBasis","text":"DimBasis(basis_vectors...) -> DimBasis\n\nCreate a dimensional basis for a number of basis_vectors (quantities, units or dimensions). A string identifier can optionally be added to each basis vector. \n\n\n\n\n\n","category":"type"},{"location":"library/#Functions-1","page":"Library","title":"Functions","text":"","category":"section"},{"location":"library/#","page":"Library","title":"Library","text":"dim_matrix","category":"page"},{"location":"library/#Dimensionless.dim_matrix","page":"Library","title":"Dimensionless.dim_matrix","text":"dim_matrix(basis_dims, all_values...)\n\nReturn the dimensional matrix for a set of basis dimensions basis_dims and all_values, a set of quantities, units or dimensions.\n\n\n\n\n\n","category":"function"},{"location":"library/#","page":"Library","title":"Library","text":"number_of_dimensions","category":"page"},{"location":"library/#Dimensionless.number_of_dimensions","page":"Library","title":"Dimensionless.number_of_dimensions","text":"number_of_dimensions(all_vars...)\n\nReturn the number of unique dimensions for the given variables all_vars.\n\n\n\n\n\n","category":"function"},{"location":"library/#","page":"Library","title":"Library","text":"number_of_dimensionless","category":"page"},{"location":"library/#Dimensionless.number_of_dimensionless","page":"Library","title":"Dimensionless.number_of_dimensionless","text":"number_of_dimensionless(all_vars...)\n\nReturn the number of dimensionless numbers that can be constructed for a problem characterized by the variables all_vars.\n\n\n\n\n\n","category":"function"},{"location":"library/#","page":"Library","title":"Library","text":"dimensionless","category":"page"},{"location":"library/#Dimensionless.dimensionless","page":"Library","title":"Dimensionless.dimensionless","text":"dimensionless(quantity, basis)\n\nMake a quantity dimensionless using a dimensional basis.\n\n\n\n\n\n","category":"function"},{"location":"library/#","page":"Library","title":"Library","text":"dimensionful","category":"page"},{"location":"library/#Dimensionless.dimensionful","page":"Library","title":"Dimensionless.dimensionful","text":"dimensionful(value, unit, basis)\n\nRestore the units of a dimensionless value using a dimensional basis.\n\n\n\n\n\n","category":"function"},{"location":"library/#","page":"Library","title":"Library","text":"change_basis","category":"page"},{"location":"library/#Dimensionless.change_basis","page":"Library","title":"Dimensionless.change_basis","text":"change_basis(var, basis, new_basis)\n\nTransform the specified quantity or unit var from a current basis to a new_basis.\n\n\n\n\n\n","category":"function"},{"location":"library/#","page":"Library","title":"Library","text":"print_dimensionless","category":"page"},{"location":"library/#Dimensionless.print_dimensionless","page":"Library","title":"Dimensionless.print_dimensionless","text":"print_dimensionless([io, ]named_dim_vector, basis)\n\nPrint the dimensionless number that can be constructed using named_dim_vector in the specified dimensional basis.\n\n\n\n\n\n","category":"function"},{"location":"#Introduction-1","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This package provides tools for dimensional analysis and similitude problems. Dimensionless.jl uses units and dimensions from Unitful.jl.","category":"page"},{"location":"#Installation-1","page":"Home","title":"Installation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Install the package using Julia's package manager:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"pkg> add Dimensionless","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Currently, Unitful's functionality is not exported from Dimensionless. To use Dimensionless, Unitful should be imported as well:","category":"page"},{"location":"#","page":"Home","title":"Home","text":"julia> using Dimensionless, Unitful","category":"page"},{"location":"example/#Example-1","page":"Example","title":"Example","text":"","category":"section"},{"location":"example/#","page":"Example","title":"Example","text":"The following example is taken from the Similitude Wikipedia article:","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"Consider a submarine modeled at 1/40th scale. The application operates in sea water at 0.5 °C, moving at 5 m/s. The model will be tested in fresh water at 20 °C. Find the power required for the submarine to operate at the stated speed.A free body diagram is constructed and the relevant relationships of force and velocity are formulated using techniques from continuum mechanics. The variables which describe the system are:Variable Application Scaled model Units\nL (Diameter of submarine) 1 1/40 (m)\nV (Speed) 5 Calculate (m/s)\nρ (Density) 1028 998 (kg/m^3)\nμ (Dynamic viscosity) 1.88x10^−3 1.00x10^−3 Pa·s (N s/m^2)\nF (Force) Calculate To be measured N (kg m/s^2)","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"Let's analyze this exemplary problem using Dimensionless.jl. We start by calculating the number of dimensions that characterize the problem:","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"julia> number_of_dimensions(u\"m\", u\"m/s\", u\"kg/m^3\", u\"Pa*s\", u\"N\")\n3","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"According to Buckinghams π theorem, the problem is governed by two dimensionless numbers (5 variables - 3 dimensions = 2 dimensional numbers):","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"julia> number_of_dimensionless(u\"m\", u\"m/s\", u\"kg/m^3\", u\"Pa*s\", u\"N\")\n2","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"In order to find the required power of the submarine, we need to find a dimensional basis for the problem. In Dimensionless.jl, this is done by creating an instance of type DimBasis. As calculated before, there need to be three linear independent basis vectors to span the dimensional space:","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"julia> basis = DimBasis(u\"m\", u\"m/s\", u\"kg/m^3\", u\"Pa*s\") # Invalid basis\nERROR: Invalid basis! There are 4 basis vectors of which only 3 are linear independent.\n\njulia> basis = DimBasis(u\"m/s\", u\"kg/m^3\", u\"Pa*s\") # Valid basis\nDimBasis{...}","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"A valid basis! Together with the two remaining variables v and F, the two dimensionless numbers can be determined:","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"julia> print_dimensionless(\"v\"=>u\"m/s\", basis)\nv L ρ μ^-1\n\njulia> print_dimensionless(\"F\"=>u\"N\", basis)\nF ρ μ^-2","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"The first found number is the Reynolds number mathitRe and the second found number is the product of the drag coefficient c_d and mathitRe^2. All equations describing the problem can be scaled using these numbers.","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"It is possible to add names and quantities to a DimBasis. Let's construct two bases that describe the submarine model and the real submarine:","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"julia> applications_basis = DimBasis(\"L\"=>1u\"m\", \"ρ\"=>1028u\"kg/m^3\", \"μ\"=>1.88e-3u\"Pa*s\")\nDimBasis{...}\n\njulia> model_basis = DimBasis(\"L\"=>0.025u\"m\", \"ρ\"=>998u\"kg/m^3\", \"μ\"=>1e-3u\"Pa*s\")\nDimBasis{...}","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"It is easy to transform values from one dimensional bases to the other in order to fill the missing values in the table:","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"julia> v_model = change_basis(5u\"m/s\", applications_basis, model_basis)\n109.58 m s^-1\n\njulia> change_basis(u\"N\", model_basis, applications_basis)\n3.4313\n\njulia> P_application_fac = v_model * change_basis(u\"W\", model_basis, applications_basis)\n17.156 m s^-1","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"Similarity between the model and the application is achieved with a model velocity of 109.6 m/s, which is 21.9 times higher than the application velocity. The measured drag forces acting on the model under these circumstances need to be multiplied by 3.4 to determine the drag forces acting on the full scale submarine. The necessary propulsion power equates to P_application = F_model  cdot  v_model  cdot  fracP_applicationP_model = F_model  cdot  172  fracms.","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"As a last point, it is often helpful to be able to remove and restore the dimensions of variables in a given dimensional basis:","category":"page"},{"location":"example/#","page":"Example","title":"Example","text":"julia> dim_vars = (1u\"cm/g\", 1u\"mm/s\")\n(1 cm g^-1, 1 mm s^-1)\n\njulia> dimless_vars = dimensionless.(dim_vars, model_basis)\n(6.2375, 24.95)\n\njulia> dimensionful.(dimless_vars, (u\"cm/g\", u\"mm/s\"), model_basis)\n(1.0 cm g^-1, 1.0 mm s^-1)","category":"page"}]
}
