var documenterSearchIndex = {"docs":
[{"location":"library/#Library","page":"Library","title":"Library","text":"","category":"section"},{"location":"library/#Types","page":"Library","title":"Types","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"DimBasis","category":"page"},{"location":"library/#Dimensionless.DimBasis","page":"Library","title":"Dimensionless.DimBasis","text":"DimBasis(basis_vectors...) -> DimBasis\n\nCreate a dimensional basis for a number of basis_vectors (quantities, units or dimensions). A string identifier can optionally be added to each basis vector. \n\n\n\n\n\n","category":"type"},{"location":"library/#Functions","page":"Library","title":"Functions","text":"","category":"section"},{"location":"library/","page":"Library","title":"Library","text":"dim_matrix","category":"page"},{"location":"library/#Dimensionless.dim_matrix","page":"Library","title":"Dimensionless.dim_matrix","text":"dim_matrix(basis_dims, all_values...)\n\nReturn the dimensional matrix for a set of basis dimensions basis_dims and all_values, a set of quantities, units or dimensions.\n\n\n\n\n\n","category":"function"},{"location":"library/","page":"Library","title":"Library","text":"num_of_dims","category":"page"},{"location":"library/#Dimensionless.num_of_dims","page":"Library","title":"Dimensionless.num_of_dims","text":"num_of_dims(all_vars...)\n\nReturn the number of unique dimensions for the given variables all_vars.\n\n\n\n\n\n","category":"function"},{"location":"library/","page":"Library","title":"Library","text":"num_of_dimless","category":"page"},{"location":"library/#Dimensionless.num_of_dimless","page":"Library","title":"Dimensionless.num_of_dimless","text":"num_of_dimless(all_vars...)\n\nReturn the number of dimensionless numbers that can be constructed for a problem characterized by the variables all_vars.\n\n\n\n\n\n","category":"function"},{"location":"library/","page":"Library","title":"Library","text":"fac_dimful","category":"page"},{"location":"library/#Dimensionless.fac_dimful","page":"Library","title":"Dimensionless.fac_dimful","text":"fac_dimful(unit, basis)\n\nReturn the scalar, dimensionless factor that a dimensionless value has to be multiplied with in order to translate it into the given unit in the specified basis. \n\n\n\n\n\n","category":"function"},{"location":"library/","page":"Library","title":"Library","text":"dimless","category":"page"},{"location":"library/#Dimensionless.dimless","page":"Library","title":"Dimensionless.dimless","text":"dimless(quantity, basis)\n\nMake a quantity dimensionless using a dimensional basis.\n\n\n\n\n\n","category":"function"},{"location":"library/","page":"Library","title":"Library","text":"dimful","category":"page"},{"location":"library/#Dimensionless.dimful","page":"Library","title":"Dimensionless.dimful","text":"dimful(value, unit, basis)\n\nRestore the units of a dimensionless value using a dimensional basis.\n\n\n\n\n\n","category":"function"},{"location":"library/","page":"Library","title":"Library","text":"change_basis","category":"page"},{"location":"library/#Dimensionless.change_basis","page":"Library","title":"Dimensionless.change_basis","text":"change_basis(var, basis, new_basis)\n\nTransform the specified quantity or unit var from a current basis to a new_basis.\n\n\n\n\n\n","category":"function"},{"location":"library/","page":"Library","title":"Library","text":"print_dimless","category":"page"},{"location":"library/#Dimensionless.print_dimless","page":"Library","title":"Dimensionless.print_dimless","text":"print_dimless([io, ]named_dim_vector, basis)\n\nPrint the dimensionless number that can be constructed using named_dim_vector in the specified dimensional basis.\n\n\n\n\n\n","category":"function"},{"location":"#Introduction","page":"Home","title":"Introduction","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Dimensionless is a package built on top of Unitful.jl. It contains tools to conduct dimensional analysis and solve similitude problems based on the Buckingham Π theorem.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Install the package using Julia's package manager:","category":"page"},{"location":"","page":"Home","title":"Home","text":"pkg> add Dimensionless","category":"page"},{"location":"","page":"Home","title":"Home","text":"Currently, Unitful's functionality is not exported from Dimensionless. To use Dimensionless, Unitful should be imported as well:","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Dimensionless, Unitful","category":"page"},{"location":"example/#Finding-a-dimensional-basis","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"","category":"section"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"All calculations in Dimensionless.jl revolve around the concept of dimensional bases. Such bases are constructed as a set of linearly independent dimensional variables spanning the dimensional space of the studied problem. Constructing such a dimensional basis is achieved by creating an instance of type DimBasis. For example, this is a valid basis:","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> basis = DimBasis(u\"m\", u\"kg/m^3\", u\"Pa*s\") \nDimBasis{...}","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Every variable with a unit that lies inside the space spanned by this basis can be expressed as a linear combination of the basis variables. At the same time, each dimensioned variable can be transformed into a dimensionless representation using the basis.","category":"page"},{"location":"example/#Removing-and-restoring-dimensions","page":"Finding a dimensional basis","title":"Removing and restoring dimensions","text":"","category":"section"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Arguably the two most important functions of this package are dimless and dimful. For a given dimensional basis, these two functions can be used to convert a single variables or multiple variables at once from a formulation with dimensions to a corresponding dimensionless formulation and vice versa: ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> dim_vars = (1u\"cm/g\", 1u\"mm/s\")\n(1 cm g^-1, 1 mm s^-1)\n\njulia> dimless_vars = dimless.(dim_vars, model_basis)\n(6.2375, 24.95)\n\njulia> dimful.(dimless_vars, (u\"cm/g\", u\"mm/s\"), model_basis)\n(1.0 cm g^-1, 1.0 mm s^-1)","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Such a conversion is particularly helpful to obtain easily generalizable results from complex calculations, e.g. large FEM problems. In addition, the conversion to dimensionless variables helps to normalize a problem formulations to avoid numerical instabilities. In addition to the functions mentioned above, there is also the function fac_dimful, which only returns the conversion factor from a dimensionless number to dimensionful representation in the given dimensional basis. For the first variable of the example above, this looks like this:  ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> fac_dimful(u\"cm/g\", model_basis)\n0.16032064128256512\n\njulia> fac_dimful(u\"cm/g\", model_basis) * dimless_vars[1]\n0.1","category":"page"},{"location":"example/#Similitude-example","page":"Finding a dimensional basis","title":"Similitude example","text":"","category":"section"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"A closely related use case is the systematic formulation of similitude problems. The following example is taken from the Similitude Wikipedia article:","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Consider a submarine modeled at 1/40th scale. The application operates in sea water at 0.5 °C, moving at 5 m/s. The model will be tested in fresh water at 20 °C. The variables which describe the system are:Variable Application Scaled model Units\nL (Diameter of submarine) 1.0 1/40 (m)\nV (Speed) 5.0 To be calculated (m/s)\nρ (Density) 1028.0 998.0 (kg/m^3)\nμ (Dynamic viscosity) 1.88e−3 1.0e−3 Pa·s (N s/m^2)\nF (Force) To be calculated To be measured N (kg m/s^2)Find the power required for the submarine to operate at the stated speed.","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Let's analyze this exemplary problem using Dimensionless.jl. The idea of similitude is that physical problems, despite appearing different in scale, geometry, or conditions, can exhibit similar behavior when analyzed using a dimensionless formulation of the underlying equations.   ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"For the example at hand, the mechanical power P of the submarine can generally be expressed as the product of the propulsion force F and the velocity V.  This relationship applies regardless of the specific problem parameters and the units used. However, the numerical values of the underlying variables change depending on the scaling of the problem and the employed units. For example, the numerical values and the used units in this equation formulated for the submarine in the real application (P_mathrmapp = F_mathrmapp cdot V_mathrmapp) will often be different from the values and units for the scaled submarine model (P_mathrmmodel = F_mathrmmodel cdot V_mathrmmodel). Similitude facilitates to systematically transform all variables in this physically meaningful equation to a dimensionless formulation without the constraints of specific units or scales: tildeP = tildeF cdot tildeV. ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"In the following, it will be explained how such a transformation can be achieved and how it can be used to transfer knowledge between differently scaled versions of the same physical problem.   ","category":"page"},{"location":"example/#Analyzing-the-number-of-independent-variables","page":"Finding a dimensional basis","title":"Analyzing the number of independent variables","text":"","category":"section"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"According to the Buckingham Π theorem, a problem defined by n variables with corresponding units can be transformed into an equation with n − m dimensionless quantities, where m represents the number of independent base units used. ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Translated to the exemplary submarine similitude problem stated above with n = 5 variables, Dimensionless.jl can be used to calculate that there are m = 3 independent base units in the present case:","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> m = num_of_dims(u\"m\", u\"m/s\", u\"kg/m^3\", u\"Pa*s\", u\"N\")\n3","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Correspondingly, there are n - m = 5 - 3 = 2 dimensionless variables, or Pi terms, that should be selected to effectively capture the essence of the physical phenomena under investigation: ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> num_of_dimless(u\"m\", u\"m/s\", u\"kg/m^3\", u\"Pa*s\", u\"N\")\n2","category":"page"},{"location":"example/#Finding-a-dimensional-basis-2","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"","category":"section"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"The first step to set up the similitude problem is to find a suitable set of m independent variables that serve as a dimensional basis. These variables must span the dimensional space of the studied problem.  As explained above, constructing such a dimensional basis is achieved by creating an instance of type DimBasis. While creating such a basis using the DimBasis function, it is helpful to name the used variables to allow pretty printing of the results later on:","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> basis = DimBasis(\"L\"=>u\"m\", \"ρ\"=>u\"kg/m^3\", \"μ\"=>u\"Pa*s\") \nDimBasis{...}","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"In this case, a suitable basis could be found by selecting the m = 3 variables L, ρ, and μ as basis variables. Every variable with a unit that lies inside the space spanned by this basis can be expressed as a corresponding dimensionless variable. For example, any velocity V can be made dimensionless by multiplying it with a factor L ρ / μ as can be determined using Dimensionless.jl:  ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> print_dimless(\"V\"=>u\"m/s\", basis)\nV L ρ μ^-1","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"The same is true for every force F:","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> print_dimless(\"F\"=>u\"N\", basis)\nF ρ μ^-2","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Note that the first found number is the Reynolds number mathitRe and the second number can be identified as the product of the drag coefficient c_mathrmd and mathitRe^2. In Dimensionless.jl, these numbers are mostly used under the hood to transform between formulations with and without dimensions. The actual conversions for single or multiple variables can be achieved by using the functions explained above: dimless and dimful.","category":"page"},{"location":"example/#Calculating-the-missing-entries-in-the-table","page":"Finding a dimensional basis","title":"Calculating the missing entries in the table","text":"","category":"section"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Let's use the values in the data table to construct two bases describing the submarine model and the submarine in the real application, respectively:","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> app_basis = DimBasis(\"L\"=>1u\"m\", \"ρ\"=>1028u\"kg/m^3\", \"μ\"=>1.88e-3u\"Pa*s\")\nDimBasis{...}\n\njulia> model_basis = DimBasis(\"L\"=>0.025u\"m\", \"ρ\"=>998u\"kg/m^3\", \"μ\"=>1e-3u\"Pa*s\")\nDimBasis{...}","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"It is easy to transform values from one dimensional base to the other in order to fill the missing values in the table row by row. Starting with the submarine's velocity of 5 m/s in the application case, the velocity of the scaled model should be selected to 109.6 m/s in order to behave physically similar: ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> V_model = change_basis(5u\"m/s\", app_basis, model_basis)\n109.58 m s^-1","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"This model velocity is 21.9 times higher than the velocity of the submarine in the real application:","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> change_basis(u\"m/s\", app_basis, model_basis)\n21.916172771074063","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Similarly, the factor needed to convert any given force in Newton from the case of the scaled submarine model to the real application can also be calculated: ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> change_basis(u\"N\", model_basis, app_basis)\n3.4313","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"A similar transformation can also be used to determine the power of the real submarine from experiments with the scaled submarine model. As stated above, the mechanical power of the submarine can generally be expressed as P = F cdot V. If an exemplary force of 1 N is measured in the scaled experiment with a suitably scaled velocity of 109.6 m/s as calculated above, this equates to a mechanical power of 1 N * 109.6 m/s = 109.6 W. This value can be transformed back to obtain the mechanical power of the real submarine while traveling with V = 5 m/s:  ","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"julia> P_app_fac = change_basis(109.58u\"W\", model_basis, app_basis)\n17.156144908079394 W","category":"page"},{"location":"example/","page":"Finding a dimensional basis","title":"Finding a dimensional basis","text":"Of course, at a velocity of V = 1096 mathrmmmathrms during the scaled experiment, the obtained force can easily be higher than 1 N. Generally, the mechanical power of the real submarine scales proportional to the measured force with a factor 17156144908079394 mathrmWmathrmN. ","category":"page"}]
}
