### A Pluto.jl notebook ###
# v0.17.7

using Markdown
using InteractiveUtils

# ╔═╡ d0487496-e685-46a0-b954-d93eebf2b5fc
md"""
### Problem set three worked on by Madeline Rafa (mhr77) and Kathleen Clifford (kmc363). 
"""

# ╔═╡ 2babb04c-7f14-4c30-a4f9-348ed31d4fbf
md"""
### Flux Balance Analysis of the Urea Cycle
We are interested in the level of information required for flux balance analysis calculations. In particular, does a flux calculation significantly improve when including different levels of information in the bounds and the dilution terms? To explore this question, calculate the flux distribution through the urea cycle in a population of [HL60 cells](https://www.atcc.org/products/ccl-240?matchtype=b&network=g&device=c&adposition=&keyword=hl60%20cells&gclid=EAIaIQobChMIyfXmmfiG9gIVlDY4Ch3MYQ-YEAAYAiAAEgLPfvD_BwE) growing in batch culture with a doubling time of $\tau_{d}$ = 20 hr for four cases:

* __Case 1__: Ignore dilution and bounds metabolite effects
* __Case 2__: Ignore dilution but include metabolite data in the bounds
* __Case 3__: Include dilution but ignore metabolite data in the bounds
* __Case 4__: Include both dilution and metabolite data in bounds

The teaching team found that the rate of oxygen uptake by the urea cycle is $q_{O2}$ = 0.25 $\mu$mol/gDW-s.   

__Assumptions__: 
* The $k_{cat}$'s: EC:3.5.3.1 = 249 s$^{-1}$; EC:2.1.33 = 88.1 s$^{-1}$; EC:4.3.2.1 = 34.5 s$^{-1}$; EC:6.3.4.5 = 203 s$^{-1}$ and EC:1.14.13.39 = 13.7 s$^{-1}$;
* The steady-state enzyme concentration in the pathway is uniform E $\simeq$ 0.01 μmol/gDW;
* Use [Park et al. Nat Chem Biol 12:482-9, 2016](https://pubmed.ncbi.nlm.nih.gov/27159581/) for $K_{m}$ and metabolite concentrations measurements
* All enzymes are maximally active (ignore allosteric effects).
"""

# ╔═╡ 54d23ccb-211c-4716-a80c-2c087198232d
md"""
##### Load/build stoichiometric array
"""

# ╔═╡ 47a3f3fe-425d-434c-be06-1cd36a56fe25
md"""
### C1: Ignore dilution and bounds metabolite effects
"""

# ╔═╡ deafe11f-9065-43f4-b500-c64ff41026bd
md"""
##### Computation C1
"""

# ╔═╡ 12eef416-b3ab-4f03-8ac0-a9783dcca9bd
md"""
### C2: Ignore dilution but include metabolite data in the bounds

Let's use general multisubstrate kinetics to compute flux bounds for the urea cycle problem.
Suppose the irreversible rate $v_{i}$ is dependent upon susbtrates $S_{j},j=1,2,\dots,\mathcal{S}$, then the multiple saturation kinetic form is given by:

$$v_{i} = V_{max,i}\left[\frac{\prod_{j}\frac{S_{j}}{K_{j}}}{\prod_{j}\left(1+\frac{S_{j}}{K_{j}}\right) - 1}\right]\qquad{i=1,2,\dots,\mathcal{R}}$$

where $V_{max,i}$ denote the maximum reaction rate (units: concentration/time), $S_{j}$ denotes the substrate concentration (units: concentration) and $K_{j}$ denotes the saturation constant for substrate $j$.

* [Liebermeister W, Klipp E. Bringing metabolic networks to life: convenience rate law and thermodynamic constraints. Theor Biol Med Model. 2006;3:41. Published 2006 Dec 15. doi:10.1186/1742-4682-3-41](https://www.ncbi.nlm.nih.gov/labs/pmc/articles/PMC1781438/)

We need to calculate the substrate/$K_{m}$ ratios to calculate the rates. To do that, we can use the metabolite concentration and $K_{m}$ data from:

* [Park JO, Rubin SA, Xu YF, Amador-Noguez D, Fan J, Shlomi T, Rabinowitz JD. Metabolite concentrations, fluxes and free energies imply efficient enzyme usage. Nat Chem Biol. 2016 Jul;12(7):482-9. doi: 10.1038/nchembio.2077. Epub 2016 May 2. PMID: 27159581; PMCID: PMC4912430.](https://pubmed.ncbi.nlm.nih.gov/27159581/)
"""

# ╔═╡ e2917845-d956-45f0-a379-ea826caf7d88
md"""
##### Are any of the enzymes in our system in the Park dataset?
To answer this question, let's take advantage of the features of [DataFrames](https://dataframes.juliadata.org/stable/) and do in-memory filtering of the dataset using the [filter](https://dataframes.juliadata.org/stable/lib/functions/#Base.filter) command. Let's filter on the [enzyme commission number (ec number)](https://en.wikipedia.org/wiki/Enzyme_Commission_number).
"""

# ╔═╡ 70239f9d-1ea8-4ad2-92a3-126cd99de4f0
md"""
### C3: Include dilution but ignore metabolite data in the bounds

We have (up to now) ignored the dilution to growth terms (as often done in practice). Typically, we do not have access to intracellular metabolite concentration data; thus, including the dilution terms makes the problem more difficult. In particular, when including the dilution due to growth terms, the intracellular metabolite balances are given by (at steady state):

$$\sum_{j=1}^{\mathcal{R}}\sigma_{ij}v_{j} - \mu{x}_{i} = 0\qquad{i=1,2,\dots,\mathcal{M}}$$

where $\sigma_{ij}$ is the stoichiometric coefficient of metabolite $i$ in reaction $j$, $v_{j}$ denotes the $j$ metabolic flux, $\mu$ denotes the specific growth rate (units: hr$^{-1}$) and $x_{i}$ denotes the concentration of metabolite $i$ (units: $\star$mol/gDW-hr). 


Let's use the [Park et al dataset](https://pubmed.ncbi.nlm.nih.gov/27159581/) to estimate values for $x_{i}$. However, how do we change the problem structure to account for the dilution terms and convert the concentration values to cell mass-specific units?
"""

# ╔═╡ cb7d76b8-85f9-4886-a4f3-3f9fb82e42dc
md"""
##### Propose a strategy (or strategies) for accounting for dilution terms 
"""

# ╔═╡ d3a2474e-f502-4247-bba6-1f7d5f88f12d
md"""
Fill me in using markdown and latex
"""

# ╔═╡ db005134-6e21-474a-9539-91faf6d0a4a4
md"""
##### Compute specific growth rate
"""

# ╔═╡ a7df4d41-1c0f-422e-8f21-c84b812ac3cd
md"""
##### Concentration conversion factor
Assume HL60 cells are spherical. Use [bionumbers]() to formulate a concentration conversion factor between M and $\mu$mol/gDW concentration units.
"""

# ╔═╡ b552fa38-cdc2-4f46-917f-4cac78694c86
md"""
##### Formulate the concentration array
"""

# ╔═╡ c0569dfd-784f-43c5-966d-6ea2b25b9992
md"""
##### Computation C3
"""

# ╔═╡ eaa02eea-2c86-4e05-bbe3-46174bb9c7d2
md"""
### C4: Include both dilution and metabolite data in bounds
"""

# ╔═╡ 58a35e5c-2f3c-4818-a290-7b8ae509320f
function ingredients(path::String)

    # this is from the Julia source code (evalfile in base/loading.jl)
    # but with the modification that it returns the module instead of the last object
    name = Symbol("lib")
    m = Module(name)
    Core.eval(m,
        Expr(:toplevel,
            :(eval(x) = $(Expr(:core, :eval))($name, x)),
            :(include(x) = $(Expr(:top, :include))($name, x)),
            :(include(mapexpr::Function, x) = $(Expr(:top, :include))(mapexpr, $name, x)),
            :(include($path))))
    m
end

# ╔═╡ acd56faf-5f71-435f-ba5b-10faa7cc7b61
begin

    # import some packages -
    using PlutoUI
    using PrettyTables
    using LinearAlgebra
    using Plots
    using GLPK
    using DataFrames
    using CSV
	using StatsPlots
	using Statistics
	using StatsBase
	
    # setup paths -
    const _PATH_TO_NOTEBOOK = pwd()
    const _PATH_TO_SRC = joinpath(_PATH_TO_NOTEBOOK, "src")
    const _PATH_TO_DATA = joinpath(_PATH_TO_NOTEBOOK, "data")
	const _PATH_TO_CONFIG = joinpath(_PATH_TO_NOTEBOOK, "config")
	const _PATH_TO_FIGS = joinpath(_PATH_TO_NOTEBOOK, "figs")
	
    # load the PS3 code lib -
    lib = ingredients(joinpath(_PATH_TO_SRC, "Include.jl"))

    # return -
    nothing
end

# ╔═╡ 0fe1a22d-a333-4834-88be-92180c39bbb3
PlutoUI.LocalResource(joinpath(_PATH_TO_FIGS,"Fig-Urea-cycle.png"))

# ╔═╡ 24f4f6d4-7f78-40fa-bb3e-505c0940ec91
begin

	# Load the stoichiometric matrix -

    # load/parse the network file -
    _PATH_TO_NETWORK_FILE = joinpath(_PATH_TO_CONFIG,"Network.net")
    list_of_reactions = lib.read_reaction_file(_PATH_TO_NETWORK_FILE);

    # build the stoichiometric array -
    (S, mna, rna) = lib.build_stoichiometric_matrix(list_of_reactions; expand=true);

	# get expanded reaction strings -
	expanded_reaction_array = lib.expand_reversible_reactions(list_of_reactions)

    # return -
    nothing
end

# ╔═╡ e953d072-7b5a-4310-b259-4e0cabbd6291
(ℳ,ℛ) = size(S)

# ╔═╡ 7558ee4f-87ca-4fa2-89bb-7a476c5e1252
begin

	# E - 
	Eₒ = 0.01 # units: μmol/gDW
	u = ones(6)
	E = Eₒ.*u
	
	# setup flux bounds array -
	flux_bounds_array = zeros(ℛ,2)
	flux_bounds_array[:,2] .= 100.0 # default value is 100 for flux units: μmol/gDW-s
	flux_bounds_array[1,2] = 203.0*E[1]
	flux_bounds_array[2,2] = 34.5*E[2]
	flux_bounds_array[3,2] = 249.0*E[3]
	flux_bounds_array[4,2] = 88.1*E[4]
	flux_bounds_array[5,2] = 13.7*E[5]
	flux_bounds_array[6,2] = 13.7*E[6]

	# O2 uptake -
	flux_bounds_array[15,1] = 0.25

	# setup species bounds array -
	species_bounds_array = zeros(ℳ,2)

	# setup the objective coefficient array -
	c_vector = zeros(ℛ)

	# what is the index of the urea flux?
	idx_b₄ = findall(x->x=="b4", rna)[1]
	c_vector[idx_b₄] = -1

	# compute the flux -
	result_case_1 = lib.flux(S,flux_bounds_array,species_bounds_array,c_vector);

	# show -
	nothing
end

# ╔═╡ 4dce933f-693c-4788-9c77-63784871a4c0
# check:
with_terminal() do
	ef = result_case_1.exit_flag
	sf = result_case_1.status_flag

	if (ef == 0.0 && sf == 5.0)
		println("Case 1: Optimal solution found. exit flag = $(ef) and status flag = $(sf)")
	else
		println("Case 1: Ooops! Check your problem setup. exit flag = $(ef) and status flag = $(sf)")
	end
end

# ╔═╡ 70d11054-4d8a-4dcc-a8cd-1f762c2dd9b8
let

	# get flux values from the result -
	calculated_flux_array = result_case_1.calculated_flux_array

	# build flux table -
	flux_table = Array{Any,2}(undef,ℛ,4)

	# populate -
	for i ∈ 1:ℛ
		flux_table[i,1] = i
		flux_table[i,2] = rna[i]
		flux_table[i,3] = calculated_flux_array[i]*(1) # units: μmol/gDW-hr
		flux_table[i,4] = expanded_reaction_array[i]
	end

	# setup header -
	header_row = (["i","name","flux","reaction"],["","","μmol/gDW-s",""])
	
	with_terminal() do
		pretty_table(flux_table; header=header_row, alignment=:l)
	end
end

# ╔═╡ 880ce921-308b-4e3c-8bd9-fb9d07d18bd3
begin

    # load the metabolite data into a DataFrame -
    path_to_data_file = joinpath(_PATH_TO_DATA, "Metabolite-NatChemBio-12-482-2016.csv")
    metabolite_table = CSV.read(path_to_data_file, DataFrame)
end

# ╔═╡ 85e1ac31-90cf-48da-b4d7-b6c009328084
begin

	# what are the ec numbers in our system? (add as strings)
	ec_number_array = ["6.3.4.5","4.3.2.1","3.5.3.1","2.1.3.3","1.14.13.39"]
	# ec_number_array = [];
	
	# use the filter command on the df to check: do we have our ec numbers?
	filter_col_key = Symbol("EC Number")
	df = filter(filter_col_key=>x->in(x,ec_number_array), metabolite_table)
end

# ╔═╡ 4b601722-c993-424c-9411-794c3b5b9362
begin

	# E - 
	Eₒ2 = 0.01 # units: μmol/gDW
	u2 = ones(6)
	E2 = Eₒ.*u
	
	# setup flux bounds array -
	#assuming all other metabolites are 1
	flux_bounds_array2 = zeros(ℛ,2)
	flux_bounds_array2[:,2] .= 100.0 # default value is 100 for flux units: μmol/gDW-s
	flux_bounds_array2[1,2] = 203.0*E2[1]*((11.9*96.7)/(((1+11.9)*(1+96.7))-1))
	flux_bounds_array2[2,2] = 34.5*E2[2]*((0.0915*0.0852)/(((1+0.0915)*(1+0.0852))-1))
	flux_bounds_array2[3,2] = 249.0*E2[3]*((0.165*1.39)/(((1+0.165)*(1+1.39))-1))
	flux_bounds_array2[4,2] = 88.1*E2[4]*((2.81*0.0119)/(((1+2.81)*(1+0.0119))-1))
	flux_bounds_array2[5,2] = 13.7*E2[5]*73.1/((1+73.1)-1)
	flux_bounds_array2[6,2] = 13.7*E2[6]*73.1/((1+73.1)-1)
		

	# O2 uptake -
	flux_bounds_array2[15,1] = 0.25

	# compute the flux -
	result_case_2 = lib.flux(S,flux_bounds_array2,species_bounds_array,c_vector);

	# show -
	nothing
end

# ╔═╡ 1ad8d299-7895-42ea-947a-7f8e42fe5efd
# check:
with_terminal() do
	ef2 = result_case_2.exit_flag
	sf2 = result_case_2.status_flag

	if (ef2 == 0.0 && sf2 == 4.0)
		println("Case 2: Optimal solution found. exit flag = $(ef2) and status flag = $(sf2)")
	else
		println("Case 2: Ooops! Check your problem setup. exit flag = $(ef2) and status flag = $(sf2)")
	end
end

# ╔═╡ 9af91fc4-717f-492f-bbe9-72500f26d358
let

	# get flux values from the result -
	calculated_flux_array2 = result_case_2.calculated_flux_array

	# build flux table -
	flux_table2 = Array{Any,2}(undef,ℛ,4)

	# populate -
	for i ∈ 1:ℛ
		flux_table2[i,1] = i
		flux_table2[i,2] = rna[i]
		flux_table2[i,3] = calculated_flux_array2[i]*(1) # units: μmol/gDW-hr
		flux_table2[i,4] = expanded_reaction_array[i]
	end

	# setup header -
	header_row2 = (["i","name","flux","reaction"],["","","μmol/gDW-s",""])
	
	with_terminal() do
		pretty_table(flux_table2; header=header_row2, alignment=:l)
	end
end

# ╔═╡ 34682c11-8139-4e7e-8a8b-bcfafa8ce628
begin

	# default value for μ -
	td = 72000 #units: s
	μ = log(2)/td # units: 1/s
	
	
	with_terminal() do
		println("Specific growth rate μ = $(μ) s⁻¹")
	end
end

# ╔═╡ e891e546-f275-4cc8-beb4-0447094f01b2
begin

	# Need: convert mol/L to μmol/gDW
	# mol/L * L 
	#multiply by cell weight
	#convert between grams and grams dry weight (since cells are 70% water)
	#convert to micro moles
	diameter = 10 #
	CF = 3327
	

	with_terminal() do
		println("Conversion factor (M -> μmol/gDW) CF = $(CF)")
	end
end

# ╔═╡ c1878272-dd06-46b3-84f3-6d05c459688a
begin
    #add an additional column with xi's to stoich matrix
	#add mu to rate vector
	#the two multiplied should equal mu
	
	# replacement flag -
	replacement_strategy_flag = 1 # use 0's, 1 if geometric mean
	
	# initialize -
	x_concentration_array = zeros(ℳ,2) # default: missing value = 0.0
	if (replacement_strategy_flag == 1)
		
		# compute geometric mean metabolite concentration -
		gmean_x = geomean(metabolite_table[!,:Concentration]) # units: M
		x_concentration_array = (gmean_x)*ones(ℳ,2)*CF # units: μmol/gDW
	end
	
	# build table of metabolites concentrations (M) for the metabolites in the our model -
	metabolites_in_our_model_array = [
		"AMP","ATP", "M_Carbamoyl_phosphate_c", "M_Diphosphate_c", 
		"fumarate", "M_H2O_c", "M_H_c",
		"arginine","aspartate","citrulline","ornithine",
		"M_N-(L-Arginino)succinate_c","nadph","nadp","M_Nitric_oxide_c",
		"M_Orthophosphate_c","M_Oxygen_c","M_Urea_c"
	];

	for (i,metabolite) ∈ enumerate(metabolites_in_our_model_array)

		# look up the metabolite value -
		df_metabolite = filter([:Metabolite,:Organism]=>(x,y)->
			(x == metabolite && y=="Homo sapiens"), metabolite_table)
		
		# is the df empty?
		if (isempty(df_metabolite) == false)
			x_concentration_array[i,1] = df_metabolite[1,:Concentration]
			x_concentration_array[i,2] = CF*df_metabolite[1,:Concentration]
		else

			if (replacement_strategy_flag == 1)
				x_concentration_array[i,1] = gmean_x
			end
		end
	end
end

# ╔═╡ a7fb5ffb-d3ec-4ddc-909d-b9e1a7920321
let

	# setup the state array -
	state_array = Array{Any,2}(undef,ℳ,4)

	for (i,metabolite) ∈ enumerate(metabolites_in_our_model_array)

		state_array[i,1] = i
		state_array[i,2] = metabolites_in_our_model_array[i]
		state_array[i,3] = x_concentration_array[i,1]
		state_array[i,4] = x_concentration_array[i,2]
	end

	# setup the header -
	header_row = (["index","metabolite","concentration", "concentration"],
		["","","M","μmol/gDW"])
	
	with_terminal() do
		
		# make the table -
		pretty_table(state_array; header=header_row)
	end
end

# ╔═╡ 7e11b015-ab19-4aba-80b1-a413ab13439f
begin

	# setup concentration array -
	conc_array = copy(x_concentration_array[:,2])
	
	# computation for C3 goes here ...
	S_C3 = S
	flux_bounds_array_C3 = flux_bounds_array
	
	# compute the flux -
	result_case_3 = lib.flux(S_C3,flux_bounds_array_C3,species_bounds_array,c_vector);

	# show -
	nothing
end

# ╔═╡ fa075bcb-6411-4f4d-aa67-f2c02b4381b7
# check:
with_terminal() do
	ef = result_case_3.exit_flag
	sf = result_case_3.status_flag

	if (ef == 0.0 && sf == 5.0)
		println("Case 3: Optimal solution found. exit flag = $(ef) and status flag = $(sf)")
	else
		println("Case 3: Ooops! Check your problem setup. exit flag = $(ef) and status flag = $(sf)")
	end
end

# ╔═╡ 00c37b03-821a-4bcb-b3d2-c52109dcfb13
let

	# get flux values from the result -
	calculated_flux_array = result_case_3.calculated_flux_array
	if (size(calculated_flux_array,1) == ℛ)
		calculated_flux_array = [calculated_flux_array ; 0.0]
	end
	
	# build flux table -
	flux_table = Array{Any,2}(undef,(ℛ+1),4)

	rna = [rna ; "μ"]
	expanded_reaction_array = [expanded_reaction_array ; "precursors -> cellmass"]

	# populate -
	for i ∈ 1:(ℛ+1)
		flux_table[i,1] = i
		flux_table[i,2] = rna[i]
		flux_table[i,3] = calculated_flux_array[i]*(1) # units: μmol/gDW-hr
		flux_table[i,4] = expanded_reaction_array[i]
	end

	# setup header -
	header_row = (["i","name","flux","reaction"],["","","μmol/gDW-s",""])
	
	with_terminal() do
		pretty_table(flux_table; header=header_row, alignment=:l)
	end
end

# ╔═╡ a53d7f87-376d-4721-86bc-89f22307fb7e
let

	residual_array = result_case_3.uptake_array

	# build flux table -
	state_table = Array{Any,2}(undef,ℳ,5)

	# populate -
	for i ∈ 1:ℳ
		state_table[i,1] = i
		state_table[i,2] = mna[i]
		state_table[i,3] = x_concentration_array[i,2] # units: μmol/gDW
		state_table[i,4] = conc_array[i]
		state_table[i,5] = residual_array[i]*(1) # units: μmol/gDW-s
	end

	# setup header -
	header_row = (["i","name","conc (measured)","conc (model)", "residual"],["","","μmol/gDW","μmol/gDW","μmol/gDW-s"])
	
	with_terminal() do
		pretty_table(state_table; header=header_row, alignment=:l)
	end
	
end

# ╔═╡ fca0ae0f-607a-4144-9b30-e237df7f43af
begin
	
	# background color plots -
    background_color_outside = RGB(1.0, 1.0, 1.0)
    background_color = RGB(0.99, 0.98, 0.96)
    CB_BLUE = RGB(68 / 255, 119 / 255, 170 / 255)
    CB_LBLUE = RGB(102 / 255, 204 / 255, 238 / 255)
    CB_GRAY = RGB(187 / 255, 187 / 255, 187 / 255)
    CB_RED = RGB(238 / 255, 102 / 255, 119 / 255)

	# show -
	nothing
end

# ╔═╡ 07ca2450-8a84-4e71-adcf-91b12a1544ee
begin

	# grab data from the [Met]/Km col -
	col_key = Symbol("[Met]/Km")
	length = 1000
	saturation_data_set = sort(metabolite_table[!,col_key])[1:length]
	number_of_bins = round(Int64,0.25*length)
	
	# make a histogram plot -
	stephist(saturation_data_set, bins = number_of_bins, normed = :true,
                background_color = background_color, background_color_outside = background_color_outside,
                foreground_color_minor_grid = RGB(1.0, 1.0, 1.0),
                lw = 2, c = CB_RED, foreground_color_legend = nothing, label = "N = $(length)")

	# label the axis -
	xlabel!("Metabolite saturation xᵢ/Kₘ (dimensionless)",fontsize=18)
	ylabel!("Instance count", fontsize=18)
end

# ╔═╡ f472e85e-8f51-11ec-25e8-e94287a542b6
html"""
<style>
main {
    max-width: 900px;
    width: 70%;
    margin: auto;
    font-family: "Roboto, monospace";
}

a {
    color: blue;
    text-decoration: none;
}
</style>"""

# ╔═╡ 1d16ad51-bd30-48d6-ac19-b875b1bfe530
html"""
<script>
	// initialize -
	var section = 0;
	var subsection = 0;
	var subsubsection = 0;
	var headers = document.querySelectorAll('h3, h5, h6');
	
	// main loop -
	for (var i=0; i < headers.length; i++) {
	    
		var header = headers[i];
	    var text = header.innerText;
	    var original = header.getAttribute("text-original");
	    if (original === null) {
	        
			// Save original header text
	        header.setAttribute("text-original", text);
	    } else {
	        
			// Replace with original text before adding section number
	        text = header.getAttribute("text-original");
	    }
	
	    var numbering = "";
	    switch (header.tagName) {
	        case 'H3':
	            section += 1;
	            numbering = section + ".";
	            subsection = 0;
	            break;
	        case 'H5':
	            subsection += 1;
	            numbering = section + "." + subsection;
	            break;
			case 'H6':
	            subsubsection += 1;
	            numbering = section + "." + subsection + "." + subsubsection;
	            break;
	    }
		// update the header text 
		header.innerText = numbering + " " + text;
	};
</script>"""

# ╔═╡ ce4ed3a0-f75c-4807-a3c6-3642e373e951
TableOfContents(title="📚 Table of Contents", indent=true, depth=5, aside=true)

# ╔═╡ Cell order:
# ╟─d0487496-e685-46a0-b954-d93eebf2b5fc
# ╟─2babb04c-7f14-4c30-a4f9-348ed31d4fbf
# ╠═0fe1a22d-a333-4834-88be-92180c39bbb3
# ╟─54d23ccb-211c-4716-a80c-2c087198232d
# ╠═24f4f6d4-7f78-40fa-bb3e-505c0940ec91
# ╠═e953d072-7b5a-4310-b259-4e0cabbd6291
# ╟─47a3f3fe-425d-434c-be06-1cd36a56fe25
# ╟─deafe11f-9065-43f4-b500-c64ff41026bd
# ╠═7558ee4f-87ca-4fa2-89bb-7a476c5e1252
# ╠═4dce933f-693c-4788-9c77-63784871a4c0
# ╠═70d11054-4d8a-4dcc-a8cd-1f762c2dd9b8
# ╟─12eef416-b3ab-4f03-8ac0-a9783dcca9bd
# ╠═880ce921-308b-4e3c-8bd9-fb9d07d18bd3
# ╟─07ca2450-8a84-4e71-adcf-91b12a1544ee
# ╟─e2917845-d956-45f0-a379-ea826caf7d88
# ╠═85e1ac31-90cf-48da-b4d7-b6c009328084
# ╠═4b601722-c993-424c-9411-794c3b5b9362
# ╠═1ad8d299-7895-42ea-947a-7f8e42fe5efd
# ╠═9af91fc4-717f-492f-bbe9-72500f26d358
# ╟─70239f9d-1ea8-4ad2-92a3-126cd99de4f0
# ╟─cb7d76b8-85f9-4886-a4f3-3f9fb82e42dc
# ╟─d3a2474e-f502-4247-bba6-1f7d5f88f12d
# ╟─db005134-6e21-474a-9539-91faf6d0a4a4
# ╠═34682c11-8139-4e7e-8a8b-bcfafa8ce628
# ╟─a7df4d41-1c0f-422e-8f21-c84b812ac3cd
# ╠═e891e546-f275-4cc8-beb4-0447094f01b2
# ╟─b552fa38-cdc2-4f46-917f-4cac78694c86
# ╠═c1878272-dd06-46b3-84f3-6d05c459688a
# ╟─a7fb5ffb-d3ec-4ddc-909d-b9e1a7920321
# ╟─c0569dfd-784f-43c5-966d-6ea2b25b9992
# ╟─7e11b015-ab19-4aba-80b1-a413ab13439f
# ╟─fa075bcb-6411-4f4d-aa67-f2c02b4381b7
# ╠═00c37b03-821a-4bcb-b3d2-c52109dcfb13
# ╟─a53d7f87-376d-4721-86bc-89f22307fb7e
# ╟─eaa02eea-2c86-4e05-bbe3-46174bb9c7d2
# ╠═acd56faf-5f71-435f-ba5b-10faa7cc7b61
# ╠═fca0ae0f-607a-4144-9b30-e237df7f43af
# ╠═58a35e5c-2f3c-4818-a290-7b8ae509320f
# ╠═f472e85e-8f51-11ec-25e8-e94287a542b6
# ╠═1d16ad51-bd30-48d6-ac19-b875b1bfe530
# ╠═ce4ed3a0-f75c-4807-a3c6-3642e373e951
