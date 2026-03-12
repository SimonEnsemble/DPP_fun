### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ cecf3058-bb8f-11f0-97f3-bda46249b7c9
begin
	import Pkg; Pkg.activate()
	using ShortestPathMolecularGraphKernels
	using CairoMakie, Printf, LinearAlgebra, DataFrames, CSV, 
		PlutoUI, JLD2, Test, Graphs, MakieThemes
	set_theme!(ggthemr(:pale))
end

# ╔═╡ 79adcae5-192f-470b-898d-7688e8041054
TableOfContents()

# ╔═╡ d7b41ed5-53af-4c48-bc5c-3c8e1554a2ee
md"# 👃 read in the Leffingwell + Goodscents compound olfactory perception data"

# ╔═╡ 713ba284-16d7-45b6-ab97-e200cb101863
raw_data = CSV.read(
	download(
		"https://raw.githubusercontent.com/SimonEnsemble/Leffingwell_Goodscents_combine/refs/heads/main/pyrfume.csv"
	), 
	DataFrame
)

# ╔═╡ f7a06153-1014-412d-91de-dc430ee1f5e8
md"
## filter out troublesome molecules
filter out the molecules that do not constitute connected graphs. this occurs e.g. when the molecule has non-bonded components, indicated by the period in the SMILES. some of the molecules also cannot be read by `MolecularGraphs.jl`."

# ╔═╡ 76144644-aac5-4d0e-9535-20344b360239
begin
	troublesome_molecules = String[]
	for molecule in raw_data[:, "molecule"]
		try 
			mg = MolGraph(molecule)
			if ! is_connected(mg.g)
				push!(troublesome_molecules, molecule)
			end
		catch ArgumentError
			println("ArgumentError for ", molecule)
			push!(troublesome_molecules, molecule)
		end
	end
	troublesome_molecules
end

# ╔═╡ 3c85b30a-b1d7-4d84-a206-b9b47f894125
data = filter(row -> ! (row["molecule"] in troublesome_molecules), raw_data)

# ╔═╡ 017b6abe-fadd-4697-96ba-7454e375f668
md"# 🧪🚗 construct the molecular graphs and find shortest paths
"

# ╔═╡ 31e64f9f-6048-467d-8906-64298cae8db9
@time begin
	mgs = MolGraph.(data[:, "molecule"])
	for mg in mgs
		find_shortest_paths!(mg)
	end
end

# ╔═╡ d100bcf3-4495-4599-9734-6d15c0b90d64
hist(
	[length(mg.spaths) for mg in mgs], bins=500,
	axis=(; xlabel="# shortest paths", ylabel="# molecules", 
		  limits=(0, 500, 0, nothing)
	)
)

# ╔═╡ 348ba4c3-ae1f-4737-bd2d-70306f427ebd
md"# 🐸 compute Gram matrix"

# ╔═╡ f65a43dc-5710-4cc5-8cf2-8da1e24024e8
gram_matrix_filename = "gram_matrix.jld2" # save cuz expensive

# ╔═╡ 30ff7345-d7b0-4741-bca4-6fe2a882df96
begin
	if ! isfile(gram_matrix_filename)
		println("$gram_matrix_filename DNE. computing!")
		Ks = Dict()
		for exact_seq_matching in [true, false]
			@time Ks[exact_seq_matching] = compute_Gram_matrix(
				mgs, exact_seq_matching
			)
		end
		jldsave(gram_matrix_filename; Ks)
		Ks
	else
		println("loading from $gram_matrix_filename")
		Ks = load(gram_matrix_filename)["Ks"]
	end
	Ks
end

# ╔═╡ c3e93eec-73e0-48ac-b72b-2a1b9670bddd
md"use as final Gram matrix the composite kernel, in terms of exact vs. src-dst pattern matching."

# ╔═╡ 47939e72-61e0-4377-af30-fa387d815e72
K = Ks[true] .+ Ks[false]

# ╔═╡ 6831a957-5af9-4afa-9934-7703cb22ca98
md"## look for indistinguishable molecules

warning: this takes a long time.
"

# ╔═╡ 9c6e9fa9-1b51-4df7-9040-50e2493e0f1a
md"look for duplicates? $(@bind look_for_duplicates PlutoUI.CheckBox())"

# ╔═╡ ec7c6096-17bd-4960-9864-c34d0461114b
if look_for_duplicates
	# loop over pairs of molecules
	ids_duplicates = []
	for i = 1:length(mgs)
		for j = i+1:length(mgs)
			if isapprox(K[:, i], K[:, j])
				push!(ids_duplicates, (i, j))
			end
		end
	end
	ids_duplicates
end

# ╔═╡ c7c99cf3-1e2d-485b-b4de-4ec287514f65
viz(mgs[1]) # viz duplicates here if pertinent

# ╔═╡ 2d129919-52e2-47a1-bc76-1b8ed486e452
md"indistinguishable pair browser. $(@bind duplicate_id PlutoUI.Slider(1:length(ids_indistinguisable)))"

# ╔═╡ 5f4aa5c8-39d3-4a57-8e06-ec15c6a914ba
if length(ids_indistinguisable) > 0
	viz(mgs[ids_indistinguisable[duplicate_id][1]])
end

# ╔═╡ 7be50944-19f7-4cd7-b303-fdbee3f35d39
if length(ids_indistinguisable) > 0
	viz(mgs[ids_indistinguisable[duplicate_id][2]])
end

# ╔═╡ a157eb93-7027-4dee-9f3b-07e2178253b9
md"# embedding via kernel-PCA

save embedding

center kernel matrix

normalize somehow? not sure.
"

# ╔═╡ e0682429-a157-4be2-bf7f-a11a6375b54b
viz(mgs[19])

# ╔═╡ 462f1af9-b341-4fd7-8a00-f0c08e49d279
md"# ⛵ explore"

# ╔═╡ 61183670-4199-40ff-9bec-e9b86239eae9
@bind a_smell PlutoUI.Select(sort(names(data)[2:end]))

# ╔═╡ e28cfbe3-ae20-42f8-8abf-3ca28d7a6956
md"interesting smells: vanilla, mint, fish, wine, chamomile, popcorn, currant"

# ╔═╡ 8fe20699-4f78-4a5a-ae81-ec77824b302c
begin
	fig = Figure()
	ax = Axis(fig[1, 1], aspect=DataAspect())
	
	ids_draw = data[:, a_smell] .== 1
	scatter!(
		X[.! ids_draw, 1], X[.! ids_draw, 2], 
		color="black",
		markersize=5
	)
	scatter!(
		X[ids_draw, 1], X[ids_draw, 2], 
		color="red",
		markersize=5
	)
	
	ax.title = a_smell

	fig
end

# ╔═╡ 775809e2-41c5-4ae5-8e80-e354b67b13b3
md"the data we are visualizing"

# ╔═╡ d545ad1f-b49e-48a1-ac66-891851d908b4
data[ids_filter, my_smell]

# ╔═╡ Cell order:
# ╠═cecf3058-bb8f-11f0-97f3-bda46249b7c9
# ╠═79adcae5-192f-470b-898d-7688e8041054
# ╟─d7b41ed5-53af-4c48-bc5c-3c8e1554a2ee
# ╠═713ba284-16d7-45b6-ab97-e200cb101863
# ╟─f7a06153-1014-412d-91de-dc430ee1f5e8
# ╠═76144644-aac5-4d0e-9535-20344b360239
# ╠═3c85b30a-b1d7-4d84-a206-b9b47f894125
# ╟─017b6abe-fadd-4697-96ba-7454e375f668
# ╠═31e64f9f-6048-467d-8906-64298cae8db9
# ╠═d100bcf3-4495-4599-9734-6d15c0b90d64
# ╟─348ba4c3-ae1f-4737-bd2d-70306f427ebd
# ╠═f65a43dc-5710-4cc5-8cf2-8da1e24024e8
# ╠═30ff7345-d7b0-4741-bca4-6fe2a882df96
# ╟─c3e93eec-73e0-48ac-b72b-2a1b9670bddd
# ╠═47939e72-61e0-4377-af30-fa387d815e72
# ╟─6831a957-5af9-4afa-9934-7703cb22ca98
# ╟─9c6e9fa9-1b51-4df7-9040-50e2493e0f1a
# ╠═ec7c6096-17bd-4960-9864-c34d0461114b
# ╠═c7c99cf3-1e2d-485b-b4de-4ec287514f65
# ╟─2d129919-52e2-47a1-bc76-1b8ed486e452
# ╠═5f4aa5c8-39d3-4a57-8e06-ec15c6a914ba
# ╠═7be50944-19f7-4cd7-b303-fdbee3f35d39
# ╠═a157eb93-7027-4dee-9f3b-07e2178253b9
# ╠═e0682429-a157-4be2-bf7f-a11a6375b54b
# ╟─462f1af9-b341-4fd7-8a00-f0c08e49d279
# ╠═61183670-4199-40ff-9bec-e9b86239eae9
# ╟─e28cfbe3-ae20-42f8-8abf-3ca28d7a6956
# ╠═8fe20699-4f78-4a5a-ae81-ec77824b302c
# ╟─775809e2-41c5-4ae5-8e80-e354b67b13b3
# ╠═d545ad1f-b49e-48a1-ac66-891851d908b4
