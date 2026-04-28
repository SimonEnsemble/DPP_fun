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
	using CairoMakie, Printf, LinearAlgebra, DataFrames, CSV, StatsBase,
		PlutoUI, JLD2, Test, Graphs, MakieThemes, Statistics, SparseArrays

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

# ╔═╡ 9044a157-a02f-43e1-a693-6ed5b85e6339
@assert length(unique(raw_data[:, "molecule"])) == nrow(raw_data)

# ╔═╡ f7a06153-1014-412d-91de-dc430ee1f5e8
md"
## filter out troublesome molecules
filter out the molecules that do not constitute connected graphs. this occurs e.g. when the molecule has non-bonded components, indicated by the period in the SMILES. some of the molecules also cannot be read by `MolecularGraphs.jl`."

# ╔═╡ 76144644-aac5-4d0e-9535-20344b360239
function get_troublesome_molecules(raw_data, name_mol)
	troublesome_molecules = String[]
	for molecule in raw_data[:, name_mol]
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
	return troublesome_molecules
end

# ╔═╡ e74a68b2-cda6-466c-9ccf-78d690c335fd
 troublesome_molecules = get_troublesome_molecules(raw_data, "molecule")

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
function compute_kernel(mgs, gram_matrix_filename)
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
	return Ks
end

# ╔═╡ 5a2535c2-e6ab-4871-bb7f-3e8666499241
Ks = compute_kernel(mgs, gram_matrix_filename)

# ╔═╡ c3e93eec-73e0-48ac-b72b-2a1b9670bddd
md"use as final Gram matrix the composite kernel, in terms of exact vs. src-dst pattern matching."

# ╔═╡ 64531922-0912-4b59-af4c-5273db60506b
begin
	K1 = normalize_Gram_matrix(Ks[true])
	K2 = normalize_Gram_matrix(Ks[false])
	
	K = K1 .+ K2
end

# ╔═╡ 7381ab55-333e-4a10-885c-9448b6d9c854
CSV.write("Gram_matrix.csv", DataFrame(K, :auto))

# ╔═╡ e2fe2ee6-d09a-4ed3-9dab-c0fbdfe355e7
K_times = Ks[true] .* Ks[false]

# ╔═╡ 4af28c47-df8a-4a4b-9e91-5c9c0da5a775
K_center = center_Gram_matrix(K)

# ╔═╡ cf524ed3-a21f-4864-8722-13b97f40d4a2
@assert maximum(abs.(mean(K_center, dims=1))) < 1e-10

# ╔═╡ 71e26f93-5005-4bc1-bcb1-796f32a1f15c
@assert maximum(abs.(mean(K_center, dims=2))) < 1e-10

# ╔═╡ 81081c9d-47ba-488a-ad48-cab848938f74
function eig_gram(Kc::AbstractMatrix)
    E = eigen(Symmetric(Kc)) 
    return E.values, E.vectors
end

# ╔═╡ a66605f1-128c-4f04-913d-0da5b7ca8806
vals, vecs = eig_gram(K_center)

# ╔═╡ 7830cfd8-9148-4deb-9101-25fa08c31e9a
function k_PCA(vals, vecs, n_components=2)
	vals_m = vals[end-n_components+1:end]
    vecs_m = vecs[:, end-n_components+1:end]

	n_PCA = vecs_m .* sqrt.(vals_m)'
	PCA_importance = vals_m / sum(vals)
	
    return n_PCA[:, 1:n_components], PCA_importance
end

# ╔═╡ 3a756f2c-cbad-4bc2-a031-9d888f5d7b1f
n_components = 2

# ╔═╡ 52991237-d58d-4fe2-a6a2-73d74db47e75
pca, pca_importance = k_PCA(vals, vecs, n_components)

# ╔═╡ 68cece24-e775-466f-87f1-a4c51465df53
function get_pac_data(data, pca, pca_importance, n_components)
	pca_data = copy(data)
	for (i, name) in enumerate(Symbol.("pca_" .* string.(1:n_components)))
	    pca_data[!, name] = pca[:, end-i+1]
	end
	for (i, name) in enumerate(Symbol.("importance_pca_" .* string.(1:n_components)))
	    pca_data[!, name] = Vector{Union{Missing, Float64}}(missing, nrow(pca_data))
    	pca_data[1, name] = pca_importance[end-i+1]
	end
	return pca_data
end

# ╔═╡ 38d74afa-eba1-49b1-b409-b03cac81ffbe
pca_data = get_pac_data(data, pca, pca_importance, n_components)

# ╔═╡ b673f7d1-9020-4a00-a277-beb1972cb3ae
CSV.write("pca_data.csv", pca_data)

# ╔═╡ 6831a957-5af9-4afa-9934-7703cb22ca98
md"## look for indistinguishable molecules

warning: this takes a long time.
"

# ╔═╡ 07bd1f30-c4b9-47cd-903f-4e943a7b97ce
idps_filename = "indistinguishable_pairs.jld2"

# ╔═╡ ec7c6096-17bd-4960-9864-c34d0461114b
if ! isfile(idps_filename)
	idps = indistinguishable_pairs(K)
	jldsave(idps_filename; idps)
else
	println("loading from $idps_filename")
	idps = load(idps_filename)["idps"]
	idps
end

# ╔═╡ 2d129919-52e2-47a1-bc76-1b8ed486e452
md"indistinguishable pair browser. $(@bind duplicate_id PlutoUI.Slider(1:length(idps)))"

# ╔═╡ 398f4803-8944-4a83-aece-64312ab8b925
mgs[idps[duplicate_id][1]].smiles

# ╔═╡ 7e0299eb-a494-4730-b143-ad62e2ec43ba
mgs[idps[duplicate_id][2]].smiles

# ╔═╡ 5f4aa5c8-39d3-4a57-8e06-ec15c6a914ba
if length(duplicate_id) > 0
	viz(mgs[idps[duplicate_id][1]])
end

# ╔═╡ 7be50944-19f7-4cd7-b303-fdbee3f35d39
if length(duplicate_id) > 0
	viz(mgs[idps[duplicate_id][2]])
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

# ╔═╡ d1f679ee-7c13-4bf9-a829-50b6bc43b729
md" # surfactant"

# ╔═╡ 1636288f-7167-440c-bbba-99ba23a8bfb0
raw_s_data = CSV.read(
    download("https://raw.githubusercontent.com/BigChemistry-RobotLab/SurfPro/main/data/surfpro_literature.csv"),
    DataFrame
)

# ╔═╡ aacad155-3d18-4cf8-bfb2-a35041cc086d
begin
	s_data = raw_s_data[
	    coalesce.(raw_s_data[!, :Surfactant_Type] .!= "mixture", false),
	    :
	]
	
	s_data = s_data[
	    coalesce.(s_data[!, :Temp_Celsius] .== 25.0, false),
	    :
	]
	
	s_data = select(s_data, [:SMILES, :CMC])
	s_data = dropmissing(s_data)
end

# ╔═╡ c82abfb5-1396-4340-8490-491c92a86981
@assert length(unique(s_data[:, "SMILES"])) == nrow(s_data)

# ╔═╡ 384fb2a3-f30d-4801-99ce-f291f1390b2e
begin
	mg_nc₁ = MolGraph("CC(C(=O)OCC)S(=O)(=O)[O-].[Na+]")
	find_shortest_paths!(mg_nc₁)
	viz(mg_nc₁)
end

# ╔═╡ 4363cda4-d95b-464f-952d-97204a614796
nv(mg_nc₁.g)

# ╔═╡ 26796e29-83dc-484b-af45-69fdfbbddd4d
@time begin
	s_mgs = MolGraph.(s_data[:, "SMILES"])
	for s_mg in s_mgs
		find_shortest_paths!(s_mg)
	end
end

# ╔═╡ c8ff6352-5873-443e-881b-66de42461763
s_gram_matrix_filename = "s_gram_matrix.jld2"

# ╔═╡ c77da8c2-852e-47b3-901b-1b7d2c2d2114
s_Ks = compute_kernel(s_mgs, s_gram_matrix_filename)

# ╔═╡ 378aae95-4367-43de-aca0-c48efda8ea86
begin
	# s_K1 = normalize_Gram_matrix(s_Ks[true])
	# s_K2 = normalize_Gram_matrix(s_Ks[false])
	
	# s_K = s_K1 .+ s_K2

	s_K = normalize_Gram_matrix(s_Ks[true] .+ s_Ks[false])
end

# ╔═╡ 23eaf7f6-ad0f-4511-a0f5-b84211d11525
CSV.write("s_Gram_matrix.csv", DataFrame(s_K, :auto))

# ╔═╡ 4fbaed84-93b4-46fd-b219-38bd7ddfaf39
s_K_center = center_Gram_matrix(s_K)

# ╔═╡ ef2a7361-9a00-4053-91e3-1eec495a2df3
s_vals, s_vecs = eig_gram(s_K_center)

# ╔═╡ 33503c50-951d-4ec0-a1f7-220822a4d972
s_pca, s_pca_importance = k_PCA(s_vals, s_vecs, n_components)

# ╔═╡ c00293e0-dbb2-4317-9c8c-4ea954636189
s_pca_data = get_pac_data(s_data, s_pca, s_pca_importance, n_components)

# ╔═╡ c025e3bd-1570-464b-b5db-3b48a8a11d11
CSV.write("s_pca_data.csv", s_pca_data)

# ╔═╡ 36de3d44-0e8c-4134-9a59-c54ac17ef207
ids_dpp = CSV.read("ids_molecule_mcmc_dpp_50.csv", DataFrame)

# ╔═╡ 26e55e20-cf02-4372-988e-3175060f0265
begin
	n_nv_dpp = zeros(50)
	for (n_i, id) in enumerate(ids_dpp[!, 2])
		mg_dpp = MolGraph(data[id, :molecule])
		n_nv_dpp[n_i] = nv(mg_dpp.g)
	end
end

# ╔═╡ 51bf4fb7-e8b2-4f5c-ba33-96e4e3104d66
ids_uniform = CSV.read("ids_molecule_uniform_50.csv", DataFrame)

# ╔═╡ b9ae3fd7-17b7-4328-b0ae-40d81fb7d4e4
begin
	n_nv_uniform = zeros(50)
	for (n_i, id) in enumerate(ids_uniform[!, 2])
		mg_uniform = MolGraph(data[id, :molecule])
		n_nv_uniform[n_i] = nv(mg_uniform.g)
	end
end

# ╔═╡ 364d85d7-e5a4-4100-b25a-06476bb653cd
begin
	f = Figure()
	_ax = Axis(f[1, 1])
	hist!(n_nv_dpp, label="dpp")
	hist!(n_nv_uniform, label="uniform")
	axislegend()
	f
end

# ╔═╡ Cell order:
# ╠═cecf3058-bb8f-11f0-97f3-bda46249b7c9
# ╠═79adcae5-192f-470b-898d-7688e8041054
# ╠═d7b41ed5-53af-4c48-bc5c-3c8e1554a2ee
# ╠═713ba284-16d7-45b6-ab97-e200cb101863
# ╠═9044a157-a02f-43e1-a693-6ed5b85e6339
# ╟─f7a06153-1014-412d-91de-dc430ee1f5e8
# ╠═76144644-aac5-4d0e-9535-20344b360239
# ╠═e74a68b2-cda6-466c-9ccf-78d690c335fd
# ╠═3c85b30a-b1d7-4d84-a206-b9b47f894125
# ╠═017b6abe-fadd-4697-96ba-7454e375f668
# ╠═31e64f9f-6048-467d-8906-64298cae8db9
# ╠═d100bcf3-4495-4599-9734-6d15c0b90d64
# ╟─348ba4c3-ae1f-4737-bd2d-70306f427ebd
# ╠═f65a43dc-5710-4cc5-8cf2-8da1e24024e8
# ╠═30ff7345-d7b0-4741-bca4-6fe2a882df96
# ╠═5a2535c2-e6ab-4871-bb7f-3e8666499241
# ╟─c3e93eec-73e0-48ac-b72b-2a1b9670bddd
# ╠═64531922-0912-4b59-af4c-5273db60506b
# ╠═7381ab55-333e-4a10-885c-9448b6d9c854
# ╠═e2fe2ee6-d09a-4ed3-9dab-c0fbdfe355e7
# ╠═4af28c47-df8a-4a4b-9e91-5c9c0da5a775
# ╠═cf524ed3-a21f-4864-8722-13b97f40d4a2
# ╠═71e26f93-5005-4bc1-bcb1-796f32a1f15c
# ╠═81081c9d-47ba-488a-ad48-cab848938f74
# ╠═a66605f1-128c-4f04-913d-0da5b7ca8806
# ╠═7830cfd8-9148-4deb-9101-25fa08c31e9a
# ╠═3a756f2c-cbad-4bc2-a031-9d888f5d7b1f
# ╠═52991237-d58d-4fe2-a6a2-73d74db47e75
# ╠═68cece24-e775-466f-87f1-a4c51465df53
# ╠═38d74afa-eba1-49b1-b409-b03cac81ffbe
# ╠═b673f7d1-9020-4a00-a277-beb1972cb3ae
# ╟─6831a957-5af9-4afa-9934-7703cb22ca98
# ╠═07bd1f30-c4b9-47cd-903f-4e943a7b97ce
# ╠═ec7c6096-17bd-4960-9864-c34d0461114b
# ╟─2d129919-52e2-47a1-bc76-1b8ed486e452
# ╠═398f4803-8944-4a83-aece-64312ab8b925
# ╠═7e0299eb-a494-4730-b143-ad62e2ec43ba
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
# ╟─d1f679ee-7c13-4bf9-a829-50b6bc43b729
# ╠═1636288f-7167-440c-bbba-99ba23a8bfb0
# ╠═aacad155-3d18-4cf8-bfb2-a35041cc086d
# ╠═c82abfb5-1396-4340-8490-491c92a86981
# ╠═384fb2a3-f30d-4801-99ce-f291f1390b2e
# ╠═4363cda4-d95b-464f-952d-97204a614796
# ╠═26796e29-83dc-484b-af45-69fdfbbddd4d
# ╠═c8ff6352-5873-443e-881b-66de42461763
# ╠═c77da8c2-852e-47b3-901b-1b7d2c2d2114
# ╠═378aae95-4367-43de-aca0-c48efda8ea86
# ╠═23eaf7f6-ad0f-4511-a0f5-b84211d11525
# ╠═4fbaed84-93b4-46fd-b219-38bd7ddfaf39
# ╠═ef2a7361-9a00-4053-91e3-1eec495a2df3
# ╠═33503c50-951d-4ec0-a1f7-220822a4d972
# ╠═c00293e0-dbb2-4317-9c8c-4ea954636189
# ╠═c025e3bd-1570-464b-b5db-3b48a8a11d11
# ╠═36de3d44-0e8c-4134-9a59-c54ac17ef207
# ╠═26e55e20-cf02-4372-988e-3175060f0265
# ╠═51bf4fb7-e8b2-4f5c-ba33-96e4e3104d66
# ╠═b9ae3fd7-17b7-4328-b0ae-40d81fb7d4e4
# ╠═364d85d7-e5a4-4100-b25a-06476bb653cd
