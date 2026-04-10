### A Pluto.jl notebook ###
# v0.20.23

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
		PlutoUI, JLD2, Test, Graphs, MakieThemes, Statistics,
		LinearAlgebra

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

# ╔═╡ f61ef77e-17ba-4dd6-86be-852db526bf4c
function normalize_kernel(K)
    d = sqrt.(diag(K))  
    K_norm = K ./ (d * d')  
    return K_norm
end

# ╔═╡ 64531922-0912-4b59-af4c-5273db60506b
begin
	K1 = normalize_kernel(Ks[true])
	K2 = normalize_kernel(Ks[false])
	
	K = K1 .+ K2
end

# ╔═╡ 7381ab55-333e-4a10-885c-9448b6d9c854
CSV.write("Gram_matrix.csv", DataFrame(K, :auto))

# ╔═╡ e2fe2ee6-d09a-4ed3-9dab-c0fbdfe355e7
K_times = Ks[true] .* Ks[false]

# ╔═╡ a11e1eb3-9a83-4ae5-a6ad-200c21a6ce5e
function center_gram(K::AbstractMatrix)
    @assert size(K, 2) == size(K, 1) "Gram matrix must be square"

    row_mean = mean(K, dims=2)   
    col_mean = mean(K, dims=1)   
    total_mean = mean(K)         

    return K .- row_mean .- col_mean .+ total_mean
end

# ╔═╡ 4af28c47-df8a-4a4b-9e91-5c9c0da5a775
K_center = center_gram(K)

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
function k_PCA(vals, vecs; n_components=2)
    perm = sortperm(vals, rev=true)
	vals_m = vals[perm][1:n_components]
    vecs_m = vecs[:, perm][:, 1:n_components]

	n_PCA = vecs_m .* sqrt.(vals_m)'
    return n_PCA[:, 1], n_PCA[:, 2]
end

# ╔═╡ 52991237-d58d-4fe2-a6a2-73d74db47e75
pca_1, pca_2 = k_PCA(vals, vecs)

# ╔═╡ 68cece24-e775-466f-87f1-a4c51465df53
begin
	pca_data = copy(data)
	pca_data[!, :pca_1] = pca_1
	pca_data[!, :pca_2] = pca_2
	CSV.write("pca_data.csv", pca_data)
end

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

# ╔═╡ Cell order:
# ╠═cecf3058-bb8f-11f0-97f3-bda46249b7c9
# ╠═79adcae5-192f-470b-898d-7688e8041054
# ╟─d7b41ed5-53af-4c48-bc5c-3c8e1554a2ee
# ╠═713ba284-16d7-45b6-ab97-e200cb101863
# ╠═9044a157-a02f-43e1-a693-6ed5b85e6339
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
# ╠═f61ef77e-17ba-4dd6-86be-852db526bf4c
# ╠═64531922-0912-4b59-af4c-5273db60506b
# ╠═7381ab55-333e-4a10-885c-9448b6d9c854
# ╠═e2fe2ee6-d09a-4ed3-9dab-c0fbdfe355e7
# ╠═a11e1eb3-9a83-4ae5-a6ad-200c21a6ce5e
# ╠═4af28c47-df8a-4a4b-9e91-5c9c0da5a775
# ╠═cf524ed3-a21f-4864-8722-13b97f40d4a2
# ╠═71e26f93-5005-4bc1-bcb1-796f32a1f15c
# ╠═81081c9d-47ba-488a-ad48-cab848938f74
# ╠═a66605f1-128c-4f04-913d-0da5b7ca8806
# ╠═7830cfd8-9148-4deb-9101-25fa08c31e9a
# ╠═52991237-d58d-4fe2-a6a2-73d74db47e75
# ╠═68cece24-e775-466f-87f1-a4c51465df53
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
