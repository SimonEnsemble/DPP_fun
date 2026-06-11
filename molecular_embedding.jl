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
	import Pkg; Pkg.activate() #"diffusion_env" DiffusionMap
	using ShortestPathMolecularGraphKernels, MCMCkDPP, Base.Threads,
		CairoMakie, Printf, LinearAlgebra, DataFrames, CSV, StatsBase,
		PlutoUI, JLD2, Test, Graphs, MakieThemes, Statistics, SparseArrays,
		Random

	set_theme!(ggthemr(:pale))
end

# ╔═╡ 79adcae5-192f-470b-898d-7688e8041054
TableOfContents()

# ╔═╡ 2d92e910-ee92-461b-ad49-c250787932c4
@bind dataset Select(["smells", "bees", "surfactants"])

# ╔═╡ 2cbd393f-5f89-40a1-8489-8b2bcd038856
datapath = "data_from_$dataset"

# ╔═╡ 13a4b4ce-dfdd-4125-b794-1b5667c22a00
isdir(datapath) || mkdir(datapath)

# ╔═╡ d7b41ed5-53af-4c48-bc5c-3c8e1554a2ee
md"## 👃 read in the data"

# ╔═╡ a79cd32d-cb9a-4129-ad20-2833bad1c2ab
dataset_to_link = Dict(
	"smells" => "https://raw.githubusercontent.com/SimonEnsemble/Leffingwell_Goodscents_combine/refs/heads/main/pyrfume.csv",
	"surfactants" => "https://raw.githubusercontent.com/BigChemistry-RobotLab/SurfPro/main/data/surfpro_literature.csv",
	"bees" => "https://raw.githubusercontent.com/j-adamczyk/ApisTox_dataset/refs/heads/master/outputs/dataset_final.csv"
)

# ╔═╡ 713ba284-16d7-45b6-ab97-e200cb101863
begin
	raw_data = CSV.read(download(dataset_to_link[dataset]), DataFrame)
	if dataset == "smells"
		rename!(raw_data, :molecule => :SMILES)
	end
	if dataset == "surfactants"
		filter!(row -> row["Surfactant_Type"] != "mixture", raw_data)
		filter!(row -> ! ismissing(row["Temp_Celsius"]), raw_data)
		filter!(row -> row["Temp_Celsius"] == 25.0, raw_data)
		raw_data = select(raw_data, [:SMILES, :CMC])
		raw_data = dropmissing(raw_data)
	end
end

# ╔═╡ 9044a157-a02f-43e1-a693-6ed5b85e6339
@assert length(unique(raw_data[:, "SMILES"])) == nrow(raw_data)

# ╔═╡ f7a06153-1014-412d-91de-dc430ee1f5e8
md"
### filter out troublesome molecules
filter out the molecules that do not constitute connected graphs. this occurs e.g. when the molecule has non-bonded components, indicated by the period in the SMILES. some of the molecules also cannot be read by `MolecularGraphs.jl`."

# ╔═╡ 76144644-aac5-4d0e-9535-20344b360239
function get_troublesome_molecules(smiles_list::Vector{String})
	troublesome_molecules = String[]
	for (i, smiles) in enumerate(smiles_list)
	    try
	        mol = MolGraph(smiles)
	        find_shortest_paths!(mol)
	    catch e
	        push!(troublesome_molecules, smiles)
	    end
	end
	return troublesome_molecules
end

# ╔═╡ e74a68b2-cda6-466c-9ccf-78d690c335fd
troublesome_molecules = get_troublesome_molecules(raw_data[:, "SMILES"])

# ╔═╡ 3c85b30a-b1d7-4d84-a206-b9b47f894125
begin
	_data = filter(row -> ! (row["SMILES"] in troublesome_molecules), raw_data)
	data = unique(_data, :SMILES)
end

# ╔═╡ 017b6abe-fadd-4697-96ba-7454e375f668
md"## 🚗 construct the molecular graphs and find shortest paths
"

# ╔═╡ 31e64f9f-6048-467d-8906-64298cae8db9
begin
	mgs = MolGraph.(data[:, "SMILES"])
	find_shortest_paths!.(mgs)
end

# ╔═╡ d100bcf3-4495-4599-9734-6d15c0b90d64
hist(
	[length(mg.spaths) for mg in mgs], bins=500,
	axis=(; xlabel="# shortest paths", ylabel="# molecules", 
		  limits=(0, 500, 0, nothing)
	)
)

# ╔═╡ 348ba4c3-ae1f-4737-bd2d-70306f427ebd
md"## 🐸 compute Gram matrix"

# ╔═╡ f65a43dc-5710-4cc5-8cf2-8da1e24024e8
gram_matrix_filename = joinpath(datapath, "gram_matrix.jld2")

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
K = normalize_Gram_matrix(Ks[true]) .+ normalize_Gram_matrix(Ks[false])

# ╔═╡ 50bd48b6-84e2-4907-89e9-62a74e621c39
CSV.write(joinpath(datapath, "Gram_matrix.csv"), DataFrame(K, :auto))

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
	# take the most important eigenvalues and eigenvectors
	vals_m = vals[end-n_components+1:end]
    vecs_m = vecs[:, end-n_components+1:end]

	pca = vecs_m .* sqrt.(vals_m)'
	pca_importance = vals_m / sum(vals)
	
    return pca, pca_importance
end

# ╔═╡ 3a756f2c-cbad-4bc2-a031-9d888f5d7b1f
n_components = 2

# ╔═╡ 52991237-d58d-4fe2-a6a2-73d74db47e75
pca, pca_importance = k_PCA(vals, vecs, n_components)

# ╔═╡ 91771f45-9b30-4247-9ec6-e577ba1b45a6
# ╠═╡ disabled = true
#=╠═╡
begin
	# diffusion_map
	P = Ks[true] .+ Ks[false]
	P[P .< 0] .= 0
	normalize_to_stochastic_matrix!(P)
	X = diffusion_map(P, 2; t=1)
end
  ╠═╡ =#

# ╔═╡ 68762c22-f3fb-48c2-938b-a1c57910eb4d
if dataset == "surfactants"
	CSV.write(
		joinpath(datapath, "PCA.csv"), 
		DataFrame("SMILES" => data[:, :SMILES], "CMC" => data[:, :CMC],
				  "x1" => pca[:, 1], "x2" => pca[:, 2]
				  # "d1" => X[:, 1], "d2" => X[:, 2])
	))
elseif dataset == "bees"
	CSV.write(
		joinpath(datapath, "PCA.csv"), 
		DataFrame("SMILES" => data[:, :SMILES], "toxic" => data[:, :label],
				  "x1" => pca[:, 1], "x2" => pca[:, 2] 
				  # "d1" => X[:, 1], "d2" => X[:, 2])
	))
elseif dataset == "smells"
	CSV.write(
		joinpath(datapath, "PCA.csv"), 
		DataFrame([name => data[:, name] for name in names(data)]...,
				  "x1" => pca[:, 1], "x2" => pca[:, 2]
				  # "d1" => X[:, 1], "d2" => X[:, 2])
	))
end

# ╔═╡ 8012f081-7e0f-44ab-a4e1-3d6a14632b04
CSV.write(
	joinpath(datapath, "PCA_importance.csv"),
	DataFrame("imp" => pca_importance)
)

# ╔═╡ 6831a957-5af9-4afa-9934-7703cb22ca98
md"## look for indistinguishable molecules

warning: this takes a long time.
"

# ╔═╡ 07bd1f30-c4b9-47cd-903f-4e943a7b97ce
idps_filename = joinpath(datapath, "indistinguishable_pairs.jld2")

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

# ╔═╡ bdbfdef3-b972-4c79-aa5a-fa5505e717ef
data[idps[duplicate_id][1], :]

# ╔═╡ 398f4803-8944-4a83-aece-64312ab8b925
mgs[idps[duplicate_id][1]].smiles

# ╔═╡ 532d5f28-4fcd-4f92-b2b3-aadf7ba83d95
data[idps[duplicate_id][2], :]

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

# ╔═╡ 8360ef9e-f7ca-46ff-969b-32b60b042c53
md" ## remove indistinguishable molecules"

# ╔═╡ 74fe9133-3528-4574-9b2c-664d9212e45a
md" ### remove all the indistinguishable molecules"

# ╔═╡ 7c08b0d6-5c63-4549-97f4-ce54a0ff9e6d
same_kernel_row = unique(collect(Iterators.flatten(idps)))

# ╔═╡ eaebbbee-7e2e-4cd4-a356-9a13d0e4af19
md"### randomly leave one of the indistinguishable molecules"

# ╔═╡ b3b741e3-c9f3-4a30-8f16-4529d2191ad7
begin
	ids_pairs = Vector{Vector{Int}}()
	
	for (id1, id2) in idps
	    new_group = nothing
	
	    for group in ids_pairs
	        if id1 in group || id2 in group
	            new_group = group
	            break
	        end
	    end
	
	    if new_group === nothing
	        push!(ids_pairs, [id1, id2])
	    else
	        if !(id1 in new_group)
	            push!(new_group, id1)
	        end
	
	        if !(id2 in new_group)
	            push!(new_group, id2)
	        end
	    end
	end
	ids_pairs
end

# ╔═╡ ede337f6-375f-444a-bf9b-31a8449fdc23
left_ids = [rand(ids) for ids in ids_pairs]

# ╔═╡ 5ce78e63-6c3d-49ef-bf51-a320041c381e
@bind remove_all CheckBox()

# ╔═╡ 832b7ad8-5b61-4d68-b4ba-c1dd082f2d30
if dataset != "bees"
	if remove_all
		perfect_data = data[setdiff(1:nrow(data), same_kernel_row), :]
	else
		perfect_data = data[setdiff(1:nrow(data), left_ids), :]
	end
	Random.seed!(123)  
	
	n = nrow(perfect_data)
	train_ratio = 0.8
	
	train_ids = sample(1:n, round(Int, train_ratio * n), replace=false)
	test_ids = setdiff(1:n, train_ids)
else
	train = CSV.read(joinpath(datapath, "maxmin_train.csv"), DataFrame)
	train = filter(row -> ! (row["SMILES"] in troublesome_molecules), train)
	train_ids = collect(1:nrow(train))
	test = CSV.read(joinpath(datapath, "maxmin_test.csv"), DataFrame)
	test = filter(row -> ! (row["SMILES"] in troublesome_molecules), test)
	perfect_data = vcat(train, test)
	test_ids = setdiff(collect(1:nrow(perfect_data)), train_ids)
end

# ╔═╡ be92a1df-44ea-44fd-abc9-0fd583d2d46e
gram_matrix_perfect_filename = joinpath(datapath, "gram_matrix_perfect.jld2")

# ╔═╡ f1db597b-7fdf-4973-8cc1-7dc7ac1653b6
perfect_data

# ╔═╡ f290cd95-b18e-4782-9af3-2670016ea2bb
begin
	mgs_perfect = MolGraph.(perfect_data[:, "SMILES"])
	find_shortest_paths!.(mgs_perfect)
	Ks_perfect = compute_kernel(mgs_perfect, gram_matrix_perfect_filename)
	K_perfect = normalize_Gram_matrix(
		Ks_perfect[true]) .+ normalize_Gram_matrix(Ks_perfect[false]
	)
	K_center_perfect = center_Gram_matrix(K_perfect)
	vals_perfect, vecs_perfect = eig_gram(K_center_perfect)
	pca_perfect, pca_importance_perfect = k_PCA(vals_perfect, vecs_perfect, n_components)
	P_perfect = Ks_perfect[true] .+ Ks_perfect[false]
	# P_perfect[P_perfect .< 0] .= 0
	# normalize_to_stochastic_matrix!(P_perfect)
	# X_perfect = diffusion_map(P_perfect, 2; t=1)
	CSV.write(joinpath(datapath, "Gram_matrix_perfect.csv"), DataFrame(K_perfect, :auto))
	
	if dataset == "surfactants"
	CSV.write(
		joinpath(datapath, "PCA_perfect.csv"), 
		DataFrame(
			"SMILES" => perfect_data[:, :SMILES], 
			"CMC" => perfect_data[:, :CMC],
			"x1" => pca_perfect[:, 1], "x2" => pca_perfect[:, 2], 
			# "d1" => X_perfect[:, 1], "d2" => X_perfect[:, 2])
	))
	elseif dataset == "bees"
		CSV.write(
			joinpath(datapath, "PCA_perfect.csv"), 
			DataFrame(
				"SMILES" => perfect_data[:, :SMILES], 
				"toxic" => perfect_data[:, :label],
				"x1" => pca_perfect[:, 1], "x2" => pca_perfect[:, 2], 
				# "d1" => X_perfect[:, 1], "d2" => X_perfect[:, 2])
		))
	elseif dataset == "smells"
		CSV.write(
			joinpath(datapath, "PCA_perfect.csv"), 
			DataFrame(
				[name => perfect_data[:, name] for name in names(perfect_data)]..., "x1" => pca_perfect[:, 1], "x2" => pca_perfect[:, 2], 
				# "d1" => X_perfect[:, 1], "d2" => X_perfect[:, 2])
		))
	end
	
	CSV.write(
		joinpath(datapath, "PCA_importance_perfect.csv"),
		DataFrame("imp" => pca_importance_perfect)
	)
end

# ╔═╡ 09893808-aee3-4b18-b5e8-8866f2398b31
md"## split data"

# ╔═╡ e1d7e437-1510-4131-92b8-4913390ba498
CSV.write(
		joinpath(datapath, "train_data.csv"),
		DataFrame("train" => train_ids)
)

# ╔═╡ 1d69b104-c8ff-4340-bf36-e501a47fea51
CSV.write(
		joinpath(datapath, "test_data.csv"),
		DataFrame("test" => test_ids)
)

# ╔═╡ 152425b0-a088-4790-91ce-84a1f1f0879c
md"## using DPP sample molecules"

# ╔═╡ 6f35c039-b1bc-4582-b359-b662c2036ba4
n_runs = 10

# ╔═╡ 01ab8dba-1f1b-4bf7-adac-325a83fbae6e
n_molecules = 10

# ╔═╡ e08c959b-caf2-4d79-88f0-92061c98fbbd
n_prop_swaps = floor(Int, n_molecules * log(n_molecules) * 10)

# ╔═╡ 93c08f69-2571-49b3-bb33-4a81812739de
function draw_dpp_samples(n_runs, K, n_molecules, n_prop_swaps)
	ids_dpp = [[] for n in 1:n_runs]
	@threads for n in 1:n_runs
		ids, na, np = mcmc_kdpp(K, n_molecules, n_prop_swaps=n_prop_swaps)
		@printf("%d/%d swaps accepted\n", na, np)
		ids_dpp[n] = ids
	end
	return DataFrame(ids_dpp, :auto)
end

# ╔═╡ 4b7ba037-2d1d-42ef-9f12-69899dafe860
@bind perfect_K CheckBox()

# ╔═╡ 0b7dde54-662b-48dd-86a8-e9241a3c555a
if perfect_K
	K_running = K_perfect
else
	K_running = K
end

# ╔═╡ a9871b55-fbd6-4644-93f4-a53338e935f0
ids_dpp = draw_dpp_samples(
	n_runs, K_running[train_ids, train_ids], n_molecules, n_prop_swaps
)

# ╔═╡ 08035090-1e8a-4488-838f-4aa6f1760091
CSV.write(joinpath(datapath, "ids_dpp_$(n_molecules).csv"), ids_dpp)

# ╔═╡ 40cf1543-8e17-4755-a28f-2627689bb59b
ids_uniform = [
	sample(collect(1:K_running[train_ids, train_ids].size[1]), n_molecules)
	for r = 1:n_runs
]

# ╔═╡ 72dd3d57-fa93-44e0-bdbe-dd569ac1b47d
CSV.write(
	joinpath(datapath, "ids_uniform_$(n_molecules).csv"), 
	DataFrame(ids_uniform, :auto)
)

# ╔═╡ Cell order:
# ╠═cecf3058-bb8f-11f0-97f3-bda46249b7c9
# ╠═79adcae5-192f-470b-898d-7688e8041054
# ╠═2d92e910-ee92-461b-ad49-c250787932c4
# ╠═2cbd393f-5f89-40a1-8489-8b2bcd038856
# ╠═13a4b4ce-dfdd-4125-b794-1b5667c22a00
# ╟─d7b41ed5-53af-4c48-bc5c-3c8e1554a2ee
# ╠═a79cd32d-cb9a-4129-ad20-2833bad1c2ab
# ╠═713ba284-16d7-45b6-ab97-e200cb101863
# ╠═9044a157-a02f-43e1-a693-6ed5b85e6339
# ╟─f7a06153-1014-412d-91de-dc430ee1f5e8
# ╠═76144644-aac5-4d0e-9535-20344b360239
# ╠═e74a68b2-cda6-466c-9ccf-78d690c335fd
# ╠═3c85b30a-b1d7-4d84-a206-b9b47f894125
# ╟─017b6abe-fadd-4697-96ba-7454e375f668
# ╠═31e64f9f-6048-467d-8906-64298cae8db9
# ╠═d100bcf3-4495-4599-9734-6d15c0b90d64
# ╟─348ba4c3-ae1f-4737-bd2d-70306f427ebd
# ╠═f65a43dc-5710-4cc5-8cf2-8da1e24024e8
# ╠═30ff7345-d7b0-4741-bca4-6fe2a882df96
# ╠═5a2535c2-e6ab-4871-bb7f-3e8666499241
# ╟─c3e93eec-73e0-48ac-b72b-2a1b9670bddd
# ╠═64531922-0912-4b59-af4c-5273db60506b
# ╠═50bd48b6-84e2-4907-89e9-62a74e621c39
# ╠═4af28c47-df8a-4a4b-9e91-5c9c0da5a775
# ╠═cf524ed3-a21f-4864-8722-13b97f40d4a2
# ╠═71e26f93-5005-4bc1-bcb1-796f32a1f15c
# ╠═81081c9d-47ba-488a-ad48-cab848938f74
# ╠═a66605f1-128c-4f04-913d-0da5b7ca8806
# ╠═7830cfd8-9148-4deb-9101-25fa08c31e9a
# ╠═3a756f2c-cbad-4bc2-a031-9d888f5d7b1f
# ╠═52991237-d58d-4fe2-a6a2-73d74db47e75
# ╠═91771f45-9b30-4247-9ec6-e577ba1b45a6
# ╠═68762c22-f3fb-48c2-938b-a1c57910eb4d
# ╠═8012f081-7e0f-44ab-a4e1-3d6a14632b04
# ╟─6831a957-5af9-4afa-9934-7703cb22ca98
# ╠═07bd1f30-c4b9-47cd-903f-4e943a7b97ce
# ╠═ec7c6096-17bd-4960-9864-c34d0461114b
# ╟─2d129919-52e2-47a1-bc76-1b8ed486e452
# ╠═bdbfdef3-b972-4c79-aa5a-fa5505e717ef
# ╠═398f4803-8944-4a83-aece-64312ab8b925
# ╠═532d5f28-4fcd-4f92-b2b3-aadf7ba83d95
# ╠═7e0299eb-a494-4730-b143-ad62e2ec43ba
# ╠═5f4aa5c8-39d3-4a57-8e06-ec15c6a914ba
# ╠═7be50944-19f7-4cd7-b303-fdbee3f35d39
# ╟─8360ef9e-f7ca-46ff-969b-32b60b042c53
# ╟─74fe9133-3528-4574-9b2c-664d9212e45a
# ╠═7c08b0d6-5c63-4549-97f4-ce54a0ff9e6d
# ╟─eaebbbee-7e2e-4cd4-a356-9a13d0e4af19
# ╠═b3b741e3-c9f3-4a30-8f16-4529d2191ad7
# ╠═ede337f6-375f-444a-bf9b-31a8449fdc23
# ╠═5ce78e63-6c3d-49ef-bf51-a320041c381e
# ╠═832b7ad8-5b61-4d68-b4ba-c1dd082f2d30
# ╠═be92a1df-44ea-44fd-abc9-0fd583d2d46e
# ╠═f1db597b-7fdf-4973-8cc1-7dc7ac1653b6
# ╠═f290cd95-b18e-4782-9af3-2670016ea2bb
# ╟─09893808-aee3-4b18-b5e8-8866f2398b31
# ╠═e1d7e437-1510-4131-92b8-4913390ba498
# ╠═1d69b104-c8ff-4340-bf36-e501a47fea51
# ╟─152425b0-a088-4790-91ce-84a1f1f0879c
# ╠═6f35c039-b1bc-4582-b359-b662c2036ba4
# ╠═01ab8dba-1f1b-4bf7-adac-325a83fbae6e
# ╠═e08c959b-caf2-4d79-88f0-92061c98fbbd
# ╠═93c08f69-2571-49b3-bb33-4a81812739de
# ╠═4b7ba037-2d1d-42ef-9f12-69899dafe860
# ╠═0b7dde54-662b-48dd-86a8-e9241a3c555a
# ╠═a9871b55-fbd6-4644-93f4-a53338e935f0
# ╠═08035090-1e8a-4488-838f-4aa6f1760091
# ╠═40cf1543-8e17-4755-a28f-2627689bb59b
# ╠═72dd3d57-fa93-44e0-bdbe-dd569ac1b47d
