### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ 84d15bb4-05f1-11f1-97e9-2d44dcee1c9d
begin
	import Pkg; Pkg.activate()
	using CairoMakie, LinearAlgebra, MakieThemes, StatsBase, Colors, Logging, CSV, DataFrames, Base.Threads
	using MCMCkDPP
	# modifying the plot scheme
	# see here for other themes
	#  https://makieorg.github.io/MakieThemes.jl/dev/themes/ggthemr/
	local my_theme = :pale
	
	set_theme!(ggthemr(my_theme))
	update_theme!(
		fontsize=20, linewidth=4, #backgroundcolor=:white,
		Axis=(; bottomspinevisible=false, leftspinevisible=false, 
			titlefont=:regular)
	)
	
	colors = parse.(Colorant, MakieThemes.GGThemr.ColorTheme[my_theme][:swatch])

	pepperoni_color = "crimson"
end

# ╔═╡ 83fa5b7f-1049-42ed-b774-027b48732760
md"# 🍕 settings"

# ╔═╡ 3b1213bb-4eb6-48fb-b2b4-b99d95a29fe7
pepperoni_radius = 3.0 # cm

# ╔═╡ 1eb2a72f-cabc-4a6a-89a4-5886be46868a
pizza_radius = 35.0 # cm

# ╔═╡ 16a85e49-635f-4d01-b484-069a0af21082
crust_radius = 3.25 # cm

# ╔═╡ 2ca3e543-a87e-4a59-812d-8401cb358e0c
γ = pizza_radius # cm

# ╔═╡ 179863df-2049-43c4-9c88-3c17c7f93828
md"# 🍕 pepperoni & the pepperoni kernel"

# ╔═╡ 42a405f0-b267-4f6a-a7f9-df90831176f0
begin
	struct Pepperoni
		# radius
		radius::Float64
		
		# center
		x::Float64
		y::Float64
	end
	
	Pepperoni() = Pepperoni(NaN, NaN, NaN)
end

# ╔═╡ c710b421-45c1-4b24-80fe-d9bc3ec7d4b6
function distance(pᵢ::Pepperoni, pⱼ::Pepperoni)
	return sqrt(
		(pᵢ.x - pⱼ.x) ^ 2 + (pᵢ.y - pⱼ.y) ^ 2
	)
end

# ╔═╡ 11b28263-2bba-417c-9f16-df437a9c3a8c
function k(
	pᵢ::Pepperoni, pⱼ::Pepperoni; # two pepperonis in question
	γ::Float64=γ # bandwidth (length-scale of repulsion force)
)
	# pepperoni-pepperoni distance
	d = distance(pᵢ, pⱼ)

	# SE kernel
	return exp(-(d / γ) ^ 2)
	# Matern ν=3/2 kernel
	# return (1+sqrt(3) * d / γ) * exp(-sqrt(3) * d / γ)
end

# ╔═╡ bceefbcd-9047-4430-b2ac-b1211b5d9060
function viz_kernel(pepperoni_radius::Float64, pizza_diameter::Float64)
	p_ref = Pepperoni(pepperoni_radius, 0.0, 0.0)
	
	ds = range(0.0, pizza_diameter, 200)
	ps = [Pepperoni(pepperoni_radius, d, 0.0) for d in ds]

	fig = Figure()
	ax = Axis(fig[1, 1], xlabel="||xᵢ - xⱼ||", ylabel="k(xᵢ, xⱼ)")
	hlines!(ax, 0.0, color="black", linewidth=1)
	vlines!(ax, 0.0, color="black", linewidth=1)
	lines!(
		ds,
		[k(p_ref, p) for p in ps]
	)

	vlines!(ax, 2 * pepperoni_radius, label="2r", linewidth=1, linestyle=:dot)
	xlims!(-1.0, pizza_diameter)
	fig
end

# ╔═╡ 4901684e-aa96-4896-9c67-7ef4c3478a85
viz_kernel(pepperoni_radius, pizza_radius)

# ╔═╡ 53b51914-94d6-464c-87ea-c7433c7cf4c7
md"# 🍕 a pizza"

# ╔═╡ c654512b-c669-44a8-82f9-13fee51b86b4
begin
	struct Pizza
		# overall radius
		radius::Float64 # cm
		
		# crust radius
		crust_radius::Float64
	
		# pepperoni
		pepperonis::Vector{Pepperoni}
	end
	
	function Pizza(radius::Float64, crust_radius::Float64)
		# uniform sampling
		return Pizza(radius, crust_radius, Pepperoni[])
	end
end

# ╔═╡ a6abf8c1-3f01-4bf1-ab3b-94b76e71b080
pizza = Pizza(pizza_radius, crust_radius)

# ╔═╡ 6b4f3254-6296-4401-8468-4b342f3b8950
function viz_pizza(
	pizza::Pizza; 
	candidate_pepperonis::Union{Nothing, Vector{Pepperoni}}=nothing,
	title::String=""
)
	the_ticks = -pizza.radius:10:pizza.radius
		
	fig = Figure()
	ax = Axis(
		fig[1, 1], 
		aspect=DataAspect(), 
		xticks=the_ticks, 
		yticks=the_ticks,
		title=title
	)
	# hidedecorations!(ax)
	
	# the pizza
	poly!(Circle(Point2f(0, 0), pizza.radius), color=:peru)

	# the cheese (on top of the sauce)
	poly!(ax, Circle(Point2f(0, 0), pizza.radius - pizza.crust_radius), color=:gold)

	# candidate points
	if ! isnothing(candidate_pepperonis)
		xs = [p.x for p in candidate_pepperonis]
		ys = [p.y for p in candidate_pepperonis]
		scatter!(
			xs, ys,
			marker=:circle, color="black", markersize=5
		)
	end

	# pepperonis
	for pepperoni in pizza.pepperonis
		the_circle = Circle(Point2f(pepperoni.x, pepperoni.y), pepperoni.radius)
		poly!(ax, the_circle, color=pepperoni_color, strokecolor="black", strokewidth=2)
	end
	
	return fig
end

# ╔═╡ 7cacd5d5-8ff7-4800-9d87-fa9ad67d81a3
viz_pizza(pizza, title="my pizza")

# ╔═╡ bbc44853-21d5-413f-84bc-2bc331d77fb8
md"# 🍕 candidate spots for the pepperoni"

# ╔═╡ bd4f3b71-b8b9-4c10-bca9-f98c95fe2a2d
function gen_candidates(pizza::Pizza, pepperoni_radius::Float64, n::Int=50)
	xs = range(-pizza.radius, pizza.radius, n)
	ys = range(-pizza.radius, pizza.radius, n)

	r_acceptable = pizza.radius - pizza.crust_radius - pepperoni_radius

	pepperonis = Pepperoni[]
	for i = 1:n
		for j = 1:n
			if xs[i] ^ 2 + ys[j] ^ 2 < r_acceptable ^ 2
				push!(
					pepperonis, 
					Pepperoni(pepperoni_radius, xs[i], ys[j])
				)
			end
		end
	end
	
	return pepperonis
end

# ╔═╡ ca6986a2-fe3d-4095-b1a8-f391bb93cdc6
candidate_pepperonis = gen_candidates(pizza, pepperoni_radius)

# ╔═╡ fbbf0014-a056-4b37-8b94-1089390fe728
viz_pizza(
	pizza, 
	candidate_pepperonis=candidate_pepperonis, 
	title="candidate pepperoni placements"
)

# ╔═╡ f425afb3-a4d4-4cf7-98f6-8e73be0e388d
md"# 🍕 uniform sampling"

# ╔═╡ bb880a76-cdb6-4634-8b74-1daa6b5edc1f
n_pepperoni = 24

# ╔═╡ 4b6def45-9bfd-41c2-925e-9a6ab69ea3ef
pizza_uniform = Pizza(
	pizza_radius,
	crust_radius,
	sample(candidate_pepperonis, n_pepperoni)
)

# ╔═╡ 5903739f-e31d-4c50-bde4-93c0f6675970
viz_pizza(pizza_uniform, title="uniform design")

# ╔═╡ 941f16c2-a778-4904-8d3d-9bb031b5ae9d
md"# 🍕 k-DPP"

# ╔═╡ 4fa54dbb-7f98-4190-af25-4e694b94831e
L = [k(pᵢ, pⱼ) for pᵢ in candidate_pepperonis, pⱼ in candidate_pepperonis]

# ╔═╡ eb6e3024-2364-462b-b39b-9beba3b43937
function viz_L(L::Matrix{Float64})
	fig = Figure()
	ax = Axis(
		fig[1, 1], xlabel="candidate pepepperoni", ylabel="candidate pepepperoni",
		aspect=DataAspect()
	)
	h = heatmap!(L)
	Colorbar(fig[1, 2], h, label="k(pᵢ, pⱼ)")	
	fig
end

# ╔═╡ 8efd12e2-9975-45d6-919a-f278178b22af
viz_L(L)

# ╔═╡ 19b16e44-560b-41db-9f5f-dd1f866ee296
L2 = ["($i, $j)" for i = 1:5, j = 1:5]

# ╔═╡ 9ff71c87-f48f-4a3a-bb52-1236883b889f
L2[[1, 3], [1, 3]]

# ╔═╡ c1ae2e45-e7a2-4988-a31a-3bb5a0eb83a0
function sample_new_item(
	# current IDs
	ids::Vector{Int64},
	# number of items
	n::Int
)
	id = sample(1:n)
	if id in ids
		return sample_new_item(ids, n) # try again
	end
	return id
end

# ╔═╡ a672c27c-8887-4e2e-92dd-cb8077a25f33
with_logger(ConsoleLogger(stdout, Logging.Debug)) do
	ids = rand(1:size(L)[1], 10)
	print(ids)
	mcmc_kdpp(L[ids, ids], 4; n_steps=20)
end

# ╔═╡ 3dbfdbf5-f3d7-4c11-b02b-f02f981f3901
ids_mcmc_dpp = mcmc_kdpp(L, n_pepperoni)

# ╔═╡ 764dc86c-86ab-40fd-b0fe-f5ea1a43ebd1
pizza_dpp = Pizza(
	pizza_radius,
	crust_radius,
	candidate_pepperonis[ids_mcmc_dpp]
)

# ╔═╡ 93669cfa-7cb4-4296-a676-94e158378d93
viz_pizza(pizza_dpp)

# ╔═╡ 59cf0e20-59a5-4d85-ac5c-4957097da90a
ids_greedy_dpp = greedy_kdpp(L, n_pepperoni)

# ╔═╡ 58932740-27a6-4eb4-95e9-4b5c63402c9f
viz_pizza(Pizza(
	pizza_radius,
	crust_radius,
	candidate_pepperonis[ids_greedy_dpp]
))

# ╔═╡ 826e374d-1f6e-4e72-a9f6-165a8ae0686c
md"## statistical analysis of DPP vs uniform"

# ╔═╡ 10135e58-4de0-4a4a-96ac-8477e0ff9352
function get_pairwise_pepperoni_distances(pizza::Pizza)
	n = length(pizza.pepperonis)
	return [
		distance(pizza.pepperonis[i], pizza.pepperonis[j]) 
		for i = 1:n for j = (i+1):n
	]
end

# ╔═╡ b5ea1703-3864-4c5d-96d6-e50ad00d076f
overlap(pᵢ, pⱼ) = distance(pᵢ, pⱼ) < pᵢ.radius + pⱼ.radius

# ╔═╡ eac7d700-7495-46d7-9343-8d9c1d962a5b
function nb_pairs_overlapping(pizza::Pizza)
	n = length(pizza.pepperonis)
	n_overlaps = 0
	for i = 1:n
		for j = (i+1):n
			if overlap(pizza.pepperonis[i], pizza.pepperonis[j])
				n_overlaps += 1
			end
		end
	end
	return n_overlaps
end

# ╔═╡ e8701378-ef4c-4493-b1c3-bd9078f717c1
nb_pairs_overlapping(pizza_uniform)

# ╔═╡ f6867941-5f28-4030-b9c8-0e29402b2cc6
nb_pairs_overlapping(pizza_dpp)

# ╔═╡ 1fc8d00f-4308-4fa5-93f0-602144d163f3
begin
	local fig = Figure()
	local ax = Axis(
		fig[1, 1], 
		xlabel="pairwise distance [cm]", 
		ylabel="# pepperoni pairs"
	)


	density!(
		get_pairwise_pepperoni_distances(pizza_uniform), 
		label="uniform", color=(colors[1], 0.5), boundary=(0.0, 100)
	)
	density!(
		get_pairwise_pepperoni_distances(pizza_dpp), 
		label="DPP", color=(colors[2], 0.5)
	)
	
	axislegend()
	fig
end

# ╔═╡ 156bb81a-6ccb-4222-8ee8-eea5c4804fd1
md"## using DPP sample molecules"

# ╔═╡ 473e3deb-9c5b-48f2-9cc8-a0f73941a48b
n_run = 100

# ╔═╡ 685219e6-e162-4121-87f0-31ec818bb170
Threads.nthreads()

# ╔═╡ e17f6889-b748-4b92-8d11-ea39c2c72d26
function n_run_compute(method, n_run, molecule_L, n_molecules)
	
	n_run_ids = [[] for n in 1:n_run]
	@threads for n in 1:n_run
		n_run_ids[n] = method(molecule_L, n_molecules)
	end
	return DataFrame(n_run_ids, :auto)
end

# ╔═╡ 96d7b6dd-5796-45e4-a1fd-c47bf01cc26e
function psd(K; eps=1e-10)
    Ksym = Symmetric((K + K') / 2)
    F = eigen(Ksym)
    vals = max.(F.values, eps)
    return Matrix(Symmetric(F.vectors * Diagonal(vals) * F.vectors'))
end

# ╔═╡ 82a7e6af-d8d1-4007-924d-cd480ba033f7
begin
	molecule_L = Matrix(CSV.read("Gram_matrix.csv", DataFrame))
	# molecule_L = psd(molecule_L)
end

# ╔═╡ 9bd786e6-3ca9-4d14-aa4c-3b055da6bf8b
s_molecule_L = Matrix(CSV.read("s_Gram_matrix.csv", DataFrame))

# ╔═╡ bf6f47f8-b9a8-4523-9a4f-f4d6f0f37ac2
n_molecules = 100

# ╔═╡ 6d5d1589-e8b3-4695-ac86-7ea62018904a
# ╠═╡ disabled = true
#=╠═╡
ids_molecule_mcmc_dpp = n_run_compute(mcmc_kdpp, n_run, molecule_L, n_molecules)
  ╠═╡ =#

# ╔═╡ 03475552-5b4f-48f4-86e2-0ab218ecc747
s_ids_molecule_mcmc_dpp = n_run_compute(mcmc_kdpp, n_run, s_molecule_L, n_molecules)

# ╔═╡ f436d6a5-f299-4a44-ab3e-05c55c7a77a2
CSV.write("ss_ids_molecule_mcmc_dpp_$(n_molecules).csv", s_ids_molecule_mcmc_dpp)

# ╔═╡ 6c278377-8aa2-410f-a6a8-1d349102d58c
#=╠═╡
CSV.write("ids_molecule_mcmc_dpp_$(n_molecules).csv", ids_molecule_mcmc_dpp)
  ╠═╡ =#

# ╔═╡ bb02beaa-58d1-4a84-a496-acf2733612e3
md" ### uniform sample molecules"

# ╔═╡ ef8375c5-f353-4cf4-a3f4-ab30c33fba98
function uniform(molecule_L, n_molecules)
	return sample(collect(1:molecule_L.size[1]), n_molecules)
end

# ╔═╡ 1605a807-12d2-4fd2-a45a-6aa8f5bd2f6a
# ╠═╡ disabled = true
#=╠═╡
ids_molecule_uniform = n_run_compute(uniform, n_run, molecule_L, n_molecules)
  ╠═╡ =#

# ╔═╡ b3410109-aa61-4b20-a484-4443ca20fd3b
#=╠═╡
CSV.write("ids_molecule_uniform_$(n_molecules).csv", ids_molecule_uniform)
  ╠═╡ =#

# ╔═╡ cbc0b68c-487f-4a7e-beea-ce9b6241d147
s_ids_molecule_uniform = n_run_compute(uniform, n_run, s_molecule_L, n_molecules)

# ╔═╡ 0f8811c6-abf1-4251-a9aa-0fbb96bc34bb
CSV.write("ss_ids_molecule_uniform_$(n_molecules).csv", s_ids_molecule_uniform)

# ╔═╡ Cell order:
# ╠═84d15bb4-05f1-11f1-97e9-2d44dcee1c9d
# ╟─83fa5b7f-1049-42ed-b774-027b48732760
# ╠═3b1213bb-4eb6-48fb-b2b4-b99d95a29fe7
# ╠═1eb2a72f-cabc-4a6a-89a4-5886be46868a
# ╠═16a85e49-635f-4d01-b484-069a0af21082
# ╠═2ca3e543-a87e-4a59-812d-8401cb358e0c
# ╟─179863df-2049-43c4-9c88-3c17c7f93828
# ╠═42a405f0-b267-4f6a-a7f9-df90831176f0
# ╠═c710b421-45c1-4b24-80fe-d9bc3ec7d4b6
# ╠═11b28263-2bba-417c-9f16-df437a9c3a8c
# ╠═bceefbcd-9047-4430-b2ac-b1211b5d9060
# ╠═4901684e-aa96-4896-9c67-7ef4c3478a85
# ╟─53b51914-94d6-464c-87ea-c7433c7cf4c7
# ╠═c654512b-c669-44a8-82f9-13fee51b86b4
# ╠═a6abf8c1-3f01-4bf1-ab3b-94b76e71b080
# ╠═6b4f3254-6296-4401-8468-4b342f3b8950
# ╠═7cacd5d5-8ff7-4800-9d87-fa9ad67d81a3
# ╟─bbc44853-21d5-413f-84bc-2bc331d77fb8
# ╠═bd4f3b71-b8b9-4c10-bca9-f98c95fe2a2d
# ╠═ca6986a2-fe3d-4095-b1a8-f391bb93cdc6
# ╠═fbbf0014-a056-4b37-8b94-1089390fe728
# ╟─f425afb3-a4d4-4cf7-98f6-8e73be0e388d
# ╠═bb880a76-cdb6-4634-8b74-1daa6b5edc1f
# ╠═4b6def45-9bfd-41c2-925e-9a6ab69ea3ef
# ╠═5903739f-e31d-4c50-bde4-93c0f6675970
# ╟─941f16c2-a778-4904-8d3d-9bb031b5ae9d
# ╠═4fa54dbb-7f98-4190-af25-4e694b94831e
# ╠═eb6e3024-2364-462b-b39b-9beba3b43937
# ╠═8efd12e2-9975-45d6-919a-f278178b22af
# ╠═19b16e44-560b-41db-9f5f-dd1f866ee296
# ╠═9ff71c87-f48f-4a3a-bb52-1236883b889f
# ╠═c1ae2e45-e7a2-4988-a31a-3bb5a0eb83a0
# ╠═a672c27c-8887-4e2e-92dd-cb8077a25f33
# ╠═3dbfdbf5-f3d7-4c11-b02b-f02f981f3901
# ╠═764dc86c-86ab-40fd-b0fe-f5ea1a43ebd1
# ╠═93669cfa-7cb4-4296-a676-94e158378d93
# ╠═59cf0e20-59a5-4d85-ac5c-4957097da90a
# ╠═58932740-27a6-4eb4-95e9-4b5c63402c9f
# ╠═826e374d-1f6e-4e72-a9f6-165a8ae0686c
# ╠═10135e58-4de0-4a4a-96ac-8477e0ff9352
# ╠═b5ea1703-3864-4c5d-96d6-e50ad00d076f
# ╠═eac7d700-7495-46d7-9343-8d9c1d962a5b
# ╠═e8701378-ef4c-4493-b1c3-bd9078f717c1
# ╠═f6867941-5f28-4030-b9c8-0e29402b2cc6
# ╠═1fc8d00f-4308-4fa5-93f0-602144d163f3
# ╟─156bb81a-6ccb-4222-8ee8-eea5c4804fd1
# ╠═473e3deb-9c5b-48f2-9cc8-a0f73941a48b
# ╠═685219e6-e162-4121-87f0-31ec818bb170
# ╠═e17f6889-b748-4b92-8d11-ea39c2c72d26
# ╠═96d7b6dd-5796-45e4-a1fd-c47bf01cc26e
# ╠═82a7e6af-d8d1-4007-924d-cd480ba033f7
# ╠═9bd786e6-3ca9-4d14-aa4c-3b055da6bf8b
# ╠═bf6f47f8-b9a8-4523-9a4f-f4d6f0f37ac2
# ╠═6d5d1589-e8b3-4695-ac86-7ea62018904a
# ╠═03475552-5b4f-48f4-86e2-0ab218ecc747
# ╠═f436d6a5-f299-4a44-ab3e-05c55c7a77a2
# ╠═6c278377-8aa2-410f-a6a8-1d349102d58c
# ╟─bb02beaa-58d1-4a84-a496-acf2733612e3
# ╠═ef8375c5-f353-4cf4-a3f4-ab30c33fba98
# ╠═1605a807-12d2-4fd2-a45a-6aa8f5bd2f6a
# ╠═b3410109-aa61-4b20-a484-4443ca20fd3b
# ╠═cbc0b68c-487f-4a7e-beea-ce9b6241d147
# ╠═0f8811c6-abf1-4251-a9aa-0fbb96bc34bb
