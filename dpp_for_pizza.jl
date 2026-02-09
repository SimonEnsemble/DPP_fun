### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 84d15bb4-05f1-11f1-97e9-2d44dcee1c9d
begin
	import Pkg; Pkg.activate()
	using CairoMakie, LinearAlgebra, MakieThemes, StatsBase, Colors, Logging
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
pepperoni_radius = 2.5 # cm

# ╔═╡ 1eb2a72f-cabc-4a6a-89a4-5886be46868a
pizza_radius = 35.0 # cm

# ╔═╡ 16a85e49-635f-4d01-b484-069a0af21082
crust_radius = 3.25 # cm

# ╔═╡ 2ca3e543-a87e-4a59-812d-8401cb358e0c
γ = pizza_radius / 2 # cm

# ╔═╡ 179863df-2049-43c4-9c88-3c17c7f93828
md"# 🍕 pepperoni & the pepperoni kernel"

# ╔═╡ 42a405f0-b267-4f6a-a7f9-df90831176f0
struct Pepperoni
	# radius
	radius::Float64
	
	# center
	x::Float64
	y::Float64
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

# ╔═╡ 90d2b32a-a50d-4e8c-a0a1-00b4849213f9
function mcmc_kdpp(
	# item-item similarity matrix (psd)
	L::Matrix{Float64}, 
	# number of items to select
	k::Int;
	# MCMC steps
	n_steps::Int=10000
)
	# infer the number of items
	n = size(L)[1]

	# initial sample
	ids = sample(1:n, k, replace=false)

	# current det(Lₓ)
	current_det = det(L[ids, ids])
	
	# MCMC
	for s = 1:n_steps
		@debug "step $s"
		@debug "\tids = $ids"
		@debug "\tcurrent det L = $current_det"
		###
		# propose an exchange move
		###
		# index of current item we propose to replace
		id_replace_item = sample(1:k)

		# new item outside of ids to replace it with
		id_new_item = sample_new_item(ids, n)

		# proposed new current set
		new_ids = deepcopy(ids)
		new_ids[id_replace_item] = id_new_item

		# new det
		new_det = det(L[new_ids, new_ids])
		@debug "\tpropose replacement of:\n\t\tcurrent item ID $id_replace_item\n\t\ti.e. item $(ids[id_replace_item])\n\twith item $id_new_item\n\t\tnew det = $new_det"

		###
		# accept or reject
		###
		if rand() < new_det / current_det
			@debug "\tACCEPT!"
			ids = new_ids
			current_det = new_det
		else
			@debug "\tREJECT!"
		end

		if current_det < 0.0
			error("negative det(Lₓ; L not PSD)")
		end
	end
	return ids
end

# ╔═╡ a672c27c-8887-4e2e-92dd-cb8077a25f33
with_logger(ConsoleLogger(stdout, Logging.Debug)) do
	ids = rand(1:size(L)[1], 10)
	print(ids)
	mcmc_kdpp(L[ids, ids], 4; n_steps=20)
end

# ╔═╡ 66bb6dbb-1ed1-4754-a5e3-7445f77e8dae
ids_dpp = with_logger(ConsoleLogger(stdout, Logging.Info)) do
	mcmc_kdpp(L, n_pepperoni; n_steps=150000)
end

# ╔═╡ d719e69b-d649-4b34-91ff-d50f355df1ce
pizza_dpp = Pizza(
	pizza_radius,
	crust_radius,
	candidate_pepperonis[ids_dpp]
)

# ╔═╡ 84369ce0-b19b-4b6b-b876-3872c57a8bf9
viz_pizza(pizza_dpp)

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
# ╠═90d2b32a-a50d-4e8c-a0a1-00b4849213f9
# ╠═a672c27c-8887-4e2e-92dd-cb8077a25f33
# ╠═66bb6dbb-1ed1-4754-a5e3-7445f77e8dae
# ╠═d719e69b-d649-4b34-91ff-d50f355df1ce
# ╠═84369ce0-b19b-4b6b-b876-3872c57a8bf9
