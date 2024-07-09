### A Pluto.jl notebook ###
# v0.19.43

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ d5b9cf97-d1b3-4154-b991-3c645971f795
 import Pkg; Pkg.activate()

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using RealTimeAudioDiffEq, PlutoUI, Plots, DifferentialEquations

# ╔═╡ 5c9116e8-dd1b-4e99-b36c-6b143d16a348
html"<button onclick='present()'>present</button>"

# ╔═╡ a6c03043-4255-4295-af49-fdd32b7d8bef
html""" <h1> </h1>"""

# ╔═╡ 1ee7a8d8-f233-4603-b613-f80ea1aead30
html""" <h1> Behind the Scenes 1 </h1> """

# ╔═╡ dc66de9c-5dda-49d7-a274-88b1b176fd02
md"""
- Function Definition
- Creating DESource
- Start Stop and Reset
"""

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function duffingvdp!(du,u,p,t)
	(μ,c,α,τ) = p
	du[1] = u[2]
	du[2] = -μ*u[2]+u[1]*(1.0-u[1]^2)+c*(u[3]+u[4])
	du[3] = τ*u[4]
	du[4] = τ*(α*u[4]*(1.0-u[3]^2)-u[3])
end	

# ╔═╡ b1e2f00a-3c9b-4f35-840d-f60f1abd1a3f
source = DESource(duffingvdp!, [0.1;0.1;0.1;0.1],[1.0,0.0,0.1,0.1]; channel_map = [1,2]);

# ╔═╡ 6fc81c93-2a64-4c7e-98dd-e8a5b712f1d6
html""" <h1> Behind the Scenes 2 </h1> """

# ╔═╡ 00a46792-c686-46c2-8f0d-70b784bb808b
md"""
Pluto UI Widgets
- Parameters (a,b)
- Temporal scaling and gain (ts,g)
- Channel Mapping
- Buttons (ticks/clock, start, stop, reset)
"""

# ╔═╡ 5b266fd9-b76c-411e-883c-96825f0404dc
chan_widget = @bind chan PlutoUI.combine() do Child
	md"""
	L : $(Child("L", Select([1 => "x1",2 => "y1",3 => "x2",4 => "y2"],default=1)))
	R : $(Child("R", Select([1 => "x1",2 => "y1",3 => "x2",4 => "y2"],default=2)))
	"""
end;

# ╔═╡ cd976d06-b5e6-4115-9ab0-7abb5390347a
begin
	ticks_button = @bind ticks PlutoUI.Clock(0.1);
	start_button = @bind start Button("START");
	stop_button = @bind stop Button("STOP");
	reset_button = @bind reset Button("RESET");
end;

# ╔═╡ 2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
let 
	start
	output_device = get_default_output_device();
	start_DESource(source, output_device; buffer_size=convert(UInt32,1024))
end	

# ╔═╡ 5b8f7326-6d7f-44ac-82b9-799f03cedf46
let 
	stop
	stop_DESource(source)
end	

# ╔═╡ 6fbf80dd-a138-481b-b302-d5a432bccde9
let
	reset
	reset_state!(source)
end;

# ╔═╡ 27fd31f9-6446-4f5f-bf66-815c554432c6
buttons = PlutoUI.ExperimentalLayout.Div([start_button,stop_button,reset_button,chan_widget], style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ 83bf7282-9f19-4cfa-a6f4-9a7297ba6b43
plot_widget = PlutoUI.ExperimentalLayout.Div(ticks_button, style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ 77ce8d02-27e6-4bcc-bbe4-25c6b1ddfade
html""" <h1> Behind the Scenes 3 </h1> """

# ╔═╡ 36aa0503-e724-48e4-9218-472223fb5d0e
md"""
Setting the atomic control parameters and Plotting the solution
"""

# ╔═╡ d8e7bdec-03dc-4f1c-865f-95d4a64ce110
set_channelmap!(source,[chan.L,chan.R]);

# ╔═╡ 9c75d205-b124-4719-ab75-8475490cfe23
theme(:dark)

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
	margin: 0 auto;
    max-width: 1500px;
	padding-left: max(283px, 10%);
    padding-right: max(383px, 10%); 
}
input[type*="range"] {
	width: 38%;
}
</style>
"""

# ╔═╡ 7ea8b061-4860-4696-b140-4147bedb8863
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
par_widget = @bind par PlutoUI.combine() do Child
	md"""
	$\dot{x_1} = y_1  \qquad \qquad \qquad \qquad \qquad\dot{x_2} = \tau y_2$
	$\dot{y_1} = -\mu y_1+x_1(1-x_1^2)+c(x_2+y_2) \qquad \qquad \dot{y_2} = \tau(\alpha y_2 (1-x_2^2)-x_2)$ 
	
	
	α : $(Child("α", Slider(-0.1:0.01:5.0,default=0.1;show_value=true))) $sp
	τ : $(Child("τ", Slider(0.1:0.001:1.0,default=0.1;show_value=true))) \
	μ : $(Child("μ", Slider(0.01:0.001:0.2,default=0.1;show_value=true))) $sp
	c : $(Child("c", Slider(0.0:0.01:2.0,default=-0.11;show_value=true))) 
	"""
end;

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin
	set_param!(source,1,par.μ)
	set_param!(source,2,par.c)
	set_param!(source,3,par.α)
	set_param!(source,4,par.τ)
end;

# ╔═╡ 83cccc37-a8eb-4451-8104-dca9f24a38d3
scale_widget = @bind scale PlutoUI.combine() do Child
	md"""
	g : $(Child("g", Slider(0.0:0.1:1.0,default=0.1;show_value=true))) 	$sp
	ts : $(Child("ts", Slider(100.0:10.0:3000.0,default=1500.0;show_value=true))) 
	"""
end;

# ╔═╡ 100eefa5-3cab-4b5b-a82c-0833f4992a74
begin 
	set_ts!(source,scale.ts)
	set_gain!(source,scale.g)
end;

# ╔═╡ 6ef0c91a-93dd-429e-902f-dab242a8995a
begin
	ticks
	sol = solve(ODEProblem(duffingvdp!,source.data.state.u,(source.data.state.t,source.data.state.t+0.2*scale.ts),source.data.control.p),Tsit5());
end;

# ╔═╡ ec7f79c7-2d4c-401c-86d2-32281bc03f56
plot_phase = plot(sol,idxs=(source.data.control.channel_map[1],source.data.control.channel_map[2]),c=:yellow,label="",border=:none,size=(500,300));

# ╔═╡ 239d2ce8-9d3f-445e-9414-93e0977146e4
PlutoUI.ExperimentalLayout.vbox([
	par_widget,
	scale_widget,
	buttons,
	plot_phase,
	plot_widget
])

# ╔═╡ Cell order:
# ╟─5c9116e8-dd1b-4e99-b36c-6b143d16a348
# ╟─a6c03043-4255-4295-af49-fdd32b7d8bef
# ╟─239d2ce8-9d3f-445e-9414-93e0977146e4
# ╟─1ee7a8d8-f233-4603-b613-f80ea1aead30
# ╠═dc66de9c-5dda-49d7-a274-88b1b176fd02
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═b1e2f00a-3c9b-4f35-840d-f60f1abd1a3f
# ╠═2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╠═5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═6fbf80dd-a138-481b-b302-d5a432bccde9
# ╟─6fc81c93-2a64-4c7e-98dd-e8a5b712f1d6
# ╟─00a46792-c686-46c2-8f0d-70b784bb808b
# ╠═1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╠═83cccc37-a8eb-4451-8104-dca9f24a38d3
# ╠═5b266fd9-b76c-411e-883c-96825f0404dc
# ╠═cd976d06-b5e6-4115-9ab0-7abb5390347a
# ╠═27fd31f9-6446-4f5f-bf66-815c554432c6
# ╠═83bf7282-9f19-4cfa-a6f4-9a7297ba6b43
# ╟─77ce8d02-27e6-4bcc-bbe4-25c6b1ddfade
# ╟─36aa0503-e724-48e4-9218-472223fb5d0e
# ╠═bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╠═100eefa5-3cab-4b5b-a82c-0833f4992a74
# ╠═d8e7bdec-03dc-4f1c-865f-95d4a64ce110
# ╠═6ef0c91a-93dd-429e-902f-dab242a8995a
# ╠═ec7f79c7-2d4c-401c-86d2-32281bc03f56
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═d5b9cf97-d1b3-4154-b991-3c645971f795
# ╟─9c75d205-b124-4719-ab75-8475490cfe23
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╟─7ea8b061-4860-4696-b140-4147bedb8863
