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

# ╔═╡ 2e014313-3cdc-43f8-8762-e972a8bd8c28
 import Pkg; Pkg.activate()

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using RealTimeAudioDiffEq, PlutoUI, Plots, DifferentialEquations

# ╔═╡ 71f589fb-2164-4fce-b000-a82970a90248
md"""
# Auto - tuned Reed Model

As an example of a self-oscillator, inspired by a musical instrument, that undergooes an Hopf bifurcation we choose modified version of the Rayleigh model of a clarinet reed:

$\dot{x}  =  y$
$\dot{y}  =  \mu y - y^3 - kx - \sigma k^{\frac{3}{2}} x^3$

when $\sigma=0$ it becomes the original Rayleigh model and undergoes a Supercritical Hopf bifurcation when $\mu$ crosses zero from the negatives. The angular frequency of the oscillator at the bifurcation is $\omega = \sqrt{k}$. As the value of $\mu$ increases the amplitude of the stable limit cycle increases proportional to $\sqrt{\mu}$ and the frequency decreases since a relaxation orbit take place close to the "N" nullcline $kx = \mu y - y^3$ 

When $\sigma>0$ the nullcline changes its shape and becomes "faster", then the frequency of the relaxation cycle increases. For some intermediate value of $\sigma$ this increase in frequency compensates the decrease due to to the jump from the Hopf to the relaxation cycle. 
This intermediate $\sigma$ value depends on $k$ in an approximate power law form with exponent $3/2$. 

The value $\sigma=0.18$ seems adequate for $k$ values within the range $(0.1,2)$

"""

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function nreed!(du,u,p,t)
    (μ,k,σ) = p
    du[1] = u[2] 
	du[2] = u[2]*(μ-u[2]^2)-k*u[1]*(1.0+σ*sqrt(k)*u[1]^2)
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = DESource(nreed!, [0.1, 0.], [-0.1, 1.0,0.0,0.0]; channel_map = [2,2]);

# ╔═╡ 41dabe99-2401-4c34-9341-f8f1c1627c64
chan_widget = @bind chan PlutoUI.combine() do Child
	md"""
	L : $(Child("L", Select([1 => "x",2 => "y"],default=1)))
	R : $(Child("R", Select([1 => "x",2 => "y"],default=2)))
	"""
end;

# ╔═╡ 3e3b7c18-335c-48a3-afde-039590c50095
begin
	ticks_button = @bind ticks Clock(0.1);
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

# ╔═╡ 9512cc09-48cd-40c8-99d1-9299c96fcf8e
let
	reset
	reset_state!(source)
end;

# ╔═╡ 248f0ef5-b462-427a-8a5c-09991da405e6
buttons = PlutoUI.ExperimentalLayout.Div([start_button,stop_button,reset_button,chan_widget], style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ c5916ed9-fcc5-4aec-8369-9dd74a38f966
plot_widget = PlutoUI.ExperimentalLayout.Div(ticks_button, style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ c000bc1a-3782-4924-b3f6-0f403645d95d
set_channelmap!(source,[chan.L,chan.R]);

# ╔═╡ 572007b9-8d58-4c4b-bda7-24cd734674f4
begin
	ticks
	sol = solve(ODEProblem(nreed!,source.data.state.u,(source.data.state.t,source.data.state.t+60.0),source.data.control.p),Tsit5());
end;

# ╔═╡ c06113d0-a5d7-4bee-80a3-dcc4d4ac5202
plot_phase = 
	plot(sol,idxs=(0,chan.L),c=:yellow,label="",border=:none,ylims=(-2.5,2.5),size=(600,200));

# ╔═╡ b9cad04e-ec07-40f6-b659-335bcb058ade
theme(:dark)

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    margin: 0 auto;
    max-width: 1800px;
	padding-left: max(283px, 10%);
    padding-right: max(383px, 10%); 
}
input[type*="range"] {
	width: 38%;
}
</style>
"""

# ╔═╡ 1795b892-ecb4-45c7-8ff8-7711fb518159
kk = [2.,1.5,1.,0.7,0.5,0.4,0.3,0.2,0.1]

# ╔═╡ 6140b3b6-380a-46a5-9141-86115f89cbab
ss = [0.24,0.21,0.18,0.15,0.13,0.11,0.1,0.08,0.05]

# ╔═╡ 42560b1c-8ba7-4ac8-954f-81dd58eeded1
begin
	scatter(kk,ss)
	plot!(kk,0.18*kk.^(1/2))
end	

# ╔═╡ 1c211ca9-3c35-4bc3-936b-26a366322bd9
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
par_widget = @bind par PlutoUI.combine() do Child
	md"""
	# New Reed
	μ : $(Child("μ", Slider(-0.01:0.005:3.0;show_value=true))) $sp
	k : $(Child("k", Slider(0.1:0.002:2.0;show_value=true))) \
	σ : $(Child("σ", Slider(0.0:0.01:1.0;show_value=true))) 
	"""
end;

# ╔═╡ 5d6537bb-b55e-4e88-b69a-28c2edecbffb
par.k

# ╔═╡ b84f016d-8eb8-4084-99d9-35727647594a
begin
	set_param!(source,1,par.μ)
	set_param!(source,2,par.k)
	set_param!(source,3,par.σ)
end;

# ╔═╡ 0b5f462e-861b-4e5f-8dd0-165a4025e3cc
scale_widget = @bind scale PlutoUI.combine() do Child
	md"""
	g : $(Child("g", Slider(0.0:0.1:1.0,default=0.1;show_value=true))) 	$sp
	ts : $(Child("ts", Slider(100.0:10.0:3000.0,default=1500.0;show_value=true))) 
	"""
end;

# ╔═╡ b07c4785-7068-477a-8b12-81bb9a0a3d95
PlutoUI.ExperimentalLayout.vbox([
	par_widget,
	scale_widget,
	buttons,
	plot_phase,
	plot_widget
])

# ╔═╡ ad7aadfb-d596-4bde-aa10-64cf1e0fb688
begin 
	set_ts!(source,scale.ts)
	set_gain!(source,scale.g)
end;

# ╔═╡ Cell order:
# ╠═2e014313-3cdc-43f8-8762-e972a8bd8c28
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╟─71f589fb-2164-4fce-b000-a82970a90248
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─b07c4785-7068-477a-8b12-81bb9a0a3d95
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═9512cc09-48cd-40c8-99d1-9299c96fcf8e
# ╠═1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╠═0b5f462e-861b-4e5f-8dd0-165a4025e3cc
# ╠═41dabe99-2401-4c34-9341-f8f1c1627c64
# ╠═3e3b7c18-335c-48a3-afde-039590c50095
# ╠═248f0ef5-b462-427a-8a5c-09991da405e6
# ╠═c5916ed9-fcc5-4aec-8369-9dd74a38f966
# ╠═5d6537bb-b55e-4e88-b69a-28c2edecbffb
# ╠═b84f016d-8eb8-4084-99d9-35727647594a
# ╠═ad7aadfb-d596-4bde-aa10-64cf1e0fb688
# ╠═c000bc1a-3782-4924-b3f6-0f403645d95d
# ╠═572007b9-8d58-4c4b-bda7-24cd734674f4
# ╠═c06113d0-a5d7-4bee-80a3-dcc4d4ac5202
# ╠═b9cad04e-ec07-40f6-b659-335bcb058ade
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╠═1795b892-ecb4-45c7-8ff8-7711fb518159
# ╠═6140b3b6-380a-46a5-9141-86115f89cbab
# ╠═42560b1c-8ba7-4ac8-954f-81dd58eeded1
# ╟─1c211ca9-3c35-4bc3-936b-26a366322bd9
