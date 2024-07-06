### A Pluto.jl notebook ###
# v0.19.42

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

# ╔═╡ 056742e9-fe30-4bac-b52d-b29c0456f6d5
 import Pkg; Pkg.activate()

# ╔═╡ faf6f058-ee06-4a6c-8dd0-970620229f2f
Pkg.add("PlutoPlotly")

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using RealTimeAudioDiffEq, PlutoUI, DifferentialEquations, PlutoPlotly

# ╔═╡ e69993f0-56aa-49bf-b91b-44d39989b5ff
md"""
# Double Hopf from two Coupled Bautin Oscillators

Each oscillator is controlled by four parameters:

$\dot{x} = y$

$\dot{y} = -kx + \mu y + \sigma x^2y  + \xi y (x^2+y^2)^2$

- parameter $k$ controls the pitch ($k=\omega^2$ at the bifurcation) and also the shape of the limit cycle (higher frequencies give rise to smoother cycles, lower frequencies give rise to relaxation oscillators with sharp peaks)

- parameter $\mu$ controls the onset of the oscillation and the amplitude (controls the Hopf bifurcation). The bifurcation occurs at $\mu=0$ (and for the supercritical case $\sigma<0$ this is the onset of the oscillation).

- parameter $\sigma$ controls the Hopf type. For $\sigma<0$ the bifurcation is supercritical and the stable equilibrium gives rise to a stable limit cycle as $\mu$ crosses zero from below. For $\sigma>0$ the bifurcation is subcritical and the stable equilibrium gives rise to an unstable limit cycle at $\mu=0$ that go backwards in parameter space.

- parameter $\xi$ must be always negative and provides the higher nonlinear saturation for the subcritical case and in that case it creates another stable limit cycle that anihilates the unstable limit cycle in a saddle-node bifurcation of periodics orbits for $\mu<0$. 

## Coupling

Taking the normal form of the Double Hopf in its "simplest" version we couple the oscillators through the cubic term:

$\dot{x_1} = y_1$
$\dot{y_1} = -k_1x_1 + \mu_1 y_1 + \sigma y_1x_1^2  + c_{21}y_1x_2^2 + \xi y_1 (x_1^2+y_1^2)^2$
$\dot{x_2} = y_2$
$\dot{y_2} = -k_2x_2 + \mu_2 y_2 + \sigma y_2x_2^2  + c_{12}y_2x_1^2 + \xi y_2 (x_2^2+y_2^2)^2$

Parameters $k$ and $\mu$ are different for each oscillator, but $\sigma$ and $\xi$ are in common (we fix $\xi=-0.1$). The coupling parameters $c_{12}$ ans $c_{21}$ are kept positive in order to have a fixed bifurcation diagram (two Neimark Sacker curves in the positive cuadrant of $\mu_1$ $\mu_2$)

"""

# ╔═╡ 87db0fcd-6fe3-4620-be1f-6634e51723cb
function dhopf(du,u,p,t)
	(μ1,μ2,k1,k2,σ,c12,c21,ξ) = p
	du[1] = u[2]
	du[2] = -k1*u[1] + u[2]*(μ1 + σ*u[1]^2 + c21*u[3]^2 +  + ξ*(u[1]^2+u[2]^2)^2)
	du[3] = u[4]
	du[4] = -k2*u[3] + u[4]*(μ2 + σ*u[3]^2 + c12*u[1]^2 +  + ξ*(u[3]^2+u[4]^2)^2)
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = DESource(dhopf, [0.1,0.0,0.1,0.0], [-0.1,-0.1,1.0,1.0,-0.1,0.0,0.1,-0.1]; channel_map = [1,3]);

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
begin
	ticks_button = @bind ticks Clock(0.1);
	start_button = @bind start Button("START");
	stop_button = @bind stop Button("STOP");
	reset_button = @bind reset Button("RESET");
end;

# ╔═╡ b116436d-784f-4f58-aae8-37299c4acee7
PlutoUI.ExperimentalLayout.Div([start_button,stop_button,reset_button,ticks_button], style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"))

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

# ╔═╡ 94521225-d7e4-40c3-b674-9a7097ad83b6
let
	reset
	reset_state!(source)
end;

# ╔═╡ 7cc1dc79-235f-40f8-9399-197c8d1625f0
theme(:dark)

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 38%;
}
</style>
"""

# ╔═╡ da90a3c7-6d94-44a6-9edc-265437fa5502
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
md"""
ts $(@bind ts Slider(100.0:10.0:3000.0,default=1500.0;show_value=true)) $sp
g $(@bind g Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
μ1 : $(@bind μ1 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) $sp
μ2 : $(@bind μ2 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) \
kc : $(@bind kc Slider(0.01:0.01:2.0,default=1.0;show_value=true)) $sp
c : $(@bind c Slider(0.0:0.01:1.0,default=0.1;show_value=true)) \
σ $(@bind σ Select([1 => "positive", -1 => "negative"]))  $sp $sp $sp
Coupling $(@bind cs Select([[1,1] => "coupling ++", [-1,-1] => "coupling --", [1,-1] => "coupling +-"]))  $sp $sp
"""

# ╔═╡ 6a52b54b-8d04-4c74-b06e-1b940de77173
begin
	m=-0.4:0.01:0.4
	m1=-0.4:0.01:0
	ns1 = -σ*0.01*m1/c+1/16.2*(m1/c).^2
	ns2 = -σ*0.01*m1/c+1/18.5*(m1/c).^2
	plot(m,m*0,c=:black,label="H1")
	plot(m*0,m,c=:gray,label="H2")
	#scatter([μ2],[μ1])
	plot(ns1,m1,c=:red,label="NS1")
	plot(m1,ns2,c=:blue,label="NS2",size=(600,400),xlims=(-0.4,0.4),ylims=(-0.4,0.4))
end	

# ╔═╡ f30f39a2-98f4-4888-9c9b-f2c557e49f31
begin
	ticks
	sol = solve(ODEProblem(dhopf,source.data.state.u,(source.data.state.t,source.data.state.t+0.2*ts),source.data.control.p));
end;

# ╔═╡ f0943d24-1701-47b1-809c-447d20161b06
begin
	plot(sol,idxs=(1,2),c=:yellow,label="")
	plot!(sol,idxs=(3,4),c=:green,label="",border=:none)
end	

# ╔═╡ dd0a3af1-2827-498b-89d6-346db1b498af
begin
	set_param!(source,1,μ1)
	set_param!(source,2,μ2)
end;

# ╔═╡ 666b788e-ff77-4ed3-9f72-a84ff84374e3
begin
	set_param!(source,3,0.1)
	set_param!(source,4,0.1*kc)
end;

# ╔═╡ 302fde26-1ee6-432c-a9d2-fb34181af5dd
begin
	set_param!(source,5,σ*0.2)
	set_param!(source,6,cs[1]*c)
	set_param!(source,7,cs[2]*c)
end;

# ╔═╡ d1665ea7-b382-451e-bcfa-e0eb11b11c11
begin 
	set_ts!(source,ts)
	set_gain!(source,g)
end;

# ╔═╡ Cell order:
# ╠═056742e9-fe30-4bac-b52d-b29c0456f6d5
# ╠═faf6f058-ee06-4a6c-8dd0-970620229f2f
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╟─e69993f0-56aa-49bf-b91b-44d39989b5ff
# ╠═87db0fcd-6fe3-4620-be1f-6634e51723cb
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╠═6a52b54b-8d04-4c74-b06e-1b940de77173
# ╟─b116436d-784f-4f58-aae8-37299c4acee7
# ╟─f0943d24-1701-47b1-809c-447d20161b06
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═1678c539-0f23-4638-a9ff-461ef268ad63
# ╠═f30f39a2-98f4-4888-9c9b-f2c557e49f31
# ╠═94521225-d7e4-40c3-b674-9a7097ad83b6
# ╠═dd0a3af1-2827-498b-89d6-346db1b498af
# ╠═666b788e-ff77-4ed3-9f72-a84ff84374e3
# ╠═302fde26-1ee6-432c-a9d2-fb34181af5dd
# ╠═d1665ea7-b382-451e-bcfa-e0eb11b11c11
# ╟─7cc1dc79-235f-40f8-9399-197c8d1625f0
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╟─da90a3c7-6d94-44a6-9edc-265437fa5502
