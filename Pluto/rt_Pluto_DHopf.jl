### A Pluto.jl notebook ###
# v0.19.41

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

# ╔═╡ d1af62f9-4313-42cd-a87d-3557f55b0511
using Pkg; Pkg.activate("../rtODE")

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using PortAudio, PortAudio.LibPortAudio, PlutoUI, DifferentialEquations, Plots

# ╔═╡ fa35c482-d74f-11ee-0e9f-77b332036253
include("../rtODE/rt_SDE.jl")

# ╔═╡ 20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
begin
	# Global Parameters
	sample_rate::Cdouble = 48000.0
	buffer_size::Culonglong = 1024
end;

# ╔═╡ 6f1daecf-e410-49da-a9e0-1a6b77895179
output_device = Pa_GetDefaultOutputDevice()

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

# ╔═╡ 305d0572-cd1f-467b-8689-d01e899c0dbe
function noise(du,u,p,t)
	du[1] = 0.001
	du[2] = 0
	du[3] = 0.001
	du[4] = 0
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = rt_SDESource(dhopf, noise, [0.1,0.0,0.1,0.0], [-0.1,-0.1,1.0,1.0,-0.1,0.0,0.1,-0.1], sample_rate, buffer_size,[1,3]);

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
@bind start Button("START")

# ╔═╡ 8a155287-3565-4e4c-b2e9-1a8d658d6957
@bind stop Button("STOP")

# ╔═╡ 2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
let 
	start
	rt_SDEStart(source, output_device)
end	

# ╔═╡ 5b8f7326-6d7f-44ac-82b9-799f03cedf46
let 
	stop
	rt_SDEStop(source)
end	

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
ts $(@bind ts Slider(100:10:3000,default=1600;show_value=true)) $sp
g $(@bind g Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
μ1 : $(@bind μ1 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) $sp
μ2 : $(@bind μ2 Slider(-1.0:0.01:1.0,default=0.1;show_value=true)) \
k : $(@bind k Slider(0.01:0.01:2.0,default=1.0;show_value=true)) $sp
c : $(@bind c Slider(0.0:0.01:1.0,default=0.1;show_value=true)) \
σ $(@bind σ Select([1 => "positive", -1 => "negative"]))  $sp $sp $sp
Coupling $(@bind cs Select([[1,1] => "coupling ++", [-1,-1] => "coupling --", [1,-1] => "coupling +-"]))  
"""

# ╔═╡ 6a52b54b-8d04-4c74-b06e-1b940de77173
begin
	m=-0.4:0.01:0.4
	m1=-0.4:0.01:0
	ns1 = -σ*0.01*m1/c+1/16.2*(m1/c).^2
	ns2 = -σ*0.01*m1/c+1/18.5*(m1/c).^2
	plot(m,m*0,c=:black,label="H1")
	plot!(m*0,m,c=:gray,label="H2")
	scatter!([μ2],[μ1])
	plot!(ns1,m1,c=:red,label="NS1")
	plot!(m1,ns2,c=:blue,label="NS2",size=(600,400),xlims=(-0.4,0.4),ylims=(-0.4,0.4))
end	

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin 
	k1 = k*0.1
	k2 = 0.1
	σ2 = σ*0.2
	c12 = cs[2]*c
	c21 = cs[1]*c
	@atomic source.data.control.ts = ts
	@atomic source.data.control.p = [μ1,μ2,k1,k2,σ2,c12,c21,-0.1]
	@atomic source.data.control.gain = g
end;

# ╔═╡ Cell order:
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═d1af62f9-4313-42cd-a87d-3557f55b0511
# ╠═fa35c482-d74f-11ee-0e9f-77b332036253
# ╟─20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
# ╠═6f1daecf-e410-49da-a9e0-1a6b77895179
# ╟─e69993f0-56aa-49bf-b91b-44d39989b5ff
# ╠═87db0fcd-6fe3-4620-be1f-6634e51723cb
# ╠═305d0572-cd1f-467b-8689-d01e899c0dbe
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─6a52b54b-8d04-4c74-b06e-1b940de77173
# ╟─1678c539-0f23-4638-a9ff-461ef268ad63
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─8a155287-3565-4e4c-b2e9-1a8d658d6957
# ╟─bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╟─da90a3c7-6d94-44a6-9edc-265437fa5502
