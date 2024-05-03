### A Pluto.jl notebook ###
# v0.19.40

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

# ╔═╡ 00af18bb-7300-432f-9d40-45ed232f0d29
using Pkg; Pkg.activate("../rtODE")

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using PortAudio, PortAudio.LibPortAudio, PlutoUI, DifferentialEquations, Plots, PolynomialRoots, JLD2

# ╔═╡ fa35c482-d74f-11ee-0e9f-77b332036253
include("../rtODE/rt_SDE.jl")

# ╔═╡ 6f1daecf-e410-49da-a9e0-1a6b77895179
output_device = Pa_GetDefaultOutputDevice()

# ╔═╡ 20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
begin
	# Global Parameters
	sample_rate::Cdouble = 48000.0
	buffer_size::Culonglong = 512
end;

# ╔═╡ 785363b9-1d03-47ba-b2a8-a82700824a80
md"""
# Bogdanov Takens Normal Form with Cubic Terms

This is a posible universal unfolding of the Bogdanov-Takens (BT) codimension 2 bifurcation (that in fact has 2 BT points) with a third parameter controlling the "hardness" of the oscillator and the senistivity to noise.

We use this system to explore the neighborhoods of a Saddle-Node on a Limit Cycle bifurcation (SNLC) which is a second route to oscillatory behavior in 2D system (if we consider the Hopf bifurcation as the first route). In this route, the oscillations appear with infinite period, and before the bifurcation the system behaves as **excitable**, that means that a small perturbation (like noise added to a variable) can trigger an amplified response. Also this SN branch (that in turn originates in a cusp bifurcation) can limit some region of **bistability** between a stationary state and a stable limit cycle. The border between these regions of excitable and bistable behavior is an Homoclinic Bifurcation curve. 

Since the region of interest is narrow zone around that SN bifurcation we will move along that curve. Then we will have three control parameters for the system and one for the noise:

- μ2 that controls the position along the SN curve going from the cusp point at μ2=-1/3 and crosses two Hopf bifurcations (the main region of interest is between these two Hopf crossings). Since we are restrictig to this curve the other bifurcation parameter of the system (μ1) is not free.

- δμ2 that controls the distance from the SN curve and acts as a correction of the previous μ2 control parameter. For positive values we are in the excitable (or bistable) region and for negative values we enter into the periodic regime

- δ is a third parameter that acts like a nonliean dissipation and control the "hardness" of the oscillations. For lower values there are many oscillations around the focus point and the system y more senitive to noise. For higher values the oscillations are more like peaks and is less sensitive to noise. Also this parameter controls the second Hopf curve and the Homoclinic curve. Then, for lower values the region that limits the SN is bistable and that region get narrower as δ rises and finally it gets confined only to the higher part of the curve (near the first Hopf curve) and all the rest correspond to the excitable region.

- η controls the amount of white noise added to the first equation of the system. Whithin the excitable or bistable regions a small amount of noise can trigger pulses or bursts respectively.

The deterministic part of the system is:

$\dot{x} = y$
$\dot{y} = \mu_1 + \mu_2 x - \delta xy + x^2 - x^3 - x^2y$


"""

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function bogdanov(du, u, p, t)
	(μ1,μ2,δ, η) = p 
	du[1] = u[2]
    du[2] = μ1+u[1]*(μ2-δ*u[2]+u[1]*(1-u[1]-u[2]))
end

# ╔═╡ aa83ec1c-9aec-4a71-a7ff-5c024d370fed
function noise(du,u,p,t)
	du[1] = p[end]
	du[2] = 0
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = rt_SDESource(bogdanov, noise, [0.1, 0.],[-0.2,-0.2,0.5,0.0], sample_rate, buffer_size);

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
@bind start Button("START")

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
md"""
ts $(@bind ts Slider(100:10:3000,default=1600;show_value=true)) \
g $(@bind g Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
dμ2 : $(@bind dμ2 Slider(-0.02:0.001:0.02,default=0.02;show_value=true))\
μ20 : $(@bind μ20 Slider(-0.33:0.001:0.4,default=-0.02;show_value=true)) \
δ : $(@bind δ Slider(0.0:0.01:2.0,default=0.5;show_value=true)) \
η : $(@bind η Slider(0.0:0.001:0.1,default=0.0;show_value=true)) \
varout : $(@bind varout Select([1,2]))
"""

# ╔═╡ 26b6d16f-d884-467b-9bf8-b2b7a80f7810


# ╔═╡ 8a155287-3565-4e4c-b2e9-1a8d658d6957
@bind stop Button("STOP")

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin 
	@atomic source.data.control.ts = ts
	@atomic source.data.control.p = [sn1(μ20),μ20-dμ2,δ,η]
	@atomic source.data.control.gain = g
	@atomic source.data.control.varout = varout
end;

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

# ╔═╡ 8d7fb8ae-ab5e-44bc-85f0-b6a1c5e472a3
begin
	sn1(μ2) =	-(3*μ2+2/3+sqrt(3*μ2+1.0)*(2*μ2+2/3))/9.0
	sn2(μ2) =	-(3*μ2+2/3-sqrt(3*μ2+1.0)*(2*μ2+2/3))/9.0
	h2(μ2,δ) = δ*(μ2-δ-δ^2)
	function bt(δ) 
		δ2 = 9*(δ+3*δ^2+3*δ^3)
		m2 = (real(roots([-(δ2+1),0,3.0+9*δ,-2.0])).^2 .-1 )./3
		m2[1]
	end	
end;	

# ╔═╡ 516e46ab-93ad-4928-8bc0-8f0da1b840d3
sn1(μ20)

# ╔═╡ 66434892-ed49-490c-aad7-6670917954f5
begin
	file_temp = load("Homall.jld2")
	dval = file_temp["dvalues"]
	hom = file_temp["hom"]
end;

# ╔═╡ 6799b5b5-61bf-4242-874b-3456d64dd02d
begin
	bt2 = bt(δ)
	μsn = -0.33:0.01:max(bt2,0.8)
	μh1 = -0.5:0.01:0
	μh2 = -0.5:0.01:0.8
	nhom = argmin(abs.(δ .- dval/10))
	plot(μsn,sn1.(μsn),c=:black,label="SN",xlabel="μ2",ylabel="μ1")
	plot!(μsn,sn2.(μsn),c=:black,label="SN")
	plot!(μh1,μh1*0,c=:red,label="Hopf 1")
	plot!(μh2,h2.(μh2,δ),c=:red,label="Hopf 2")
	plot!(hom[nhom][2,:],hom[nhom][1,:],c=:green,label="Homoclinic")
	scatter!([0],[0],c=:red,ms=3,label="BT1")
	scatter!([bt2],[sn2(bt2)],c=:red,ms=3,label="BT2")
	scatter!([μ20-dμ2],[sn1(μ20)],c=:blue,label="",xlims=(-0.5,0.8),ylims=(-0.5,0.1))
end	

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    max-width: 1200px;
}
input[type*="range"] {
	width: 90%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═00af18bb-7300-432f-9d40-45ed232f0d29
# ╠═fa35c482-d74f-11ee-0e9f-77b332036253
# ╠═6f1daecf-e410-49da-a9e0-1a6b77895179
# ╠═20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
# ╟─785363b9-1d03-47ba-b2a8-a82700824a80
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═aa83ec1c-9aec-4a71-a7ff-5c024d370fed
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─1678c539-0f23-4638-a9ff-461ef268ad63
# ╟─6799b5b5-61bf-4242-874b-3456d64dd02d
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╠═516e46ab-93ad-4928-8bc0-8f0da1b840d3
# ╠═26b6d16f-d884-467b-9bf8-b2b7a80f7810
# ╟─8a155287-3565-4e4c-b2e9-1a8d658d6957
# ╠═bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═8d7fb8ae-ab5e-44bc-85f0-b6a1c5e472a3
# ╠═66434892-ed49-490c-aad7-6670917954f5
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
