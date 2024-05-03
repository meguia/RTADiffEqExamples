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
using PortAudio, PortAudio.LibPortAudio, PlutoUI, DifferentialEquations, JLD2, Plots

# ╔═╡ fa35c482-d74f-11ee-0e9f-77b332036253
include("../rtODE/rt_ODE.jl")

# ╔═╡ 6f1daecf-e410-49da-a9e0-1a6b77895179
output_device = Pa_GetDefaultOutputDevice()

# ╔═╡ 1b0e9cf3-1dbb-4ae6-b6d4-3e533748e7a7
md"""
# Self - Oscilator: modified reed model

As an example of a self-oscillator, inspired by a musical instrument, that undergooes an Hopf bifurcation we choose modified version of the Rayleigh model of a clarinet reed:

$\dot{x}  =  y$
$\dot{y}  =  - kx + \mu(y-v_0)-(y-v_0)^3$

when $v_0=0$ it becomes the original Rayleigh model and undergoes a Supercritical Hopf bifurcation when $\mu$ crosses zero from the negatives. The angular frequency of the oscillator at the bifurcation is $\omega = \sqrt{k}$. As the value of $\mu$ increases the amplitude of the stable limit cycle increases proportional to $\sqrt{\mu}$ and the frequency slightly decreases.

In the extended parameter space $(v_0,\mu)$ the Hopf bifurcation follow the curve $\mu=3v_0^2$. Crossing the bifurcation at values of $v_0$ differente from zero give rise to more asymetrical limit cycles ith greater harmonic richness and which have a much steeper drop in pitch immediately after crossing the bifurcation.



"""

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function freed!(du,u,p,t)
    (μ,K,v0) = p
	v = (u[1]-v0)
    du[2] = u[1] 
	du[1] = v*(μ-v*v)-K*u[2]
end

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
@bind start Button("START")

# ╔═╡ 8a155287-3565-4e4c-b2e9-1a8d658d6957
@bind stop Button("STOP")

# ╔═╡ 225ba501-3a46-443a-91ea-c22dc7cbdad5
begin
	saved_values = load("freed.jld2")
	mu_v = saved_values["mu"]
	v0_v = saved_values["v0"]
	knotes = saved_values["knotes"]
	ampl_v = saved_values["ampl"]
	freq_v = saved_values["freq"]
end;	

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
md"""
dμ : $(@bind dμ Slider(-0.1:0.001:0.3,default=-0.1;show_value=true)) \
v0 : $(@bind v0 Slider(0.0:0.01:1.5,default=0.3;show_value=true)) \
K : $(@bind K Slider(knotes,default=1.0;show_value=true)) \
ts $(@bind ts Slider(100:10:3000,default=1600;show_value=true)) \
g $(@bind g Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
"""

# ╔═╡ 68c7bb4f-5345-4f15-a974-bdecdc7c0209
begin
	ik = findmin(abs.(knotes .- K))[2]
	#contourf(v0_v,mu_v,freq_v[ik,:,:]',levels=levels=freq_v[ik,1,1]*2 .^(-12/12:1/12:0),colorbar=false)
	plot!(v0_v,3*v0_v.^2,c=:black,lw=3,label="")
	scatter!([v0],[3*v0^2+dμ],c=:blue,label="",xlabel="V0",ylabel="μ")
end	

# ╔═╡ 20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
begin
	# Global Parameters
	sample_rate::Cdouble = 48000.0
	buffer_size::Culonglong = 512
end;

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = rt_ODESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.], sample_rate, buffer_size)

# ╔═╡ 2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
let 
	start
	rt_ODEStart(source, output_device)
end	

# ╔═╡ 5b8f7326-6d7f-44ac-82b9-799f03cedf46
let 
	stop
	rt_ODEStop(source)
end	

# ╔═╡ 20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
begin
	# Global Parameters
	sample_rate::Cdouble = 48000.0
	buffer_size::Culonglong = 512
end;

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = rt_ODESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.], sample_rate, buffer_size)

# ╔═╡ 2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
let 
	start
	rt_ODEStart(source, output_device)
end	

# ╔═╡ 5b8f7326-6d7f-44ac-82b9-799f03cedf46
let 
	stop
	rt_ODEStop(source)
end	

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin 
	@atomic source.data.control.ts = ts
	@atomic source.data.control.p = [3*v0^2+dμ, K, v0]
	@atomic source.data.control.gain = g
end;

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 90%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═00af18bb-7300-432f-9d40-45ed232f0d29
# ╟─fa35c482-d74f-11ee-0e9f-77b332036253
# ╠═6f1daecf-e410-49da-a9e0-1a6b77895179
# ╟─1b0e9cf3-1dbb-4ae6-b6d4-3e533748e7a7
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─1678c539-0f23-4638-a9ff-461ef268ad63
# ╠═68c7bb4f-5345-4f15-a974-bdecdc7c0209
# ╠═1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─8a155287-3565-4e4c-b2e9-1a8d658d6957
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═225ba501-3a46-443a-91ea-c22dc7cbdad5
# ╠═bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╠═20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
# ╠═b0744443-8d19-41dc-abe8-9ba90ca91ca7
