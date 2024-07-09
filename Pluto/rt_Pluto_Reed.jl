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

# ╔═╡ b6b4f91c-24e8-43f6-a9c4-d636c445ea65
 import Pkg; Pkg.activate()

# ╔═╡ 5833dc44-fa81-4197-bd29-0b10d182d34c
Pkg.add("JLD2")

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using PortAudio.LibPortAudio, PlutoUI, Plots, DifferentialEquations, JLD2

# ╔═╡ df188812-0e85-4b96-b172-b2bf010d475b
include("C:/Users/Camilo/RealTimeAudioDiffEq.jl/src/rt_DE.jl")

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
    (μ,k,v0) = p
	v = (u[2]-v0)
    du[1] = u[2] 
	du[2] = v*(μ-v*v)-k*u[1]
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = DESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.]; channel_map = [1,2])

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
md"""
dμ : $(@bind dμ Slider(-0.1:0.001:0.8,default=-0.1;show_value=true)) \
v0 : $(@bind v0 Slider(0.0:0.01:0.5,default=0.3;show_value=true)) \
k : $(@bind k Slider(0.1:0.01:1.0,default=1.0;show_value=true)) \
ts $(@bind ts Slider(100.0:10.0:3000.0,default=1600.0;show_value=true)) \
g $(@bind g Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
"""

# ╔═╡ 3e3b7c18-335c-48a3-afde-039590c50095
begin
	start_button = @bind start Button("START");
	stop_button = @bind stop Button("STOP");
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

# ╔═╡ 248f0ef5-b462-427a-8a5c-09991da405e6
buttons = PlutoUI.ExperimentalLayout.Div([start_button,stop_button], style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ 09d28e82-ac78-4b7b-9ed5-2fa8f8ceaeef
buttons

# ╔═╡ 225ba501-3a46-443a-91ea-c22dc7cbdad5
begin
	saved_values = load("freed.jld2")
	mu_v = saved_values["mu"]
	v0_v = saved_values["v0"]
	knotes = saved_values["knotes"]
	ampl_v = saved_values["ampl"]
	freq_v = saved_values["freq"]
end;	

# ╔═╡ 68c7bb4f-5345-4f15-a974-bdecdc7c0209
begin
	ik = findmin(abs.(knotes .- k))[2]
	contourf(v0_v,mu_v,freq_v[ik,:,:]',levels=levels=freq_v[ik,1,1]*2 .^(-12/12:1/12:0),colorbar=false)
	plot!(v0_v,3*v0_v.^2,c=:black,lw=3,label="")
	scatter!([v0],[3*v0^2+dμ],c=:blue,label="",xlabel="V0",ylabel="μ")
end	

# ╔═╡ b84f016d-8eb8-4084-99d9-35727647594a
begin
	set_param!(source,1,3*v0^2+dμ)
	set_param!(source,2,k)
	set_param!(source,3,v0)
end;

# ╔═╡ ad7aadfb-d596-4bde-aa10-64cf1e0fb688
begin 
	set_ts!(source,ts)
	set_gain!(source,g)
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
# ╠═5833dc44-fa81-4197-bd29-0b10d182d34c
# ╠═b6b4f91c-24e8-43f6-a9c4-d636c445ea65
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═df188812-0e85-4b96-b172-b2bf010d475b
# ╠═1b0e9cf3-1dbb-4ae6-b6d4-3e533748e7a7
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─68c7bb4f-5345-4f15-a974-bdecdc7c0209
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─09d28e82-ac78-4b7b-9ed5-2fa8f8ceaeef
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═3e3b7c18-335c-48a3-afde-039590c50095
# ╠═248f0ef5-b462-427a-8a5c-09991da405e6
# ╟─225ba501-3a46-443a-91ea-c22dc7cbdad5
# ╠═b84f016d-8eb8-4084-99d9-35727647594a
# ╠═ad7aadfb-d596-4bde-aa10-64cf1e0fb688
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
