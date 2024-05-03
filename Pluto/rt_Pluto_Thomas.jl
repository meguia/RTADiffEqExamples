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
using PortAudio, PortAudio.LibPortAudio, PlutoUI, DifferentialEquations

# ╔═╡ fa35c482-d74f-11ee-0e9f-77b332036253
include("../rtODE/rt_ODE.jl")

# ╔═╡ 20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
begin
	# Global Parameters
	sample_rate::Cdouble = 48000.0
	buffer_size::Culonglong = 512
end;

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function freed!(du,u,p,t)
    (μ,K,v0) = p
	v = (u[2]-v0)
    du[1] = u[2] 
	du[2] = v*(μ-v*v)-K*u[1]
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = rt_ODESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.], sample_rate, buffer_size)

# ╔═╡ 6f1daecf-e410-49da-a9e0-1a6b77895179
output_device = Pa_GetDefaultOutputDevice()

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
@bind start Button("START")

# ╔═╡ 2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
let 
	start
	rt_ODEStart(source, output_device)
end	

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
md"""
ts $(@bind ts Slider(100:10:3000,default=1600;show_value=true)) \
g $(@bind g Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
μ : $(@bind μ Slider(-0.5:0.005:6.0,default=-0.1;show_value=true))\
v0 : $(@bind v0 Slider(0.0:0.01:1.0,default=0.0;show_value=true)) \
K : $(@bind K Slider(0.1:0.1:5.0,default=2.0;show_value=true)) 
"""

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin 
	@atomic source.data.control.ts = ts
	@atomic source.data.control.p = [μ, K, v0]
	@atomic source.data.control.gain = g
end	

# ╔═╡ 8a155287-3565-4e4c-b2e9-1a8d658d6957
@bind stop Button("STOP")

# ╔═╡ 5b8f7326-6d7f-44ac-82b9-799f03cedf46
let 
	stop
	rt_ODEStop(source)
end	

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
# ╠═20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╠═6f1daecf-e410-49da-a9e0-1a6b77895179
# ╟─1678c539-0f23-4638-a9ff-461ef268ad63
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╠═bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─8a155287-3565-4e4c-b2e9-1a8d658d6957
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═b0744443-8d19-41dc-abe8-9ba90ca91ca7
