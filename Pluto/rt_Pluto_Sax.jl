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

# ╔═╡ 00af18bb-7300-432f-9d40-45ed232f0d29
using Pkg; Pkg.activate("../PortAudioODE/rtODE")

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using PortAudio, PortAudio.LibPortAudio, PlutoUI, DifferentialEquations

# ╔═╡ fa35c482-d74f-11ee-0e9f-77b332036253
include("../../PortAudioODE/rtODE/rt_ODE.jl")

# ╔═╡ 20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
begin
	# Global Parameters
	sample_rate::Cdouble = 48000.0
	buffer_size::Culonglong = 512
end;

# ╔═╡ 6f1daecf-e410-49da-a9e0-1a6b77895179
output_device = Pa_GetDefaultOutputDevice()

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function saxRN!(dx,x,p,t)
    γ,ζ = p
	ω0 = 4.224
	ω = [1.195,2.483,3.727,4.405,5.153,6.177,6.749,7.987]
	α = [0.0176,0.0355,0.0653,0.2693,0.0703,0.166,0.0945,0.1165]
	C = [0.1761,0.4705,0.6494,0.328,0.541,0.2249,0.3822,0.4099]
    P = sum(x[3:2:end])
    Fc = 100.0*min(real(x[1])+1,0)^2*(1-x[2])
    u = ζ*max(real(x[1])+1,0)*sign(γ-P)*sqrt(abs(γ-P))
    dx[1] = x[2]
    dx[2] = -ω0*x[2]+ω0^2*(P-γ-x[1]+Fc)
    for n=3:2:length(x)
        m = (n-1)÷2
        dx[n] = -α[m]*x[n]-2*ω[m]*x[n+1]+2*C[m]*u
        dx[n+1] = -α[m]*x[n+1]+0.5*ω[m]*x[n]
    end
    dx
end

# ╔═╡ 08508733-0b35-4f45-a0ae-59a787637fe8
u0 = vcat([-0.3, 0],repeat([0.0,0.001],8));

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = rt_ODESource(saxRN!, u0, [0.3, 0.6], sample_rate, buffer_size)

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
@bind start Button("START")

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
md"""
ts $(@bind ts Slider(100:10:3000,default=1600;show_value=true)) \
g $(@bind g Slider(0.0:0.1:1.0,default=0.1;show_value=true)) \
γ : $(@bind γ Slider(0.0:0.005:1.0,default=0.1;show_value=true))\
ζ : $(@bind ζ Slider(0.0:0.01:1.0,default=0.6;show_value=true)) \
"""

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin
	@atomic source.data.control.ts = ts
	@atomic source.data.control.p = [γ, ζ]
	@atomic source.data.control.gain = g
end	

# ╔═╡ 8a155287-3565-4e4c-b2e9-1a8d658d6957
@bind stop Button("STOP")

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
# ╠═fa35c482-d74f-11ee-0e9f-77b332036253
# ╟─20f8ad8d-eb47-49ca-b7f2-262dbc8cc707
# ╠═6f1daecf-e410-49da-a9e0-1a6b77895179
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═08508733-0b35-4f45-a0ae-59a787637fe8
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─1678c539-0f23-4638-a9ff-461ef268ad63
# ╠═bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─8a155287-3565-4e4c-b2e9-1a8d658d6957
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
