### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using Atomix,PortAudio.LibPortAudio, PlutoUI, DifferentialEquations, Plots

# ╔═╡ fa35c482-d74f-11ee-0e9f-77b332036253
include("../../PortAudioODE/rtODE/rt_ODE.jl")

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function fduff!(du,u,p,t)
    du[1] = u[2]
    du[2] = -p[1]*u[2]+u[1]*(p[2]-u[1]*u[1])+p[3]*cos(p[4]*t)
end	

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = DESource(fduff!, [0.1, 0.], [0.15, 1.0, 0.,0.]; channel_map = [1,2]);

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
begin
	ticks_button = @bind ticks Clock(0.1,true);
	start_button = @bind start Button("START");
	stop_button = @bind stop Button("STOP");
end;

# ╔═╡ 2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
let 
	start
	output_device = get_default_output_device()
	start_DESource(source, output_device; buffer_size=convert(UInt32,1024))
end;

# ╔═╡ 5b8f7326-6d7f-44ac-82b9-799f03cedf46
let 
	stop
	stop_DESource(source)
end	

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
par_widget = @bind par PlutoUI.combine() do Child
	md"""
	# Forced Duffing
	μ : $(Child("μ", Slider(0.0:0.005:1.0,default=0.1;show_value=true))) 
	β : $(Child("β", Slider(-1.0:0.01:1.0,default=0.5;show_value=true))) \
	A : $(Child("A", Slider(0.0:0.1:5.0,default=0.0;show_value=true))) 
	ω : $(Child("ω", Slider(0.0:0.01:1.0,default=0.5;show_value=true))) \
	ts : $(Child("ts", Slider(100:10:3000,default=1600;show_value=true))) 
	g : $(Child("g", Slider(0.0:0.1:1.0,default=0.1;show_value=true))) \
	tail : $(Child("tail", Slider(10:10:300,default=100;show_value=true)))
	"""
end;

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin 
	@atomic source.data.control.ts = par.ts
	@atomic source.data.control.p = [par.μ, par.β, par.A, par.ω]
	@atomic source.data.control.gain = par.g
end;

# ╔═╡ 18bcdad5-90a4-43c5-a978-94d7746c2887
begin
	ticks
	sol = solve(ODEProblem(fduff!,source.data.state.u,
		(source.data.state.t,source.data.state.t+par.tail),source.data.control.p),Tsit5())
end;

# ╔═╡ 224ecd89-2aa1-41f2-8a6b-65a279681a05
plot_phase = plot(sol,idxs=(1,2),c=:yellow,label="",border=:none,size=(600,400));

# ╔═╡ a81bba1b-efc0-43b6-92f5-0d0f2ff2d17f
PlutoUI.ExperimentalLayout.vbox([par_widget,start_button,stop_button,plot_phase,ticks_button])

# ╔═╡ 8635d63b-f995-45d4-966c-40790ee2c12a
theme(:dark)

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 40%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═fa35c482-d74f-11ee-0e9f-77b332036253
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─a81bba1b-efc0-43b6-92f5-0d0f2ff2d17f
# ╟─2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╟─224ecd89-2aa1-41f2-8a6b-65a279681a05
# ╟─1678c539-0f23-4638-a9ff-461ef268ad63
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╟─18bcdad5-90a4-43c5-a978-94d7746c2887
# ╟─8635d63b-f995-45d4-966c-40790ee2c12a
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
