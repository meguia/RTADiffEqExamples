### A Pluto.jl notebook ###
# v0.20.0

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

# ╔═╡ 2a36d245-8692-4a38-ac8b-5eef2fc5d727
 import Pkg; Pkg.activate()

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using RealTimeAudioDiffEq, PlutoUI, DifferentialEquations, JLD2, Plots

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function freed!(du,u,p,t)
    (μ,k,v0) = p
	v = (u[2]-v0)
    du[1] = u[2] 
	du[2] = v*(μ-v*v)-k*u[1]
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
begin 
	source1 = DESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.]; channel_map = [0,1])
	source2 = DESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.]; channel_map = [1,1])
	source3 = DESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.]; channel_map = [1,0])
end;

# ╔═╡ 6f1daecf-e410-49da-a9e0-1a6b77895179
output_device = get_default_output_device();

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    max-width: 1200px;
}
input[type*="range"] {
	width: 20%;
}
</style>
"""

# ╔═╡ efc7afd3-c113-4910-949f-2ced2bec0646
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ 2afa73bf-92e3-4ee7-867f-0135b5a73656
R1_widget = @bind R1 PlutoUI.combine() do Child
	md"""
	ON : $(Child("on", CheckBox())) 
	L : $(Child("L", Select([0 => "",1 => "x",2 => "y"],default=1))) 
	R : $(Child("R", Select([0 => "",1 => "x",2 => "y"],default=1))) |
	dμ : $(Child("dμ", Slider(-0.1:0.001:0.3,default=-0.1;show_value=true))) $sp
	v0 : $(Child("v0", Slider(0.0:0.01:0.5,default=0.3;show_value=true))) $sp
	k : $(Child("k", Slider(0.1:0.01:2.0,default=1.0;show_value=true))) \
	"""
end;

# ╔═╡ df2bdf00-c0a3-499c-89b9-f880bc91e66b
if R1.on
	start_DESource(source1, output_device)
else
	if isactive(source1)
		stop_DESource(source1)
	end	
end

# ╔═╡ 36a5578f-429f-47a0-a702-3d298928a8c2
R2_widget = @bind R2 PlutoUI.combine() do Child
	md"""
	ON : $(Child("on", CheckBox())) 
	L : $(Child("L", Select([0 => "",1 => "x",2 => "y"],default=1))) 
	R : $(Child("R", Select([0 => "",1 => "x",2 => "y"],default=1))) |
	dμ : $(Child("dμ", Slider(-0.1:0.001:0.3,default=-0.1;show_value=true))) $sp
	v0 : $(Child("v0", Slider(0.0:0.01:0.5,default=0.3;show_value=true))) $sp
	k : $(Child("k", Slider(0.1:0.01:2.0,default=1.0;show_value=true))) \
	"""
end;

# ╔═╡ d5efe08e-543a-42d6-8eec-2f50d73f6dee
if R2.on
	start_DESource(source2, output_device)
else	
	if isactive(source2)
		stop_DESource(source2)
	end	
end

# ╔═╡ 08318199-8b82-431f-84ae-502a2b4124e8
R3_widget = @bind R3 PlutoUI.combine() do Child
	md"""
	ON : $(Child("on", CheckBox())) 
	L : $(Child("L", Select([0 => "",1 => "x",2 => "y"],default=1))) 
	R : $(Child("R", Select([0 => "",1 => "x",2 => "y"],default=1))) |
	dμ : $(Child("dμ", Slider(-0.1:0.001:0.3,default=-0.1;show_value=true))) $sp
	v0 : $(Child("v0", Slider(0.0:0.01:0.5,default=0.3;show_value=true))) $sp
	k : $(Child("k", Slider(0.1:0.01:2.0,default=1.0;show_value=true))) \
	"""
end;

# ╔═╡ f3f2860a-6f55-4c96-9d21-9f86dd926268
if R3.on
	start_DESource(source3, output_device)
else	
	if isactive(source3)
		stop_DESource(source3)
	end	
end

# ╔═╡ d63e789a-8351-4410-a824-e4ed0d5b9a72
begin
	set_param!(source1,1,3.0*R1.v0^2+R1.dμ)
	set_param!(source1,2,R1.k)
	set_param!(source1,3,R1.v0)
	set_param!(source2,1,3.0*R2.v0^2+R2.dμ)
	set_param!(source2,2,R2.k)
	set_param!(source2,3,R2.v0)
	set_param!(source3,1,3.0*R3.v0^2+R3.dμ)
	set_param!(source3,2,R3.k)
	set_param!(source3,3,R3.v0)
end;

# ╔═╡ f0fdff97-4ed2-4513-a8a9-62847ddabd7f
begin
	set_channelmap!(source1,[R1.L,R1.R]);
	set_channelmap!(source2,[R2.L,R2.R]);
	set_channelmap!(source3,[R3.L,R3.R]);
end;

# ╔═╡ e71b8de3-1cbb-42a1-88ce-5ac453a7e9ca
scale_widget = @bind scale PlutoUI.combine() do Child
	md"""
	g : $(Child("g", Slider(0.0:0.1:1.0,default=0.1;show_value=true))) 	$sp
	ts : $(Child("ts", Slider(100.0:10.0:3000.0,default=1500.0;show_value=true))) 
	"""
end;

# ╔═╡ 71523504-077d-47db-af0c-d13f8b207638
PlutoUI.ExperimentalLayout.vbox([
	R1_widget,
	R2_widget,
	R3_widget,
	scale_widget
])

# ╔═╡ 602a76ac-7007-4b0d-b7d2-50416f0216ea
begin 
	set_ts!(source1,scale.ts)
	set_gain!(source1,scale.g)
	set_ts!(source2,scale.ts)
	set_gain!(source2,scale.g)
	set_ts!(source3,scale.ts)
	set_gain!(source3,scale.g)
end;

# ╔═╡ Cell order:
# ╠═2a36d245-8692-4a38-ac8b-5eef2fc5d727
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╠═71523504-077d-47db-af0c-d13f8b207638
# ╠═2afa73bf-92e3-4ee7-867f-0135b5a73656
# ╠═36a5578f-429f-47a0-a702-3d298928a8c2
# ╠═08318199-8b82-431f-84ae-502a2b4124e8
# ╠═e71b8de3-1cbb-42a1-88ce-5ac453a7e9ca
# ╠═df2bdf00-c0a3-499c-89b9-f880bc91e66b
# ╠═d5efe08e-543a-42d6-8eec-2f50d73f6dee
# ╠═f3f2860a-6f55-4c96-9d21-9f86dd926268
# ╠═6f1daecf-e410-49da-a9e0-1a6b77895179
# ╠═d63e789a-8351-4410-a824-e4ed0d5b9a72
# ╠═602a76ac-7007-4b0d-b7d2-50416f0216ea
# ╠═f0fdff97-4ed2-4513-a8a9-62847ddabd7f
# ╠═b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╟─efc7afd3-c113-4910-949f-2ced2bec0646
