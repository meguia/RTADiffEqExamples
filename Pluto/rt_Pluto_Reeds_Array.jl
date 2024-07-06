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

# ╔═╡ 87c3d8bb-1c0b-460c-8acf-d03992de9aa5
 import Pkg; Pkg.activate()

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using RealTimeAudioDiffEq, PlutoUI, DifferentialEquations

# ╔═╡ 8aa9a79e-1442-44f0-afef-06c935a6ea05
using MarkdownLiteral: @mdx

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function freed!(du,u,p,t)
    (μ,k,v0) = p
	v = (u[2]-v0)
    du[1] = u[2] 
	du[2] = v*(μ-v*v)-k*u[1]
end

# ╔═╡ dc40a3e4-ed4c-4e04-ac6f-8bb2be28de86
Nreeds = 7

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
sources = [DESource(freed!, [0.1, 0.], [-0.1, 1.0, 0.]; channel_map = [mod(n,2),mod(n+1,2)]) for n=1:Nreeds];

# ╔═╡ 8870af1e-c3fe-4d24-a69b-f2ab210efc7f
osc_widget = @bind osc PlutoUI.combine() do Child
        @mdx(""" 
        $([[" **R$(n)**", Child("R$(n)", CheckBox())]
			 for n = 1:Nreeds]) 
        """)
	end;

# ╔═╡ 20b1d7b0-8eec-47eb-a3cc-f46b8479683b
pars_widget = @bind pars PlutoUI.combine() do Child
        @mdx("""
        $([["**R$(n)** : dμ$(n)", Child("dμ$(n)", Slider(-0.1:0.001:0.3)), 
			"V$(n)", Child("V$(n)", Slider(0.0:0.01:0.5)),
			"  k$(n)", Child("k$(n)", Slider(0.1:0.01:2.0)),
		br] for n = 1:Nreeds]) 
        """)
	end;

# ╔═╡ c0c29b3e-5f27-480e-958b-d49d97cda89d
scale_widget = @bind scale PlutoUI.combine() do Child
	md"""
	ALL : 
	g : $(Child("g", Slider(0.0:0.1:1.0,default=0.1;show_value=true))) 
	ts : $(Child("ts", Slider(100.0:10.0:3000.0,default=1500.0;show_value=true))) \
	"""
end;

# ╔═╡ 9b7e8539-1478-4378-b3a6-435fe7b88e94
PlutoUI.ExperimentalLayout.vbox([osc_widget,pars_widget,scale_widget])

# ╔═╡ e0671032-3ff4-47f7-8c52-fba441dca4ce
for n=1:Nreeds
	μ = 3*pars[3*n-1]^2+pars[3*n-2]
 	set_param!(sources[n],1,μ)
	set_param!(sources[n],2,pars[3*n])
	set_param!(sources[n],3,pars[3*n-1])
	set_ts!(sources[n],scale.ts)
	set_gain!(sources[n],scale.g)
end	

# ╔═╡ 6f1daecf-e410-49da-a9e0-1a6b77895179
output_device = get_default_output_device();

# ╔═╡ df2bdf00-c0a3-499c-89b9-f880bc91e66b
for n=1:Nreeds
	if osc[n]
		start_DESource(sources[n], output_device)
	else
		if isactive(sources[n])
			stop_DESource(sources[n])
		end	
	end
end	

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 25%;
}
</style>
"""

# ╔═╡ Cell order:
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═dc40a3e4-ed4c-4e04-ac6f-8bb2be28de86
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╠═9b7e8539-1478-4378-b3a6-435fe7b88e94
# ╠═df2bdf00-c0a3-499c-89b9-f880bc91e66b
# ╟─8870af1e-c3fe-4d24-a69b-f2ab210efc7f
# ╟─20b1d7b0-8eec-47eb-a3cc-f46b8479683b
# ╟─c0c29b3e-5f27-480e-958b-d49d97cda89d
# ╟─e0671032-3ff4-47f7-8c52-fba441dca4ce
# ╟─6f1daecf-e410-49da-a9e0-1a6b77895179
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╠═87c3d8bb-1c0b-460c-8acf-d03992de9aa5
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═8aa9a79e-1442-44f0-afef-06c935a6ea05
