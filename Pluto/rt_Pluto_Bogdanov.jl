### A Pluto.jl notebook ###
# v0.20.16

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

# ╔═╡ 8635f9ab-060e-4091-ab4a-f01d6dbb479e
 import Pkg; Pkg.activate()

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using RealTimeAudioDiffEq, PlutoUI, Plots, DifferentialEquations, PolynomialRoots, JLD2

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
function bogdanov!(du, u, p, t)
	(μ1,μ2,δ, η) = p 
	du[1] = u[2]
    du[2] = μ1+u[1]*(μ2-δ*u[2]+u[1]*(1-u[1]-u[2]))
end

# ╔═╡ aa83ec1c-9aec-4a71-a7ff-5c024d370fed
function noise!(du,u,p,t)
	du[1] = p[end]
	du[2] = 0
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = DESource(bogdanov!, noise!, [0.1, 0.],[-0.2,-0.2,0.5,0.0]; channel_map = [2,2]);

# ╔═╡ 3b5889be-879b-4526-8b96-993ef07891b9
begin
	ticks_button = @bind ticks PlutoUI.Clock(0.1);
	start_button = @bind start Button("START");
	stop_button = @bind stop Button("STOP");
end;

# ╔═╡ 8a155287-3565-4e4c-b2e9-1a8d658d6957
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

# ╔═╡ d42d2c83-a7f0-45b0-80ec-3034e86d1045
buttons = PlutoUI.ExperimentalLayout.Div([start_button,stop_button], style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

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

# ╔═╡ 66434892-ed49-490c-aad7-6670917954f5
begin
	file_temp = load("Homall.jld2")
	dval = file_temp["dvalues"]
	hom = file_temp["hom"]
end;

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
	margin: 0 auto;
    max-width: 1500px;
	padding-left: max(283px, 10%);
    padding-right: max(383px, 10%); 
}
input[type*="range"] {
	width: 30%;
}
</style>
"""

# ╔═╡ 434c0ab9-0ab4-4d2b-8ae4-91a4dd6005c0
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ 1b21621d-ddc2-42dc-945f-60f4809d7ba3
par_widget = @bind par PlutoUI.combine() do Child
	md"""
	## Bogdanov_Takens normal form with cubic terms and additive noise
	$\dot{x} = y  + \eta(t)$
	$\dot{y} = \mu_1 + \mu_2 x - \delta xy + x^2 - x^3 - x^2y$ 
	
	μ2 : $(Child("μ2", Slider(-0.33:0.001:0.4,default=-0.02;show_value=true))) $sp
	dμ2 : $(Child("dμ2", Slider(-0.02:0.001:0.02,default=0.02;show_value=true))) \
	δ : $(Child("δ", Slider(0.0:0.1:2.0,default=0.5;show_value=true))) $sp
	η : $(Child("η", Slider(0.0:0.001:0.1,default=0.0;show_value=true))) 
	"""
end;

# ╔═╡ 6799b5b5-61bf-4242-874b-3456d64dd02d
begin
	bt2 = bt(par.δ)
	μsn = -0.33:0.01:max(bt2,0.8)
	μh1 = -0.5:0.01:0
	μh2 = -0.5:0.01:0.8
	nhom = argmin(abs.(par.δ .- dval/10))
	plot_phase = plot(μsn,sn1.(μsn),c=:black,label="SN",xlabel="μ2",ylabel="μ1")
	plot!(plot_phase,μsn,sn2.(μsn),c=:black,label="SN")
	plot!(plot_phase,μh1,μh1*0,c=:red,label="Hopf 1")
	plot!(plot_phase,μh2,h2.(μh2,par.δ),c=:red,label="Hopf 2")
	plot!(plot_phase,hom[nhom][2,:],hom[nhom][1,:],c=:green,label="Homoclinic")
	scatter!(plot_phase,[0],[0],c=:red,ms=3,label="BT1")
	scatter!(plot_phase,[bt2],[sn2(bt2)],c=:red,ms=3,label="BT2")
	scatter!(plot_phase,[par.μ2-par.dμ2],[sn1(par.μ2)],c=:blue,label="",xlims=(-0.5,0.8),ylims=(-0.5,0.1),size=(600,400))
end;

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin 
	set_param!(source,1,sn1(par.μ2))
	set_param!(source,2,par.μ2-par.dμ2)
	set_param!(source,3,par.δ)
	set_param!(source,4,par.η)
end;

# ╔═╡ db1c4406-490f-44a7-9bd6-822d383f8c06
scale_widget = @bind scale PlutoUI.combine() do Child
	md"""
	g : $(Child("g", Slider(0.0:0.1:1.0,default=0.1;show_value=true))) 	$sp
	ts : $(Child("ts", Slider(100.0:10.0:3000.0,default=1500.0;show_value=true))) 
	"""
end;

# ╔═╡ 75747068-fa44-4a25-b494-408de1a162e8
PlutoUI.ExperimentalLayout.vbox([
	par_widget,
	scale_widget,
	buttons,
	plot_phase
])

# ╔═╡ 30f6fd14-bee3-469d-a9c5-0364c4e7a150
begin 
	set_ts!(source,scale.ts)
	set_gain!(source,scale.g)
end;

# ╔═╡ Cell order:
# ╠═8635f9ab-060e-4091-ab4a-f01d6dbb479e
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╟─785363b9-1d03-47ba-b2a8-a82700824a80
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═aa83ec1c-9aec-4a71-a7ff-5c024d370fed
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─75747068-fa44-4a25-b494-408de1a162e8
# ╟─8a155287-3565-4e4c-b2e9-1a8d658d6957
# ╟─5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╟─6799b5b5-61bf-4242-874b-3456d64dd02d
# ╟─1b21621d-ddc2-42dc-945f-60f4809d7ba3
# ╟─db1c4406-490f-44a7-9bd6-822d383f8c06
# ╟─3b5889be-879b-4526-8b96-993ef07891b9
# ╟─d42d2c83-a7f0-45b0-80ec-3034e86d1045
# ╟─bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╟─30f6fd14-bee3-469d-a9c5-0364c4e7a150
# ╟─8d7fb8ae-ab5e-44bc-85f0-b6a1c5e472a3
# ╟─66434892-ed49-490c-aad7-6670917954f5
# ╠═b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╟─434c0ab9-0ab4-4d2b-8ae4-91a4dd6005c0
