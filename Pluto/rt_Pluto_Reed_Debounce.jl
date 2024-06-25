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

# ╔═╡ 2e014313-3cdc-43f8-8762-e972a8bd8c28
 import Pkg; Pkg.activate()

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using RealTimeAudioDiffEq, PlutoUI, Plots, DifferentialEquations,HypertextLiteral

# ╔═╡ 71f589fb-2164-4fce-b000-a82970a90248
md"""
# Auto - tuned Reed Model

As an example of a self-oscillator, inspired by a musical instrument, that undergooes an Hopf bifurcation we choose modified version of the Rayleigh model of a clarinet reed:

$\dot{x}  =  y$
$\dot{y}  =  \mu y - y^3 - kx - \sigma k^{\frac{3}{2}} x^3$

when $\sigma=0$ it becomes the original Rayleigh model and undergoes a Supercritical Hopf bifurcation when $\mu$ crosses zero from the negatives. The angular frequency of the oscillator at the bifurcation is $\omega = \sqrt{k}$. As the value of $\mu$ increases the amplitude of the stable limit cycle increases proportional to $\sqrt{\mu}$ and the frequency decreases since a relaxation orbit take place close to the "N" nullcline $kx = \mu y - y^3$ 

When $\sigma>0$ the nullcline changes its shape and becomes "faster", then the frequency of the relaxation cycle increases. For some intermediate value of $\sigma$ this increase in frequency compensates the decrease due to to the jump from the Hopf to the relaxation cycle. 
This intermediate $\sigma$ value depends on $k$ in an approximate power law form with exponent $3/2$. 

The value $\sigma=0.18$ seems adequate for $k$ values within the range $(0.1,2)$

"""

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function nreed!(du,u,p,t)
    (μ,k) = p
    du[1] = u[2] 
	du[2] = u[2]*(μ-u[2]^2)-k*u[1]*(1.0+0.18*sqrt(k)*u[1]^2)
end

# ╔═╡ 46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
source = DESource(nreed!, [0.1, 0.], [0.1, 1.0]; channel_map = [2,2]);

# ╔═╡ 8f260354-0a78-40c0-ad7b-f3daead3a9f7
begin
	μ_array = -0.1:0.01:3.0
	k_array = 0.1:0.01:2.0
end	

# ╔═╡ b9cad04e-ec07-40f6-b659-335bcb058ade
theme(:dark)

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    margin: 0 auto;
    max-width: 1800px;
	padding-left: max(283px, 10%);
    padding-right: max(383px, 10%); 
}
input[type*="range"] {
	width: 38%;
}
</style>
"""

# ╔═╡ 1c211ca9-3c35-4bc3-936b-26a366322bd9
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ b5b950d4-6852-42a9-895d-a9ce2eb6ca08
function renderlodash(el, wait_ms::Real, method::String, options)
	@htl("""
	<span>
	
	$(el)

	<script>
		
	const span = currentScript.parentElement
	
	const el = span.firstElementChild
	
	console.log(el)

		const val = { current: el.value }
		
		const relay = _[$(method)](() => {
			span.value = val.current
			span.dispatchEvent(new CustomEvent("input"))
		}, $(wait_ms), $(options))
		
		
		;(async () => {
			for (const value of Generators.input(el)) {
				val.current = await value
				relay()
			  }
		
		})();
		
		
		el.addEventListener("input", (e) => {
			e.stopPropagation()
		})

	</script>
	</span>
	""")
end

# ╔═╡ 61d2bc84-7830-471e-93a2-18317f19440b
begin
	struct ThrottleDebounced{T}
		x::T
		wait_ms::Real
		method::String
		options
	end
	
	
	
	Base.get(t::ThrottleDebounced{T}) where T = if hasmethod(Base.get, Tuple{T})
		Base.get(t.x)
	else
		missing
	end
	
	function Base.show(io::IO, m::MIME"text/html", t::ThrottleDebounced)
		
		
		Base.show(io, m, renderlodash(t.x, t.wait_ms, t.method, t.options))
	end
	
	
	
end

# ╔═╡ 0b5f462e-861b-4e5f-8dd0-165a4025e3cc
scale_widget = @bind scale PlutoUI.combine() do Child
	md"""
	g : $(Child("g", Slider(0.0:0.1:1.0,default=0.1;show_value=true))) 	$sp
	ts : $(Child("ts", Slider(100.0:10.0:3000.0,default=1500.0;show_value=true))) 
	"""
end;

# ╔═╡ ad7aadfb-d596-4bde-aa10-64cf1e0fb688
begin 
	set_ts!(source,scale.ts)
	set_gain!(source,scale.g)
end;

# ╔═╡ 41dabe99-2401-4c34-9341-f8f1c1627c64
chan_widget = @bind chan PlutoUI.combine() do Child
	md"""
	L : $(Child("L", Select([1 => "x",2 => "y"],default=1)))
	R : $(Child("R", Select([1 => "x",2 => "y"],default=2)))
	"""
end;

# ╔═╡ c000bc1a-3782-4924-b3f6-0f403645d95d
set_channelmap!(source,[chan.L,chan.R]);

# ╔═╡ 3e3b7c18-335c-48a3-afde-039590c50095
begin
	ticks_button = @bind ticks Clock(0.1);
	start_button = @bind start Button("START");
	stop_button = @bind stop Button("STOP");
	reset_button = @bind reset Button("RESET");
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

# ╔═╡ 9512cc09-48cd-40c8-99d1-9299c96fcf8e
let
	reset
	reset_state!(source)
end;

# ╔═╡ 248f0ef5-b462-427a-8a5c-09991da405e6
buttons = PlutoUI.ExperimentalLayout.Div([start_button,stop_button,reset_button,chan_widget], style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ c5916ed9-fcc5-4aec-8369-9dd74a38f966
plot_widget = PlutoUI.ExperimentalLayout.Div(ticks_button, style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ 572007b9-8d58-4c4b-bda7-24cd734674f4
begin
	ticks
	sol = solve(ODEProblem(nreed!,source.data.state.u,(source.data.state.t,source.data.state.t+60.0),source.data.control.p),Tsit5());
end;

# ╔═╡ c06113d0-a5d7-4bee-80a3-dcc4d4ac5202
plot_phase = 
	plot(sol,idxs=(0,chan.L),c=:yellow,label="",border=:none,ylims=(-2.5,2.5),size=(600,200));

# ╔═╡ b07c4785-7068-477a-8b12-81bb9a0a3d95
PlutoUI.ExperimentalLayout.vbox([
	scale_widget,
	buttons,
	plot_phase,
	plot_widget
])

# ╔═╡ d3115ee2-72d0-42d2-ae10-880d32d0c77e
function debounced(el, wait=0; leading::Bool=false, maxwait::Union{Nothing,Real}=nothing, trailing::Bool=true)
	ThrottleDebounced(
		el,
		wait * 1000,
		"debounce",
		Dict(
			:leading=>leading,
			:maxwait=>maxwait,
			:trailing=>trailing,
		)
	)
end

# ╔═╡ d5d21944-f61a-4a71-9753-289ced41f7ef
@bind μidx debounced(Slider(-0.1:0.01:3.0),1.0)

# ╔═╡ 161cad06-1aca-46ed-89d0-a0a02b7eb123
@bind kidx debounced(Slider(0.1:0.01:2.0),1.0)

# ╔═╡ b84f016d-8eb8-4084-99d9-35727647594a
begin
	set_param!(source,1,μ_array[μidx])
	set_param!(source,2,k_array[kidx])
end;

# ╔═╡ Cell order:
# ╠═2e014313-3cdc-43f8-8762-e972a8bd8c28
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╟─71f589fb-2164-4fce-b000-a82970a90248
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─d5d21944-f61a-4a71-9753-289ced41f7ef
# ╟─161cad06-1aca-46ed-89d0-a0a02b7eb123
# ╟─b07c4785-7068-477a-8b12-81bb9a0a3d95
# ╠═2b6e2f6a-2a89-43ca-b75e-e6a28f34737d
# ╠═5b8f7326-6d7f-44ac-82b9-799f03cedf46
# ╠═9512cc09-48cd-40c8-99d1-9299c96fcf8e
# ╠═0b5f462e-861b-4e5f-8dd0-165a4025e3cc
# ╠═41dabe99-2401-4c34-9341-f8f1c1627c64
# ╠═3e3b7c18-335c-48a3-afde-039590c50095
# ╠═248f0ef5-b462-427a-8a5c-09991da405e6
# ╠═c5916ed9-fcc5-4aec-8369-9dd74a38f966
# ╠═b84f016d-8eb8-4084-99d9-35727647594a
# ╠═ad7aadfb-d596-4bde-aa10-64cf1e0fb688
# ╠═c000bc1a-3782-4924-b3f6-0f403645d95d
# ╠═572007b9-8d58-4c4b-bda7-24cd734674f4
# ╠═c06113d0-a5d7-4bee-80a3-dcc4d4ac5202
# ╠═8f260354-0a78-40c0-ad7b-f3daead3a9f7
# ╠═b9cad04e-ec07-40f6-b659-335bcb058ade
# ╠═b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╟─1c211ca9-3c35-4bc3-936b-26a366322bd9
# ╠═d3115ee2-72d0-42d2-ae10-880d32d0c77e
# ╠═61d2bc84-7830-471e-93a2-18317f19440b
# ╠═b5b950d4-6852-42a9-895d-a9ce2eb6ca08
