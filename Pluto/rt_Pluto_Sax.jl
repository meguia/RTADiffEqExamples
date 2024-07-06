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

# ╔═╡ c7744a6e-da88-45c5-b2af-50c5d85b7a7b
 import Pkg; Pkg.activate()

# ╔═╡ 9922fbbc-b68a-4ce1-a790-7c6c03c894ec
using  RealTimeAudioDiffEq, PlutoUI, Plots, DifferentialEquations

# ╔═╡ 1d42bd2f-0518-4392-8abd-b14afb0f1b59
function saxRN!(dx,x,p,t)
    (γ,ζ) = p
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
source = DESource(saxRN!, u0, [0.3, 0.6] ; channel_map = [1,2]);

# ╔═╡ bbb5fd9a-7b93-4032-9cfd-60a098c35c93
var_list = [1 => "x",2 => "v",3 => "p1",4 => "p2",5 => "p3",6 => "p4",7 => "p5",8 => "p6",9 => "p7",10 => "p8"];

# ╔═╡ 05a5baa1-4373-42c8-88be-bfed98bcd9d7
chan_widget = @bind chan PlutoUI.combine() do Child
	md"""
	L : $(Child("L", Select( var_list ,default=1)))
	R : $(Child("R", Select( var_list ,default=2)))
	"""
end;

# ╔═╡ bce16403-6dac-4b30-9327-0fd17f04d2a9
begin
	ticks_button = @bind ticks Clock(0.1);
	start_button = @bind start Button("START");
	stop_button = @bind stop Button("STOP");
	reset_button = @bind reset Button("RESET");
end;

# ╔═╡ 1678c539-0f23-4638-a9ff-461ef268ad63
let 
	start
	output_device = get_default_output_device();
	start_DESource(source, output_device; buffer_size=convert(UInt32,2048))
end	

# ╔═╡ 5fd39e58-16db-4137-863c-a50e593b6e58
let 
	stop
	stop_DESource(source)
end	

# ╔═╡ dda4f990-6911-47a0-b94a-a5fc1b660759
let
	reset
	reset_state!(source)
end;

# ╔═╡ 37232245-ec07-4b70-a47f-903c6fb2a75f
buttons = PlutoUI.ExperimentalLayout.Div([start_button,stop_button,reset_button,chan_widget], style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ 514493d0-e0a5-48b5-af95-f4f76f1887c6
plot_widget = PlutoUI.ExperimentalLayout.Div(ticks_button, style=Dict(	"display" => "flex","flex-direction" => "row","background" => "gray"));

# ╔═╡ e2e9a567-86bf-42ea-8ccd-d9b5d3929cdc
set_channelmap!(source,[chan.L,chan.R]);

# ╔═╡ e8454fdd-b322-4e67-bfa1-1c66be424837
begin
	ticks
	sol = solve(ODEProblem(saxRN!,source.data.state.u,(source.data.state.t,source.data.state.t+30.0),source.data.control.p),Tsit5());
end;

# ╔═╡ 7820d89a-1e54-4800-9ca8-5f04ac08b7de
begin 
	plot_phase = plot(sol,idxs=(0,source.data.control.channel_map[1]),c=:blue,label="");
	plot!(plot_phase,sol,idxs=(0,source.data.control.channel_map[2]),c=:red,label="",border=:none,size=(600,200));
end;

# ╔═╡ b0744443-8d19-41dc-abe8-9ba90ca91ca7
html"""
<style>
main {
    max-width: 1000px;
}
input[type*="range"] {
	width: 38%;
}
</style>
"""

# ╔═╡ ead99107-492e-447c-bb5f-8407d1909658
sp = html"&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp&nbsp";

# ╔═╡ 3c72c446-0afa-41d2-b0ec-635f7d7ab695
par_widget = @bind par PlutoUI.combine() do Child
	md"""
	# Saxophone
	γ : $(Child("γ", Slider(0.0:0.005:1.0,default=0.1;show_value=true))) $sp
	ζ : $(Child("ζ", Slider(0.0:0.01:1.0,default=0.6;show_value=true))) \
	""" 
end;

# ╔═╡ 38ed81aa-ac15-40f7-9066-efd070f84fa4
begin
	set_param!(source,1,par.γ)
	set_param!(source,2,par.ζ)
end;

# ╔═╡ fa2fc75a-41b8-4197-92f6-71febfa64f30
scale_widget = @bind scale PlutoUI.combine() do Child
	md"""
	g : $(Child("g", Slider(0.0:0.01:0.2,default=0.1;show_value=true))) 	$sp
	ts : $(Child("ts", Slider(100.0:10.0:3000.0,default=1600.0;show_value=true))) 
	"""
end;

# ╔═╡ fdbfbbc9-9887-40b3-acdb-08fcc9f8d98e
begin 
	dash_cell =  PlutoRunner.currently_running_cell_id[] |> string
	PlutoUI.ExperimentalLayout.vbox([
		par_widget,
		buttons,
		plot_phase,
		plot_widget,
		scale_widget
	])
end

# ╔═╡ 3d419dab-b634-480d-93eb-3fb92c623a0a
begin 
	set_ts!(source,scale.ts)
	set_gain!(source,scale.g)
end;

# ╔═╡ ba3f15e0-6321-42db-aede-ccc20c0a73b0
notebook = PlutoRunner.notebook_id[] |> string

# ╔═╡ 3eda34b2-ddf6-4466-bd1f-539270a0b77e
dash_url = "http://localhost:1234/edit?id=" * notebook * "&isolated_cell_id=" * dash_cell * "&"

# ╔═╡ 21c23174-807d-4367-8d6a-b40a25c2c9d4
Markdown.parse("[Open Saxophone Dashboard]($dash_url)")

# ╔═╡ Cell order:
# ╠═c7744a6e-da88-45c5-b2af-50c5d85b7a7b
# ╠═9922fbbc-b68a-4ce1-a790-7c6c03c894ec
# ╠═1d42bd2f-0518-4392-8abd-b14afb0f1b59
# ╠═08508733-0b35-4f45-a0ae-59a787637fe8
# ╠═46a395c6-2cd0-42c8-b4e5-d0d0f71638e3
# ╟─fdbfbbc9-9887-40b3-acdb-08fcc9f8d98e
# ╟─21c23174-807d-4367-8d6a-b40a25c2c9d4
# ╟─1678c539-0f23-4638-a9ff-461ef268ad63
# ╟─5fd39e58-16db-4137-863c-a50e593b6e58
# ╟─dda4f990-6911-47a0-b94a-a5fc1b660759
# ╟─3c72c446-0afa-41d2-b0ec-635f7d7ab695
# ╟─fa2fc75a-41b8-4197-92f6-71febfa64f30
# ╟─bbb5fd9a-7b93-4032-9cfd-60a098c35c93
# ╟─05a5baa1-4373-42c8-88be-bfed98bcd9d7
# ╟─bce16403-6dac-4b30-9327-0fd17f04d2a9
# ╟─37232245-ec07-4b70-a47f-903c6fb2a75f
# ╟─514493d0-e0a5-48b5-af95-f4f76f1887c6
# ╟─38ed81aa-ac15-40f7-9066-efd070f84fa4
# ╟─3d419dab-b634-480d-93eb-3fb92c623a0a
# ╟─e2e9a567-86bf-42ea-8ccd-d9b5d3929cdc
# ╠═e8454fdd-b322-4e67-bfa1-1c66be424837
# ╠═7820d89a-1e54-4800-9ca8-5f04ac08b7de
# ╟─b0744443-8d19-41dc-abe8-9ba90ca91ca7
# ╟─ead99107-492e-447c-bb5f-8407d1909658
# ╟─ba3f15e0-6321-42db-aede-ccc20c0a73b0
# ╟─3eda34b2-ddf6-4466-bd1f-539270a0b77e
