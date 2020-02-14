# Loading main functions
include("funcs2D.jl")

# setting the computation mode (1 modeling, 2 illumination)
comp = parse(Int32, ARGS[2])

# Loading models
if ARGS[1] == "1"
    @printf( "\nLoading simple model.\n")
    data = load_data("synthmodel2D", "..")
else
    @printf( "\nLoading SEG/EAGE model.\n")
    data = load_data("SEGsaltmodel", "..")
end

# model dimensions (nz is depth and nx is the horizontal axis)
(nz, nx) = size(data.V)
@printf( "Model size: z %d, x %d\n\n", nz, nx )

# peak source frequency
freq = 15.f0

# Setting number of time steps 
if ARGS[1] == "1"
    # simple model
    nt = 750
    # ns number of shots, z location then x location
    ns = 1
else
    # EAGE/SEG
    nt = 1500
    # ns number of shots, z location then x location
    ns = 1
end


SourceVec = zeros(Int32, (ns, 2))
# Setting shots (aka sources)
if ARGS[1] == "1"
    # sources for simple model
    SourceVec = [ 
		#=Int32(5) Int32(40)
		Int32(5) Int32(70)=#
		Int32(5) Int32(100)
		#Int32(5) Int32(140)
		#=Int32(5) Int32(160)=#
             	] 
else
    # sources for SEG/EAGE model
    SourceVec = [ 
		#=Int32(5) Int32(60) 
		Int32(5) Int32(120) 
		Int32(5) Int32(180) 
		Int32(5) Int32(240) 
		Int32(5) Int32(300) =#
		Int32(5) Int32(360) 
		#=Int32(5) Int32(400) 
		Int32(5) Int32(460) 
		Int32(5) Int32(520)
		Int32(5) Int32(580)  
		Int32(5) Int32(640)=#  
		] 
end

if comp == 1
    # Computing propagation (be careful with your RAM!)
    waves = zeros(Float32, (nx, nz, nt, ns))
    traces = zeros(Float32, (nt, nx, ns))
    (waves, traces ) = modeling( data, SourceVec, freq, nt );

    for i = 1:ns
	savefile( traces[:,:,i], "gather_shot$i.bin" )
	savefile( waves[:,:,:,i], "wavefields_shot$i.bin" )
    end

    @printf( "\nDone modeling.\n")
end

if comp == 2
    # Computing illumination
    illmap = zeros(Float32, (nx, nz))
    illmap = illuminating( data, SourceVec, freq, nt );

    savefile( illmap, "source_ill_shots.bin" )

    @printf( "\nDone illuminating.\n")
end
