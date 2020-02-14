
# Loading main functions
include("funcs2D_imaging.jl")

# Loading models
if ARGS[1] == "1"
    @printf( "\nLoading simple model.")
    data = load_data("synthmodel2D", "..")
    nt = 1000
else
    @printf( "\nLoading SEG/EAGE model.")
    data = load_data("SEGsaltmodel", "..")
    nt = 2000
end

# setting the computation mode ( # 1 imaging and receiver illumination)
comp = parse(Int32, ARGS[2])

# model dimensions (nz is depth and nx is the horizontal axis)
(nz, nx) = size(data.V)
@printf( "\nModel size: z %d, x %d\n\n", nz, nx )

# peak source frequency
freq = 15.f0

# number of shots, offsets and skip steps
ns = 3 
offset = 20
skip = 100
# shot locations
(shotx) = linspace(offset, nx-offset, ns)
SourceVec = zeros(Int32, (ns, 2))
SourceVec[:,1] = Int32(5)
for i = 1:ns
    SourceVec[i,2] = round(shotx[i])
end

# Computing image
if comp == 1 
    images = zeros(Float32, (nz, nx, ns) )
    illumR = zeros(Float32, (nz, nx, ns) )
    #wss = zeros(Float32, (nz, nx, Int64(round(nt/skip)), ns) )
    #wsr = zeros(Float32, (nz, nx, Int64(round(nt/skip)), ns) )
    #(wss, wsr, images, illumR) = migrating( data, SourceVec, freq, nt );
    (images, illumR) = migrating( data, SourceVec, freq, nt, skip );
    img = zeros( Float32, (nz, nx) )
    ill = zeros( Float32, (nz, nx) )
    for i = 1:ns
	img[:,:] += images[:,:,i]
	ill[:,:] += illumR[:,:,i]
    end
    savefile( img, "image_shots.bin" )
    savefile( ill, "recv_ill_shots.bin" )

    @printf( "\nDone migrating and illuminating from the receivers.\n")
end
