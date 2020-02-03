#=
=
=    Authors: Mauricio Araya Polo
=    Date: 09/2016 - present
=
=#

global const bcoeffs = [ -205.f0/72.f0, 8.f0/5.f0, 
                         -1.f0/5.f0, 8.f0/315.f0, 
                         -1.f0/560.1f0 ]
global const STENCIL = Int32(4)

type Struct
    V::Array{Float32, 2}
    dx::Float32
    dz::Float32
    # model dependent image filtering or blanking down to index iz
    zero2iz::Int
    rampup2iz::Int
end

function illuminating( data, SourceVec::Array{Int32,2}, freq, nt )
    local nz::Int32, nx::Int32

    (nz, nx) = size(data.V)
    local const num_shots = countnz(SourceVec[:,2]) 
    @show num_shots

    # need the params and coeffs
    coeffs = Array(Float32, (16))
    abcoeffs = Array(Float32, (3))
    params = Array(Int32, (14))
    ctop = zeros(Float32, (8, Int64(nx) ) )
    cbottom = zeros(Float32, (8, Int64(nx) ) )
    cleft = zeros(Float32, (Int64(nz), 8) )
    cright = zeros(Float32, (Int64(nz), 8) )
    preparation!( data, freq, nt, coeffs, ctop, 
                  cbottom, cleft, cright, abcoeffs, params )
    source = zeros(Float32, (nt))
    source = convert( Array{Float32, 1}, ricker( freq, coeffs[16], nt ) )

    ill = zeros(Float32, (Int64(nz), Int64(nx), num_shots) )
    local iSource::Int32
    tic()
    for iSource = 1:num_shots
        @printf( "shot %d, loc x %d, loc z %d\n", iSource, 
                 SourceVec[iSource,2], SourceVec[iSource,1] )

	ill[:,:,iSource] = prop( data, freq, nt, source, coeffs, ctop, 
                                 cbottom, cleft, cright, abcoeffs, params, 
                                 SourceVec[iSource,:], ~, 1 )

    end
    toc()
    return ill
end

function prop( data, freq, nt, source, coeffs, ctop, 
               cbottom, cleft, cright, abcoeffs, params, 
               sloc::Array{Int32, 2}, rdata, mode )
    (nz, nx) = size(data.V)
    local xmax = nx + STENCIL
    local zmax = nz + STENCIL
    if mode == 1
	waves = zeros(Float32, (nz, nx) )
    else
	waves = zeros(Float32, (nz, nx, Int64(nt)) )
	if mode == 2
	    traces = zeros(Float32, (Int64(nt), nx) )
	end
    end

    pin = zeros(Float32, (nz + Int32(2*STENCIL), nx + Int32(2*STENCIL)) )
    pout = zeros(Float32, (nz + Int32(2*STENCIL), nx + Int32(2*STENCIL)) )
    tmp = Ptr{Float32}
    top = zeros(Float32, (4, Int64(nx) ) )
    bottom = zeros(Float32, (4, Int64(nx) ) )
    left = zeros(Float32, (Int64(nz), 4) )
    right = zeros(Float32, (Int64(nz), 4) )

    println("Calling 2D propagator")
    for it = 1:nt

	#if (it-1)%(nt/10) == 0
	#    @printf( "\n\t\tStep %d %f %f", it-1, pout[14, 342], source[it] )
	#end

	# needed for the ABC
        saving_abc2D!( pin, top, bottom, left, right )

	# PDE solver
        ccall( (:fwd_step2D, "./libprops.so"), Void, (Ptr{Void}, Ptr{Void},
                                                      Ptr{Cfloat}, Ptr{Cint},
                                                      Ptr{Cfloat}),
                                                      pin, pout, data.V, params,
                                                      coeffs )

	# source/receiver intro
        if mode == 3
	    rec_injection2D!( pin, data.V, rdata, coeffs[16], Int32(nt-it+1) )
	else
	    shot_injection2D!( pin, data.V, source[it], sloc, coeffs[16], Int32(it) )
	end

	# absorbing boundary condition
	#pin[1:STENCIL, 1:STENCIL] = rand(4,4)
	#pout[1:STENCIL, 1:STENCIL] = rand(4,4)
	higdon_abc!( pin, pout, data.V, top, bottom, left, right, params,
		     coeffs, abcoeffs, ctop, cbottom, cleft, cright )

	if mode == 1
	    waves[:,:] = waves[:,:] + pin[5:zmax, 5:xmax] .* pin[5:zmax, 5:xmax]
	else
	    waves[:,:,it] = pin[5:zmax, 5:xmax]
	    if mode == 2
		traces[it,:] = pin[5,5:xmax]
	    end
	end

        tmp = pin
        pin = pout
        pout = tmp
    end

    if mode == 2 
	return waves, traces
    else
	return waves
    end
end

function preparation!( data, freq, nt, coeffs, ctop, cbottom, 
		       cleft, cright, abcoeffs, params )

    vmax = maximum(data.V[:,:])
    dt = 8.0f0 * (0.09f0 * data.dz./vmax/convert(Float32,sqrt(2)))
    @show dt
    coeffs[ 14 ] = data.dx
    coeffs[ 15 ] = data.dz
    coeffs[ 16 ] = dt

    if check_num( data.dx, data.dz, dt, data.V, freq )
        quit()
    end

    compute_coeffs!( data.dz, data.dx, coeffs )
    (nz, nx) = size(data.V)
    compute_abcoeffs!( data.dz, data.dx, dt, nx, nz, data.V, 
                       ctop, cbottom, cleft, cright, abcoeffs )

    params[1] = STENCIL
    params[2] = nz
    params[3] = nx
    plane = (nz + Int32(2*STENCIL)) * (nx + Int32(2*STENCIL))
    params[5] = plane
    params[6] = Int32(2*plane)
    params[7] = Int32(3*plane)
    params[8] = Int32(4*plane)
    pln = nz * nx
    params[9] = pln 
    col = nz + Int32(2*STENCIL)
    params[10] = col  
    params[11] = Int32(2*col)
    params[12] = Int32(3*col)
    params[13] = Int32(4*col)
    aux = STENCIL + STENCIL*nz
    params[14] = aux
end

function compute_coeffs!( dz::Float32, dx::Float32, coeffs )
    dxi2 = 1.f0 / dx
    dxi2 *= dxi2
    b0x = bcoeffs[1] * dxi2
    coeffs[2] = bcoeffs[2] * dxi2
    coeffs[3] = bcoeffs[3] * dxi2
    coeffs[4] = bcoeffs[4] * dxi2
    coeffs[5] = bcoeffs[5] * dxi2

    dzi2 = 1.f0 / dz
    dzi2 *= dzi2
    b0z = bcoeffs[1] * dzi2
    coeffs[6] = bcoeffs[2] * dzi2
    coeffs[7] = bcoeffs[3] * dzi2
    coeffs[8] = bcoeffs[4] * dzi2
    coeffs[9] = bcoeffs[5] * dzi2

    coeffs[1] = b0x + b0z
end

function shot_injection2D!( p0::Array{Float32,2}, vel::Array{Float32,2},
                           data::Float32, source::Array{Int32,2},
                           dt::Float32, t::Int32 )
    const local twostencil::Int32 = Int32(2*STENCIL)
    const local nz::Int32 = Int32(size(p0,1))
    const local nx::Int32 = Int32(size(p0,2))
    const local jmin::Int32 = source[1,2] + STENCIL
    const local jmax::Int32 = min(jmin + twostencil -1 , nx - STENCIL)
    const local imin::Int32 = source[1,1] + STENCIL
    const local imax::Int32 = min(imin + twostencil -1, nz - STENCIL)
    local dist::Float32; local tmp::Float32; local tmp2::Float32
    local j::Int32; local i::Int32
    const local dt2::Float32 = dt * dt
    if abs(data) > 0.f0
        for j = jmin:jmax
            for i = imin:imax
                dist = Float32((j - jmin) * (j - jmin) +
                               (i - imin) * (i - imin) )
                dist = Float32(exp(-dist))
                tmp2 = vel[i-STENCIL+1,j-STENCIL+1]
                tmp = data * dist * dt2 * tmp2 * tmp2
                #if t == 2
                #    @printf(" %d %d %f %f %f %f %f %2.20f\n", 
                #    i, j, tmp, dist, p0[i,j], data, tmp2, p0[i,j] + tmp * dist)
                #end
                p0[i,j] = p0[i,j] + tmp * dist

            end
        end
    end
end

function check_num( dx, dz, dt, V, freq )
    # CFL check
    v_max = maximum(V)
    dt_max = min(dx, dz) * sqrt(3/8) / v_max
    if dt_max < dt
        println("WARNING: chosen dt violates CFL condition (dt_max = %f). The simulation can be unstable.\n", dt_max)
        return true
    end

    # Wavelength check
    v_min = minimum(V)
    np_min = v_min / freq / max(dx, dz)
    if np_min < 8
        println("WARNING: only %f grid points per wavelength in at least one dimension. This can lead to dispersion.\n", np_min)
        return true
    end
    return false
end

function ricker( freq, dt, nt )
    t_peak = Float32(1 / freq)
    t = linspace(0, (nt - 1) * dt, nt) - t_peak
    t = t.^2 * freq^2 * pi^2
    amplitudes = (1 - 2*t) .* exp(-t)
end

function saving_abc2D!( p0::Array{Float32, 2}, top::Array{Float32, 2},
                        bottom::Array{Float32, 2}, left::Array{Float32, 2},
                        right::Array{Float32, 2} )
    const local nx::Int32 = Int32(size(p0,2))
    const local nz::Int32 = Int32(size(p0,1))
    local cmax::Int32 = nx-STENCIL-STENCIL

    #@sync @parallel for k = 1:cmax
    @inbounds @simd for k = 1:cmax
        top[1, k] = p0[1, STENCIL+k]
        top[2, k] = p0[2, STENCIL+k]
        top[3, k] = p0[3, STENCIL+k]
        top[4, k] = p0[4, STENCIL+k]

        bottom[1, k] = p0[nz  , STENCIL+k]
        bottom[2, k] = p0[nz-1, STENCIL+k]
        bottom[3, k] = p0[nz-2, STENCIL+k]
        bottom[4, k] = p0[nz-3, STENCIL+k]
    end
    cmax = nz-STENCIL-STENCIL
    @inbounds @simd for k = 1:cmax
        left[k, 1] = p0[STENCIL+k, 1]
        left[k, 2] = p0[STENCIL+k, 2]
        left[k, 3] = p0[STENCIL+k, 3]
        left[k, 4] = p0[STENCIL+k, 4]

        right[k, 1] = p0[STENCIL+k, nx  ]
        right[k, 2] = p0[STENCIL+k, nx-1]
        right[k, 3] = p0[STENCIL+k, nx-2]
        right[k, 4] = p0[STENCIL+k, nx-3]
    end
end

function higdon_abc!( p0::Array{Float32, 2}, p1::Array{Float32, 2},
                      vel::Array{Float32, 2}, top::Array{Float32, 2},
                      bottom::Array{Float32, 2}, left::Array{Float32, 2},
                      right::Array{Float32, 2}, params, coeffs, abcoeffs,
                      ctop::Array{Float32, 2}, cbottom::Array{Float32, 2},
                      cleft::Array{Float32, 2}, cright::Array{Float32, 2} )
    local nx::Int32 = Int32(size(p0,2))
    local nz::Int32 = Int32(size(p0,1))
    local lmin::Int32 = STENCIL+1
    local lmax::Int32 = nx-STENCIL
    local cmax::Int32 = nx-Int32(2*STENCIL)
    local j::Int32, k::Int32, tmp::Float32, tmp2::Float32

    # top & bottom
    ccall( (:bc2Dz, "./libprops.so"), Void, (Ptr{Void}, Ptr{Void}, Ptr{Cfloat}, Ptr{Cfloat},
                                             Ptr{Cfloat}, Ptr{Cint}, Ptr{Cfloat}, Ptr{Cfloat}, 
                                             Ptr{Cfloat}, Ptr{Cfloat}),
                                             p0, p1, vel, top, 
                                             bottom, params, coeffs, abcoeffs, 
                                             ctop, cbottom )

    # left & right
    ccall( (:bc2Dx, "./libprops.so"), Void, (Ptr{Void}, Ptr{Void}, Ptr{Cfloat}, Ptr{Cfloat},
                                             Ptr{Cfloat}, Ptr{Cint}, Ptr{Cfloat}, Ptr{Cfloat}, 
                                             Ptr{Cfloat}, Ptr{Cfloat}),
                                             p0, p1, vel, left, 
                                             right, params, coeffs, abcoeffs, 
                                             cleft, cright )

end

function compute_abcoeffs!( dz, dx, dt, nx, nz, vel, 
                            ctop, cbottom, cleft, cright, abcoeffs )
    local beta2::Float32, beta1::Float32, b::Float32, qxt::Float32, qxt2::Float32

    beta2 = cos((pi/180.f0)*60.f0)
    beta1 = 1.f0
    b = 0.25f0
    qxt = b/(b - 1.f0)
    qxt2 = 2.f0 * qxt

    abcoeffs[2] = 1.f0 / (dz * dz)
    abcoeffs[3] = 1.f0 / (dx * dx)
    abcoeffs[1] = -2.f0 * (abcoeffs[2] + abcoeffs[3])

    # top
    local c::Array{Float32,2} = vel[1,:] .* (dt/dz)

    local tmp1::Array{Float32,2} = 1.f0 ./ ((c .+ beta1) .* (1.f0 - b))
    local qx::Array{Float32,2} = (b*beta1 .+ (b .* c) .- c) .* tmp1
    local qt::Array{Float32,2} = (b*beta1 .+ (b .* c) .- beta1) .* tmp1
    tmp1 = 1.f0 ./ ((c .+ beta2) .* (1.f0 - b))
    local rx::Array{Float32,2} = (b*beta2 .+ (b .* c) .- c) .* tmp1
    local rt::Array{Float32,2} = (b*beta2 .+ (b .* c) .- beta2) .* tmp1

    ctop[ 1, : ] = qx .+ rx
    ctop[ 2, : ] = qx .* rx
    ctop[ 3, : ] = qt .+ rt
    ctop[ 4, : ] = qx .* rt .+ qt .* rx .+ qxt2
    ctop[ 5, : ] = qxt .* (qx .+ rx)
    ctop[ 6, : ] = qt .* rt
    ctop[ 7, : ] = qxt .* (qt .+ rt)
    ctop[ 8, : ] = qxt .* qxt

    # bottom
    c = vel[nz,:] .* (dt/dz)

    tmp1 = 1.f0 ./ ((c .+ beta1) .* (1.f0 - b))
    qx = (b*beta1 .+ (b .* c) .- c) .* tmp1
    qt = (b*beta1 .+ (b .* c) .- beta1) .* tmp1
    tmp1 = 1.f0 ./ ((c .+ beta2) .* (1.f0 - b))
    rx = (b*beta2 .+ (b .* c) .- c) .* tmp1
    rt = (b*beta2 .+ (b .* c) .- beta2) .* tmp1

    cbottom[ 1, : ] = qx .+ rx
    cbottom[ 2, : ] = qx .* rx
    cbottom[ 3, : ] = qt .+ rt
    cbottom[ 4, : ] = qx .* rt .+ qt .* rx .+ qxt2
    cbottom[ 5, : ] = qxt .* (qx .+ rx)
    cbottom[ 6, : ] = qt .* rt
    cbottom[ 7, : ] = qxt .* (qt .+ rt)
    cbottom[ 8, : ] = qxt .* qxt

    # left & right
    for j = 1:nz
        # left
        local tmp::Float32 = vel[j,1] * (dt/dx)

        tmp12 = 1.f0 / ((tmp + beta1) * (1.f0 - b))
        qx2 = (b*(beta1 + tmp) - tmp) * tmp12
        qt2 = (b*(beta1 + tmp) - beta1) * tmp12
        tmp12 = 1.f0 / ((tmp + beta2) * (1.f0 - b))
        rx2 = (b*(beta2 + tmp) - tmp) * tmp12
        rt2 = (b*(beta2 + tmp) - beta2) * tmp12

        cleft[ j, 1 ] = qx2 + rx2
        cleft[ j, 2 ] = qx2 * rx2
        cleft[ j, 3 ] = qt2 + rt2
        cleft[ j, 4 ] = qx2 * rt2 + qt2 * rx2 + qxt2
        cleft[ j, 5 ] = qxt * (qx2 + rx2)
        cleft[ j, 6 ] = qt2 * rt2
        cleft[ j, 7 ] = qxt * (qt2 + rt2)
        cleft[ j, 8 ] = qxt * qxt

        # right
        tmp = vel[j,nx] * (dt/dx)

        tmp12 = 1.f0 / ((tmp + beta1) * (1.f0 - b))
        qx2 = (b*beta1 + (b * tmp) - tmp) * tmp12
        qt2 = (b*beta1 + (b * tmp) - beta1) * tmp12
        tmp12 = 1.f0 / ((tmp + beta2) * (1.f0 - b))
        rx2 = (b*beta2 + (b * tmp) - tmp) * tmp12
        rt2 = (b*beta2 + (b * tmp) - beta2) * tmp12

        cright[ j, 1 ] = qx2 + rx2
        cright[ j, 2 ] = qx2 * rx2
        cright[ j, 3 ] = qt2 + rt2
        cright[ j, 4 ] = qx2 * rt2 + qt2 * rx2 + qxt2
        cright[ j, 5 ] = qxt * (qx2 + rx2)
        cright[ j, 6 ] = qt2 * rt2
        cright[ j, 7 ] = qxt * (qt2 + rt2)
        cright[ j, 8 ] = qxt * qxt
    end
end

function load_data( name, data_dir )
    if name == "SEGsaltmodel"
        nz = 201
        nx = 676
        vel = zeros( Float32, (nz, nx) )
        fpath = "./$name.bin"
        f = open(fpath, "r")
        read!(f, vel)
        close(f)
        zero2iz = 40
        rampup2iz = 85
        data = Struct(vel, 10.f0, 10.f0, zero2iz, rampup2iz )
        return data
    end
    if name == "synthmodel2D"
        nz = 200
        nx = 200
        vel = zeros( Float32, (nz, nx) )
        fpath = "./$name.bin"
        f = open(fpath, "r")
        read!(f, vel)
        close(f)
        zero2iz = 50
        rampup2iz = 70
        data = Struct(vel, 10.f0, 10.f0, zero2iz, rampup2iz )
        return data
    end
end

function filtering2D( img )
    (mi, ni) = size(img)
    flt=zeros( Float32, (3,3) ) 
    flt[2,2] = 1
    (mf, nf) = size(flt)
    OUT = zeros( Float32, (mi+mf-1,ni+nf-1,8) ) 
    FLT=zeros(3,3,8)

    for i = 1:8
        flt1 = flt
        flt1[i] = -1
        FLT[:,:,i] = flt1
        out = conv2(flt1,img)
        OUT[:,:,i] = out
    end

    return OUT
end

# Reverse time migration (RTM)
function migrating( data, SourceVec::Array{Int32,2}, freq, nt, skip )
    local nz::Int32, nx::Int32

    (nz, nx) = size(data.V)
    local const num_shots = size(SourceVec, 1)
    topv = mean(data.V[1,:])

    # need the params and coeffs
    coeffs = Array(Float32, (16))
    abcoeffs = Array(Float32, (3))
    params = Array(Int32, (14))
    ctop = zeros(Float32, (8, Int64(nx) ) )
    cbottom = zeros(Float32, (8, Int64(nx) ) )
    cleft = zeros(Float32, (Int64(nz), 8) )
    cright = zeros(Float32, (Int64(nz), 8) )
    preparation!( data, freq, nt, coeffs, ctop, 
                  cbottom, cleft, cright, abcoeffs, params )
    source = zeros(Float32, (nt))
    source = convert( Array{Float32, 1}, ricker( freq, coeffs[16], nt ) )

    img = zeros(Float32, (Int64(nz), Int64(nx), num_shots ) )
    illR = zeros(Float32, (Int64(nz), Int64(nx), num_shots ) )
    #wss = zeros(Float32, (Int64(nz), Int64(nx), Int64(round(nt/skip)), num_shots ) )
    #wsr = zeros(Float32, (Int64(nz), Int64(nx), Int64(round(nt/skip)), num_shots ) )
    local iSource::Int32
    #tic()
    for iSource = 1:num_shots
        @printf( "shot %d, loc x %d, loc z %d\n", iSource, 
                 SourceVec[iSource,2], SourceVec[iSource,1] )

        rec_data = zeros(Float32, (Int64(nt), Int64(nx)))

        #tic()
        (wv_src, rec_data) = prop( data, freq, nt, source, coeffs, ctop, 
                                   cbottom, cleft, cright, abcoeffs, params, 
                                   SourceVec[iSource,:], rec_data, 2 )
        #toc()

        #tic()
        mask = direct_mask_sourcespec( nx, nt, data.dz, data.dx, coeffs[16], 
                                       SourceVec[iSource,:], topv, freq, 1, 20);
        rec_data = rec_data .* mask;
        #toc()

        #tic()
        wv_rec = prop( data, freq, nt, source, coeffs, ctop, 
                       cbottom, cleft, cright, abcoeffs, params, 
                       Int32[0 0], rec_data, 3 )
        #toc()

        #tic()
        img[:,:,iSource], illR[:,:,iSource] = imaging_condition( wv_src, wv_rec )
	#for j = 1:Int32(round(nt/skip))
	#    wss[:,:,j,iSource] = wv_src[:,:,j*100]
	#    wsr[:,:,j,iSource] = wv_rec[:,:,j*100]
	#end

        #toc()
    end
    #toc()

    return img, illR
    #return wss, wsr, img, illR
end

function direct_mask_sourcespec( nx, nt, dz, dx, dt, iSource, v, f,
				 buffer_numwaves, transition_timesteps )
    mask = zeros(Float32, (Int64(nt), Int64(nx)))
    delay = 1/f * v * 2;               # peak of ricker wavelet
    buffer = v / f * buffer_numwaves # one wavelength
    transition_zone = transition_timesteps * dt * v

    # iSource is assumed to be a tuple [iSourceZ, iSourceX]
    iSourceZ = iSource[1,1]
    iSourceX = iSource[1,2]

    for ix = 1:nx
        for it = 1:nt
            dst_travelled = (it - 1) * dt * v
            dst_source = sqrt((dx*(iSourceX - ix))^2 + (dz*(iSourceZ - 1))^2)
            if dst_travelled <= dst_source + buffer + delay
                mask[it, ix] = 0.f0
            elseif dst_travelled <= dst_source + buffer + delay + transition_zone
                mask[it, ix] = exp(((dst_travelled - dst_source - buffer - delay) / transition_zone)^2) / exp(1)
            else
                mask[it, ix] = 1.f0
            end
        end
    end
    return mask
end

################################################################
# divides by S tapered 1:0 top to bottom
################################################################
function imaging_condition( S, R )
    (nz, nx, nt) = size(S) # assumed same as for R
    img = zeros(Float32, (nz, nx) )
    illumS = zeros(Float32, (nz, nx) )
    illumR = zeros(Float32, (nz, nx) )
    for it = 1:nt
        img = img .+ S[:,:,it] .* R[:,:,nt - it + 1]
	illumS += S[:,:,it] .* S[:,:,it]
	illumR += R[:,:,it] .* R[:,:,it]
    end
    # simultaneously taper mean(abs(S)) 0:1
    scale = mean(abs(illumS))
    for i in 1:nz
            # taper = ( 1.f0 - cos((rampup2iz-i)/(rampup2iz-zero2iz)*pi) ) / 2.f0
            taper = (nz-i)/(nz-1)
            # taper down 1:0 actual illumination
            illumS[i,:] *= taper
            # taper up 0:1 average illumination scale
            illumS[i,:] += scale * (1.f0 - taper)
    end
    
    img = img ./ illumS
    return img, illumR 
end

function rec_injection2D!( p0::Array{Float32,2}, vel::Array{Float32,2},
                           source::Array{Float32,2}, dt::Float32, t::Int32 )
    const local nz::Int32 = Int32(size(p0,1))
    const local nx::Int32 = Int32(size(p0,2))
    const local jmin::Int32 = STENCIL + 1
    const local jmax::Int32 = nx - STENCIL
    local j::Int32; local i::Int32
    const local dt2::Float32 = dt * dt

    for j = jmin:jmax
	tmp2 = vel[1,j-STENCIL]
	tmp = dt2 * tmp2 * tmp2
	p0[jmin,j] = p0[jmin,j] + tmp * source[t,j-STENCIL]
    end
end

function savefile( data, filename )

    f = open( filename, "w" )
    write(f, data)
    close(f)

end
