# Siesmic-Migration
Run the code in the Oracel VM VirtualBox	
Main steps and notes:
=====================
su root
1. Compile the propagation lib. gcc compiler is required
> make 

2. Run the provided Julia script, modeling2d.jl
> ./julia-ae26b25d43/bin/julia modeling2d.jl MODEL ACTION
where MODEL can be 1 for simple and 2 for SEG/EAGE slide
where ACTION can be 1 for modeling and 2 for illumination

3. Modify modeling2d.jl with different location of shots  and run to get the illumination

4. Modify imaging2D.jl according to with different shot location and run tp get the image
> ./julia-ae26b25d43/bin/julia imaging2D.jl

5. Graphics is an issue with Julia, use the python scripts:
> gather.py
> wavefield.py
> ill.py
> image.py
