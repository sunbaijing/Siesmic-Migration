	ESCI 570 - fall 2019
	Mauricio Araya, PhD

Main steps and notes:
=====================

1. Compile the propagation lib. gcc compiler is required
> make 

2. Run the provided Julia script, modeling2d.jl
> ./julia-ae26b25d43/bin/julia modeling2d.jl MODEL ACTION
where MODEL can be 1 for simple and 2 for SEG/EAGE slide
where ACTION can be 1 for modeling and 2 for illumination

3. Modify modeling2d.jl according to Plan.txt and then re-run

4. Modify imaging2D.jl according to Plan.txt and run
> ./julia-ae26b25d43/bin/julia imaging2D.jl

5. Graphics is an issue with Julia, use the python scripts:
> gather.py
> wavefield.py
> ill.py
> image.py
