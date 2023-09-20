# SIG code

The "SIG.jl" script provides the SIG recipo to reproduce the magnetic field mapped in the El Gordo galaxy cluster (see  Hu et al. 2023, arXiv:2306.10011). 

"ElGordo_hres_sync_blanked.fits" is the MeerKAT 1.28 GHz image of the El Gordo galaxy cluster. The SIG.jl works also for other clusters, by changing the input images.

#### 1. System requirements

Julia language (version > 1.1, see https://julialang.org/)

#### 2. Installaiton guide 

To run SIG.jl, one needs to install some necessary packages in Julia:

```julia
using Pkg
Pkg.add(["PyCall","PyPlot","LsqFit","Images","LinearAlgebra","StatsBase","FITSIO","Statistics","ImageFiltering]")
```

Then we need to include the two external packages "LazCore.jl" and "LazType.jl". In the root directory of Julia, which is typically "C:\Users\your_user_name\\.julia", creat a new directory "config" under the "./julia" directory and then put the "starup.jl" in the "config" direcotry. 

In the starup.jl, please modify the "MODULE="F:\\julia" to the location where you save "LazCore.jl" and "LazType.jl" files.

```julia
println("Hello World!")
MODULE="F:\\julia";
push!(LOAD_PATH,MODULE)
```

Restart the Julia terminal, you should see "Hello World" at the beginning.

#### 3. Demo

A demo of SIG is provided in the SIG.jl file. The running time depends on how fast the cpu is. Typically it takes several hours and will produce the Fig. 6 in Hu et al. (2023).

#### 4. Instructions for use

To reproduce other figures 1, 2, 4 , 5 in Hu et al. (2023), one can change the input image (ElGordo_hres_sync_blanked.fits) to the corresponding images of that cluster.