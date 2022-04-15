# Caustics Engineering

This repo is for generating 3D surface meshes that project caustic images. It is written in Julia.

See the write-up [here](https://mattferraro.dev/posts/caustics-engineering)!

# To install
If you don't have the packages, run in this repo `julia`, and:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile() 
```

# To Run

To run the cat example from the blogpost, run the file without an image path.

```
julia ./run.jl [path-to-image.png]
```
