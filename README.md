# Caustics Engineering

This repo is for generating 3D surface meshes that project caustic images. It is written in Julia.

See the write-up [here](https://mattferraro.dev/posts/caustics-engineering)!

# To Run from source

To run the cat example from the blogpost, run line by line from ` src/scratchpad.jl`.

_OR_ from the command line.

```
julia ./run.jl
```
The image file is currently hard-coded.

# To Run in docker

In this case you don't have install Julia locally - instead you will be able to build and run the project in a docker container.
To do that, follow the guide:

1. Put your source image to `data/input.jpg`. This is the file the project will process.
1. Clean your build (you can omit this step if you are running the code for the first time).
It is required if you modify your source code only (changing input image doesn't require to execute this step):

    ```sh
    make clean 
    ```

1. Build the project:

    ```sh
    make build
    ```

1. Run the project:

    ```sh
    make run
    ```

1. Find your output in `data/output.obj`.
