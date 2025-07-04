# Gravity1D
Exact N-body integrator for 1D gravity.

## Aim of the code

This repository aims to compute the exact evolution of a one-dimensional self-gravitating system of N particles. It is a heap-based, event-driven algorithm which computes the collision times of each particles, as described by [Noullez & al. (2001)](https://ui.adsabs.harvard.edu/abs/2001cond.mat..1336N/abstract). Following the observations of [Schulz & al. (2013)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.431...49S/abstract), we use a 106-bit mantissa (via the `DoubleFloats.jl` package) to suppress the accumulation of roundoff errors during long simulations.

## Installation

Install Julia by following the instructions at [https://julialang.org/downloads/platform/](https://julialang.org/downloads/platform/).

To invoke Julia in the Terminal, you need to make sure that the julia command-line program is in your `PATH`. 

On MacOS, we may have to create a link in `/usr/local/bin`.

## Julia packages

Open the terminal in the folder `packages` and type

```
$ julia Install_pkg.jl
```

to install the necessary packages.

## Local Julia environment

One may create a local Julia environment by following the instructions at [https://jkrumbiegel.com/pages/2022-08-26-pkg-introduction/](https://jkrumbiegel.com/pages/2022-08-26-pkg-introduction/).

One can create a local environment called `LocalEnv` by following the next instructions. First, open `Julia` on the console.

```
$ module load julia/1.11.3
$ julia
```

Then, open the Pkg mode and create the local environment

```
julia> ]
(@v1.11) pkg> generate LocalEnv
```


One may import packages within that local environment by following the instructions at [https://pkgdocs.julialang.org/v1/environments/](https://pkgdocs.julialang.org/v1/environments/).

First, go to the folder containing `LocalEnv/` and open Julia. Then, add the packages (here, `SpecialFunctions`) in Pkg mode after activating the local environment.

```
julia> ]
(@v1.11) pkg> activate LocalEnv
(LocalEnv) pkg> add SpecialFunctions
```

One may also deactivate the local environment by writing in Pkg mode

```
(LocalEnv) pkg> activate
```

Alternatively, one may activate the local environment without having to switch to Pkg mode by writing 

```
julia> using Pkg; Pkg.activate("LocalEnv")
```

Similarly, deactivating the local environment can by done through the line

```
julia> using Pkg; Pkg.activate()
```


## Integration

One may run a simulation by going to the `src/` folder and running the bash command

```
$ julia Main.jl
```

We refer to the file `Args.jl` for a list of the arguments of the code. As an illustration, the command

```
$ julia Main.jl --N 100 --tmax 100.0 --save_freq 10 --model plummer --seed 0 --output output --save_final_state true
```

runs a simulation of 100 particles, randomly sampled from a Plummer distribution with random seed 0, until `t=100` time units. It saves a snapshot every 10 time units using the name `output_` in a `data/` folder created for this purpose. Finally, it saves the exact final state of the system at `t=100` in a `restart_data_output_seed_0.jld2` binary file in a `data/restart/` folder. In that case, the run may be restarted from that exact final state via the command

```
$ julia Main.jl --N 100 --tmax 200.0 --save_freq 10 --model plummer --seed 0 --output output --save_final_state true --restart restart_data_output_seed_0.jld2
```

starting at `t=100` and ending at `t=200`.

The models that can be used as initial conditions are:

- `plummer` : A one-dimensional analog to the Plummer model (see [Roule & Fouvry, 2022](https://ui.adsabs.harvard.edu/abs/2022PhRvE.106d4118R/abstract)),
- `harmonic` : The one-dimensional harmonic system,
- `cold` : A one-dimension, cold homogeneous system (see [Schulz & al., 2013](https://ui.adsabs.harvard.edu/abs/2013MNRAS.431...49S/abstract)).


## Post-treatment

We provide a series of Julia scripts to treat the data saved during the simulation in the `post_treatment` folder.


## Running a job array on a cluster with SLURM

We provide a bash script using SLURM to launch a job array on a computing cluster in the folder `slurm`. It may be used via the command

```
$ sbatch job_run.sh
```