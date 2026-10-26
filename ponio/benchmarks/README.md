# Benchmarks

The programs in this directory are profiling workloads. The examples write files and run for a few milliseconds, which is too short for a sampling profiler and spends most of the time in I/O. A benchmark writes nothing, takes its size and its method from the command line, and prints a throughput so that runs of different sizes or methods can be compared.

| Program | Problem | What it measures |
|---------|---------|------------------|
| `heat_bench` | heat equation 1D, `std::valarray` state | cost of ponio's stage combinations against the right hand side, for explicit and stabilized methods |

## Build

```
cmake . -G Ninja -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARKS=ON
cmake --build build
```

or `pixi run build_bench`, which keeps the default build type `RelWithDebInfo` (`-O2`). The benchmarks are compiled with `-g -fno-omit-frame-pointer` on top of the flags of the build type: `-g` gives the profiler source lines and inline frames without changing the generated code, frame pointers give it cheap call stacks.

## Run

```
./build/ponio/benchmarks/heat_bench --method rk_44 --nx 1000000 --steps 2000
```

```
method rk_44, nx 1000000, 2000 steps, 8000 rhs evaluations
time 4.01982 s, 497.534 Mcell.step/s, 1990.14 Mcell.rhs/s
checksum 1
```

`Mcell.step/s` is the number of cells advanced by one time step per second. `Mcell.rhs/s` divides by the number of right hand side evaluations instead, which is the figure to compare between methods with different numbers of stages. The program exits with an error when the time step is outside the stability domain of the method.

Choose `--steps` so that a run lasts a few seconds: the sampling error of a profiler falls as the square root of the number of samples.

## Profile with nunatak

[nunatak](https://github.com/hpc-maths/nunatak) runs `perf` on the program and attributes the time to source lines:

```
nunatak run --name heat-rk44 -- ./build/ponio/benchmarks/heat_bench --method rk_44 --nx 1000000 --steps 2000
```

Before comparing two runs with `nunatak compare`, profile the same binary twice and compare those two runs. The difference between them is the noise of the machine, and a smaller difference between two versions of the code means nothing.

## Memory placement

With a large `--nx`, a method with many stages reads many arrays at once: `rkc_202` combines up to 21 arrays of 8 MB for `--nx 1000000`. glibc allocates each of them with `mmap`, so they all start at the same offset in a page and compete for the same cache sets. On an aarch64 machine, `rkc_202` runs 1.6 times faster when glibc places these arrays on the heap instead:

```
GLIBC_TUNABLES=glibc.malloc.mmap_threshold=33554432 ./build/ponio/benchmarks/heat_bench --method rkc_202 --nx 1000000 --steps 50
```

A change of ponio that allocates its stages at another time can therefore change the throughput of these methods without changing the code that runs. Run both versions under the same tunable before concluding.
