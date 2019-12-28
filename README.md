# lbm3d-kokkos
## What is it?
lbm3d-kokkos is a C++ weakly compressible fluid flow solver written for CPUs or GPUs. The app is built on top of Sandia National Lab's performance portable Kokkos library. The majority of the app was written while attending a Kokkos workshop hosted by ORNL in July of 2018. Boundary conditions are specific to the 3D lid driven cavity. If you want something else, feel free to modify/extend the code.
## Build
```bash
mkdir proj
cd proj
git clone https://github.com/vasylivy/lbm3d-kokkos.git
```
Edit `lbmd3d-kokkos/lbm3d/Makefile` and set the device, architecture, and path for Kokkos. For example, setting `KOKKOS_DEVICES = OpenMP` and `KOKKOS_ARCH = SNB` will compile for Sandy Bridge Intel CPUs + OpenMP. Seting `KOKKOS_DEVICES = Cuda` and `KOKKOS_ARCH = SNB, Maxwell52`, will compile for a Sandy Bridge Intel CPU +  Nvidia GTX Titan X. A list of archiecture flags, as well as other Kokkos compilation flags can be found [here](https://github.com/kokkos/kokkos/wiki/Compiling). The default Kokkos path is `KOKKOS_PATH = ../../Kokkos/kokkos`. While still in the proj directory, execute the following cmds to download Kokkos
```bash
mkdir Kokkos
cd Kokkos
git clone https://github.com/kokkos/kokkos.git
```
To build the lbm3d executable
```bash
cd ../lbm3d-kokkos/lbm3d
mkdir release
make -j
```
## lattice Boltzmann method
The LBM implemention uses the BGK collision operator, D3Q19 quadrature for velocity space, and half-way bounce back for the boundary conditions. For collision-streaming, an A-B pattern is used. That is 19 distriubtions are read from array A, collided, and then streamed to array B. After applying boundary conditions, pointers to the arrays are then swapped. This collision-streaming pattern requires 2x the memory storage but is easy to implement. Details and further references regarding LBM can be found in Kruger's et al. text [[1]](https://link.springer.com/content/pdf/10.1007/978-3-319-44649-3.pdf).

## Running the 3D lid driven cavity
For an OpenMP build, the app is `lbm.host`. Similarly, `lbm.cuda` for a Cuda build. Navigate to `lbm3d-kokkos/cavity3d` folder and edit the bash script to your choosing. Executing the bash script `./run.sh` will output macroscopic variables under  `cavity3d/output/X.bin` according to the specified output rate until the steady state is reached or the number of maximum steps is exceeded.

## Validation
As the name suggests, the lid driven cavity is an enclosed flow driven by the movement of a lid. Here, the lid velocity is set to only have the x-component. Below a critical Reynolds number, the flow is steady. Figures below compare the horizontal component of velocity at the vertical centerline for different Reynolds numbers (`Re = 100, 400, 1000`) against Jiang et al. (1994) least-square FEM results. The LBM simulations use a resolution of `126x126x126`.

<img src="https://github.com/vasylivy/lbm3d-kokkos/blob/master/cavity3d/re100_128x128x128.png" alt="Re100" width="250"/>
<img src="https://github.com/vasylivy/lbm3d-kokkos/blob/master/cavity3d/re400_128x128x128.png" alt="Re400" width="250"/>
<img src="https://github.com/vasylivy/lbm3d-kokkos/blob/master/cavity3d/re1000_128x128x128.png" alt="Re1000" width="250"/>

## Performance
