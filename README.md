# FISHPACK in C

# Quickstart
Currently, only macos is supported, but if you have raylib and glfw locally installed and can modify the Makefile to point to the correct thirdparty library you can make it work on any platform.

## Build on macOS
Run `./setup_macos.sh` to download local copies of `raylib-5.5_macos` and `glfw-3.4_macos` into `thirdparty/`.

After having downloaded the required third party libraries, you can build and run the code with `make`.

![viridis_bump](https://github.com/user-attachments/assets/a413a415-d7e9-45b8-b29c-5c9a06ee240d)

![hot_cold_sinusoidal](https://github.com/user-attachments/assets/bbe26f84-5800-47d7-9fab-389d0947278b)

![grayscale](https://github.com/user-attachments/assets/80c790c0-b0a7-427b-94cd-4baf98a36307)

## Options in the rendering
When the rendiring is running you have the following options
- Press `Enter` to toggle orbiting.
- Press the `0` key to toggle the wirefarme of the mesh.

### TODO
- [x] implement `plot_surface` with `GenMeshHeightmap(Image heightmap, Vector3 size)`
- [x] after `plot_surface` is availabe, see if `solve_poisson2d` works as intended
- [x] implement inhom Dirichlet boundary for 2d Poisson
- [x] implement 9 point stencil for 2d Poisson
- [x] add shaders to modify the colormap of the surface plot
- [x] allow the viewer to look around in the surface plot
- [x] implement viridis colorplot
- [x] Move out the parameters to a params struct, its getting crowded for the function signature
- [x] try adding lighting to the surface plotting scene to see if it makes the curves more pronounced
- [x] implement plasma colormap, and some other colormaps that look nice
- [x] indicate what colormaps are available
- [x] make the cmap a settable parameter of the plotter
- [x] refactor the colormap code to have only one fragment shader and read in the colormap as a uniform array of vec3
- [x] move the colormaps into separate headers
- [ ] make the lighting smoother on the surface plot by interpolating the surface normals on the triangle faces
