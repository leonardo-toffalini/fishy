# FISHPACK in C

![viridis_bump](https://github.com/user-attachments/assets/a413a415-d7e9-45b8-b29c-5c9a06ee240d)

![hot_cold_sinusoidal](https://github.com/user-attachments/assets/bbe26f84-5800-47d7-9fab-389d0947278b)

![grayscale](https://github.com/user-attachments/assets/80c790c0-b0a7-427b-94cd-4baf98a36307)

### TODO
- [x] implement `plot_surface` with `GenMeshHeightmap(Image heightmap, Vector3 size)`
- [x] after `plot_surface` is availabe, see if `solve_poisson2d` works as intended
- [x] implement inhom Dirichlet boundary for 2d Poisson
- [x] implement 9 point stencil for 2d Poisson
- [x] add shaders to modify the colormap of the surface plot
- [x] allow the viewer to look around in the surface plot
- [x] implement viridis colorplot
- [ ] Move out the parameters to a params struct, its getting crowded for the function signature
- [x] try adding lighting to the surface plotting scene to see if it makes the curves more pronounced
- [x] implement plasma colormap, and some other colormaps that look nice
- [x] indicate what colormaps are available
- [x] make the cmap a settable parameter of the plotter
