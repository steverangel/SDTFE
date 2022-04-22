#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <iostream>

extern "C" {
#include "utils.h"
#include "triangulate.h"
#include "dtfe.h"
}

namespace py = pybind11;


py::array_t<double> compute_density_py(
    py::array_t<float, py::array::f_style | py::array::forcecast> particle_pos,
    std::array<double, 3> center,
    unsigned grid_dim,
    float sample_factor,
    double box_width,
    double box_depth,
    double particle_mass,
    double mc_box_width,
    unsigned mc_sample_count,
    std::optional<std::array<double, 3>> rotation_angle
) {
    // convert to qhdata
    // std::cout << "Particle position dim: " << particle_pos.ndim() << std::endl;
    // std::cout << "Particle position shape: " << particle_pos.shape(0) << ", " << particle_pos.shape(1) << std::endl;
    // std::cout << "Particle position flags: " << particle_pos.flags() << std::endl;
    if(particle_pos.shape(1) != 3) {
        throw py::value_error("incorrect shape of particle position array");
    }
    size_t n_particles = particle_pos.shape(0);
    double *particle_data = convert_to_qhdata(particle_pos.data(), n_particles);

    // shuffle and rotate
    size_t n_shuffle = n_particles*sample_factor;
    if(n_shuffle < 5) {
        throw py::value_error("not enough particles");
    }
    shuffle_buff(particle_data, n_particles, n_shuffle);
    // adjust particle mass for subsampling
    particle_mass /= sample_factor;

    if(rotation_angle) {
        rotate3d(particle_data, n_shuffle, rotation_angle->data(), center.data());
    }

    // std::cout << "Calling QHULL with " << n_shuffle << " particles (originally " << n_particles << ")" << std::endl;
    // triangulate
    size_t n_tetra;
    int *tetra_data = triangulate(particle_data, n_shuffle, &n_tetra, "qhull d Qz Qt Qbb");

    // calculate density
    py::array_t<double> rho({grid_dim, grid_dim});
    std::fill(rho.mutable_data(), rho.mutable_data() + grid_dim*grid_dim, 0.);
    compute_density(particle_data, n_shuffle, tetra_data, n_tetra, grid_dim, box_width, box_depth, rho.mutable_data(), particle_mass, center[0], center[1], center[2], mc_sample_count, mc_box_width);

    // cleanup
    free(tetra_data);
    free(particle_data);

    return rho;
}

py::array_t<double> compute_3d_density_py(
    py::array_t<float, py::array::f_style | py::array::forcecast> particle_pos,
    std::array<double, 3> center,
    unsigned grid_dim,
    float sample_factor,
    double box_length,
    double particle_mass,
    unsigned supersampling
) {
    // convert to qhdata
    // std::cout << "Particle position dim: " << particle_pos.ndim() << std::endl;
    // std::cout << "Particle position shape: " << particle_pos.shape(0) << ", " << particle_pos.shape(1) << std::endl;
    // std::cout << "Particle position flags: " << particle_pos.flags() << std::endl;
    if(particle_pos.shape(1) != 3) {
        throw py::value_error("incorrect shape of particle position array");
    }
    size_t n_particles = particle_pos.shape(0);
    double *particle_data = convert_to_qhdata(particle_pos.data(), n_particles);

    // shuffle and rotate
    size_t n_shuffle = n_particles*sample_factor;
    if(n_shuffle < 5) {
        throw py::value_error("not enough particles");
    }
    shuffle_buff(particle_data, n_particles, n_shuffle);
    // adjust particle mass for subsampling
    particle_mass /= sample_factor;

    // std::cout << "Calling QHULL with " << n_shuffle << " particles (originally " << n_particles << ")" << std::endl;
    // triangulate
    size_t n_tetra;
    int *tetra_data = triangulate(particle_data, n_shuffle, &n_tetra, "qhull d Qz Qt Qbb");

    // calculate density
    py::array_t<double> rho({grid_dim, grid_dim, grid_dim});
    std::fill(rho.mutable_data(), rho.mutable_data() + grid_dim*grid_dim*grid_dim, 0.);
    compute_3d_density(particle_data, n_shuffle, tetra_data, n_tetra, grid_dim, box_length,
    rho.mutable_data(), particle_mass, center[0], center[1], center[2], supersampling);

    // cleanup
    free(tetra_data);
    free(particle_data);

    return rho;
}

py::array_t<double> interpolate_field2grid(
    py::array_t<float, py::array::f_style | py::array::forcecast> particle_pos,
    py::array_t<float> particle_value,
    std::array<double, 3> center,
    unsigned grid_dim,
    float sample_factor,
    double box_length,
    unsigned supersampling
) {
    // convert to qhdata
    // std::cout << "Particle position dim: " << particle_pos.ndim() << std::endl;
    // std::cout << "Particle position shape: " << particle_pos.shape(0) << ", " << particle_pos.shape(1) << std::endl;
    // std::cout << "Particle position flags: " << particle_pos.flags() << std::endl;
    if(particle_pos.shape(1) != 3) {
        throw py::value_error("incorrect shape of particle position array");
    }
    size_t n_particles = particle_pos.shape(0);
    if(n_particles != particle_value.shape(0)) {
        throw py::value_error("different length of particle position and value arrays");
    }

    double *particle_data = convert_to_qhdata(particle_pos.data(), n_particles);

    // shuffle and rotate
    // TODO: if we want to enable this functionality, we need to make sure we shuffle the
    // particle_value array accordingly
    if(sample_factor != 1.0f) {
        throw py::value_error("`sample_factor != 1.0f` not supported at the moment");
    }
    /*
    size_t n_shuffle = n_particles*sample_factor;
    if(n_shuffle < 5) {
        throw py::value_error("not enough particles");
    }
    shuffle_buff(particle_data, n_particles, n_shuffle);
    */
    size_t n_shuffle = n_particles;

    // std::cout << "Calling QHULL with " << n_shuffle << " particles (originally " << n_particles << ")" << std::endl;
    // triangulate
    size_t n_tetra;
    int *tetra_data = triangulate(particle_data, n_shuffle, &n_tetra, "qhull d Qz Qt Qbb");
    for(size_t i=0; i<n_particles; ++i) {
        particle_data[i*4 + 3] = particle_value.data()[i];
    }

    // calculate density
    py::array_t<double> grid({grid_dim, grid_dim, grid_dim});
    std::fill(grid.mutable_data(), grid.mutable_data() + grid_dim*grid_dim*grid_dim, 0.);
    interpolate_to_3d_grid(particle_data, n_shuffle, tetra_data, n_tetra, grid_dim, box_length,
    grid.mutable_data(), center[0], center[1], center[2], supersampling);

    // cleanup
    free(tetra_data);
    free(particle_data);

    return grid;
}

py::array_t<float> compute_particle_density_py(
    py::array_t<float, py::array::f_style | py::array::forcecast> particle_pos,
    double particle_mass
) {
    // convert to qhdata
    if(particle_pos.shape(1) != 3) {
        throw py::value_error("incorrect shape of particle position array");
    }
    size_t n_particles = particle_pos.shape(0);
    double *particle_data = convert_to_qhdata(particle_pos.data(), n_particles);

    // triangulate
    size_t n_tetra;
    int *tetra_data = triangulate(particle_data, n_particles, &n_tetra, "qhull d Qz Qt Qbb");

    // calculate density at particle location
    pre_compute_vol(&particle_data, n_particles, tetra_data, n_tetra, particle_mass);

    // particle_data is an array of size [n_particles, 4]
    // the 4th column should contain the densities
    py::array_t<float> rho(n_particles);
    for(size_t i=0; i<n_particles; ++i) {
        rho.mutable_data()[i] = particle_data[4*i + 3];
    }

    // cleanup
    free(tetra_data);
    free(particle_data);

    return rho;
}

PYBIND11_MODULE(pydtfe, m) {
    m.doc() = "python interface to SDTFE";
    m.def("compute_surface_density", &compute_density_py, R"Delim(
computes the surface density using the DTFE tesselation of a particle
distribution

Parameters
----------
particle_pos: np.ndarray
    3d positions of the particles (shape ``(N, 3)``)

center: List[float]
    position of the halo center (length 3)

grid_dim: int
    resolution of the surface density map (will be square)

sample_factor: float
    subsampling factor between [0, 1]. If smaller than 1.0, the code will only
    use a random subset of the particles, corresponding to the fraction
    specified

box_width: float
    physical size of the width (and height) of the map

box_depth: float
    physical depth of the surface density map (density will be integrated
    between [center_z - box_depth/2, center_z + box_depth/2])

particle_mass: float
    mass of a single particle

mc_box_width: float
    width of monte-carlo sampling of rays in each surface cell

mc_sample_count: int
    number of monte-carlo samples to be drawn per cell

rotation_angle: Optional[List[float]]

Returns
-------
rho: np.ndarray
    the surface density map, in units of ``mass / length^2``
    (shape: ``(grid_dim, grid_dim)``)
)Delim",
        py::arg("particle_pos"),
        py::arg("center"),
        py::arg("grid_dim"),
        py::arg("sample_factor"),
        py::arg("box_width"),
        py::arg("box_depth"),
        py::arg("particle_mass"),
        py::arg("mc_box_width"),
        py::arg("mc_sample_count"),
        py::arg("rotation_angle") = nullptr
    );

    m.def("compute_3d_density", &compute_3d_density_py, R"Delim(
computes the volumetric density using the DTFE tesselation of a particle
distribution

Parameters
----------
particle_pos: np.ndarray
    3d positions of the particles (shape ``(N, 3)``)

center: List[float]
    position of the halo center (length 3)

grid_dim: int
    resolution of the density map (will be a cube)

sample_factor: float
    subsampling factor between [0, 1]. If smaller than 1.0, the code will only
    use a random subset of the particles, corresponding to the fraction
    specified

box_depth: float
    physical size of the density cube

particle_mass: float
    mass of a single particle

supersampling: int
    number of sampled points in each volume cell. The density will be averaged
    over a ``supersampling ** 3`` grid within the volume cell. For example, if
    set to 1, the density will be evaluated once at the center of each cell.
    For 2, it will be evaulated on a regular 2x2x2 mesh and averaged.

Returns
-------
rho: np.ndarray
    the volumetric density map, in units of ``mass / length^3``
    (shape: ``(grid_dim, grid_dim, grid_dim)``)

)Delim",
        py::arg("particle_pos"),
        py::arg("center"),
        py::arg("grid_dim"),
        py::arg("sample_factor"),
        py::arg("box_length"),
        py::arg("particle_mass"),
        py::arg("supersampling")=1u
    );

    m.def("interpolate_field2grid", &interpolate_field2grid, R"Delim(
interplates a scalar field from particle locations to a grid using the DTFE

Parameters
----------
particle_pos: np.ndarray
    3d positions of the particles (shape ``(N, 3)``)

particle_value: np.ndarray
    value of the scalar field at the particle position (shape ``(N,)``)

center: List[float]
    position of the halo center (length 3)

grid_dim: int
    resolution of the density map (will be a cube)

sample_factor: float
    subsampling factor between [0, 1]. If smaller than 1.0, the code will only
    use a random subset of the particles, corresponding to the fraction
    specified

box_depth: float
    physical size of the density cube

supersampling: int
    number of sampled points in each volume cell. The density will be averaged
    over a ``supersampling ** 3`` grid within the volume cell. For example, if
    set to 1, the density will be evaluated once at the center of each cell.
    For 2, it will be evaulated on a regular 2x2x2 mesh and averaged.

Returns
-------
grid: np.ndarray
    the scalar field interpolated to the grid points
    (shape: ``(grid_dim, grid_dim, grid_dim)``)


Note
----
currently, only `sampling_factor=1.0` is allowed. Adding support for subsampling
requires some changes to the code.

)Delim",
        py::arg("particle_pos"),
        py::arg("particle_value"),
        py::arg("center"),
        py::arg("grid_dim"),
        py::arg("sample_factor"),
        py::arg("box_length"),
        py::arg("supersampling")=1u
    );

    m.def("compute_particle_density", &compute_particle_density_py, R"Delim(
computes the DTFE density at the particle locations

Parameters
----------
particle_pos: np.ndarray
    3d positions of the particles (shape ``(N, 3)``)

particle_mass: float
    mass of a single particle

Returns
-------
rho: np.ndarray
    the density at the particle locations, in units of ``mass / length^3``
    (shape: ``(n_particles)``)

)Delim",
        py::arg("particle_pos"),
        py::arg("particle_mass")
    );
}