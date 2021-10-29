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
    std::cout << "Particle position dim: " << particle_pos.ndim() << std::endl;
    std::cout << "Particle position shape: " << particle_pos.shape(0) << ", " << particle_pos.shape(1) << std::endl;
    std::cout << "Particle position flags: " << particle_pos.flags() << std::endl;
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
    if(rotation_angle) {
        rotate3d(particle_data, n_shuffle, rotation_angle->data(), center.data());
    }

    std::cout << "Calling QHULL with " << n_shuffle << " particles (originally " << n_particles << ")" << std::endl;
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

PYBIND11_MODULE(pydtfe, m) {
    m.doc() = "python interface to SDTFE";
    m.def("compute_surface_density", &compute_density_py,
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
}