#include <array>
#include <chrono>
#include <cmath>
#include <fstream>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <valarray>
#include <vector>

#include "algebra_types.h"
#include "algebra_utils.h"
#include "detector_setup.h"
#include "geometry_types.h"
#include "histogram.h"
#include "sampled_distribution.h"
#include "utilities.h"
#include "simulation_algorithms.h"

constexpr int g_verbosity { 3 };

auto main() -> int
{
    std::cout << "detector MC simulator\n";

    // set up random number generator
    std::random_device rd; // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    // the following parameters should be adopted to the individual situation
    constexpr std::size_t nr_events { 100'000 }; //<! the total number of tracks to be simulated
    constexpr double theta_max { toRad(90.) }; //<! the maximum theta angle taken into account
    constexpr double theta_step { toRad(1.) }; //<! the desired granularity of the simulated angular distributions
    // define rotation axis, here x-axis
    const Vector detector_rotation_axis { R3::Base::X };
    // define the rotation angle
    // set to 0, if setup should not be rotated
    constexpr double detector_rotation_angle { toRad(0.) };
    
    constexpr std::size_t nr_bins { static_cast<int>(theta_max / theta_step) + 1 };
    std::cout << "nr of bins: " << nr_bins << "\n";

    // vector of 2d polygon vertices defining the shape (in x-y-plane) of the detector.
    // note, that the points have to be in geometrical sequential order in
    // counter-clockwise orientation (i.e. a sequence of points defining the detector outline
    // going in ccw direction). The first point is considered to close the polygon together
    // with the last element in the vector

    // definition of the large double-paddle detector in IIPI-JLU lab
    const std::vector<Point> large_paddle_points_upper {
        { -150., -87.5 },
        { 150., -87.5 },
        { 150., 87.5 },
        { -150., 87.5 }
    };
    const std::vector<Point> large_paddle_points_lower {
        { -150., -100. },
        { 150., -100. },
        { 150., 100. },
        { -150., 100. }
    };

    // definition of the MuonPi standard-size (octagon) detector
    const std::vector<Point> octagon_points {
        { -126.5, -20. },
        { -91.5, -62.5 },
        { 91.5, -62.5 },
        { 126.5, -20. },
        { 126.5, 20. },
        { 91.5, 62.5 },
        { -91.5, 62.5 },
        { -126.5, 20. }
    };

    // definition of the MuonPi half-size detector
    const std::vector<Point> half_size_detector_points {
        { -130., -20. },
        { -60., -62.5 },
        { 0., -62.5 },
        { 0., 62.5 },
        { -60., 62.5 },
        { -130., 20. }
    };

    // definition of the MuonPi hexagon (small-size) detector
    constexpr double hex_length_a { 34.64 };
    constexpr double hex_length_b { 30.0 };
    const std::vector<Point> hexagon_detector_points {
        { -hex_length_a, 0. },
        { -hex_length_a / 2, -hex_length_b },
        { hex_length_a / 2., -hex_length_b },
        { hex_length_a, 0. },
        { hex_length_a / 2, hex_length_b },
        { -hex_length_a / 2, hex_length_b }
    };

    // definition of the large detector bars for JLU cosmic detector array
    const std::vector<Point> large_bar_points {
        { -500., -50. },
        { 500., -50. },
        { 500., 50. },
        { -500., 50. },
    };

    const std::vector<Point> ref_volume_points {
        { -500., -500. },
        { 500., -500. },
        { 500., 500. },
        { -500., 500. },
    };

    // create 3d objects of type ExtrudedObject defined by the 2d outline,
    // a global position offset and a thickness
    ExtrudedObject ref_volume { ref_volume_points, { 0., 0., -200. }, 400. };
    //ref_volume.set_position( R3::Origin );
    ExtrudedObject detector1 { large_paddle_points_lower, { 0., 0., -100. }, 7. };
    ExtrudedObject detector2 { large_paddle_points_upper, { 0., 0., 100. }, 7. };

    // create 3d objects of type ExtrudedObject but using the constructor for generation of a
    // circular shape specified by a global position offset, radius, thickness and an optional
    // number of vertex points to generate the circle
    ExtrudedObject round_detector1 { { 0., 0., 0. }, 50., 10. };
    ExtrudedObject round_detector2 { { 0., 0., 100. }, 50., 10. };

    // construct a detector setup with two detectors
    DetectorSetup setup { { /*ref_volume,*/ detector1, detector2 } };

    // simulate the effective area (geometric aperture) at theta=0 of the detector system
    [[maybe_unused]] const double effective_area_sqm { simulate_geometric_aperture(setup, gen, nr_events) };

    // add a rotation to the system
    setup.rotate(detector_rotation_axis, detector_rotation_angle);

    // uncomment the following block to calculate the double differential acceptance
    // as function of phi and theta
    /*
    [[maybe_unused]] const auto acceptance_phi_theta = theta_phi_scan<361, 46>(setup, gen, nr_events, 0., theta_max, -pi(), pi());
*/

    // initialize the histogram vector
    std::vector<Histogram> histos {};

    // run a scan over theta angle (uniformly distributed)
    // to record the detector acceptance, if required
//    Append(histos, theta_scan(setup, gen, nr_events, 0., theta_max, nr_bins));

    // now, run the full simulation and append the resulting histograms
    // to the already existing histogram vector
    //cosmic_simulation(setup, gen, nr_events * nr_bins, &histos, nr_bins, theta_max);

    auto acceptance_dataseries { cosmic_simulation_detector_sweep(setup, gen, nr_events, detector_rotation_axis, toRad(-90.), toRad(90.), 90) };
    acceptance_dataseries.export_file("detector_sweep_acceptances.dat");

    // export each histogram into a separate file (human readable ASCII format)
    for (auto histo : histos) {
        histo.export_file(histo.getName() + ".hist");
    }

    exit(0);
}
