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

constexpr int g_verbosity { 3 };
constexpr double toDeg(double x) { return x * 180 / pi(); }
constexpr double toRad(double x) { return x * pi() / 180; }

/** @brief The cumulative distribution function (CDF) for the cos^2(x) PDF distribution
 * This CDF is used for the calculation of the Probability Density Function (PDF)
 * of the cos^2(x) distribution for generating random values following the angular distribution
 * of the muon tracks
 */
auto cos2cdf = [](double x) {
    //return cdf to following pdf: cos^2(x)
    return (2 / pi()) * (x / 2. + sin(2. * x) / 4.) + 0.5; //from Wolfram Alpha
};

template <typename T>
void Append(std::vector<T>& a, const std::vector<T>& b)
{
    a.reserve(a.size() + b.size());
    a.insert(a.end(), b.begin(), b.end());
}

/** @brief theta_scan
 * This function does Monte-Carlo simulated evaluation of particle hits for a given detector setup iterating over the specified theta angular range. Tracks are distributed as follows:
 * - phi: uniformly distributed in -pi..pi
 * - theta: iterating from theta_min to theta_max in nr_bins steps
 * @param setup DetectorSetup object containing the detector definitions
 * @param gen random number generator (Mersenne-Twister engine of type std::mt19937)
 * @param nr_events The number of tracks to be generated per theta bin. This is not identical to the number of MC events (=events having a hit in the reference detector)
 * @param theta_min the minimum zenith angle theta from which events are generated
 * @param theta_max the maximum zenith angle theta up to which events are generated
 * @param nr_bins The number of bins in the resulting acceptance histogram
 * @return A std::vector of Histogram objects containing the resulting distributions. Currently, only the angular acceptance (detection efficiency vs. theta) is provided
 * @note In case of an error, the returned histogram vector is empty
 * @note The setup object must have the ref_detector iterator set to any valid detector it contains.
*/
std::vector<Histogram> theta_scan(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_min, double theta_max, std::size_t nr_bins)
{
    std::vector<Histogram> histos {};

    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr << "no reference detector defined in DetectorSetup!\n";
        return histos;
    }

    Histogram acc_hist("acceptance_scan_theta",
        nr_bins,
        0., theta_max);

    const auto bounds { setup.ref_detector()->bounding_box() };

    std::uniform_real_distribution<> distro_x {
        bounds.first[0],
        bounds.second[0],
    };
    std::uniform_real_distribution<> distro_y {
        bounds.first[1],
        bounds.second[1],
    };
    std::uniform_real_distribution<> distro_z {
        bounds.first[2],
        bounds.second[2]
    };
    std::uniform_real_distribution<> distro_phi(-pi(), pi());

    const double theta_step { (theta_max - theta_min) / (nr_bins - 1) };
    std::cout << "theta step: " << theta_step <<" rad = " << toDeg(theta_step) << " deg\n";
    std::cout << "#theta acceptance acceptance_error\n";
    double theta { theta_min };
    for (size_t bin { 0 }; bin < nr_bins; ++bin) {
        std::size_t mc_events { 0 };
        std::size_t detector_events { 0 };
        for (std::size_t n = 0; n < nr_events; ++n) {
            double phi { distro_phi(gen) };
            Line line { Line::generate(
                { distro_x(gen), distro_y(gen), distro_z(gen) }, theta, phi) };

            bool coincidence { false };
            LineSegment refdet_path { setup.ref_detector()->intersection(line) };
            if (refdet_path.length() > DEFAULT_EPSILON) {
                mc_events++;
                coincidence = true;
            }
            for (auto detector { setup.detectors().cbegin() };
                 detector != setup.detectors().end();
                 ++detector) {
                if (detector == setup.ref_detector())
                    continue;
                LineSegment det_path { detector->intersection(line) };
                if (det_path.length() < DEFAULT_EPSILON) {
                    coincidence = false;
                }
            }
            if (coincidence) {
                //std::cout << n << " " << std::setw(2) << toDeg(theta) << " " << toDeg(phi) << " " << det1_path.length() << " " << det2_path.length() << "\n";
                detector_events++;
            }
        }
        std::cout << toDeg(theta) << " " << static_cast<double>(detector_events) / mc_events << " " << std::sqrt(detector_events) / mc_events << "\n";
        //        acc_hist.fill(theta, detector_events);
        acc_hist.fill(theta, static_cast<double>(detector_events) / mc_events);
        std::cout << std::flush;
        theta += theta_step;
    }
    histos.push_back(acc_hist);
    return histos;
}

double simulate_geometric_aperture(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta = 0.)
{
    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr << "no reference detector defined in DetectorSetup!\n";
        return {};
    }

    auto bounds { setup.ref_detector()->bounding_box() };
    Point dimensions { bounds.second - bounds.first };
    std::cout << "detector bounds: min=" << bounds.first << " max=" << bounds.second << "\n";
    std::cout << "detector dimensions=" << dimensions << "\n";

    bounds.first -= dimensions * 5;
    bounds.second += dimensions * 5;

    std::cout << "simulation bounds: min=" << bounds.first << " max=" << bounds.second << "\n";

    std::uniform_real_distribution<> distro_x {
        bounds.first[0],
        bounds.second[0]
    };
    std::uniform_real_distribution<> distro_y {
        bounds.first[1],
        bounds.second[1]
    };
    const double simulation_plane_z_pos { setup.ref_detector()->bounding_box().first[2] };
    std::uniform_real_distribution<> distro_z {
        setup.ref_detector()->bounding_box().first[2],
        setup.ref_detector()->bounding_box().second[2]
    };
    std::uniform_real_distribution<> distro_phi(-pi(), pi());

    std::size_t mc_events { 0 };
    std::size_t detector_events { 0 };
    for (std::size_t n = 0; n < nr_events; ++n) {
        const double phi { (inEpsilon(theta))?0.:distro_phi(gen) };
        Line line { Line::generate( { distro_x(gen), distro_y(gen), simulation_plane_z_pos }, theta, phi) };
        bool coincidence { true };
        LineSegment refdet_path { setup.ref_detector()->intersection(line) };
        mc_events++;
        for (auto detector { setup.detectors().cbegin() };
             detector != setup.detectors().end();
             ++detector)
        {
            LineSegment det_path { detector->intersection(line) };
            if (det_path.length() < DEFAULT_EPSILON) {
                coincidence = false;
            }
        }
        if (coincidence) {
/*
            std::cout << "coincidence detected n="<< n << " " << std::setw(2) << toDeg(theta) << " " << toDeg(phi) << "\n";
*/
            detector_events++;
        }
    }
    double acceptance { static_cast<double>(detector_events) / mc_events };
    std::cout << "events simulated:"<<mc_events<<"  events detected:"<<detector_events<<" acceptance:" << static_cast<double>(detector_events) / mc_events << " acceptance error: " << std::sqrt(detector_events) / mc_events << "\n";
    
    dimensions = { (bounds.second - bounds.first) };
    const double simulation_area { 1e-6 * dimensions[0] * dimensions[1] };
    double effective_area { acceptance * simulation_area };
    std::cout << "effective area: " << effective_area << " +-" << std::sqrt(detector_events) / mc_events * simulation_area << " m^2\n";
    return effective_area;
}

template <int PHI_BINS = 256, int THETA_BINS = 256>
std::array<std::array<double, THETA_BINS>, PHI_BINS> theta_phi_scan(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_min, double theta_max, double phi_min, double phi_max)
{
    std::array<std::array<double, THETA_BINS>, PHI_BINS> phi_theta_acceptance {};

    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr << "no reference detector defined in DetectorSetup!\n";
        return phi_theta_acceptance;
    }

    const auto bounds { setup.ref_detector()->bounding_box() };
    std::uniform_real_distribution<> distro_x {
        bounds.first[0],
        bounds.second[0],
    };
    std::uniform_real_distribution<> distro_y {
        bounds.first[1],
        bounds.second[1],
    };
    std::uniform_real_distribution<> distro_z {
        bounds.first[2],
        bounds.second[2]
    };

    const double phi_step { (phi_max - phi_min) / (PHI_BINS - 1) };
    const double theta_step { (theta_max - theta_min) / (THETA_BINS - 1) };
    std::cout << "#phi theta acceptance acceptance_error\n";

    for (std::size_t phi_bin { 0 }; phi_bin < PHI_BINS; phi_bin++) {
        double phi { phi_min + phi_bin * phi_step };
        for (std::size_t theta_bin { 0 }; theta_bin < THETA_BINS; theta_bin++) {
            double theta { theta_min + theta_bin * theta_step };
            std::size_t mc_events { 0 };
            std::size_t detector_events { 0 };
            for (std::size_t n = 0; n < nr_events; ++n) {
                Line line { Line::generate(
                    { distro_x(gen), distro_y(gen), distro_z(gen) }, theta, phi) };

                bool coincidence { false };
                LineSegment refdet_path { setup.ref_detector()->intersection(line) };
                if (refdet_path.length() > DEFAULT_EPSILON) {
                    mc_events++;
                    coincidence = true;
                }
                for (auto detector { setup.detectors().cbegin() };
                     detector != setup.detectors().end();
                     ++detector) {
                    if (detector == setup.ref_detector())
                        continue;
                    LineSegment det_path { detector->intersection(line) };
                    if (det_path.length() < DEFAULT_EPSILON) {
                        coincidence = false;
                    }
                }
                if (coincidence) {
                    //std::cout << n << " " << std::setw(2) << toDeg(theta) << " " << toDeg(phi) << " " << det1_path.length() << " " << det2_path.length() << "\n";
                    detector_events++;
                }
            }
            double acceptance { static_cast<double>(detector_events) / mc_events };
            std::cout << toDeg(phi) << " " << toDeg(theta) << " " << acceptance << " " << std::sqrt(detector_events) / mc_events << "\n";
            phi_theta_acceptance[phi_bin][theta_bin] = acceptance;
            std::cout << std::flush;
        }
    }
    return phi_theta_acceptance;
}

/** @brief cosmic_simulation
 * This function does a full Monte-Carlo simulated evaluation of particle hits for a given detector setup. Tracks are distributed as follows:
 * - phi: uniformly distributed in -pi..pi
 * - theta: distributed in 0..pi/2 according to the standard model for muons: ~cos^n(theta) with n=2
 * @param setup DetectorSetup object containing the detector definitions
 * @param gen random number generator (Mersenne-Twister engine of type std::mt19937)
 * @param nr_events The total number of tracks to be generated. This is not identical to the number of MC events (=events having a hit in the reference detector)
 * @param theta_max the maximum zenith angle theta up to which the result histograms are scaled. Note, that the tracks are generated within the full theta range (according to the underlaying distribution) independent of this parameter.
 * @param nr_bins The number of bins in the resulting histograms
 * @return A std::vector of Histogram objects containing the resulting distributions for theta ("theta_distribution"), accepted theta ("accepted_theta_distribution"), phi ("phi_distribution") and ("accepted_phi_distribution").
 * @note In case of an error, the returned histogram vector is empty
 * @note The setup object must have the ref_detector iterator set to any valid detector it contains.
*/
std::vector<Histogram> cosmic_simulation(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_max, std::size_t nr_bins, int coinc_level = -1)
{
    std::vector<Histogram> histos {};
    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr << "no reference detector defined in DetectorSetup!\n";
        return histos;
    }

    if (coinc_level < 0) coinc_level = setup.detectors().size();

    Histogram theta_hist("theta_distribution", nr_bins, 0., theta_max);
    Histogram phi_hist("phi_distribution", nr_bins, -pi(), pi());
    Histogram theta_acc_hist("accepted_theta_distribution", nr_bins, 0., theta_max);
    Histogram phi_acc_hist("accepted_phi_distribution", nr_bins, -pi(), pi());

    const auto bounds { setup.ref_detector()->bounding_box() };
    std::uniform_real_distribution<> distro_x {
        bounds.first[0],
        bounds.second[0],
    };
    std::cout << "x-bounds: min=" << bounds.first[0] << " max=" << bounds.second[0] << "\n";

    std::uniform_real_distribution<> distro_y {
        bounds.first[1],
        bounds.second[1],
    };
    std::cout << "y-bounds: min=" << bounds.first[1] << " max=" << bounds.second[1] << "\n";

    std::uniform_real_distribution<> distro_z {
        bounds.first[2],
        bounds.second[2]
    };
    std::cout << "z-bounds: min=" << bounds.first[2] << " max=" << bounds.second[2] << "\n";

    std::uniform_real_distribution<> distro_phi(-pi(), pi());
    SampledDistribution distro_theta(cos2cdf, 0.0, pi() / 2, 65536);

    std::size_t mc_events { 0 };
    std::size_t detector_events { 0 };
    for (std::size_t n = 0; n < nr_events; ++n) {
        const double theta { distro_theta(gen) };
        const double phi { distro_phi(gen) };
        Line line { Line::generate(
            { distro_x(gen),
                distro_y(gen),
                /*bounds.first[2]*/
                distro_z(gen) },
            theta, phi) };

        unsigned int coincidence { 0 };
        LineSegment refdet_path { setup.ref_detector()->intersection(line) };
        if (refdet_path.length() > 0.) {
            theta_hist.fill(theta);
            phi_hist.fill(phi);
            mc_events++;
            if (mc_events % 100'000 == 0)
                std::cout << mc_events / 1000UL << "k MC events\n";
            coincidence++;
        }
        for (auto detector { setup.detectors().cbegin() };
             detector != setup.detectors().end();
             ++detector) {
            if (detector == setup.ref_detector())
                continue;
            LineSegment det_path { detector->intersection(line) };
            if (det_path.length() > 0.) {
                coincidence++;
            }
        }
        if (coincidence >= coinc_level) {
            //std::cout << n << " " << std::setw(2) << toDeg(theta) << " " << toDeg(phi) << " " << det1_path.length() << " " << det2_path.length() << "\n";
            theta_acc_hist.fill(theta);
            phi_acc_hist.fill(phi);
            detector_events++;
        }
    }
    std::cout << "MC events: " << mc_events << " detected events: " << detector_events << " acceptance: " << static_cast<double>(detector_events) / mc_events << " err(acceptance): " << std::sqrt(detector_events) / mc_events << "\n";
    histos.push_back(std::move(theta_hist));
    histos.push_back(std::move(phi_hist));
    histos.push_back(std::move(theta_acc_hist));
    histos.push_back(std::move(phi_acc_hist));
    return histos;
}

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
        { -126.5, -20. },
        { -91.5, -62.5 },
        { 0., -62.5 },
        { 0., 62.5 },
        { -91.5, 62.5 },
        { -126.5, 20. }
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
        {  500., -50. },
        {  500.,  50. },
        { -500.,  50. },
    };

    // create 3d objects of type ExtrudedObject defined by the 2d outline,
    // a global position offset and a thickness
    ExtrudedObject detector1 { large_paddle_points_lower, { 0., 0., 0. }, 8. };
    ExtrudedObject detector2 { large_paddle_points_upper, { 0., 0., 200. }, 8. };

    // create 3d objects of type ExtrudedObject but using the constructor for generation of a 
    // circular shape specified by a global position offset, radius, thickness and an optional
    // number of vertex points to generate the circle
    ExtrudedObject round_detector1 { { 0., 0., 0. }, 50., 10. };
    ExtrudedObject round_detector2 { { 0., 0., 100. }, 50., 10. };

    // construct a detector setup with the two detectors
    DetectorSetup setup { { detector1, detector2 } };

    // simulate the effective area (geometric aperture) at theta=0 of the detector system
    [[maybe_unused]] const double effective_area_sqm { simulate_geometric_aperture(setup, gen, nr_events) };

    // uncomment the following block to calculate the double differential acceptance
    // as function of phi and theta
/*
    [[maybe_unused]] const auto acceptance_phi_theta = theta_phi_scan<361, 46>(setup, gen, nr_events, 0., theta_max, -pi(), pi());
*/

    // initialize the histogram vector
    std::vector<Histogram> histos {};

    // run a scan over theta angle (uniformly distributed)
    // to record the detector acceptance, if required
    Append(histos, theta_scan(setup, gen, nr_events, 0., theta_max, nr_bins));

    // now, run the full simulation and append the resulting histograms
    // to the already existing histogram vector
    Append(histos,
        cosmic_simulation(setup, gen, nr_events * nr_bins, theta_max, nr_bins));

    // export each histogram into a separate file (human readable ASCII format)
    for (auto histo : histos) {
        histo.export_file(histo.getName() + ".hist");
    }

    exit(0);
}
