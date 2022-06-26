#include <chrono>
#include <cmath>
#include <functional>
#include <future>
#include <iostream>
#include <string>
#include <thread>
#include <valarray>
#include <vector>
#include <random>
#include <iomanip>
#include <fstream>
#include <numeric>

#include "sampled_distribution.h"
#include "histogram.h"
#include "algebra_utils.h"
#include "algebra_types.h"
#include "geometry_types.h"
#include "detector_setup.h"


constexpr int g_verbosity{3};
constexpr double toDeg(double x) { return x*180/pi(); }
constexpr double toRad(double x) { return x*pi()/180; }

/** @brief The cumulative distribution function (CDF) for the cos^2(x) PDF distribution
 * This CDF is used for the calculation of the Probability Density Function (PDF)
 * of the cos^2(x) distribution for generating random values following the angular distribution
 * of the muon tracks
 */ 
auto cos2cdf = [](double x) { 
    //return cdf to following pdf: cos^2(x)
    return (2/pi())*(x/2. + sin(2.*x)/4.) + 0.5; //from Wolfram Alpha
};

template <typename T>
void Append(std::vector<T>& a, const std::vector<T>& b)
{
    a.reserve(a.size() + b.size());
    a.insert(a.end(), b.begin(), b.end());
}


std::vector<Histogram> theta_scan(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_min, double theta_max, std::size_t nr_bins)
{
    std::vector<Histogram> histos {};

    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr<<"no reference detector defined in DetectorSetup!\n";
        return histos;
    }

    Histogram acc_hist("acceptance_scan_theta",
                       nr_bins,
                       0., theta_max);

    std::uniform_real_distribution<> distro_x{
        setup.ref_detector()->bounding_box().first[0],
        setup.ref_detector()->bounding_box().second[0],
    };
    std::uniform_real_distribution<> distro_y{
        setup.ref_detector()->bounding_box().first[1],
        setup.ref_detector()->bounding_box().second[1],
    };
    std::uniform_real_distribution<> distro_z{
        setup.ref_detector()->bounding_box().first[2],
        setup.ref_detector()->bounding_box().second[2]
    };
    std::uniform_real_distribution<> distro_phi(-pi(), pi());

    const double theta_step { (theta_max-theta_min)/(nr_bins-1) };
    std::cout << "#theta acceptance acceptance_error\n"; 
    for ( double theta { theta_min }; theta < theta_max; theta+=theta_step) {
        std::size_t mc_events {0};
        std::size_t detector_events {0};
        for (std::size_t n = 0; n < nr_events; ++n) {
            double phi { distro_phi(gen) };
            Line line { Line::generate(
                { distro_x(gen), distro_y(gen), distro_z(gen) }, theta, phi) 
            };

            bool coincidence { false };
            LineSegment refdet_path { setup.ref_detector()->intersection(line) };
            if (refdet_path.length() > DEFAULT_EPSILON) { 
                mc_events++;
                coincidence = true;
            }
            for ( auto detector { setup.detectors().cbegin() };
                 detector != setup.detectors().end();
                ++detector)
            {
                if (detector == setup.ref_detector()) continue;
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
        std::cout<<std::flush;
    }
    histos.push_back(acc_hist);
    return histos;
}

/** @brief cosmic_simulation
 * This function does a full Monte-Carlo simulated evaluation of particle hits for a given detector setup. Tracks are distributed as follows:
 * - phi: uniformly distributed in -pi..pi
 * - theta: distributed in 0..pi/2 according to the standard model for muons: ~cos^n(theta) with n=2
 * @arg setup DetectorSetup object containing the detector definitions
 * @arg gen random number generator (Mersenne-Twister engine of type std::mt19937)
 * @arg nr_events The total number of tracks to be generated. This is not identical to the number of MC events (=events having a hit in the reference detector)
 * @arg theta_max the maximum zenith angle theta up to which the result histograms are scaled. Note, that the tracks are generated within the full theta range (according to the underlaying distribution) independent of this parameter.
 * @arg nr_bins The number of bins in the resulting histograms
 * @return A std::vector of Histogram objects containing the resulting distributions for theta ("theta_distribution"), accepted theta ("accepted_theta_distribution"), phi ("phi_distribution") and ("accepted_phi_distribution").
 * @note In case of an error, the returned histogram vector is empty
 * @note The setup object must have the ref_detector iterator set to any valid detector it contains.
*/
std::vector<Histogram> cosmic_simulation(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_max, std::size_t nr_bins)
{
    std::vector<Histogram> histos {};
    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr<<"no reference detector defined in DetectorSetup!\n";
        return histos;
    }

    Histogram theta_hist("theta_distribution",nr_bins, 0., theta_max);
    Histogram phi_hist("phi_distribution",nr_bins, -pi(), pi());
    Histogram theta_acc_hist("accepted_theta_distribution",nr_bins, 0., theta_max);
    Histogram phi_acc_hist("accepted_phi_distribution",nr_bins, -pi(), pi());
    
    std::uniform_real_distribution<> distro_x{
        setup.ref_detector()->bounding_box().first[0],
        setup.ref_detector()->bounding_box().second[0],
    };
    std::cout<<"x-bounds: min="<<setup.ref_detector()->bounding_box().first[0]<<" max="<<setup.ref_detector()->bounding_box().second[0]<<"\n";

    std::uniform_real_distribution<> distro_y{
        setup.ref_detector()->bounding_box().first[1],
        setup.ref_detector()->bounding_box().second[1],
    };
    std::cout<<"y-bounds: min="<<setup.ref_detector()->bounding_box().first[1]<<" max="<<setup.ref_detector()->bounding_box().second[1]<<"\n";

    std::uniform_real_distribution<> distro_z{
        setup.ref_detector()->bounding_box().first[2],
        setup.ref_detector()->bounding_box().second[2]
    };
    std::cout<<"z-bounds: min="<<setup.ref_detector()->bounding_box().first[2]<<" max="<<setup.ref_detector()->bounding_box().second[2]<<"\n";

    std::uniform_real_distribution<> distro_phi(-pi(), pi());
    SampledDistribution distro_theta(cos2cdf, 0.0, pi()/2, 65536);

    std::size_t mc_events {0};
    std::size_t detector_events {0};
    for (std::size_t n = 0; n < nr_events; ++n) {
            double theta { distro_theta(gen) };
            double phi { distro_phi(gen) };
            Line line { Line::generate(
                { distro_x(gen), distro_y(gen), distro_z(gen) }, theta, phi) 
            };
                        
            bool coincidence { false };
            LineSegment refdet_path { setup.ref_detector()->intersection(line) };
            if (refdet_path.length() > DEFAULT_EPSILON) { 
                theta_hist.fill(theta);
                phi_hist.fill(phi);
                mc_events++;
                if (mc_events%100000 == 0) std::cout<<mc_events<<" MC events\n";
                coincidence = true;
            }
            for ( auto detector { setup.detectors().cbegin() };
                 detector != setup.detectors().end();
                ++detector)
            {
                if (detector == setup.ref_detector()) continue;
                LineSegment det_path { detector->intersection(line) };
                if (det_path.length() < DEFAULT_EPSILON) {
                    coincidence = false;
                }
            }
            if (coincidence) {
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

auto main() -> int {
    std::cout << "detector MC simulator\n";
    
    // set up random number generator
    std::random_device rd;  // Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

    // the following parameters should be adopted to the individual situation
    constexpr std::size_t nr_events { 100'000 }; //<! the total number of tracks to be simulated
    constexpr double theta_max { toRad(62.) }; //<! the maximum theta angle taken into account
    constexpr double theta_step { toRad(1.) }; //<! the desired granularity of the simulated angular distributions

    constexpr std::size_t nr_bins { static_cast<int>(theta_max/theta_step)+1 };
    std::cout<<"nr of bins: "<<nr_bins<<"\n";

    // vector of 2d polygon vertices defining the shape (in x-y-plane) of the detector.
    // note, that the points have to be in geometrical sequential order in
    // counter-clockwise orientation (i.e. a sequence of points defining the detector outline
    // going in ccw direction). The first point is considered to close the polygon together
    // with the last element in the vector

    // definition of the large double-paddle detector in IIPI-JLU lab

    const std::vector<Point> large_paddle_points_upper{
        {-150., -87.5},
        {150., -87.5},
        {150., 87.5},
        {-150., 87.5}
    };
    const std::vector<Point> large_paddle_points_lower{
        {-150., -100.},
        {150., -100.},
        {150., 100.},
        {-150., 100.}
    };

    // definition of the MuonPi standard-size (octagon) detector
    const std::vector<Point> octagon_points{
        {-126.5, -20.},
        {-91.5, -62.5},
        {91.5, -62.5},
        {126.5, -20.},
        {126.5, 20.},
        {91.5, 62.5},
        {-91.5, 62.5},
        {-126.5, 20.}
    };

    // definition of the MuonPi half-size detector
    const std::vector<Point> half_size_detector_points{
        {-126.5, -20.},
        {-91.5, -62.5},
        {0., -62.5},
        {0., 62.5},
        {-91.5, 62.5},
        {-126.5, 20.}
    };

    // create 3d objects of type ExtrudedObject defined by the 2d outline,
    // a global position offset and a thickness
    ExtrudedObject detector1{large_paddle_points_lower, {0.,0.,0.}, 9.5};
    ExtrudedObject detector2{large_paddle_points_upper, {0.,0.,200.}, 8.};

    // construct a detector setup with the two detectors
    DetectorSetup setup { { detector1, detector2 } };

    // first, run a scan over theta angle (uniformly distributed)
    // to record the detector acceptance
    std::vector<Histogram> histos { theta_scan(setup, gen, nr_events, 0., theta_max, nr_bins) };

    // now, run the full simulation and append the resulting histograms
    // to the already existing histogram vector
    Append(histos,
           cosmic_simulation(setup, gen, nr_events*nr_bins, theta_max, nr_bins)
    );

    // export each histogram into a separate file (human readable ASCII format)
    for ( auto histo: histos ) {
        histo.export_file(histo.getName()+".hist");
    }

    exit(0);
}
