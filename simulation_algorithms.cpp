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
#include "simulation_algorithms.h"
#include "utilities.h"

void theta_scan(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_min, double theta_max, std::size_t nr_bins, std::vector<Histogram>* histos)
{
    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr << "no reference detector defined in DetectorSetup!\n";
        return;
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
    std::cout << "theta step: " << theta_step << " rad = " << toDeg(theta_step) << " deg\n";
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
    if (histos != nullptr)
        histos->push_back(acc_hist);
}

double simulate_geometric_aperture(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta)
{
    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr << "no reference detector defined in DetectorSetup!\n";
        return {};
    }

    auto bounds { setup.ref_detector()->bounding_box() };
    Point dimensions { bounds.second - bounds.first };
    std::cout << "detector bounds: min=" << bounds.first << " max=" << bounds.second << "\n";
    std::cout << "detector dimensions=" << dimensions << "\n";

    //bounds.first -= dimensions * 5;
    //bounds.second += dimensions * 5;

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
        const double phi { (inEpsilon(theta)) ? 0. : distro_phi(gen) };
        Line line { Line::generate({ distro_x(gen), distro_y(gen), simulation_plane_z_pos }, theta, phi) };
        bool coincidence { true };
        LineSegment refdet_path { setup.ref_detector()->intersection(line) };
        mc_events++;
        for (auto detector { setup.detectors().cbegin() };
             detector != setup.detectors().end();
             ++detector) {
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
    std::cout << "events simulated:" << mc_events << "  events detected:" << detector_events << " acceptance:" << static_cast<double>(detector_events) / mc_events << " acceptance error: " << std::sqrt(detector_events) / mc_events << "\n";

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

DataItem<double> cosmic_simulation(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, std::vector<Histogram>* histos, std::size_t nr_bins, double theta_max, int coinc_level)
{
    DataItem<double> data_item {};
    if (setup.ref_detector() == setup.detectors().end()) {
        std::cerr << "no reference detector defined in DetectorSetup!\n";
        return DataItem<double> {};
    }

    if (coinc_level < 0)
        coinc_level = setup.detectors().size();

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
        Line line {
            Line::generate({ distro_x(gen), distro_y(gen), distro_z(gen) }, theta, phi)
        };

        theta_hist.fill(theta);
        phi_hist.fill(phi);

        unsigned int coincidence { 0 };
        LineSegment refdet_path { setup.ref_detector()->intersection(line) };
        if (refdet_path.length() > 0.) {
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
    if (histos != nullptr) {
        histos->push_back(std::move(theta_hist));
        histos->push_back(std::move(phi_hist));
        histos->push_back(std::move(theta_acc_hist));
        histos->push_back(std::move(phi_acc_hist));
    }
    return { static_cast<double>(detector_events) / mc_events, std::sqrt(detector_events) / mc_events };
}

MeasurementVector<double, double> cosmic_simulation_detector_sweep(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, const Vector& detector_rotation_axis, double detector_min_angle, double detector_max_angle, std::size_t nr_angles, int coinc_level)
{
    MeasurementVector<double, double> data_series {};
    //std::cout<<"rot matrix of orig. setup:\n"<<setup.ref_detector()->get_rotation_matrix();
    auto rotated_setup { setup };
    //std::cout<<"rot matrix of copied setup:\n"<<rotated_setup.ref_detector()->get_rotation_matrix();
    const double dtheta { (detector_max_angle - detector_min_angle) / std::max<std::size_t>(1, (nr_angles - 1)) };
    std::cout << "min angle=" << detector_min_angle << ", dtheta=" << dtheta << "\n";
    double angle { detector_min_angle };
    rotated_setup.rotate(detector_rotation_axis, detector_min_angle);
    for (std::size_t i = 0; i < nr_angles; ++i) {
        std::cout << "current angle=" << toDeg(angle) << "deg\n";
        DataItem<double> item { cosmic_simulation(rotated_setup, gen, nr_events, nullptr, 90, toRad(90.), coinc_level) };
        data_series.emplace_back(DataItem<double>({ angle, dtheta }), std::move(item));
        angle += dtheta;
        rotated_setup.rotate(detector_rotation_axis, dtheta);
        //std::cout<<"rot matrix of rotated setup:\n"<<rotated_setup.ref_detector()->get_rotation_matrix();
    }
    return data_series;
}
