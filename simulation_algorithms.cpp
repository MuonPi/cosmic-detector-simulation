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

constexpr double MEAN_ELOSS { 0.18 }; ///!< mean energy loss in MeV/mm
constexpr double ELOSS_WIDTH { 0.2 * MEAN_ELOSS }; ///!< width of e loss distribution


/** @brief The cumulative distribution function (CDF) for the cos^2(x) PDF distribution
 * This CDF is used for the calculation of the Probability Density Function (PDF)
 * of the cos^2(x) distribution for generating random values following the angular distribution
 * of the muon tracks
 */
auto cos2cdf = [](double x, const std::vector<double>& params = std::vector<double> {}) {
    //return cdf to following pdf: cos^2(x)
    return (2 / pi()) * (x / 2. + sin(2. * x) / 4.) + 0.5; //from Wolfram Alpha
};

/** @brief The cumulative distribution function (CDF) for the Moyal PDF distribution
 * as an approximation of the Landau-Gaus distribution describing the energy loss
 * of charged particles in finite absorbers.
 * This CDF is used for the calculation of the Probability Density Function (PDF)
 * of the energy loss distribution for generating random values of energy deposit
 * due to straggling.
 * Note: erfc(e^(-(x - μ)/(2 σ))/sqrt(2))≈erfc(exp(-(0.5 (x - μ))/σ))/sqrt(2))
 */
auto moyal_cdf = [](double x, const std::vector<double>& params) {
    assert(params.size() == 2);
    //return cdf to following pdf: Moayal distribution
    return std::erfc(std::exp(-0.5 * (x - params.at(0)) / params.at(1)) / c_sqrt2);
};

void theta_scan(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_min, double theta_max, std::size_t nr_bins, std::vector<Histogram>* histos)
{

    if (setup.ref_volume() == ExtrudedObject::invalid_volume()) {
        std::cerr << "no reference volume defined in DetectorSetup!\n";
        return;
    }

    Histogram acc_hist("acceptance_scan_theta",
        nr_bins,
        0., theta_max);

    const auto bounds { setup.ref_volume().bounding_box() };

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
            LineSegment refdet_path { setup.ref_volume().intersection(line) };
            if (refdet_path.length() > DEFAULT_EPSILON) {
                mc_events++;
                coincidence = true;
            }
            for (auto detector { setup.detectors().cbegin() };
                 detector != setup.detectors().end();
                 ++detector) {
                if (detector == setup.detectors().cbegin())
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

double simulate_geometric_aperture(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta, std::vector<Histogram>* histos)
{
    if (setup.ref_volume() == ExtrudedObject::invalid_volume()) {
        std::cerr << "no reference volume defined in DetectorSetup!\n";
        return std::numeric_limits<double>::quiet_NaN();
    }

    auto bounds { setup.ref_volume().bounding_box() };
    Point dimensions { bounds.second - bounds.first };
    //std::cout << "ref volume bounds: min=" << bounds.first << " max=" << bounds.second << "\n";
    //std::cout << "ref volume dimensions=" << dimensions << "\n";

    const std::size_t n_detectors { setup.detectors().size() };
    std::map<std::size_t, double> pathlength_values {};
    std::vector<Histogram> pathlength_histos {};
    std::vector<Histogram> eloss_histos {};
    std::size_t det_index { 0 };
    for (auto det : setup.detectors()) {
        auto bounds { det.bounding_box() };
        double max_pl = norm(bounds.second - bounds.first);
        //std::cout << "det " << std::to_string(det_index+1) << ": max dim=" << max_pl << std::endl;
        Histogram pl_hist("pathlength_fixtheta_det" + std::to_string(det_index + 1), 500U, 0., max_pl * 1.0);
        Histogram eloss_hist("eloss_fixtheta_det" + std::to_string(det_index + 1), 500U, 0., MEAN_ELOSS * max_pl * 2.0);
        pathlength_histos.emplace_back(std::move(pl_hist));
        eloss_histos.emplace_back(std::move(eloss_hist));
        det_index++;
    }

    std::uniform_real_distribution<> distro_x {
        bounds.first[0],
        bounds.second[0]
    };
    std::uniform_real_distribution<> distro_y {
        bounds.first[1],
        bounds.second[1]
    };
    const double simulation_plane_z_pos { bounds.first[2] };
    std::uniform_real_distribution<> distro_z {
        bounds.first[2],
        bounds.second[2]
    };
    std::uniform_real_distribution<> distro_phi(-pi(), pi());
    SampledDistribution distro_eloss(moyal_cdf, 0., 20., { MEAN_ELOSS, ELOSS_WIDTH });

    std::size_t mc_events { 0 };
    std::size_t detector_events { 0 };
    for (std::size_t n = 0; n < nr_events; ++n) {
        const double phi { (inEpsilon(theta)) ? 0. : distro_phi(gen) };
        double eloss { distro_eloss(gen) };
        Line line { Line::generate({ distro_x(gen), distro_y(gen), distro_z(gen) }, theta, phi) };
        unsigned int coincidence { 0 };
        LineSegment refdet_path { setup.ref_volume().intersection(line) };

        mc_events++;
        std::valarray<bool> hitvector(false, setup.detectors().size());
        std::size_t det_index { 0 };
        for (auto detector { setup.detectors().cbegin() };
             detector != setup.detectors().end();
             ++detector, ++det_index) {
            LineSegment det_path { detector->intersection(line) };
            if (det_path.length() > DEFAULT_EPSILON) {
                pathlength_values[det_index] = det_path.length();
                coincidence++;
                hitvector[det_index] = true;
            }
        }

        bool trigger { setup.isTrigger(hitvector) };
        if (trigger) {
            for (auto [detindex, pathlength_value] : pathlength_values) {
                pathlength_histos.at(detindex).fill(pathlength_value);
                eloss_histos.at(detindex).fill(pathlength_value * eloss);
                //std::cout << detindex << " " << std::setw(2) << toDeg(theta) << " " << toDeg(phi) << " " << pathlength_value << "\n";
            }
            detector_events++;
        }
    }

    if (histos != nullptr) {
        if (!pathlength_histos.empty()) {
            histos->insert(histos->end(), pathlength_histos.begin(), pathlength_histos.end());
        }
        if (!eloss_histos.empty()) {
            histos->insert(histos->end(), eloss_histos.begin(), eloss_histos.end());
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

    if (setup.ref_volume() == ExtrudedObject::invalid_volume()) {
        std::cerr << "no reference volume defined in DetectorSetup!\n";
        return phi_theta_acceptance;
    }

    const auto bounds { setup.ref_volume().bounding_box() };
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
                LineSegment refdet_path { setup.ref_volume().intersection(line) };
                if (refdet_path.length() > DEFAULT_EPSILON) {
                    mc_events++;
                    coincidence = true;
                }
                for (auto detector { setup.detectors().cbegin() };
                     detector != setup.detectors().end();
                     ++detector) {
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

DataItem<double> cosmic_simulation(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, std::vector<Histogram>* histos, std::size_t nr_bins, double theta_max)
{
    DataItem<double> data_item {};

    if (setup.ref_volume() == ExtrudedObject::invalid_volume()) {
        std::cerr << "no reference volume defined in DetectorSetup!\n";
        return DataItem<double> {};
    }

    const std::size_t n_detectors { setup.detectors().size() };
    std::map<std::size_t, double> pathlength_values {};
    std::vector<Histogram> pathlength_histos {};
    std::vector<Histogram> eloss_histos {};
    std::size_t det_index { 0 };
    for (auto det : setup.detectors()) {
        auto bounds { det.bounding_box() };
        double max_pl = norm(bounds.second - bounds.first);
        //std::cout << "det " << std::to_string(det_index+1) << ": max dim=" << max_pl << std::endl;
        Histogram pl_hist("pathlength_det" + std::to_string(det_index + 1), 500U, 0., max_pl * 1.0);
        Histogram eloss_hist("eloss_det" + std::to_string(det_index + 1), 500U, 0., MEAN_ELOSS * max_pl * 2.0);
        pathlength_histos.emplace_back(std::move(pl_hist));
        eloss_histos.emplace_back(std::move(eloss_hist));
        det_index++;
    }

    Histogram theta_hist("theta_distribution", nr_bins, 0., theta_max);
    Histogram phi_hist("phi_distribution", nr_bins, -pi(), pi());
    Histogram theta_acc_hist("accepted_theta_distribution", nr_bins, 0., theta_max);
    Histogram phi_acc_hist("accepted_phi_distribution", nr_bins, -pi(), pi());

    const auto bounds { setup.ref_volume().bounding_box() };
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
    SampledDistribution distro_theta(cos2cdf, 0., pi() / 2., std::vector<double> {});
    SampledDistribution distro_eloss(moyal_cdf, 0., 20., { MEAN_ELOSS, ELOSS_WIDTH });

    std::size_t detector_events { 0 };
    std::size_t coinc_events { 0 };
    for (std::size_t n = 0; n < nr_events; ++n) {
        const double theta { distro_theta(gen) };
        const double phi { distro_phi(gen) };
        double eloss { distro_eloss(gen) };
        //std::cout << "eloss=" << eloss << std::endl;
        Line line {
            Line::generate({ distro_x(gen), distro_y(gen), distro_z(gen) }, theta, phi)
        };

        theta_hist.fill(theta);
        phi_hist.fill(phi);

        if (n % 100'000 == 0 && n > 0)
            std::cout << n / 1000UL << "k MC events\n";

        det_index = 0;
        bool any_detector_hit { false };
        unsigned int coincidence { 0 };
        std::valarray<bool> hitvector(false, setup.detectors().size());
        for (auto detector { setup.detectors().cbegin() };
             detector != setup.detectors().end();
             ++detector, ++det_index) {
            LineSegment det_path { detector->intersection(line) };
            if (det_path.length() > 0.) {
                pathlength_values[det_index] = det_path.length();
                coincidence++;
                hitvector[det_index] = true;
                if (!any_detector_hit) {
                    any_detector_hit = true;
                    detector_events++;
                }
            }
        }
        bool trigger { setup.isTrigger(hitvector) };
        if (trigger) {
            theta_acc_hist.fill(theta);
            phi_acc_hist.fill(phi);
            for (auto [detindex, pathlength_value] : pathlength_values) {
                pathlength_histos.at(detindex).fill(pathlength_value);
                eloss_histos.at(detindex).fill(pathlength_value * eloss);
                //std::cout << detindex << " " << std::setw(2) << toDeg(theta) << " " << toDeg(phi) << " " << pathlength_value << "\n";
            }
            coinc_events++;
        }
    }
    std::cout << "MC events: " << nr_events << " detector hit events: " << detector_events << " coinc events: " << coinc_events << " acceptance: " << static_cast<double>(coinc_events) / nr_events << " err(acceptance): " << std::sqrt(coinc_events) / nr_events << "\n";
    if (histos != nullptr) {
        histos->push_back(std::move(theta_hist));
        histos->push_back(std::move(phi_hist));
        histos->push_back(std::move(theta_acc_hist));
        histos->push_back(std::move(phi_acc_hist));
        if (!pathlength_histos.empty()) {
            histos->insert(histos->end(), pathlength_histos.begin(), pathlength_histos.end());
        }
        if (!eloss_histos.empty()) {
            histos->insert(histos->end(), eloss_histos.begin(), eloss_histos.end());
        }
    }
    return { static_cast<double>(coinc_events) / nr_events, std::sqrt(coinc_events) / nr_events };
}

MeasurementVector<double, double> cosmic_simulation_detector_sweep(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, const Vector& detector_rotation_axis, double detector_min_angle, double detector_max_angle, std::size_t nr_angles)
{
    MeasurementVector<double, double> data_series {};
    auto rotated_setup { setup };
    const double dtheta { (detector_max_angle - detector_min_angle) / std::max<std::size_t>(1, (nr_angles - 1)) };
    std::cout << "min angle=" << detector_min_angle << ", dtheta=" << dtheta << "\n";
    double angle { detector_min_angle };
    rotated_setup.rotate(detector_rotation_axis, detector_min_angle);
    for (std::size_t i = 0; i < nr_angles; ++i) {
        std::cout << "current angle=" << toDeg(angle) << "deg\n";
        DataItem<double> item { cosmic_simulation(rotated_setup, gen, nr_events, nullptr, 90, toRad(90.)) };
        data_series.emplace_back(DataItem<double>({ angle, dtheta }), std::move(item));
        angle += dtheta;
        rotated_setup.rotate(detector_rotation_axis, dtheta);
        //        std::cout<<"rot matrix of rotated setup:\n"<<rotated_setup.ref_volume()->get_rotation_matrix();
    }
    return data_series;
}
