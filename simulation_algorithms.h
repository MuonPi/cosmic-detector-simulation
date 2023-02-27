#pragma once

#include <fstream>
#include <random>
#include <vector>

#include "algebra_types.h"
#include "algebra_utils.h"
#include "detector_setup.h"
#include "geometry_types.h"
#include "histogram.h"
#include "sampled_distribution.h"
#include "utilities.h"

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
 * @param histos A pointer to std::vector of Histogram objects containing the resulting distributions. Currently, only the angular acceptance (detection efficiency vs. theta) is provided
 * @return nothing
 * @note In case of an error, the returned histogram vector is empty
 * @note The setup object must have the ref_detector iterator set to any valid detector it contains.
*/
void theta_scan(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_min, double theta_max, std::size_t nr_bins, std::vector<Histogram>* histos = nullptr);

/** @brief simulate_geometric_aperture
 * This function does Monte-Carlo simulated evaluation of the effective detector area.
 * @param setup DetectorSetup object containing the detector definitions
 * @param gen random number generator (Mersenne-Twister engine of type std::mt19937)
 * @param nr_events The number of tracks to be generated.
 * @param theta the theta angle of the generated tracks (default=0)
 * @return the absolute effective area in m^2
 * @note The generated tracks are distributed over an area at least 5 times larger
 * than the boundaries of the ref detector. The number of tracks intersecting all detectors
 * to the number of tracks generated is taken as a measure of the area provided that
 * the scale of the detector dimensions is in millimeters.
*/
double simulate_geometric_aperture(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta = 0.);

/** @brief theta_scan
 * This function ...
*/
template <int PHI_BINS = 256, int THETA_BINS = 256>
std::array<std::array<double, THETA_BINS>, PHI_BINS> theta_phi_scan(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, double theta_min, double theta_max, double phi_min, double phi_max);

/** @brief cosmic_simulation
 * This function does a full Monte-Carlo simulated evaluation of particle hits for a given detector setup. Tracks are distributed as follows:
 * - phi: uniformly distributed in -pi..pi
 * - theta: distributed in 0..pi/2 according to the standard model for muons: ~cos^n(theta) with n=2
 * @param setup DetectorSetup object containing the detector definitions
 * @param gen random number generator (Mersenne-Twister engine of type std::mt19937)
 * @param nr_events The total number of tracks to be generated. This is not identical to the number of MC events (=events having a hit in the reference detector)
 * @param histos A pointer to std::vector of Histogram objects containing the resulting distributions for theta ("theta_distribution"), accepted theta ("accepted_theta_distribution"), phi ("phi_distribution") and ("accepted_phi_distribution").
 * @param nr_bins The number of bins in the resulting histograms
 * @param theta_max the maximum zenith angle theta up to which the result histograms are scaled. Note, that the tracks are generated within the full theta range (according to the underlaying distribution) independent of this parameter.
 * @param coinc_level the number of detectors which have to be crossed by a track in order to assert a coincidence event
 * @return A DataItem<double> object containing the total acceptance (ratio of detected hits to generated tracks) and it's statistical error.
 * @note In case of an error, the histogram vector is not filled and the returned DataItem object is default initialized.
 * @note The setup object must have the ref_detector iterator set to any valid detector it contains.
*/
DataItem<double> cosmic_simulation(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, std::vector<Histogram>* histos = nullptr, std::size_t nr_bins = 90, double theta_max = pi() / 2, int coinc_level = -1);

/** @brief cosmic_simulation_detector_sweep
 * This function does a full Monte-Carlo simulated evaluation of particle hits for a given detector setup for a scan over a given range of detector rotations around a given axis. Tracks are distributed as follows:
 * - phi: uniformly distributed in -pi..pi
 * - theta: distributed in 0..pi/2 according to the standard model for muons: ~cos^n(theta) with n=2
 * @param setup DetectorSetup object containing the detector definitions
 * @param gen random number generator (Mersenne-Twister engine of type std::mt19937)
 * @param nr_events The total number of tracks to be generated. This is not identical to the number of MC events (=events having a hit in the reference detector)
 * @param detector_rotation_axis A vector specifying the alignment of the rotation axis
 * @param detector_min_angle the smallest rotation angle
 * @param detector_min_angle the largest rotation angle
 * @param nr_angles the number of steps in which the detector will be rotated between detector_min_angle and detector_max_angle
 * @param coinc_level the number of detectors which have to be crossed by a track in order to assert a coincidence event
 * @return A MeasurementVector<double> (aka std::vector<std::pair<DataItem<double>, DataItem<double>>) object containing the total acceptance (ratio of detected hits to generated tracks) and it's statistical error.
 * @note In case of an error, the histogram vector is not filled and the returned MeasurementVector object is default initialized and empty.
 * @note The setup object must have the ref_detector iterator set to any valid detector it contains.
*/
MeasurementVector<double, double> cosmic_simulation_detector_sweep(const DetectorSetup& setup, std::mt19937& gen, std::size_t nr_events, const Vector& detector_rotation_axis, double detector_min_angle, double detector_max_angle, std::size_t nr_angles, int coinc_level = -1);
