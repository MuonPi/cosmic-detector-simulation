#include "detector_setup.h"
#include "algebra_types.h"
#include "algebra_utils.h"
#include "geometry_types.h"

#include <iostream>

DetectorSetup::DetectorSetup(const std::vector<ExtrudedObject>& detectorlist)
    : m_detectors(detectorlist)
{
    if (!m_detectors.empty())
        m_ref_detector = m_detectors.begin();
}

DetectorSetup::DetectorSetup(const DetectorSetup& other)
    : m_detectors(other.m_detectors)
    , m_name(other.m_name)
{
    //!TODO: also copy the iterator pointing to the reference detector
    if (!m_detectors.empty())
        m_ref_detector = m_detectors.begin();
}

void DetectorSetup::rotate(const Vector& rot_axis, double rot_angle)
{
    if (inEpsilon(rot_angle))
        return;
    for (auto& detector : m_detectors) {
        Point pos { detector.position() };
        pos = ::rotate(pos, rot_axis, rot_angle);
        detector.set_position(pos);
        detector.add_rotation(rot_axis, rot_angle);
    }
}
