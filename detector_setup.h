#pragma once

#include <vector>

#include "geometry_types.h"

/** @brief DetectorSetup - class for managing geometric objects of type ExtrudedObject
 * This class stores a list of ExtrudedObject detector objects for convenience.
 * An arbitrary number of detector objects with individual alignments can be supplied to the
 * constructor. Additionaly, a global rotation of the whole setup can be set.
 * @note The first detector of the vector which is supplied to the constructor is set as
 * reference detector. The assignment of the ref detector can be changed with the
 * DetectorSetup::ref_detector() method
*/
class DetectorSetup {
public:
    DetectorSetup() = delete;
    DetectorSetup(const std::vector<ExtrudedObject>& detectorlist);

    auto detectors() -> std::vector<ExtrudedObject>& { return m_detectors; }
    auto detectors() const -> const std::vector<ExtrudedObject>& { return m_detectors; }
    std::vector<ExtrudedObject>::iterator& ref_detector() { return m_ref_detector; }
    const std::vector<ExtrudedObject>::iterator& ref_detector() const { return m_ref_detector; }
    void rotate(const Vector& rot_axis, double rot_angle);

private:
    std::vector<ExtrudedObject> m_detectors {};
    std::vector<ExtrudedObject>::iterator m_ref_detector {};
    std::string m_name {};
};
