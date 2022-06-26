#pragma once

#include <vector>

#include "geometry_types.h"

class DetectorSetup {
public:
    DetectorSetup() = delete;
    DetectorSetup(const std::vector<ExtrudedObject>& detectorlist)
    : m_detectors(detectorlist)
    {
        if (!m_detectors.empty()) m_ref_detector = m_detectors.begin();
    }
    auto detectors() -> std::vector<ExtrudedObject>& { return m_detectors; }
    auto detectors() const -> const std::vector<ExtrudedObject>& { return m_detectors; }
    std::vector<ExtrudedObject>::iterator& ref_detector() { return m_ref_detector; }
    const std::vector<ExtrudedObject>::iterator& ref_detector() const { return m_ref_detector; }
private:
    std::vector<ExtrudedObject> m_detectors {};
    std::vector<ExtrudedObject>::iterator m_ref_detector {};
    std::string m_name {};
};
