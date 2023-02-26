#include <vector>
#include <cmath>
#include "utilities.h"

void export_file(const MeasurementVector<double,double>& data, const std::string& filename)
{
    std::ofstream ofs(filename, std::ofstream::out);
    ofs << "# Data series\n";
    ofs << "# " << data.size() << " entries\n";
    ofs << "# x err(x) y err(y)\n";
    for (const auto& [xval,yval]: data) {
        ofs << xval.value << " " << xval.error << " " << yval.value << " " << yval.error << "\n";
    }
    ofs.close();
}
