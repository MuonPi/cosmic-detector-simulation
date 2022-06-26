#ifndef HISTOGRAM_H
#define HISTOGRAM_H

#include <map>
#include <string>

class Histogram {
public:
    Histogram() = default;
    Histogram(const std::string& name, int nrBins, double min, double max) noexcept;
    ~Histogram();
    void clear();
    void setName(const std::string& name);
    void setUnit(const std::string& unit);
    void setNrBins(int bins);
    int getNrBins() const;
    void setMin(double val);
    double getMin() const;
    void setMax(double val);
    double getMax() const;
    double getRange() const;
    double getCenter() const;
    double getBinCenter(int bin) const;
    double getBinWidth() const;
    int getLowestOccupiedBin() const;
    int getHighestOccupiedBin() const;
    void fill(double x, double mult = 1.);
    void setBinContent(int bin, double value);
    double getBinContent(int bin) const;
    double getMean();
    double getMedian();
    double getMpv();
    double getRMS();
    double getUnderflow() const;
    double getOverflow() const;
    double getEntries() const;
    double getSum(double from, double to) const;
    void rescale(double center, double width);
    void rescale(double center);
    void rescale();
    void setAutoscale(bool autoscale = true);
    void export_file(const std::string& filename);

    const std::string& getName() const { return fName; }
    const std::string& getUnit() const { return fUnit; }

protected:
    int value2Bin(double value) const;
    double bin2Value(int bin) const;

    std::string fName { "defaultHisto" };
    std::string fUnit { "A.U." };
    int fNrBins { 100 };
    double fMin { 0.0 };
    double fMax { 1.0 };
    double fOverflow { 0 };
    double fUnderflow { 0 };
    std::map<int, double> fHistogramMap {};
    bool fAutoscale { false };
};

#endif // HISTOGRAM_H
