#pragma once

#include <vector>

#include <networkit/Globals.hpp>

namespace NetworKit {

class SPA
{
public:
    SPA(size_t size);

    ~SPA() = default;

    void accumulate(double const value, size_t const pos);

    size_t output(std::vector<double>& values, std::vector<index>& columnIdx, size_t const nzcur);

    size_t nnz() const;

    void reset();

    void print();

private:
    std::vector<bool> occupied_b;
    std::vector<double> values_w;
    std::vector<size_t> non_zeros_ls;
};

} /* namespace NetworKit */
