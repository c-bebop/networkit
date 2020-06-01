#pragma once

#include <algorithm>
#include <vector>
#include <utility>
#include <list>
#include <stdexcept>
#include <string>

#include <networkit/Globals.hpp>

namespace NetworKit {

class SPA
{
public:
    SPA(size_t size)
    : occupied_b(size, false)
    , values_w(size, 0.)
    {

    }

    ~SPA() = default;

    void accumulate(size_t const pos, double const value)
    {
        if (pos >= occupied_b.size())
        {
            throw std::out_of_range("Position " + std::to_string(pos) + " does not exist!");
        }

        if (occupied_b[pos])
        {
            values_w[pos] += value;
        } 
        else
        {
            occupied_b[pos] = true;
            values_w[pos] = value;
            non_zeros_ls.push_back(pos);
        }
    }

    size_t output(std::vector<index>& columnIdx, std::vector<double>& values, size_t const nzcur)
    {
        size_t non_zero_count = 0;

        for (auto& e : non_zeros_ls) 
        {
            columnIdx[nzcur + non_zero_count] = e;
            values[nzcur + non_zero_count] = values_w[e];
            ++non_zero_count;
        }

        return non_zero_count;
    }

    size_t nnz() const 
    {
        return non_zeros_ls.size();
    }

    void reset()
    {
        std::for_each(std::begin(non_zeros_ls), std::end(non_zeros_ls), [this](size_t i)
        { 
            occupied_b[i] = false; 
        });

        non_zeros_ls.clear();
    }

private:
    std::vector<bool> occupied_b;
    std::vector<double> values_w;
    std::vector<size_t> non_zeros_ls;
};

} /* namespace NetworKit */
