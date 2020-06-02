#pragma once

#include <algorithm>
#include <vector>
#include <utility>
#include <list>
#include <stdexcept>
#include <string>
#include <iostream>

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

    void accumulate(double const value, size_t const pos)
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

    size_t output(std::vector<double>& values, std::vector<index>& columnIdx, size_t const nzcur)
    {
        size_t non_zero_count = 0;

        for (auto& e : non_zeros_ls) 
        {
            size_t const next_i = nzcur + non_zero_count; 

            if (next_i >= columnIdx.size())
            {
                columnIdx.push_back(e);
            }
            else
            {
                columnIdx[next_i] = e;
            }

            if (next_i >= values.size())
            {
                values.push_back(values_w[e]);
            }
            else
            {
                values[next_i] = values_w[e];
            }

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

    void print()
    {
        std::cout << "values_w: ";

        for (auto e : values_w)
        {
            std::cout << "[" << e << "] ";
        }

        std::cout << std::endl;
    }

private:
    std::vector<bool> occupied_b;
    std::vector<double> values_w;
    std::vector<size_t> non_zeros_ls;
};

} /* namespace NetworKit */
