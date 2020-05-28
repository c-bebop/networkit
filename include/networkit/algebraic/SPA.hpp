#pragma once

#include <algorithm>
#include <vector>
#include <utility>

#include <networkit/Globals.hpp>

namespace NetworKit {

class SPAGala
{
public:
    SPAGala(size_t size)
    : b(size, false)
    {

    }

    ~SPAGala() = default;

    void accumulate(size_t const pos, double const value)
    {
        if (pos >= b.size())
        {
            return;
        }

        auto it = std::find(ls.begin(), ls.end(), pos);

        if (b[pos])
        {
            auto index = std::distance(ls.begin(), it);
            w[index] += value;
        } 
        else
        {
            b[pos] = true;

            if (it == std::end(ls))
            {
                ls.push_back(pos);
                w.push_back(value);
            }
            else
            {
                auto index = std::distance(ls.begin(), it);
                w[index] = value;
            }
        }
    }

    std::vector<std::pair<size_t, double>> output()
    {
        std::vector<std::pair<size_t, double>> out;

        for (size_t i = 0; i < ls.size(); ++i)
        {
            std::pair<size_t, double> entry;
            entry.first = ls[i];
            entry.second = w[i];

            out.emplace_back(std::move(entry));
        }

        return out;
    }

    size_t nnz() const 
    {
        return ls.size();
    }

    void reset()
    {
        std::for_each(std::begin(b), std::end(b), [](std::vector<bool>::reference x) { x = false; } );
        std::for_each(std::begin(w), std::end(w), [](double& x) { x = 0.; } ); // problem with output when resetting
        ls.clear();
    }

private:
    std::vector<bool> b;
    std::vector<size_t> ls;
    std::vector<double> w;
};

} /* namespace NetworKit */
