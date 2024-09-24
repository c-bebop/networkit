/*
 * BSuitorMatcher.hpp
 *
 *  Created on: 07.08.2023
 *      Author: Fabian Brandt-Tumescheit
 *              Frieda Gerharz
 */

#ifndef NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_

#include <algorithm>
#include <string_view>

#include <networkit/graph/Graph.hpp>
#include <networkit/matching/BMatcher.hpp>

namespace NetworKit {

/**
 * @ingroup matching
 * B-Suitor matching finding algorithm.
 */
class BSuitorMatcher : public BMatcher {
protected:
    struct MatchingNode {
        node id;
        edgeweight weight;

        MatchingNode() = default;
        MatchingNode(node n, edgeweight w) : id(n), weight(w) {}

        bool operator==(const MatchingNode &other) const {
            return id == other.id && weight == other.weight;
        }
        bool operator!=(const MatchingNode &other) const {
            return id != other.id || weight != other.weight;
        }
    };

    struct MatchingNodeInfo {
        std::vector<MatchingNode> partners;
        MatchingNode min; // (none, 0) if partners still has free capacity
        count max_size;

        MatchingNodeInfo() = default;

        MatchingNodeInfo(count b) {
            partners.reserve(b);
            min = MatchingNode(none, 0);
            max_size = b;
        }

        bool hasPartner(node u) {
            return std::find_if(partners.begin(), partners.end(),
                                [u](const MatchingNode &v) { return v.id == u; })
                   != partners.end();
        }

        MatchingNode popMinIfFull() {
            if (partners.size() < max_size) {
                return {none, 0};
            } else {
                MatchingNode min_copy = min;
                remove(min.id);
                return min_copy;
            }
        }

        void insert(const MatchingNode &u) {
            assert(partners.size() < max_size);
            partners.emplace_back(u);
            if (partners.size() == max_size && !partners.empty()) {
                min = *std::min_element(partners.begin(), partners.end(),
                                        [](const MatchingNode &x, const MatchingNode &y) {
                                            if (x.weight == y.weight) {
                                                return x.id > y.id;
                                            }
                                            return x.weight < y.weight;
                                        });
            }
        }

        void remove(node u) {
            partners.erase(std::remove_if(partners.begin(), partners.end(),
                                          [u](const MatchingNode &v) { return v.id == u; }),
                           partners.end());
            min = MatchingNode(none, 0);
        }

        void sort() {
            std::sort(partners.begin(), partners.end(),
                      [](const MatchingNode &u, const MatchingNode &v) {
                          return (u.weight > v.weight || (u.weight == v.weight && u.id < v.id));
                      });
        }
    };

public:
    /**
     * Computes a 1/2-approximate maximum weight b-matching of an undirected weighted Graph @c G
     * using the sequential b-Suitor algorithm published by Khan et al. in "Efficient Approximation
     * Algorithms For Weighted B-Matching", SIAM Journal on Scientific Computing, Vol. 38, Iss. 5
     * (2016).
     *
     * @param G An undirected graph.
     * @param b A vector of @a b values that represents the max number of edges per vertex @a v in
     * the b-Matching (b.at(v)).
     */
    BSuitorMatcher(const Graph &G, const std::vector<count> &b);

    /**
     * @param G An undirected graph.
     * @param b A value @a b that represents the max number of edges per vertex in the b-Matching.
     * Defaults to the ordinary 1-Matching.
     */
    BSuitorMatcher(const Graph &G, count b = 1);

    /**
     * @param G  An undirected graph.
     * @param path  A path to a file containing @a b values that represents the max number of edges
     * per vertex in the b-Matching.
     */
    BSuitorMatcher(const Graph &G, std::string_view &path);

    ~BSuitorMatcher() override = default;

    /**
     * Runs the algorithm.
     */
    void run() override;

protected:
    std::vector<std::unique_ptr<MatchingNodeInfo>> suitors;
    std::vector<std::unique_ptr<MatchingNodeInfo>> proposed;
    const std::vector<count> b;

    /**
     * Reads values from a file at @a path into the vector of b-values.
     *
     * @param size
     * @param path
     * @return std::vector<count>
     */
    std::vector<count> readBValuesFromFile(count size, std::string_view &path) const;

    /**
     * Iterates up to @a b times over the heaviest neighbors of node @a u and makes
     * them to suitors if eligible.
     *
     * @param u
     */
    void findSuitors(node u);

    /**
     * Finds the heaviest unmatched neighbor that @a u has not yet proposed to
     * if it exists. For equally weighted edges w(u, t), w(u, v) and t < v, w(u, t) is considered
     * smaller than w(u, v) to break ties.
     *
     * @param y
     * @return Node
     */
    MatchingNode findPreferred(node u);

    /**
     * Makes @a v a suitor of @a u and recursively calls itself for previous worse
     * suitors of @a u that got replaced with their new best match.
     *
     * @param u
     * @param w
     * @param v
     */
    void makeSuitor(node u, edgeweight w, node v);

    /**
     * Checks the symmetry of pairs of nodes. It must hold that v is in suitors(u) iff u is
     * in suitors(v).
     *
     */
    bool isSymmetrical() const;

private:

    /**
     * Creates the b-matching for given graph G. Function run() has to be called before calling
     * buildMatching. After invoking buildBMatching(), use getBMatching() to retrieve the resulting
     * b-matching.
     */
    void buildBMatching();

};
} // namespace NetworKit

#endif // NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
