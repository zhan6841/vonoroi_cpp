// Created by Anlan Zhang, 10/26/2022
#ifndef __VORONOI_ADAPT_HPP__
#define __VORONOI_ADAPT_HPP__

#include <vector>

class PartitionDecision {
public:
    std::vector<std::vector<std::vector<float>>> pb_set;
    std::vector<std::vector<int>> neighbors;
    std::vector<int> index_set;

    PartitionDecision(int oxts_length);
    ~PartitionDecision();
};

class VoronoiAdapt {
public:
    const int dimension = 4;

    VoronoiAdapt();
    ~VoronoiAdapt();

    PartitionDecision voronoi_basic(std::vector<std::vector<float>> &oxtsSet);
};

#endif