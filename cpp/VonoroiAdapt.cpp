#include "VonoroiAdapt.hpp"
#include "Vonoroi.hpp"

#include <stdio.h>
#include <math.h>

// PartitionDecision
PartitionDecision::PartitionDecision(int oxts_length) {
    this->pb_set = std::vector<std::vector<std::vector<float>>> (oxts_length, std::vector<std::vector<float>>());
    this->neighbors = std::vector<std::vector<int>> (oxts_length, std::vector<int>());
}

PartitionDecision::~PartitionDecision() {}

// static functions
static std::vector<float> getMinMaxOxtsSetXY(std::vector<std::vector<float>> &oxtsSet) {
    std::vector<float> result = {0.0f, 0.0f, 0.0f, 0.0f};

    if (oxtsSet.size() < 1) {
        return result;
    }

    std::vector<float> oxts0 = oxtsSet[0];
    float minX = oxts0[0];
    float maxX = oxts0[0];
    float minY = oxts0[1];
    float maxY = oxts0[1];

    for (size_t i = 1; i < oxtsSet.size(); i++) {
        std::vector<float> tmpOxts = oxtsSet[i];
        minX = (tmpOxts[0] < minX) ? tmpOxts[0] : minX;
        maxX = (tmpOxts[0] > maxX) ? tmpOxts[0] : maxX;
        minY = (tmpOxts[1] < minY) ? tmpOxts[1] : minY;
        maxY = (tmpOxts[1] > maxY) ? tmpOxts[1] : maxY;
    }

    result[0] = minX;
    result[1] = maxX;
    result[2] = minY;
    result[3] = maxY;

    return result;
}

static std::vector<float> transformation(std::vector<float> &oxts1, std::vector<float> &oxts2) {
    // transformation matrix - translation (to the perspective of oxts1)
    float da = oxts2[0] - oxts1[0];  // south --> north
    float db = oxts2[1] - oxts1[1];  // east --> west
    float dx = da * (float)std::cos(oxts1[5]) + db * (float)std::sin(oxts1[5]);
    float dy = da * (-(float)std::sin(oxts1[5])) + db * (float)std::cos(oxts1[5]);
    // float dz = oxts2[2] - oxts1[2];
    // float[] translation = {dx, dy, dz, 0.0f};
    std::vector<float> translation = {dx, dy};

    return translation;
}

static std::vector<float> calculatePb(std::vector<float> &p1, std::vector<float> &p2, bool normalize){
        float x1 = p1[0];
        float y1 = p1[1];
        float x2 = p2[0];
        float y2 = p2[1];

        std::vector<float> result = {-65536.0f, -65536.0f, -65536.0f};
        float a = -65536.0f;
        float b = -65536.0f;
        float c = -65536.0f;

        if (x1 == x2 && y1 == y2) {
            fprintf(stderr, "calculatePb: Two same coordinates! (%f,%f) (%f,%f)\n", x1, y1, x2, y2);
            return result;
        }
        else if (x1 == x2) {
            a = 0.0f;
            b = 1.0f;
        }
        else if (y1 == y2) {
            a = 1.0f;
            b = 0.0f;
        }
        else {
            a = x2 - x1;
            b = y2 - y1;
        }

        if (normalize == true) {
            float r = (float)std::sqrt((double)a * (double)a + (double)b * (double)b);
            a = a / r;
            b = b / r;
        }
        c = -(x1+x2)/2*a - (y1+y2)/2*b;

        if (x1*a + y1*b + c < 0.0f) {
            a = -a;
            b = -b;
            c = -c;
        }

        result[0] = a;
        result[1] = b;
        result[2] = c;
        return result;
    }

// VoronoiAdapt
VoronoiAdapt::VoronoiAdapt() {}

VoronoiAdapt::~VoronoiAdapt() {}

PartitionDecision VoronoiAdapt::voronoi_basic(std::vector<std::vector<float>> &oxtsSet) {
    PartitionDecision decision(oxtsSet.size());
    std::vector<std::vector<float>> vehSet;

    for (size_t i = 0; i < oxtsSet.size(); i++) {
        decision.index_set.push_back(i);
        vehSet.push_back(oxtsSet[i]);
    }
    int numVehs = decision.index_set.size();

    if (numVehs == 1) {
        return decision;
    } else if (numVehs == 2) {
        decision.neighbors[0].push_back(1);
        decision.neighbors[1].push_back(0);
    } else {
        // Apply Voronoi; Prepare Data
        double *xValuesIn = new double[oxtsSet.size()];
        double *yValuesIn = new double[oxtsSet.size()];
        for (size_t i = 0; i < oxtsSet.size(); i++) {
            std::vector<float> tempOxts = oxtsSet[i];
            xValuesIn[i] = (double)tempOxts[0];
            yValuesIn[i] = (double)tempOxts[1];
        }
        std::vector<float> minMaxOxtsSetXY = getMinMaxOxtsSetXY(oxtsSet);
        double minX = (double)minMaxOxtsSetXY[0] - 1.0;
        double maxX = (double)minMaxOxtsSetXY[1] + 1.0;
        double minY = (double)minMaxOxtsSetXY[2] - 1.0;
        double maxY = (double)minMaxOxtsSetXY[3] + 1.0;

        // get all neighbor pairs
        Voronoi myVoronoi(0.0);
        myVoronoi.generate_voronoi(xValuesIn, yValuesIn, oxtsSet.size(), minX, maxX, minY, maxY);
        std::vector<std::vector<int>> allNeighborPairs = myVoronoi.get_all_neighbor_pairs();
        for (size_t i = 0; i < allNeighborPairs.size(); i++) {
            std::vector<int> pair = allNeighborPairs[i];
            decision.neighbors[pair[0]].push_back(pair[1]);
            decision.neighbors[pair[1]].push_back(pair[0]);
        }
        delete [] xValuesIn;
        delete [] yValuesIn;
    }

    std::vector<std::vector<std::vector<float>>> vehs(numVehs, std::vector<std::vector<float>>());

    for (int i = 0; i < numVehs; i++) {
        for (size_t tmp = 0; tmp < decision.neighbors[i].size(); tmp++) {
            int j = decision.neighbors[i][tmp];
            vehs[i].push_back(transformation(vehSet[i], vehSet[j]));
        }
        std::vector<float> zerors = {0.0f, 0.0f};
        for (size_t x = 0; x < decision.neighbors[i].size(); x++) {
            decision.pb_set[i].push_back(calculatePb(zerors, vehs[i][x], true));
        }
    }
    return decision;
}