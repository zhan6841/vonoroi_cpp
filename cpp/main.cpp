#include "VonoroiAdapt.hpp"

#include <stdio.h>
#include <vector>

int main(int argc, char *argv[]) {
    printf("Test VoronoiBasic ...\n");

    std::vector<std::vector<float>> oxtsSet = {
        {-471.842, -830.177, 29.4617, 0.0280928, 0.00442842, 3.14157},
        {-500.038, -837.272, 29.4785, 0.00509157, -0.00424474, 0.52358},
        {-489.501, -835.234, 29.463, 0.00195096, 0.00334952, -2.61813}
    };

    VoronoiAdapt myVoronoiAdapt;
    PartitionDecision decision = myVoronoiAdapt.voronoi_basic(oxtsSet);

    for (size_t i = 0; i < decision.pb_set.size(); i++) {
        for (size_t j = 0; j < decision.pb_set[i].size(); j++) {
            std::vector<float> tempPb = decision.pb_set[i][j];
            printf("pbSet[%d,%d]: %f %f %f\n", (int)i, (int)j, tempPb[0], tempPb[1], tempPb[2]);
        }
    }
    // printf("neighbors: %d" + decision.neighbors);
    // printf("indexSet: %d" + decision.index_set);

    return 0;
}