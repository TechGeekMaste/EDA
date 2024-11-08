#include <cmath>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <iostream>
#include "problem.h"
#include "bobyqa.h"

// 模拟退火参数结构体
struct SimulatedAnnealingParams {
    double initial_temp = 1000.0;  // 初始温度
    double final_temp = 1.0;       // 结束温度
    double cooling_rate = 0.95;    // 降温率
    int max_iterations = 100;      // 每轮温度的最大迭代次数
    double search_radius = 10.0;   // 搜索范围半径
};
// 计算候选缓冲区位置的成本
double calculateCost(const ProblemData::Coordinate& candidate_location,
                     const std::vector<int>& layer1_indices,
                     const std::vector<Cluster>& Clusters_Layer1) 
{
    double max_diff = 0;
    std:: vector<double> maxDiss(layer1_indices.size(), 0);
    std::vector<double> manhattanDistance_S(layer1_indices.size(), 0);

    for (size_t i = 0; i < layer1_indices.size(); ++i) {
        int layer1_idx = layer1_indices[i];
        const Cluster& cluster1 = Clusters_Layer1[layer1_idx];

        double max_Distance = 0;
        for (const auto &ff : cluster1.ffs) {
            double manhattanDistance = std::pow(std::abs(ff.ff_x - cluster1.GetCX()) + std::abs(ff.ff_y - cluster1.GetCY()), 2);
            max_Distance = std::max(manhattanDistance , max_Distance);
        }
        maxDiss[i] = max_Distance ;

        manhattanDistance_S[i] = std::pow(std::abs(candidate_location.x - cluster1.GetCX()) + 
                                              std::abs(candidate_location.y - cluster1.GetCY()), 2);

    }
    //计算目标成本
    for (size_t i = 0; i < layer1_indices.size(); ++i){
      for (size_t j = 0; j < layer1_indices.size(); ++j){
        if (i != j){
          double diff = std::abs((maxDiss[i] + manhattanDistance_S[i]) - (maxDiss[j] + manhattanDistance_S[j]));
          max_diff = std::max(max_diff, diff);
        }
      }
    }
    return max_diff;
}

#ifdef BOBYQA_IMPL
struct BOBYQA_DATA {
  const std::vector<int>& var1;
  const std::vector<Cluster>& var2;
};

REAL bobyqa_obj(const INTEGER n, const REAL* x, void* data) {
  BOBYQA_DATA* params = (BOBYQA_DATA*)data;
  const std::vector<int>& var1 = params->var1;
  const std::vector<Cluster>& var2 = params->var2;
  return calculateCost(ProblemData::Coordinate{x[0], x[1]} , var1, var2);
}
#endif

// 模拟退火算法优化缓冲区位置
ProblemData::Coordinate optimizeBufferLocation_SA(
    const Cluster& cluster2,
    const std::vector<int>& layer1_indices,
    const std::vector<Cluster>& Clusters_Layer1,
    const SimulatedAnnealingParams& params)
{
  static double cost_sum = 0;
    ProblemData::Coordinate best_location = {cluster2.GetCX(), cluster2.GetCY()};
    double min_cost = calculateCost(best_location, layer1_indices, Clusters_Layer1);

    // 计算 cluster2 到其对应的第一层 buffer 的最远曼哈顿距离
    double max_manhattan_distance_x = 0.0;
    double max_manhattan_distance_y = 0;
    for (int layer1_idx : layer1_indices)
    {
        const Cluster& cluster1 = Clusters_Layer1[layer1_idx];
        double manhattan_distance_x = std::abs(cluster1.GetCX() - cluster2.GetCX()) + std::abs(cluster1.GetCY() - cluster2.GetCY() );
        double manhattan_distance_y = std::abs(cluster1.GetCY() - cluster2.GetCY()) + std::abs(cluster1.GetCX() - cluster2.GetCX() );
        max_manhattan_distance_x = std::max(max_manhattan_distance_x, manhattan_distance_x);
        max_manhattan_distance_y = std::max(max_manhattan_distance_y, manhattan_distance_y);

    }
    double search_radius_x = max_manhattan_distance_x ; // 更新 search_radius
    double search_radius_y = max_manhattan_distance_y ;
    std::cout << "Updated search_radius: " << search_radius_x << "  "  << search_radius_y << std::endl;

#ifndef BOBYQA_IMPL
    double temperature = params.initial_temp;

    srand(time(0)); // 初始化随机数种子
    while (temperature > params.final_temp) {
        for (int i = 0; i < params.max_iterations; ++i) {
            // 生成候选位置（在当前缓冲区位置附近随机扰动）
            ProblemData::Coordinate candidate_location = {
                best_location.x + (rand() % int(search_radius_x * 2 + 1) - search_radius_x),
                best_location.y + (rand() % int(search_radius_y * 2 + 1) - search_radius_y)
            };
            candidate_location.x = min(max(candidate_location.x, cluster2.GetCX() + search_radius_x), cluster2.GetCX() + search_radius_x);
            candidate_location.y = min(max(candidate_location.y, cluster2.GetCY() + search_radius_y), cluster2.GetCY() + search_radius_y);

            double candidate_cost = calculateCost(candidate_location, layer1_indices, Clusters_Layer1);

            // 判断是否接受新位置
            if (candidate_cost < min_cost || 
                exp((min_cost - candidate_cost) / temperature) > (rand() / double(RAND_MAX))) 
            {
                best_location = candidate_location;
                min_cost = candidate_cost;
            }
        }

        // 逐步降温
        temperature *= params.cooling_rate;
    }
#else
    double* workspace = new double[10000];
    REAL x[2] = {cluster2.GetCX(), cluster2.GetCY()};
    //std::cout << x[0] << " " << x[1] << "\n"; 
    REAL xu[2] = {cluster2.GetCX() + search_radius_y, cluster2.GetCY() + search_radius_y};
    REAL xl[2] = {cluster2.GetCX() - search_radius_x, cluster2.GetCY() - search_radius_x};  //使用更新的search_radius
    BOBYQA_DATA data = {layer1_indices, Clusters_Layer1};
    bobyqa(2, 5,
    bobyqa_obj, (void*)&data,
    x, xl, xu,
    (search_radius_x ) * 0.7, 1 ,
    0, 1000000,  workspace);
    //std::cout << x[0] << " " << x[1] << "\n"; 
    best_location = ProblemData::Coordinate{x[0], x[1]};
    delete workspace;

#endif
    cost_sum += calculateCost(best_location, layer1_indices, Clusters_Layer1);
    std::cout << "cost_sum: " << cost_sum << "\n";
    return best_location;
}

// 优化所有第二层聚类团的缓冲区位置
std::vector<ProblemData::Coordinate> optimizeAllBufferLocations_SA(
    const std::vector<Cluster>& Clusters_Layer2,
    const std::vector<Cluster>& Clusters_Layer1,
    const std::vector<std::vector<int>>& cluster_layer_mapping,
    const SimulatedAnnealingParams& params)
{
    std::vector<ProblemData::Coordinate> best_buffer_locations;

    for (size_t cluster2_idx = 0; cluster2_idx < Clusters_Layer2.size(); ++cluster2_idx) {
        const Cluster& cluster2 = Clusters_Layer2[cluster2_idx];
        const std::vector<int>& layer1_indices = cluster_layer_mapping[cluster2_idx];

        // 使用模拟退火算法为当前第二层聚类团优化缓冲区位置
        ProblemData::Coordinate best_location = optimizeBufferLocation_SA(cluster2, layer1_indices, Clusters_Layer1, params);

        best_buffer_locations.push_back(best_location); // 将最优位置添加到结果列表
    }

    return best_buffer_locations;
}
