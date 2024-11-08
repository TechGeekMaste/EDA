#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <iomanip>
#include <cstdlib>  //包含命令行参数处理库
#include "constrain.h"
#include "problem.h"
#include <stdio.h>
#include <chrono> // 包含时间库
using namespace std;
#include "kdtree.h"
#include "Cluster.h"
#include "GSR.h"
#include "GSR1.h"

void printClusterInfo(const std::vector<Cluster> *clusters, int current_buf_num, std::ostream& out)
{
    for (size_t j = 1; j < clusters->size(); ++j)
    {
        double cx = clusters->at(j).GetCX();
        double cy = clusters->at(j).GetCY();
        out << "- BUF" << j - 1 + current_buf_num << " BUF ( " << cx << " " << cy << " ) ;" << std::endl;
    }
}

void printNetInfo(const std::vector<Cluster> *clusters, int lastlayer_buf_num, int current_buf_num, bool isFirstLayer, std::ostream& out)
{
    for (size_t j = 1; j < clusters->size(); ++j)
    {
        out << "- net_" << j - 1 + current_buf_num << " ( BUF" << j - 1 + current_buf_num << " ) ( ";

        const std::vector<int> &children = clusters->at(j).GetChildIndex();
        for (size_t k = 0; k < children.size(); ++k)
        {
            out << (isFirstLayer ? "FF" : "BUF") << (isFirstLayer ?children[k] + lastlayer_buf_num : children[k] + lastlayer_buf_num-1);
            if (k < children.size() - 1)
            {
                out << " ";
            }
        }
        out << " ) ;" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // 检查命令行参数
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <constraints_file> <problem_file>" << std::endl;
        std::cerr << "       Use -h for help." << std::endl;
        return 1;
    }

    // 检查是否需要帮助信息
    if (std::string(argv[1]) == "-h") {
        cout << "This program requires two input files:" << std::endl;
        cout << "1. <constraints_file>: Path to the constraints file." << std::endl;
        cout << "2. <problem_file>: Path to the problem definition file." << std::endl;
        return 0;
    }

    // 从命令行参数中读取文件路径
    string constraints_file = argv[1];
    string problem_file = argv[2];
    auto start = chrono::high_resolution_clock::now(); // 记录开始时间

        // 创建文件输出流对象并打开文件
    std::ofstream output_file("solution.def");
    if (!output_file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    //设置输出格式
    output_file << fixed << setprecision(0);
    ConstraintsData C_data;
    ProblemData P_data;

    // cout << "Constrain Data:" << "\n";
    readConstraints(constraints_file, C_data);
    C_data.display();

    // cout << "Problem Data:" << "\n";
    readProblemDef(problem_file, P_data);
    //P_data.display();

    output_file << "UNITS DISTANCE MICRONS " <<  P_data.units_distance << " ;" << endl;
    output_file << "DIEAREA " ;
    for (const auto& coord : P_data.die_area){
      output_file << "( " << coord.x << " "<< coord.y << " ) " ;
    }
    output_file << ";" << endl;
    output_file << "FF " << "( " << P_data.ff_size.x << " " << P_data.ff_size.y << " ) ;" << endl;
    output_file << "BUF " << "( " << P_data.buf_size.x << " " << P_data.buf_size.y << " ) ;" << endl;
    output_file << "CLK " << "( " << P_data.clk_position.x << " " << P_data.clk_position.y << " ) ;" << endl;

    int layer = 1;

    // 初始化 AllClusters 和 points
    std::vector<std::vector<Cluster> *> AllClusters; // 定义并初始化 AllClusters

    std::vector<Point> points;                 // kd树数组
    //std::vector<Point> points_ff;                       // kd树数组
    //std::vector<Point> points_buf;
    //std::vector<std::vector<int>> cluster_layer_mapping; // 声明 cluster_layer_mapping

    // 第一层聚类
    GSR1(P_data, C_data, AllClusters, points, layer);

    std::vector<Cluster> *clusters = AllClusters[AllClusters.size() - 1]; // 获取指向 vector<Cluster> 的指针
    double rc_sum = std::numeric_limits<double>::infinity(); //设为无穷大
    size_t previous_cluster_size = clusters->size();
    while (rc_sum > C_data.max_net_rc )
    {
        layer++;
        GSR(P_data, C_data, AllClusters, points, layer );
        clusters = AllClusters[AllClusters.size() - 1]; // 更新指向 vector<Cluster> 的指针
        rc_sum = 0;
        for (const auto& cluster : *clusters){
          double distance = CountDistance(cluster.GetCX(), cluster.GetCY(), P_data.clk_position.x, P_data.clk_position.y);
          double RC = CountRC(distance, C_data.net_unit_r, C_data.net_unit_c);
          rc_sum += RC;   //更新rc_sum
        }
        cout << "rc_sum:"<< rc_sum <<endl;
        //clusters->at(1).Print();

        //检查是否聚类到极限
        if (clusters->size() == previous_cluster_size){
          break;
        }
        previous_cluster_size = clusters->size();

    }
   //clusters = AllClusters[2];
  // clusters->at(32).Print();
    

    // 格式输出
    int lastlayer_buf_num = 0; // 第一层至上一层buf总数
    int current_buf_num = 0;
    
    for (size_t i = 0; i < AllClusters.size(); ++i) {
        clusters = AllClusters[i];
        current_buf_num += clusters->size() - 1;
    }
    //输出ff和buffer总数
    output_file << "COMPONENTS " << current_buf_num + P_data.ff_count  << " ;" << endl;
    //输出ff坐标
    int index = 1;
    for (const auto& pos : P_data.ff_positions){
      output_file << "- FF" << index << " FF " << "( " << pos.x - P_data.ff_size.x/2 << " " << pos.y - P_data.ff_size.y/2 << " ) ;" << endl;
      index++;
    }

    lastlayer_buf_num = 0;
    current_buf_num = 0;
    for (size_t i = 0; i < AllClusters.size(); ++i)
    {
        clusters = AllClusters[i];
        printClusterInfo(clusters, current_buf_num, output_file);
        current_buf_num += clusters->size() - 1;
        lastlayer_buf_num = clusters->size() - 1; // 更新 lastlayer_buf_num
                                                  // output_file << std::endl; // 添加换行以便于阅读
    }
 
    output_file << "END COMPONENTS ;" << std::endl;
    output_file << "NETS " << current_buf_num + 1 << " ;" << std::endl;

    lastlayer_buf_num = 0;
    current_buf_num = 0;
    
    for (size_t i = 0; i < AllClusters.size(); ++i)
    {

        if (!i)
        {
            clusters = AllClusters[i];
            printNetInfo(clusters, 0, current_buf_num, i == 0, output_file);
            current_buf_num += clusters->size() - 1; // 更新 current_buf_num
        }
        else
        {
            double tempt = clusters->size() - 1;
            clusters = AllClusters[i];
            printNetInfo(clusters, lastlayer_buf_num, current_buf_num, i == 0, output_file);
            current_buf_num += clusters->size() - 1; // 更新 current_buf_num
            lastlayer_buf_num += tempt;
        }
    }

    output_file << "- net_" << current_buf_num << " ( CLK ) ( ";
    const auto &finalClusters = AllClusters.back(); // 使用最后一层聚类的指针
    for (size_t k = 1; k < finalClusters->size(); ++k)
    {
        output_file << "BUF" << finalClusters->at(k).GetCIndex() + lastlayer_buf_num - 1;
        if (k < finalClusters->size() - 1)
        {
            output_file << " ";
        }
    }
    output_file << " ) ;" << std::endl;
    output_file << "END NETS ;" << std::endl;

    auto end = chrono::high_resolution_clock::now(); // 记录结束时间
    chrono::duration<double> duration = end - start; // 计算用时
    cout << "Time taken: " << duration.count() << " seconds" << std::endl;

    return 0;
}
