// Course_Example_Computation_Geometry_Shandong_University.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include"course1_convexhull.h"

int main()
{
    std::vector<point> smallPointSet = generate_pointset_demo();
    test_convexhull(smallPointSet);
    
    std::vector<point> largePointSet = generate_pointset_large(2000);
    test_convexhull(largePointSet);

    test_convexhull(largePointSet,false);
    return 0;
}
