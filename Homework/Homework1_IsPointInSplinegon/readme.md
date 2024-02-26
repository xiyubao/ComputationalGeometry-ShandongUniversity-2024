实验1：判断点是否在splinegon内

提示：
1.了解判断点是否在凹多边形的算法
2.了解splinegon的定义
3.根据前两个知识，设计算法
4.用文件夹中的测试用例测试算法的正确性

测试：
输入格式：
第一行： s,n,m
其中s是多边形的个数,n是多边形的边数，m是测试点的个数
    以下重复s次
    n行： Sx,Sy,Ex,Ey,Cx,Cy
    其中二次bezier曲线的起点S=(Sx,Sy),终点E=(Ex,Ey),控制点C=(Cx,Cy)
    m行：Px,Py
    其中测试点P=(Px,Py)

输出格式
b11,b12,...,b1m
b21,b22,...,b2m
...
bs1,bs2,...,bsm
如果Pij在splinegon_i内,bij=true,否则bij=false.

1.上交源代码
2.上交输出文件

要求(优先满足前两条，第三条稍复杂，供大家思考)：
1.所有测试用例全部输出正确
2.算法的时间复杂度是O(n)
3.算法中只进行一次与曲线相关的运算