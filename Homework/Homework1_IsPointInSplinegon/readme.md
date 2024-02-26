实验1：判断点是否在splinegon内

提示：
1.了解判断点是否在凹多边形的算法 (射线法：https://www.cnblogs.com/muyefeiwu/p/11260366.html 【该文章思路正确，代码实现有问题，判断浮点数相等要用abs(a-b)<1e-6而不是a==b;同理,判断a>b要用a>b+1e-6】)
2.了解splinegon的定义
3.根据前两个知识，设计算法

文件夹中提供了测试用的splinegon库：big_random_point/XX/XXXXX.spilnegon
输入格式：

    第一行： n
    
    其中n是多边形的边数
   
    第2~n+1行： Sx,Sy,Ex,Ey,Cx,Cy
        
    其中二次bezier曲线的起点S=(Sx,Sy),终点E=(Ex,Ey),控制点C=(Cx,Cy)
    
要求(优先第一条，第二条稍复杂，供大家思考)：
1.算法的时间复杂度是O(n)
2.算法中只进行一次与曲线相关的运算
