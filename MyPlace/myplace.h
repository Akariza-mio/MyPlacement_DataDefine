#ifndef MYPLACE_H
#define MYPLACE_H

#include "placedata.h"
#include "common.h"
#include <cstddef>

// CoreRegion：整个布局区域的外包矩形；当前随机布局会使用它。
// 坐标单位与 Module::center、Module::width、Module::height 一致。
// 例如 CoreRegion{0, 0, 100, 80} 表示左下角 (0,0)、右上角 (100,80)。
// 区域宽度 = xmax - xmin，区域高度 = ymax - ymin。
struct CoreRegion
{
    float xmin; // 左边界的横坐标；不是单元中心的最小横坐标。
    float ymin; // 下边界的纵坐标。
    float xmax; // 右边界的横坐标；不是区域宽度。
    float ymax; // 上边界的纵坐标；不是区域高度。
};

// 前向声明：先告诉编译器存在这些类，完整定义在下面。
class Bin_2D;
class MyPlacer;

// Bin_2D：把整个区域切成二维小网格后，其中一个网格的数据。
// 可以把它理解为统计局部拥挤程度的小格子，不是电路中的一个单元。
// 这些字段为任务5的密度优化预留；当前 randomPlace() 不使用它们。
// 几何关系如下是应保持的约定，声明字段本身不会自动计算这些值。
// 本类未给 float 字段设置默认值，后续使用前应显式初始化。
class Bin_2D
{
public:
    POS_2D center; // 网格中心 (x,y)，通常为 (ll + ur) / 2；不是单元中心。
    POS_2D ll;     // lower left：网格左下角坐标 (x,y)。
    POS_2D ur;     // upper right：网格右上角坐标 (x,y)。
    float width;   // 网格宽度，等于 ur.x - ll.x。
    float height;  // 网格高度，等于 ur.y - ll.y。
    float area;    // 网格面积，等于 width * height。

    // 以下含义参照任务5《全局布局设计》3.2节；当前代码尚未计算这些值。
    // 虽然名字含 Density，资料中的用法是先累计面积贡献，汇总后再除以
    // 网格面积得到密度，因此不能直接把单个字段当成 0~1 的占用比例。
    float nodeDensity;     // 可移动标准单元和宏单元在本网格中的面积贡献。
    float fillerDensity;   // 虚拟填充单元的面积贡献；filler 用于调节密度分布。
    float terminalDensity; // 固定终端/固定单元的面积贡献，按模型进行密度缩放。
    float DarkDensity;     // 不可放置区域的面积贡献，如 SiteRow 覆盖范围外的空隙。
    float phi;             // 静电密度模型中本网格的电势值，用于后续求电场和密度梯度。
};

// MyPlacer：布局算法的执行者，通过 db 访问单元并更新位置。
// 当前提供随机布局、聚类驱动布局和结果导出。
// 网格、面积统计和梯度字段是后续算法的预留数据。
class MyPlacer
{
public:
    // _db：外部准备好的布局数据库地址，包含单元、网络、布局行等。
    // 这里只保存指针，不复制数据库；外部需保证 db 在布局器使用期间有效。
     MyPlacer(PlaceData *_db)
    {
        db = _db;
    }

    // 随机初始化布局。
    // region：允许放置单元的矩形区域；const 表示不修改它，& 表示按引用传入。
    // seed：随机种子，控制随机数序列；省略时使用默认值 42。
    //       同一运行环境、单元顺序和种子可以复现相同结果，方便调试。
    // 返回值为 void：位置直接写回 db 中的单元，不另外返回一份坐标数组。
    // 修改 Nodes 和 Terminals 中非固定单元的 center，允许重叠。
    // width/height 必须表示当前方向下的实际占用尺寸；方向保持不变。
    // 区域无效、单元尺寸无效或放不下时抛出 std::invalid_argument 异常。
    // 两个列表应保存各单元的唯一指针，坐标修改前会完成输入检查。
    void randomPlace(const CoreRegion& region, unsigned int seed = 42);

    // 先按接线分组，再安排组内位置。seed 只影响位置。
    // 只写回可移动单元的 center，固定单元保持原位，初始位置允许重叠。
    // 如提供 groupIds，成功时按 Nodes、Terminals 中可移动单元的顺序返回组号。
    // 输入检查失败时不改变原坐标或 groupIds。
    void clusterPlace(const CoreRegion& region, unsigned int seed = 42,
                      std::vector<std::size_t>* groupIds = nullptr);

    // 将当前布局导出为 Bookshelf .pl 文件，只读取数据，不重新随机布局。
    // filename：输出文件路径；相对路径以程序运行目录为基准，父目录必须已存在。
    //           如 "test_data/random.pl"；文件存在时覆盖，数据检查在打开文件前完成。
    // 输出 Nodes 和 Terminals 中的全部单元，并按 isFixed 添加 /FIXED 标记。
    // 内部 center 是中心坐标，输出的是左下角：x = center.x - width/2，y 同理。
    // 当前方向约定：模板默认值 orientation == 0 对应 N（无旋转）。
    // 其他方向尚无统一编码，暂时拒绝导出，待任务3提供方向约定后再扩展映射。
    // 返回 void；结尾 const 表示不修改布局器成员，实现中也不修改任何单元。
    // 数据或方向无效时抛出 std::invalid_argument，文件操作失败时抛出 std::runtime_error。
    void writePl(const std::string& filename) const;

    // 二维网格容器：外层 vector 装一组内层 vector，内层保存 Bin_2D 指针。
    // bins[i][j] 可访问某个网格；两个下标与 x/y 的对应由后续网格初始化约定。
    // 当前为空，randomPlace 不创建也不读取这些网格。
    vector<vector<Bin_2D *>> bins;

    float StdCellArea; // 预留：标准单元的总面积，供密度计算使用；当前构造函数未赋值。
    float MacroArea;   // 预留：宏单元（较大的电路模块）的总面积；当前构造函数未赋值。
        
    PlaceData *db; // 指向 PlaceData 数据库；例如 db->Nodes 是单元指针列表。

    // 梯度表示目标函数对坐标变化的敏感程度；沿负梯度方向走，通常使目标减小。
    // VECTOR_3D 有 x/y/z 三个分量，二维布局通常只使用 x/y。
    // 以下容器当前为空，尚未计算；大小及下标约定需在任务5中建立。
    // 按任务5资料：线长梯度对应原始单元，密度/总梯度对应原始单元+filler，
    // fillerGradient 对应填充单元；不能假定四个数组长度相同。
    vector<VECTOR_3D> wirelengthGradient; // 线长目标对单元位置的梯度，供缩短连线使用。
    vector<VECTOR_3D> densityGradient;    // 密度惩罚对位置的梯度，供缓解局部拥挤使用。
    vector<VECTOR_3D> totalGradient;      // 按权重组合各目标后的梯度，供优化器更新位置。
    vector<VECTOR_3D> fillerGradient;     // 虚拟填充单元的梯度，供更新填充单元位置使用。
};
#endif
