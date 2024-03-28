## 关于semd
可压缩Navier-Stokes方程GPU解算器

## windows下的依赖安装
### Visual Studio 2019 安装
参照vistual studio官方[安装教程](https://visualstudio.microsoft.com/zh-hans/vs/)
### cuda 安装
对于windows平台，需要[cuda10.2](https://developer.nvidia.com/cuda-10.2-download-archive?target_os=Windows&target_arch=x86_64&target_version=10&target_type=exelocal)或以上版本

cuda其他版本[下载](https://developer.nvidia.com/cuda-toolkit-archive)
### VTKm安装
zaran依赖[VTKm](https://m.vtk.org/)库，VTKm的代码[仓库地址](https://gitlab.kitware.com/vtk/vtk-m)，

关于VTKm的介绍，可下载这个[PPT](https://m-old.vtk.org/images/f/f3/VTKm_Tutorial_VIS19.pptx)，

更多信息，参考官方[用户指南](https://gitlab.kitware.com/vtk/vtk-m-user-guide/-/wikis/home)，
[zanran](https://gitee.com/feaSYtech/zaran)使用了VTKm的最新版本为v1.6

##### 下载vtkm
```
git clone https://gitlab.kitware.com/vtk/vtk-m/
git checkout v1.9.0 
```
##### 编译安装

在vtkm文件夹上右键，用Vistual Studio打开

<img src="docs/images/open-vtkm.png" />

VS2017之后，内置了对cmake的支持，省去使用cmake-gui生成sln的步骤，更方便使用。
打开vtkm源码之后，vistual studio搜索到vtkm文件夹下的存在CMakeLists.txt，调用内置的cmake构建器
可使用管理配置进行编译配置，如配置Release或Debug编译，以及编译器选择。编译器可选择msvc或clang，
甚至远程服务器上的编译环境。
<img src="docs/images/vs-config.png" />

编译前对VTKm进行配置
<img src="docs/images/vtkm-options.png" />
其中：

+ VTKm_ENABLE_CUDA：       是否打开GPU编译
+ VTKm_ENABLE_RENDERING：  是否打开内置渲染模块，打开
+ VTKm_USER_DOUBLE_PRECISION：是否使用双精度，关闭，桌面显卡的双精度性能很差 

其他选项默认即可。
配置好选项之后，保存，vistual studio自动使用内置cmake生成，之后进行编译，安装即可。

<img src="docs/images/compile_build.png" />

编译过程会使用并行编译，根据机器性能，cpu版本约10分钟，gpu版本约30分钟。安装的路径在`vtkm/out/install/*`下，
根据不同配置名称生成相应的文件夹，方便第三方库切换使用。

## semd的编译


```
git clone https://gitee.com/sciai_cq/semd.git
```
同样，右键semd文件夹，使用vistual studio打开，配置

<img src="docs/images/zaran_config.png" />

其中`VTKm_DIR`选择安装好的vtkm路径，semd会根据VTKm的配置是否启用cuda，编译gpu版本

## semd的运行
编译之后，在`semd\out\build\x64-release-gpu\framwork`为可执行文件路径
<img src="docs/images/zaran.png" />
其中`semd.exe`为可执行文件，使用cmd或powershell即可运行
```
semd.exe -i xxx.i
```
zaran需要一个配置文件，设置运行的参数。
#### config文件
配置文件为 `.i` 作为后缀
```
[Mesh]
 dims = '200 200'
 x_range = '0 4'
 y_range = '0 2'
 z_range = '0 1'
[]


[Executioner]
	steps = 100000
	dt = 4e-04
[]
```
*配置文件的结构和参数会随着程序开发改变*

2DRiemann问题的运行结果：

<img src="docs/images/riemann2d.png" />


## 开发路线图
- [x] 加入fmt
- [x] 加入gtest，放在contrib目录下
- [x] 加入getpot，放在contrib目录下
- [x] 文件目录结构重构
- [ ] zaran作为命名空间

