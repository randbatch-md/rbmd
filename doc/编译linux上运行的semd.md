# 编译 linux 上运行的 semd 步骤

编译 linux 上可运行的 semd 可执行程序，依赖于预先生成好的 docker 镜像 semd-builder，该镜像会预先添加生成 semd 所需的所有依赖(包括 g++、cmake、vtkm 等)。
如果没有生成 semd-builder 则需要先生成该[镜像](#编译推送docker镜像)，如果生成过了，则直接[编译](#编译生成可执行程序semd)。

以下步骤依赖于[docker](https://www.docker.com/)环境，没有安装，可[下载](https://www.docker.com/products/docker-desktop/)安装。

### 编译推送 semd-builder docker 镜像

- 登录 harbor

  我们本地的 harbor 镜像仓库部署在`192.168.30.35:18082`，构建镜像前需登录(登录过则可以忽略):

  ```sh
  docker login 192.168.30.35:18082 -u admin
  ```

- 添加 Vtkm_Release

  编译时需要将 Vtkm_Release(这个文件并没有放到 git 仓库)放到根目录下的 docker/vendor 目录下:

  ```
    semd/
    |-- docker/
    |   |-- vendor/
    |   |   |-- Vtkm_Release/
  ```

- 构建并推送镜像

  在项目根目录执行以下命令(需指定版本号):

  ```sh
  #用法: sh ./docker/build-image.sh <version>
  sh ./docker/build-image.sh 1.0.0
  ```

  构建成功后会输出相关信息:

  ```
  Congratulations! The semd-builder image:192.168.30.35:18082/library/semd-builder:1.0.0 has been built successfully!!!
  ```

### 编译生成可执行程序 semd

通过 docker 镜像直接编译:

```
#用法: docker run -v <local semd repo path>:/workspace/semd  192.168.30.35:18082/library/semd-builder:<version>
docker run -v xx/semd:/workspace/semd  192.168.30.35:18082/library/semd-builder:1.0.0
```

编译成功后会输出如下信息:

```
Congratulations on the successful build of SEMD! The output path is out/tests/test_semd!
```

此时，在项目根目录的 out/tests 文件夹下可以找到名为 test_semd 的可执行程序。
