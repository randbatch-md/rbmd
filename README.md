<!-- Improved compatibility of back to top link: See: https://github.com/othneildrew/Best-README-Template/pull/73 -->
<a name="readme-top"></a>
<!--
*** Thanks for checking out the Best-README-Template. If you have a suggestion
*** that would make this better, please fork the repo and create a pull request
*** or simply open an issue with the tag "enhancement".
*** Don't forget to give the project a star!
*** Thanks again! Now go create something AMAZING! :D
-->



<!-- PROJECT SHIELDS -->
<!--
*** I'm using markdown "reference style" links for readability.
*** Reference links are enclosed in brackets [ ] instead of parentheses ( ).
*** See the bottom of this document for the declaration of the reference variables
*** for contributors-url, forks-url, etc. This is an optional, concise syntax you may use.
*** https://www.markdownguide.org/basic-syntax/#reference-style-links
-->
<!--
[![Contributors][contributors-shield]][contributors-url]
[![Forks][forks-shield]][forks-url]
[![Stargazers][stars-shield]][stars-url]
[![Issues][issues-shield]][issues-url]
[![MIT License][license-shield]][license-url]
[![LinkedIn][linkedin-shield]][linkedin-url]
-->

<!-- TABLE OF CONTENTS -->

<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>
    </li>
    <li>
      <a href="#getting-started">Installation</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
          <ul>
              <li><a href="#onLinux">On Linux</a></li>
          </ul>
        <li><a href="#installation">Installation</a></li>
          <ul>
              <li><a href="#remark">Remark</a></li>
          </ul>
        <li><a href="#DockerImage">Docker Image</a></li>
      </ul>
    </li>
    <li><a href="#gettingStarted">Getting Started</a></li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#contributing">Contributing</a></li>
    <li><a href="#license">License</a></li>
    <li><a href="#contact">Contact</a></li>
    <li><a href="#acknowledgments">Acknowledgments</a></li>
  </ol>
</details>





<!-- ABOUT THE PROJECT -->
## About The Project 
<!--
[![Product Name Screen Shot][product-screenshot]](https://example.com)
-->

RBMD (Random Batch Molecular Dynamics) is a GPU-CPU heterogeneous accelerating molecular dynamics software baed on random batch series algorithms.


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- GETTING STARTED -->

## Installation 

### Prerequisites
#### On Linux
The third-party library vtk-m needs to be installed. First, clone the VTK-m from github on `your_folder`
```
git clone https://gitlab.kitware.com/vtk/vtk-m.git
```
Then switch to the branch: v1.9.0 (`cat version.txt` to check the version of vtkm for install)
```
git checkout -b 1.9.0 v1.9.0
```
Finally, compile VTK-m with cmake (CUDA support required)  , CMake/3.20, GCC/9.0, and CUDA/11.0 have been confirmed to be necessary for compiling VTK-m. Furthermore, the compilation process is related to the GPU driver, and if deploying on a server, it needs to be compiled on the compute nodes.
```
cmake .. -DBUILD_SHARED_LIBS=OFF -DCMAKE_BUILD_TYPE=Release -DVTKm_ENABLE_CUDA=ON -DVTKm_USE_64BIT_IDS=OFF -DVTKm_ENABLE_TESTING=OFF -DVTKm_ENABLE_RENDERING=OFF -DBUILD_TESTING=OFF -DCMAKE_INSTALL_PREFIX=your_folder/vtkm_out
```
And install in the `vtkm-out`
```
make -j$(nproc)
make install
```
For rbmd compilation,

1. Install the source code on `your_folder` and switch to the `dev` branch
   ```
   git clone https://github.com/randbatch-md/rbmd.git
   git checkout dev
   ```
2. cmake RBMD based on the same packages :CMake/3.22, GCC/9.3, and CUDA/11.8 in the `build` folder
   ```
   cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_LINKER=your_folder/rbmd/framework/tools/Linux/Release/jsoncpp/lib64
   ```
3. compile the executable file `rbmd` and then `rbmd` is on the `framework` folder 
   ```
   make -j$(nproc)
   ```
   

#### Remark

1. If the step for `make` rbmd is wrong in the first time, the `Linux` folder needs to be deleted
```
rm -rf your_folder/rbmd/framework/tools/Linux
```
and complie again.

2. If in the end of  `make`, the error `cannot find -ljsoncpp` is occured, the environment variables needs to be specified and the step for make is changed
```
export LD_LIBRARY_PATH=your_folder/rbmd/framework/tools/Linux/Release/jsoncpp/lib64:$LD_LIBRARY_PATH
export LIBRARY_PATH=your_folder/rbmd/framework/tools/Linux/Release/jsoncpp/lib64:$LIBRARY_PATH
cmake  ..  -DCMAKE_BUILD_TYPE=Release -DCMAKE_LINKER=your_folder/rbmd/framework/tools/Linux/Release/jsoncpp/lib64
```

### Docker Image

You can also pull the docker image and run docker dircetly

```
docker pull ghcr.io/randbatch-md/rbmd:1.0.0
docker run --rm --gpus all -it -v $PWD:/app rbmd  /bin/bash -c "rbmd -j rbmd.json"
```


<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Getting started

run the example json file `rbmd.json` in `framework`, then rbmd is executed
```
./rbmd -j rbmd.json
```


<!-- USAGE EXAMPLES -->
## Usage


_For more examples, please refer to the [https://www.randbatch.com/guide/](https://www.randbatch.com/guide/)_

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTRIBUTING -->
## Contributing

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- LICENSE -->
## License

Distributed under the GPL-3.0 License. See `LICENSE.txt` for more information.

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- CONTACT -->
## Contact

<!--
Your Name - [@your_twitter](https://twitter.com/your_username) - email@example.com
-->
Weizhu scientific computing platform Link: [https://github.com/randbatch-md](https://github.com/randbatch-md)

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

<!--
Use this space to list resources you find helpful and would like to give credit to. I've included a few of my favorites to kick things off!
-->

<p align="right">(<a href="#readme-top">back to top</a>)</p>



