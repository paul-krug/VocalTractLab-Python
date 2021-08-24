# VocalTractLabBackend-dev
This repo contains the VocalTractLab backend source code and the C/C++ API for on-going development work. *This is not the place for official stable releases!* 

You can find the bleeding edge builds here including new and experimental features, some of which may break your old projects. For official stable releases, please check [the VocalTractLab website](https://www.vocaltractlab.de).

Please feel free to fork this repo and make your own contributions, though! 

The main branch of this repo is reviewed on a semi-regular basis for inclusion into the official release.

This repo may be included in other repos as a submodule wherever the backend source code or the C/C++ API is needed.

## Getting started
- Clone the current main branch:
```
git clone https://github.com/TUD-STKS/VocalTractLabBackend-dev
```
### Build using CMake (Windows, Linux, macOS)
- Get the latet release of [CMake for your platform](https://cmake.org/)
- Create a folder ``out``inside the cloned repository folder
- Open a shell/command prompt and navigate to ``out``
- Configure the project and generate a build system:
```
cmake ..
```
- Build the library (still from within the folder ``out``)
```
cmake --build .
```

### Build using Visual Studio 2019 (Windows)
- Open ``VocalTractLabApi.sln``
- Build the project ``VocalTractLabApi``

