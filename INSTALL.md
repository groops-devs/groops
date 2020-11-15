# Installation

This guide provides step-by-step instructions for compiling and running GROOPS from scratch.

- [Get the GROOPS Source Code](#get-the-groops-source-code)
- [Microsoft Windows](#microsoft-windows)
- [Linux](#linux)
    - [Ubuntu](#ubuntu)
    - [OpenSUSE](#opensuse)

## Overview

While GROOPS is intended to be a standalone software package, some functionality depends on external libraries.
Hard dependencies are:

 - [the Expat XML parser](https://libexpat.github.io)
 - an implementation of the Linear Algebra Package (LAPACK), for example:
    - [OpenBLAS](https://github.com/xianyi/OpenBLAS)
    - [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html)

Additional libraries extend the feature set of GROOPS and can be optionally enabled at compile time.
At the moment, these include:

- [NetCDF](https://www.unidata.ucar.edu/software/netcdf) for reading and writing NetCDF files
- [zlib](https://zlib.net) for reading and writing compressed files
- the Essential Routines for Fundamental Astronomy ([liberfa](https://github.com/liberfa/erfa)) for high-precision
  Earth rotation

Another optional dependency is an implementation of the Message Passing Interface standard (MPI,
for example [MPICH](https://www.mpich.org/) or [Microsoft MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi)).
Resource intensive tasks and algorithms are designed and implemented to be optionally run in
parallel on distributed systems.

To visualize data sets, GROOPS requires an installation of the [Generic Mapping Tools (GMT)](https://www.generic-mapping-tools.org/).

## Get the GROOPS Source Code

You can download the source code of a specific version on the
[Releases](https://github.com/groops-devs/groops/releases) page, or
clone the repository to always get the latest updates:

```
git clone https://github.com/groops-devs/groops.git
```

## Microsoft Windows

GROOPS under Windows requires CMake, and 64bit C++14 and Fortan compilers.
A convenient way to install all required tools is through [MSYS2](https://www.msys2.org).
This installation guide assumes that the GROOPS source code is located in `C:\groops`.

1. Download the MSYS2 installer and follow the [installation guide](https://www.msys2.org/#installation).

2. **Important**: After successful installation, close the MSYS2 terminal from step 1 and open the **MSYS2 MinGW 64-bit terminal**
    through `Start Menu > MSYS2 64-bit > MSYS2 MinGW 64-bit`.

    The command prompt in the terminal window should now read `username@hostname MINGW64`.

3. Install compilers, cmake, expat, OpenBLAS, and LAPACK:
    ```
    pacman -S mingw-w64-x86_64-toolchain mingw-w64-x86_64-cmake expat mingw64/mingw-w64-x86_64-openblas
    ```

4. *(Optional)* Install the NetCDF library:
    ```
    pacman -S mingw-w64-x86_64-netcdf
    ```

5. *(Optional)* Download and install liberfa:

    5.1. Install the `tar` utility and required build tool:
    ```
    pacman -S tar make
    ```

    5.2. Download and build the ERFA library:
    ```
    mkdir -p /c/groops/lib && cd /c/groops/lib
    wget https://github.com/liberfa/erfa/releases/download/v1.7.0/erfa-1.7.0.tar.gz
    tar -xvf erfa-1.7.0.tar.gz
    cd erfa-1.7.0
    ./configure
    make
    make install
    ```

6. *(Optional)* Install Microsoft MPI:

    6.1 Download and install the [Microsoft MPI](https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi).

    6.2 Install the MSYS2 `msmpi` package:
    ```
    pacman -S mingw-w64-x86_64-msmpi
    ```

7. Create the build directory and compile GROOPS:
    ```
    mkdir /c/groops/source/build && cd /c/groops/source/build
    cmake.exe .. -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release
    mingw32-make.exe
    ```

8. Make sure to also read the [post-installation steps](#windows-post-installation-steps)

### Graphical User Interface (GUI)

The GROOPS GUI is based on [Qt](https://www.qt.io/) and is tested with Qt version 5.15.1.
Precompiled Windows binaries of the GROOPS GUI are provided with each [release](https://github.com/groops-devs/groops/releases).
To compile the GUI yourself, you need to:

1. [Download and install](https://www.qt.io/download-qt-installer) Qt (registration required).

2. When prompted to choose which Qt components to install, select `Select Package Categories > LTS` and
   then select `Qt > Qt 5.15.1` or a newer version. Under `Developer and Designer Tools`, `Qt Creator Debugger Support`,
   `Debugging Tools for Windows`, `cmake`, and `Ninja` should be selected automatically.

3. Open the project file `C:\groops\gui\groopsGui.pro` in Qt Creator and build the project.

### Generic Mapping Tools (GMT)

The Generic Mapping Tools (GMT) are an optional dependency of GROOPS and enable the generation of high-quality
figures.
GMT provides [Windows binaries](https://github.com/GenericMappingTools/gmt/releases) which can be easily installed.
The current GROOPS release is tested against GMT version 6.0.0.

### Windows post-installation steps

After the installation of GROOPS and GROOPS GUI, we recommend some post-installation configuration steps to make
working with GROOPS easier.

1. To use the GROOPS and GROOPS GUI binaries without directory prefix, you have to add the required
   directories to the system path.

    1.1. Open the Control Panel through the Windows Start Menu: `Windows System > Control Panel`.

    1.2. In the Control Panel window, go to `User Accounts > User Accounts`.

    1.3. There you should click on `Change my environment variables`, which will open a new window.

    1.4. In the environment variable window, select `Path` and click `Edit...`. A pop-up window will appear
    where you can add the following directories to your system path:
    ```
    "C:\groops\bin"
    "C:\msys64\mingw64\bin"
    "C:\Program Files\Microsoft MPI\Bin"
    ```
    In case you manually compiled the GUI, additionally add the directory:
    ```
    "C:\Qt\5.15.1\mingw81_64\bin"
    ```

2. *(Optional)* Set the environment variable `OPENBLAS_NUM_THREADS` to the number of threads to use for matrix operations.

## Linux

Most Linux distributions provide all GROOPS dependencies through their package managers.
We provide a detailed installation guide for Ubuntu and OpenSUSE, the installation steps
are however very similar for other distributions.

### Ubuntu

The installation procedure for Ubuntu is representative for all Debian based distributions,
however the individual package names may differ.
Check your distribution's documentation to find the correct packages.
This installation guide assumes that the GROOPS source code is located in `$HOME/groops`.

1. First, make sure your system is up to date:
    ```
    sudo apt update && sudo apt upgrade
    ```

2. Install dependencies and build tools:
    ```
    sudo apt-get install g++ gfortran cmake libexpat1-dev libopenblas-dev
    ```

3. *(Optional)* Install the NetCDF development package:
    ```
    sudo apt-get install libnetcdf-dev
    ```

4. *(Optional)* Install liberfa development packages:
    ```
    sudo apt-get install liberfa-dev
    ```

5. *(Optional)* Install MPI development packages:
    ```
    sudo apt-get install mpi-default-dev
    ```

6. Create the build directory and compile GROOPS:
    ```
    mkdir $HOME/groops/source/build && cd $HOME/groops/source/build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make
    ```

7. Make sure to also read the [post-installation steps](#linux-post-installation-steps).

#### Graphical User Interface (GUI)

The GROOPS GUI depends on Qt packages.
To install the required packages, run:
```
sudo apt-get install qtbase5-dev
```
Then, change into the `gui` directory and compile the source code:
```
cd $HOME/groops/gui
qmake
make
```

#### Generic Mapping Tools (GMT)

Ubuntu provides packages for the Generic Mapping Tools:
```
sudo apt-get install gmt gmt-gshhg
```

### OpenSUSE

1. First, make sure your system is up to date:
    ```
    sudo zypper up
    ```
2. Install dependencies and build tools:
    ```
    sudo zypper install gcc-c++ gcc-fortran cmake libexpat-devel openblas-devel
    ```

3. *(Optional)* Install the NetCDF development package:
    ```
    sudo zypper install netcdf-devel
    ```

4. *(Optional)* Install liberfa development packages:

    4.1 Add the OpenSUSE Science Repository (change OpenSUSE release version if necessary):
    ```
    sudo zypper addrepo http://download.opensuse.org/repositories/science/openSUSE_Leap_15.2/ science
    ```

    4.2 Install the required package:
    ```
    sudo zypper install erfa-devel
    ```

5. *(Optional)* Install MPI development packages:
    ```
    sudo zypper install mpich-devel
    ```

6. Create the build directory and compile GROOPS:
    ```
    mkdir $HOME/groops/source/build && cd $HOME/groops/source/build
    cmake ..
    make
    ```

7. Make sure to also read the [post-installation steps](#linux-post-installation-steps).

#### Graphical User Interface (GUI)

The GROOPS GUI depends on Qt packages.
To install the required packages, run:
```
sudo zypper install libqt5-qtbase-devel
```
Then, change into the `gui` directory and compile the source code:
```
cd $HOME/groops/gui
qmake-qt5
make
```

#### Generic Mapping Tools (GMT)

The OpenSUSE packages for the Generic Mapping Tools are available in the `GEO` repository
(change OpenSUSE release version if necessary):
```
sudo zypper addrepo http://download.opensuse.org/repositories/Application:/Geo/openSUSE_Leap_15.2/ GEO
```
Then install the packages:
```
sudo zypper install GMT GMT-doc GMT-coastlines
```

### Arch Linux

GROOPS is packaged for the [Arch User Repository](https://wiki.archlinux.org/index.php/Arch_User_Repository).
You can install the [groops-git](https://aur.archlinux.org/packages/groops-git/) package providing the core GROOPS executables,
and the [groopsgui-git](https://aur.archlinux.org/packages/groopsgui-git/) package providing the GUI and documentation.

The easiest way to do this is through an [AUR helper](https://wiki.archlinux.org/index.php/AUR_helpers). If you are using `yay`,
for example, you can install GROOPS and the GUI by executing
```
yay -S groops-git groopsgui-git
```

If you want to develop for GROOPS, a manual installation is preferrable:

1. First, make sure your system is up to date:
    ```
    sudo pacman -Syu
    ```
2. Install dependencies and build tools:
    ```
    sudo pacman -S git cmake gcc gcc-gfortran bash expat lapack zlib
    ```
3. *(Optional)* Install the NetCDF development package:
    ```
    sudo pacman -S netcdf-cxx
    ```
4. *(Optional)* Install liberfa development packages. liberfa is available as an [AUR package](https://aur.archlinux.org/packages/erfa/).

5. *(Optional)* Install an MPI development package, eg. `openmpi`:
    ```
    sudo pacman -S openmpi
    ```
6. Fetch the GROOPS git repository:
    ```
    git clone https://github.com/groops-devs/groops.git
    ```
7. Create the build directory and compile GROOPS, Adjust the paths if you checked out the GROOPS git repository to a different directory.
    ```
    mkdir $HOME/groops/source/build && cd $HOME/groops/source/build
    cmake ..
    make
    ```
8. Make sure to also read the [post-installation steps](#linux-post-installation-steps).

```
#### Graphical User Interface (GUI)

The GROOPS GUI depends on Qt packages.
To install the required packages, run:
```
sudo pacman  -S qt5-base
```
Then, change into the `gui` directory and compile the source code:
```
cd $HOME/groops/gui
qmake
make
```

#### Generic Mapping Tools (GMT)

The Generic Mapping Tools are available from the [Arch User Repository](https://wiki.archlinux.org/index.php/Arch_User_Repository).
Install the [gmt6](https://aur.archlinux.org/packages/gmt6) and [gmt-coast](https://aur.archlinux.org/packages/gmt-coast) packages.


### Linux post-installation steps

After the installation of GROOPS and GROOPS GUI, we recommend some post-installation configuration steps to make
working with GROOPS easier.

1. To use the GROOPS and GROOPS GUI binaries without directory prefix, you have to add the required
   directories to the system path:

   ```
   echo "export PATH=$PATH:$HOME/groops/bin" >> $HOME/.bashrc
   source $HOME/.bashrc
   ```

2. *(Optional)* Set the environment variable `OPENBLAS_NUM_THREADS` to the number of threads to use for matrix operations.
