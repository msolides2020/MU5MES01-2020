# How to install FEniCS

## On your laptop

Instructions on how to install FEniCS are available https://fenicsproject.org/download/.
You also need a working `python3` installation with several additional packages.

**Except for Ubuntu distributions (see below)**, we suggest at first to try to use Anaconda, which has an advanced package manager that can provide you without effort most of the software that you need. Install the Anaconda version based on *python3*.

### On ubuntu (preferred installation)

Install basic Python 3 packages

```
sudo apt update
sudo apt upgrade
sudo apt install black cython3 jupyter-core jupyter-client jupyter-notebook pylint python3-autopep8 python3-numpy python3-flake8 python3-h5py python3-jupyter-console python3-matplotlib python3-mpi4py python3-notebook python3-petsc4py python3-pytest python3-requests python3-scipy python3-slepc4py python3-skimage python3-sphinx python3-sympy spyder pelican jupyter-nbconvert jupyter-qtconsole python3-qtconsole
```

Then install the FEniCS specific packages

```
sudo apt-get install software-properties-common
sudo add-apt-repository ppa:fenics-packages/fenics
sudo apt-get update
sudo apt-get install --no-install-recommends fenics
```

### Anaconda install (Mac and Linux)

This procedure will work on Mac and any Linux distribution.

First install Anaconda, following the instruction here  `https://docs.continuum.io/anaconda/install/`.

Then run following commands in your terminal:

```
conda config --add channels conda-forge
conda config --set channel_priority strict
conda create -n fenicsproject
source activate fenicsproject
conda install fenics mshr
```

The above commands will create a virtual environment called `fenicsproject`, in which the FEniCS library will be installed. This avoids conflicts between various packages. This means that  **you need to activate the `fenicsproject` environment** when you install additional packages or perform FEniCS computations. If the environment is activated, your terminal prompt will display `(fenicsproject)`, like so:

```
(fenicsproject) maurini@maurinis-MacBook-Pro:~$
```

We suggest also to install other useful packages

```
conda install gmsh pip jupyter spyder matplotlib scipy vtkplotter meshio
```

If you need to install further packages, just type `conda install packagename` (preferred) or `pip3 install packagename`.

### On windows

FEniCS is not distributed for Windows boxes. For Windows 10, the preferred option is the [Windows subsystem for linux (WSL)](https://docs.microsoft.com/en-us/windows/wsl/install-win10): install the Ubuntu distribution, then refer to the section above.

If you cannot install the WSL, you will need to create Ubuntu virtual machine. Get in touch with your instructors in that case.

## Cloud-based installations
You can run python programs jupyter notebooks and FEniCS on online servers. The basic service is free and can be a solution if all other installation systems fail.

* *advantages:* You do not need to install anything on your machine and you do note use the ressources of your machine

* *disadvantages:* It can be slow (especially Azure, colab seems better). You can use only jupyter notebooks, you share your date with google or microsoft, you do not a full control of the system, you need to be online with a good network connection.

### Google colab
* You need a google account
* You can use this example to start https://colab.research.google.com/drive/1yhw-E21eEa4-ej27abiJiHj064iRaetS#scrollTo=iF68ZCvDTvo5
* You can save the notebooks and your working environment on your google drive

### Microsoft Azure
Another possibility is using azure notebbok (you need a Microsoft account)
https://fenics-cmaurini.notebooks.azure.com/nb/notebooks/FEniCS%20demo.ipynb


### Using docker virtual machines (suggested method if other fail).:
* See https://fenicsproject.org/download/ on the docker section and https://fenics-containers.readthedocs.io/en/latest/index.html

## How to test the installation

Run following commands in your terminal:

```
wget https://fenicsproject.org/docs/dolfin/latest/python/_downloads/598330b504d63e359baad030e1010987/demo_poisson.py
python3 demo_poisson.py
```
