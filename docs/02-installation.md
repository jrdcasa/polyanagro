# Installation guide

**``Polyanagro``** program is implemented in Python with some parts in C. This library have some dependencies from other python packages . The full list of dependencies can be followed from both **[requirements.txt](../requirements.txt)** and **[setup.py](../setup.py)** files.

These dependencies are handled by the **setup.py** script and they should be installed automatically following these steps. Anyway, few software need to be installed manually (see Pre-installation section)

-------------------

## Pre-installation

-------------------
Previous to install the program  some software need to be installed in your system. Installation of the software using the Ubuntu package repository is given between parentheses, for other distributions please check  the documentation of your distribution.

* `git` (sudo apt-get install git)
* `python3-venv` (sudo apt update; sudo apt-get install python3-venv)
* `python3-dev` (sudo apt-get install python3-dev python-dev)
* `cmake` (sudo apt-get install cmake)
* `libgraphviz-dev` (sudo apt install libgraphviz-dev)
* `swig` (sudo apt-get install swig)

-------------------

## Create a virtualenv environment

-------------------
It is highly recommended to create a virtual environment. This is a Python environment such that the Python interpreter, libraries and scripts installed into it are isolated from those installed in other virtual environments, and (by default) any libraries installed in a “system” Python.

First, we create an isolated Python environment to install the required packages (see dependencies below). Then, activate the virtual environment.

```bash
python3 -m venv <name_of_env>
source <name_of_env>/bin/activate
```

**``WARNING:``** This virtual environment **must** be activate in order to use the program.

Example:

```bash
ubuntu@ubuntu2004:~$ python3 -m venv sandbox_polyanagro
ubuntu@ubuntu2004:~$ source sandbox_polyanagro/bin/activate
(sandbox_topology) ubuntu@ubuntu2004:~$ pip install --upgrade pip
(sandbox_topology) ubuntu@ubuntu2004:~$ pip list
Package       Version
------------- -------
pip           22.0.4
pkg_resources 0.0.0
setuptools    44.0.0

(sandbox_topology) ubuntu@ubuntu2004:~$ deactivate
ubuntu@ubuntu2004:~$ 
```

Note that the python environment is activated, (sandbox_polyanagro) before the linux prompt. The environment can be deactivate using the **deactivate** command.

-------------------

## Clone the github repository and install Polyangro

-------------------
With the python virtual environment activated, clone and install requested python libraries

```bash
git clone https://github.com/jrdcasa/polyanagro.git
cd polyanagro
python -m pip install wheel
python -m pip install pygraphviz
python setup.py install 
```

If installation is correct you shouldubuntuu see the following message:

```test

(...)
Installed /home/ubuntu/Desktop/sandbox_polyangro/lib/python3.8/site-packages/topology-1.1-py3.8-linux-x86_64.egg
Processing dependencies for topology==1.1
Finished processing dependencies for topology==1.1

```

An **install.log** file has been generated with information about the installation process

-------------------

## Check topology

-------------------

```text
(sandbox_topology) ubuntu@ubuntu2004:~/Desktop/topology$ vi install.log 
(sandbox_topology) ubuntu@ubuntu2004:~/Desktop/topology$ 
(sandbox_topology) ubuntu@ubuntu2004:~/Desktop/topology$ python
Python 3.8.10 (default, Nov 26 2021, 20:14:08) 
[GCC 9.3.0] on linux
Type "help", "copyright", "credits" or "license" for more information.
>>> import topology
>>> 
```
