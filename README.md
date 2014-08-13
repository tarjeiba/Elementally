# Elementally
The finite element solver ELementally was merely written to get a better grasp
of Python, it is not meant to solve any problems not already handled far more
throroughly elsewhere.

## Short list of useful Git commands.

`git add <FILENAME>` -- please use this to add each file specifically to your
current branch, as this encertains that no extraneous files are added to the
repository.

`git checkout <BRANCHNAME> <FILENAME>` -- this fetches *<FILENAME>* from *<BRANCHNAME>* to your current branch.

The following command is from [GitHub Help](https://help.github.com/articles/pushing-to-a-remote "GitHub"), and hasn't been tested here yet:

`git push <REMOTENAME> <BRANCHNAME>` -- to push your local branch to a remote repository.

## Status
As of 13 August 2014, solely the 1D Poisson equation has been fully realized.

## Software packages
Elementally leans heavily on the following packages:
* [Meshpy](mathema.tician.de/software/meshpy "MeshPy")
* [Scipy](www.scipy.org "Scipy")
* [Numpy](www.numpy.org "Numpy")

## Guidelines
Since this software is written as a means to learning Python, emphasis is laid
on restraining all code to follow the [PEP 8](legacy.python.org/dev/peps/
pep-0008 "Style Guide for Python Code")

## Theory
This brief theory section constitutes mainly of [Numerical Models for
Differential Problems](http://www.amazon.co.uk/Numerical-Models-Differential-
Problems-MS-ebook/dp/B00JY5G3T4/ref=sr_1_2?ie=UTF8&qid=1407928272&sr=8-2&
keywords=quarteroni "Amazon.co.uk"). Page references are made where applicable.

### The Poisson Equation
*(Note that Github's markdown parser does not support equations, as such are
considered an unsafe feature due to their execution of Javascripts.)*


#### Implementation in 1D

HOLD
