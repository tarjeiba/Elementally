# Elementally
The finite element solver Elementally was written only to get a better grasp
of Python, it is not meant to solve any problems not already handled far more
throroughly elsewhere.

## Short list of useful Git commands.

`git add <FILENAME>` -- please use this to add each file specifically to your
current branch, as this encertains that no extraneous files are added to the
repository.

`git checkout <BRANCHNAME> <FILENAME>` -- this fetches *<FILENAME>* from
*<BRANCHNAME>* to your current branch.

The following command is from [GitHub Help](https://help.github.com/articles/
pushing-to-a-remote "GitHub"), and hasn't been tested here yet:

`git push <REMOTENAME> <BRANCHNAME>` -- to push your local branch to a remote
repository.

## Status
As of 13 August 2014, solely the 1D Poisson equation has been fully realized.
As of October 2014, a solver for the 2D Poisson equation has also been implemented.

## Software packages
Elementally leans heavily on the following packages:
* [Meshpy](http://mathema.tician.de/software/meshpy "MeshPy")
* [Scipy](http://www.scipy.org "Scipy")
* [Numpy](http://www.numpy.org "Numpy")

## Guidelines
Since this software is written as a means to learning Python, emphasis is laid
on restraining all code to follow the [PEP 8](legacy.python.org/dev/peps/
pep-0008 "Style Guide for Python Code")

## Theory

As the theory behind this is heavy on mathematical notation, we have separated
it into its own iPython Notebooks. The notebooks can be viewed using [nbviewer](http://nbviewer.ipython.org/github/tarjeiba/Elementally/tree/master/docs/).
