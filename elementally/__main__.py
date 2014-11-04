import sys

import solvers

def main(args=None):
    """The main routine."""
    if args is None:
        args = sys.argv[1:]

    if len(sys.argv) <= 1:
        print "Please specify whether you wan't to run the 1d or 2d Poisson"
        print "solver by. Choose solver by"
        print "    python elementally poisson_1d|poisson_2d"
        return 1


    if sys.argv[1] is not None and sys.argv[1] == 'poisson_1d':
        import solvers.poisson_1d
        return 0
    elif sys.argv[1] is not None and sys.argv[1] == 'poisson_2d':
        import solvers.poisson_2d
        return 0
    else:
        print "Choose solver by"
        print "    python elementally poisson_1d|poisson_2d"
        return 1
    

    # Do argument parsing here (e.g., with argparse) and anything else
    # you want the project to do.

if __name__ == "__main__":
    main()
