# F2PY augmented Machine Learning (FML)

FML is a Python library to efficiently carry out machine learning calculation. Everything is accessed via a simple Python interface, while efficient, parallel code (written in F90 with OpenMP) is executed behind the scenes.


## 1) Installation

By default FML compiles with GCC's gfortran and your system's standard BLAS+LAPACK libraries (-lblas -llapack). I recommend you switch to MKL for the math libraries, but the difference between ifort and gfortran shouldn't be to significant. BLAS and LAPACK implementations are the only libraries required for FML. Additionally you need a functional Python 2.x interpreter and NumPy (with F2PY) installed. Most Linux systems can install BLAS and LAPACK like this:

    sudo apt-get install libblas-dev liblapack-dev

Ok, on to the installation instructions:

1.1) First you clone this repository: 

    git clone https://github.com/andersx/fml.git

1.2) Then you simply compile by typing make in the fml folder:

    make

1.3) To make everything accessible to your Python export the fml root-folder to your PYTHONPATH. E.g. this is what I export:

    export PYTHONPATH=/home/andersx/dev/fml:$PYTHONPATH

## 2) How to use:


## 2.1) Generate descriptors from XYZ files

FML is very low-level - here is a short example to generate descriptors for your molecules via the `fml.Molecule` class. It is convenient to store a list of `fml.Molecule`s in `cPickle` format.

    import fml
    import cPickle

    ...
   
    # Read some properties 
    properties = reader("properties.csv")

    # Get list of filenames
    filenames = sorted(os.listdir(xyz_dir))


    # Max number of atoms in molecules
    nmax = 23

    mols = []

    # Generate filenames
    for filename in filenames:

        # Initialize the molecule
        mol = Molecule()
        
        # Read atoms from XYZ
        mol.read_xyz(xyz_dir + filename)
        
        # Generate and store an ARAD descriptor
        mol.generate_arad_descriptor(size=nmax)
        
        # Store some associated property with the molecule
        mol.properties = properties[filename]
        
        mols.append(mol)

    # Save the molecules for later use
    with open("mols.cpickle", "wb") as f:
        cPickle.dump(mols, f, protocol=2)
