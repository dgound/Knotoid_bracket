# Knotoid invariants
This script computes the Jones polynomial for knotoids and its extensions, as were first introduced in [1]. Depending on whether the knotoid is on S^2
or in R^2, the script can output the Jones polynomial, the extended bracket, the Turaev loop bracket and the extended Turaev loop bracket polynomials.

## Dependencies
The script has been tested on Python 3.7 on OSX 10.14.5 and on Ubuntu Linux 16.04.6. The following python packages are required:

1. NetworkX
2. SymPy

To see the installed packages of your python distribution type <code>pip list</code>. If you use conda type <code>conda list</code>.
To install NetworkX using Pip type <code> pip install networkx </code>. If you use conda type <code> conda install networkx </code>.
To install SymPy using Pip type <code> pip install sympy </code>. If you use conda type <code> conda install sympy</code>.


## Usage
The script uses (extended) oriented Gauss codes in order to encode a knotoid diagram. It can evaluate either a single knotoid diagram or a multiple diagrams contained in a file.

To evaluate a single diagram type `python bracket.py -g GAUSSCODE`. The extended oriented Gauss Code consists of two or three parts, depending on 
whether the knotoid is considered in S^2 or in the plane respectively. The first and third parts are comma separated while the second part is just a sequence of signs.
For more on Gauss codes of knotoids please refer to [2].

Example for the evaluation of the Turaev loop bracket for a planar knotoid: `python bracket.py -g "1,-2,-1,2 ++ 0,2,3"`

Example for the evaluation of the Jones polynomial for a knotoid in S^2: `python bracket.py -g "1,-2,3,-1,2,-3 +++"`

To use the script in order to evaluate a file containing multiple diagrams we use the following command: `python bracket.py -f INPUTFILENAME`. Note that
the file must contain a single Gauss code in each line and that the Gauss code must follow the form mentioned above. That is, the first and, if the knotoid diagram is planar, the third part
must be comma separated. The results are saved to a file that, by default, is named `results.txt`. By adding the option `-o` in the command line the user can specify an output filename.
For example, `python bracket.py -f INPUTFILENAME -o OUTPUTFILENAME`

By adding the option `-e`, the script can compute the *extended bracket polynomial* for knotoids in S^2 and the *extended Turaev loop bracket polynomial* for planar knotoids.

Example for single knotoid evaluation: `python bracket.py -g GAUSSCODE-e`.

Eample for multiple knotoids evaluation: `python bracket.py -f INPUTFILENAME -o OUTPUTFILENAME -e`.


## Credits

### Contributors
Dimos Goundaroulis

### Institutions
Center for Integrative Genomics, University of Lausanne.

SIB Swiss Institute of Bioinformatics.


## References
[1] V. Turaev, *Knotoids*. Osaka J. Math. **1** no 49 (2012), pp. 195-223. [Link to paper](https://projecteuclid.org/download/pdf_1/euclid.ojm/1332337244)
 
[2] D. Goundaroulis, J. Dorier and A. Stasiak, *A systematic classification of knotoids on the plane and on the sphere*. arXiv preprint [arXiv:1902.07277](https://arxiv.org/abs/1902.07277) (2019).
