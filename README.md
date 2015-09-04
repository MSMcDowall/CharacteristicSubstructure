#Characteristic Substructure Application

Created by Mary McDowall as part of MSc IT degree.

##Description:

    Given a file (.txt or .smi) containing the SMILES strings of a set of molecules, the characteristic substructure of
    the molecules will be found. The characteristic substructure is a composite structure containing subgraphs which
    appear frequently in the molecules. The method is defined in 'Finding Characteristic Substructures for 
    Metabolite Classes' by Marcus Ludwig et al. German Conference on Bioinformatics 2012.
     
    It is also possible to retrieve a list of the subgraphs which appear most frequently within all the molecules 
    (called representative substructures). 
    
    The frequency threshold below sets the percentage of molecules that a subgraph must appear in for it be considered
    representative. The default is set at 80% of the molecules. If the thresehol is set lower then structures which 
    appear in a lower number of molecules will be included in the characteristic substructure or in the list of 
    representative structures.
    
##Input

    The file containing the SMILES strings (either .txt or .smi) is specified in the command line arguments.

    If the user only wants the characteristic substructure then the "-c" flag can be used.
    If the user only wants the representative substructures then the "-r" flag can be used.
    If the user wishes to view the command line argument help statements then the "-h" flag can be used.
    If the user wishes to alter the frequency threshold then enter the new threshold as a decimal after the file name.
    Example: python source\algorithm.py -c data\folder\FileName.txt 0.6
          or python source\algorithm.py -r data\folder\FileName.txt
          or python source\algorithm.py -h

##Output

    The results file names will consist of the input file (without the extension), the relative frequency threshold
    and their results type and be located in the same directory as the input file.
    Example: FileName0.6_CharacteristicSubstructure.txt
          or FileName0.8_RepresentativeStructures.txt

    The results in the CharacteristicSubstructure file will include( all graphs in adjacency dictionary format):

        Characteristic substructure
        1. Structure contained in characteristic substructure
        2. Structure contained in characteristic substructure
        ...
        SMILES for input molecules  :   Membership array for structures in molecules 
        
    For the membership array if a 0 is present at position 2 then Structure 2 is not present in the molecule.
    If a 1 is present then the structure can be found in the molecule. 
    
    The results in the RepresentativeStructures file will include( all graphs in adjacency dictionary format):

        1. Representative structure
        2. Representative structure
        ...
        SMILES for input molecules  :   Membership array for structures in molecules
    
    
    
##Dependencies

    NetworkX is required to run the algorithm.py file as it provides the subgraph isomorphism algorithm
    MatplotLib is required if the user wishes to draw the molecule though it is not required to run the algorithm.py
    
##Contents 

###source Folder:

    algorithm.py contains an implementation of the method set down in 'Finding Characteristic Substructures for 
    Metabolite Classes' by Marcus Ludwig et al. German Conference on Bioinformatics 2012
    
    graph.py contains a graph data structure implemented as an adjacency dictionary with Graph, Vertex and Edge classes
    
    molecule.py contains Molecule, Atom and Bond classes which inherits from graph.py classes 
        
    smiles_parser.py contains a basic SMILES parser which tokenizes a SMILES string and creates a Molecule object
    
    draw_molecule.py uses the NetworkX package and Matplotlib to draw the molecules including vertex and edge labels
    
###tests Folder:

    Unit Tests:
        unit tests on the graph.py and molecule.py files as well as the basic functions in algorithm.py
        
    Feature Tests:
        feature tests that ensure smiles_parser.py and algorithm.py produce the correct results
        

