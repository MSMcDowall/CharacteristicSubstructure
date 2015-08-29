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
    representative. The default is set at 80% of the molecules. 
    
##Input/Output

    The SMILES data can either be copied into the SMILES.txt file within the directory (the default input file)
    or an alternative can be specified as a command line argument.
    
    The results file names will consist of the input file (without the extension) and their results type.
    Example: SMILES_CharacteristicSubstructure.txt or SMILES_RepresentativeStructures.txt 

###To find both the characteristic substructure and the representative substructures ([] optional):

    python source/algorithm.py [SMILES file name] [Frequency threshold in [0,1], default:0.8]
    
###To find just the characteristic substructure ([] optional):

    python source/algorithm.py [SMILES_file_name] 0 [Frequency threshold in [0,1], default:0.8]
    
###To find just the representative substructures ([] optional):

    python source/algorithm.py [SMILES_file_name] 1 [Frequency threshold in [0,1], default:0.8]
    
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
        

