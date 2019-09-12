//
// Atomic Configuration
//

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "structure.h"
#include "logger.h"

const double ANGSTROM_TO_BOHR = 1.88973;

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */
AtomicStructure::AtomicStructure(): isAtom(false), isCell(false), numberOfAtoms(0), 
    listOfAtoms(NULL), tableOfDistances(NULL) {}

AtomicStructure::~AtomicStructure() 
{ 
    // free memory for list of atoms
    if (listOfAtoms!=NULL) {
        for (int i=0; i<numberOfAtoms; i++)
            delete listOfAtoms[i];
        delete[] listOfAtoms; 
    }

    // free memory for table of distances
    if (tableOfDistances!=NULL) {
        for (int i=0; i<numberOfAtoms; i++)
            delete[] tableOfDistances[i];
        Log(INFO) << "terminated!";
        delete[] tableOfDistances;
    }    
}

int AtomicStructure::getNumberOfAtoms() 
{ 
    return numberOfAtoms; 
}

int AtomicStructure::getNumberOfAtomsForElement(const char * element) 
{ 
    return getListOfAtomIndexForElement(element).size(); 
}

std::stringstream readLineToStringStream(std::ifstream& inFile) 
{
    std::string line;
    std::getline(inFile, line);
    std::stringstream ss(line);
    return ss;
}

void AtomicStructure::readFileFormatRuNNer(const char * filename)
{
    std::ifstream inFile(filename);
    if (!inFile) 
        throw std::runtime_error( (Log(ERROR) << "Unable to open file " << filename).toString() );

    std::string line, keyword;

    // pre-reading the structure file
    while ( std::getline(inFile, line) ) 
    {
        std::stringstream ss(line);
        ss >> keyword;

        if (keyword == "atom")
           numberOfAtoms++;

        else if (keyword == "end")
            break; // read only the first frame
    }
    
    // allocate memory for list of atoms
    listOfAtoms = new Atom*[numberOfAtoms];

    // go back to begining of the file
    inFile.clear();
    inFile.seekg(0, std::ios::beg);

    // read structure file
    int atomIndex = 0;
    int cellIndex = 0;
    while ( std::getline(inFile, line) ) 
    {
        std::stringstream ss(line);
        ss >> keyword;

        if (keyword == "lattice") 
        {
            if (cellIndex == 9)
                    throw std::runtime_error(
                        (Log(ERROR) << "Unexpected number of data for cell").toString()
                    );
            // read cell parameters
            for (int i=0; i<3; i++)
                ss >> cell[cellIndex++];
        }
        else if (keyword == "atom") 
        {
            double position[3], force[3], ddummy;
            std::string element;
            ss >> position[0] >> position[1] >> position[2] >> element >> ddummy 
                >> ddummy >> force[0] >> force[1] >> force[2];
            listOfAtoms[atomIndex] =  new Atom(atomIndex, element.c_str(), position, force);
            atomIndex++;
        }
        else if (keyword == "end")
            break; // read only the first frame
    }
    inFile.close();

    // set atomic and cell data are available
    isAtom = true;
    isCell = true;

    // report number of atoms
    Log(INFO) << "Read " << filename << " (" << getNumberOfAtoms() << " atoms)";
    Log(INFO) << "Cell (PBC)";
}

void AtomicStructure::readFileFormatXYZ(const char * filename)
{
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error( (Log(ERROR) << "Unable to open file " << filename).toString() );

    // read number of atoms
    readLineToStringStream(inFile) >> numberOfAtoms;

    // allocate memory for list of atoms
    listOfAtoms = new Atom*[numberOfAtoms];

    // skip the second line
    readLineToStringStream(inFile);

    // read atomic names and coordinates
    int atomIndex = 0;
    for (int nLine=0; nLine<numberOfAtoms; nLine++) 
    {
        double position[3];
        std::string element;
        readLineToStringStream(inFile) >> element >> position[0] >> position[1] >> position[2];
        for (int d=0; d<3; d++)
            position[d]*=ANGSTROM_TO_BOHR;
        listOfAtoms[atomIndex] = new Atom(atomIndex, element.c_str(), position);
        atomIndex++;
    }
    inFile.close();

    // set atomic data is available
    isAtom = true;

    // report number of atoms
    Log(INFO) << "Read " << filename << " (" << getNumberOfAtoms() << " atoms)";
}

void AtomicStructure::setCell(const double cell[9])
{
    for(int d=0; d<9; d++)
        this->cell[d] = cell[d]*ANGSTROM_TO_BOHR;

    // set cell data is available
    isCell = true;

    // report number of atoms
    Log(INFO) << "Set cell sizes";
}

void AtomicStructure::applyPBC(double &dx, double &dy, double &dz)
{
    // TODO: extend it to non-orthogonal box
    const double lx = cell[0];
    const double ly = cell[4];
    const double lz = cell[8];

    if ( isCell )
    {
        if ( dx > lx*0.5 ) dx -= lx;
        else if ( dx < -lx*0.5 ) dx += lx;  

        if ( dy > ly*0.5 ) dy -= ly;
        else if ( dy < -ly*0.5 ) dy += ly;

        if ( dz > lz*0.5 ) dz -= lz;
        else if ( dz < -lz*0.5 ) dz += lz;   
    }
}

double AtomicStructure::distance(Atom &atom_i, Atom &atom_j)
{
    double xij = atom_i.x - atom_j.x;
    double yij = atom_i.y - atom_j.y;
    double zij = atom_i.z - atom_j.z;
    applyPBC(xij, yij, zij);
    return  sqrt( xij*xij + yij*yij + zij*zij );
}

double AtomicStructure::distance(Atom &atom_i, Atom &atom_j, double drij[3])
{
    double xij = atom_i.x - atom_j.x;
    double yij = atom_i.y - atom_j.y;
    double zij = atom_i.z - atom_j.z;;
    applyPBC(xij, yij, zij);
    drij[0] = xij; drij[1] = yij; drij[2] = zij; 
    return  sqrt( xij*xij + yij*yij + zij*zij );
}

std::vector<int> AtomicStructure::getListOfAtomIndexForElement(const char * element)
{
    std::vector<int> listOfAtomindexForElement;
    for (int i=0; i<numberOfAtoms; i++) {
        if( listOfAtoms[i]->isElement(element) )
            listOfAtomindexForElement.push_back( listOfAtoms[i]->index );
    }

    if ( listOfAtomindexForElement.size() == 0 )
        throw std::runtime_error( (Log(ERROR) << "Cannot find the element in list of atoms").toString());

    return listOfAtomindexForElement;
}

std::vector<int> AtomicStructure::getListOfAtomIndex()
{
    std::vector<int> listOfAtomIndex;
    for (int i=0; i<numberOfAtoms; i++)
            listOfAtomIndex.push_back( listOfAtoms[i]->index );

    return listOfAtomIndex;
}

// const Atom& Atoms::operator [] (unsigned int i) const { return listOfAtoms[i]; }

bool AtomicStructure::isPBC() 
{ 
    return isCell; 
} 

void AtomicStructure::readFileFormatRuNNer() 
{ 
    const char * filename = "input.data";
    readFileFormatRuNNer(filename); 
    Log(WARN) << "Read " << filename << " (as default RuNNer structure file)";
}

void AtomicStructure::calculateTableOfDistances()
{
    if ( !isAtom )
        throw std::runtime_error( (Log(ERROR) << "No atomic structure available for table of distances").toString() );

    if ( !isCell )
        throw std::runtime_error( (Log(WARN) << "No PBC is applied for table of distances").toString() );

    // get number of atoms in atomic structure
    const int size = numberOfAtoms;
    if ( size==0 )
        throw std::runtime_error( 
            (Log(ERROR) << "Unexpected size for table of distances (" << size << ")").toString() 
            );

    // allocate matrix
    tableOfDistances = new Distance*[size];
    for (int i = 0; i < numberOfAtoms; ++i)
        tableOfDistances[i] = new Distance[size];
    
    // TODO: optimize memory usage for matrix
    // loop over atom i and j
    for (int i=0; i<size; i++) 
    {
        Atom& atom_i = getAtom(i);
        
        for (int j=0; j<i; j++) 
        {
            Atom& atom_j = getAtom(j);

            // calculate distance between atoms i and j from atomic structure
            double drij[3];
            double rij = distance(atom_i, atom_j, drij);

            // fill table of distances
            tableOfDistances[i][j].set(rij, drij);
            tableOfDistances[j][i].set(rij, drij, -1.0);
        }
    }
}

