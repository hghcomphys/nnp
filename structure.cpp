//
// Atomic Configuration
//

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "structure.h"
#include "log.h"

const double ANGSTROM_TO_BOHR = 1.88973;

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */
AtomicStructure::AtomicStructure(): isAtom(false), isCell(false), atomIndex(0) {}

AtomicStructure::~AtomicStructure() { 
    for (auto atom: listOfAtoms)
        delete atom;
    listOfAtoms.clear(); 
}

void AtomicStructure::addAtom(Atom* atom) {
    listOfAtoms.push_back( atom );
    atomIndex++;
}

// std::vector<Atom*> AtomicStructure::getListOfAtoms() { return listOfAtoms; }

int AtomicStructure::getNumberOfAtoms() { return listOfAtoms.size(); }

int AtomicStructure::getNumberOfAtomsForElement(const std::string& element) { return getListOfAtomIndexForElement(element).size(); }

std::stringstream readLineToStringStream(std::ifstream& inFile) {
    std::string line;
    std::getline(inFile, line);
    std::stringstream ss(line);
    return ss;
}

void AtomicStructure::readFileFormatXYZ(const std::string& filename)
{
    std::ifstream inFile(filename);
    if (!inFile)
        throw std::runtime_error( (Log(ERROR) << "Unable to open file " + filename).toString() );

    // TODO: need improvement
    // read number of atoms
    int nAtoms;
    readLineToStringStream(inFile) >> nAtoms;
    // skip the second line
    readLineToStringStream(inFile);
    // read atomic names and coordinates
    for( int nLine=0; nLine<nAtoms; nLine++)     {
        double x, y, z;
        std::string element;
        readLineToStringStream(inFile) >> element >> x >> y >> z;
        addAtom( new Atom(x*ANGSTROM_TO_BOHR, y*ANGSTROM_TO_BOHR, z*ANGSTROM_TO_BOHR, element, atomIndex) );
    }
    inFile.close();

    // set atomic data is available
    isAtom = true;

    // report number of atoms
    Log(INFO) << "Read " + filename << " (" << getNumberOfAtoms() << " atoms)";
}

void AtomicStructure::setCell(double cell[9])
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

std::vector<int> AtomicStructure::getListOfAtomIndexForElement(const std::string &element)
{
    std::vector<int> listOfAtomindexForElement;
    for (auto it: listOfAtoms) {
        if( it->element == element )
            listOfAtomindexForElement.push_back( it->index );
    }

    if ( listOfAtomindexForElement.size() == 0 )
        throw std::runtime_error( (Log(ERROR) << "Cannot find the element in list of atoms").toString());

    return listOfAtomindexForElement;
}

std::vector<int> AtomicStructure::getListOfAtomIndex()
{
    std::vector<int> listOfAtomIndex;
    for (auto it: listOfAtoms)
            listOfAtomIndex.push_back( it->index );

    return listOfAtomIndex;
}

// const Atom& Atoms::operator [] (unsigned int i) const { return listOfAtoms[i]; }

bool AtomicStructure::isPBC() const { return isCell; } 

void AtomicStructure::readFileFormatRuNNer(const std::string& filename)
{
    std::ifstream inFile(filename);
    if (!inFile) 
        throw std::runtime_error( (Log(ERROR) << "Unable to open file " + filename).toString() );

    std::string line, keyword;
    int cellIndex = 0;
    while ( std::getline(inFile, line) ) {
        std::stringstream ss(line);
        ss >> keyword;

        if (keyword == "lattice") {
            if (cellIndex == 9)
                    throw std::runtime_error(
                        (Log(ERROR) << "Unexpected number of data for cell").toString()
                    );
            
            for (int i=0; i<3; i++)
                ss >> cell[cellIndex++];
        }
        else if (keyword == "atom") {
            double x, y, z;
            std::string element;
            double ddummy, fx, fy, fz;
            ss >> x >> y >> z >> element >> ddummy >> ddummy >> fx >> fy >> fz;
            addAtom( new Atom(x, y, z, element, atomIndex, fx, fy, fz) );
        }
        else if (keyword == "end")
            break; // read only the first frame
    }
    inFile.close();

    // set atomic and cell data are available
    isAtom = true;
    isCell = true;

    // report number of atoms
    Log(INFO) << "Read " + filename << " (" << getNumberOfAtoms() << " atoms)";
    Log(INFO) << "Cell (PBC)";
}

void AtomicStructure::readFileFormatRuNNer() { 
    const std::string filename = "input.data";
    readFileFormatRuNNer(filename); 
    Log(WARN) << "Read " + filename + " (as default RuNNer structure file)";
}