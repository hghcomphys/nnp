//
// Atomic Configuration
//

#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include "atoms.h"
#include "constants.h"

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(double x, double y, double z, std::string element, int index): x(x), y(y), z(z), element(element), index(index) {}

double Atom::getX() { return x; }

double Atom::getY() { return y; }

double Atom::getZ() { return z; }

double Atom::getIndex() { return index; }

std::string Atom::getElement() { return element; }

int Atom::getAtomicNumber(const std::string& element) 
{
    // TODO: extend to all elements
    if(element == "H")
        return 1;
    else if (element == "C")
        return 6; 
    else if (element == "O")
        return 8;
}

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */
Atoms::Atoms(): isAtom(false), isCell(false), atomIndex(0) {}

Atoms::~Atoms() { listOfAtoms.clear(); }

void Atoms::addAtom(const Atom& atom) {
    listOfAtoms.push_back( atom );
    atomIndex++;
}

std::vector<Atom>& Atoms::getListOfAtoms() { return listOfAtoms; }

int Atoms::getNumberOfAtoms() { return listOfAtoms.size(); }

int Atoms::getNumberOfAtomsForElement(const std::string& element) { return getListOfIndexForElement(element).size(); }

std::stringstream readLineToStringStream(std::ifstream& inFile) {
    std::string line;
    std::getline(inFile, line);
    std::stringstream ss(line);
    return ss;
}

void Atoms::readFileFormatXYZ(const std::string& filename)
{
    std::ifstream inFile(filename);
    if (!inFile) {
        std::string ss; 
        ss = "Unable to open file " + filename;
        throw std::runtime_error(ss);
    }
    // TODO: need improvement
    // read number of atoms
    int nAtoms;
    readLineToStringStream(inFile) >> nAtoms;
    // skip the second line
    readLineToStringStream(inFile);
    // read atomic names and coordinates
    for( int nLine=0; nLine<nAtoms; nLine++)
    {
        double x, y, z;
        std::string element;
        readLineToStringStream(inFile) >> element >> x >> y >> z;
        addAtom( Atom(x*ANGSTROM_TO_BOHR, y*ANGSTROM_TO_BOHR, z*ANGSTROM_TO_BOHR, element, atomIndex) );
    }
    inFile.close();

    // set atomic data is available
    isAtom = true;
}

void Atoms::setCell(double cell[])
{
    for(int d=0; d<9; d++)
        this->cell[d] = cell[d]*ANGSTROM_TO_BOHR;

    // set cell data is available
    isCell = true;
}

void Atoms::applyPBC(double &dx, double &dy, double &dz)
{
    // TODO: extend it to non-orthogonal box
    const double lx = cell[1];
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

double Atoms::distance(Atom &atom_i, Atom &atom_j)
{
    double xij = atom_i.getX() - atom_j.getX();
    double yij = atom_i.getY() - atom_j.getY();
    double zij = atom_i.getZ() - atom_j.getZ();
    applyPBC(xij, yij, zij);
    return  sqrt( xij*xij + yij*yij + zij*zij );
}

double Atoms::distance(Atom &atom_i, Atom &atom_j, double drij[3])
{
    double xij = atom_i.getX() - atom_j.getX();
    double yij = atom_i.getY() - atom_j.getY();
    double zij = atom_i.getZ() - atom_j.getZ();
    applyPBC(xij, yij, zij);
    drij[0] = xij; drij[1] = yij; drij[2] = zij; 
    return  sqrt( xij*xij + yij*yij + zij*zij );
}

std::vector<int> Atoms::getListOfIndexForElement(const std::string &element)
{
    std::vector<int> listOfindexForElement;
    for (Atom &atom: listOfAtoms) {
        if( atom.getElement() == element )
            listOfindexForElement.push_back( atom.getIndex() );
    }

    if ( listOfindexForElement.size() == 0 )
        throw std::runtime_error("Cannot find the element in list of atoms");

    return listOfindexForElement;
}

// const Atom& Atoms::operator [] (unsigned int i) const { return listOfAtoms[i]; }

bool Atoms::isPBC() { return isCell; } 

void Atoms::readFileFormatRuNNer(const std::string& filename)
{
    std::ifstream inFile(filename);
    if (!inFile) {
        std::string ss; 
        ss = "Unable to open file " + filename;
        throw std::runtime_error(ss);
    }

    std::string line, keyword;
    int cellIndex = 0;
    while ( std::getline(inFile, line) ) {
        std::stringstream ss(line);
        ss >> keyword;

        if (keyword == "lattice") {
            if (cellIndex == 9)
                    throw std::runtime_error("Unexpected number of data for cell");
            for (int i=0; i<3; i++)
                ss >> cell[cellIndex++];
        }
        else if (keyword == "atom") {
            double x, y, z;
            std::string element;
            ss >> x >> y >> z >> element;
            addAtom( Atom(x, y, z, element, atomIndex) );
        }
        else if (keyword == "end")
            break; // read only the first frame
    }
    inFile.close();

    // set atomic data is available
    isAtom = true;
    isCell = true;
}

void Atoms::readFileFormatRuNNer() { readFileFormatRuNNer("input.data"); }