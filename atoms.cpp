//
// Atomic Configuration
//

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include "atoms.h"
#include "constants.h"

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(double x, double y, double z, std::string element, int index): x(x), y(y), z(z), element(element), index(index) {}

double Atom::getX() { return x; }

double Atom::getY() { return y; }

double Atom::getZ() { return y; }

double Atom::getIndex() { return index; }

std::string Atom::getElement() { return element; }

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */
Atoms::Atoms(): isAtom(false), isCell(false) {}

Atoms::~Atoms() { atoms.clear(); }

std::vector<Atom>& Atoms::getListOfAtoms() { return atoms; }

int Atoms::getNumberOfAtoms() { return atoms.size(); }

std::stringstream readLineToStringStream(std::ifstream& inFile) {
    std::string line;
    std::getline(inFile, line);
    std::stringstream sline(line);
    return sline;
}

void Atoms::readXYZ(std::string filename)
{
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "ERROR: Unable to open file " << filename << "\n";
        exit(1);   // call system to stop
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
        atoms.push_back( Atom(x*ANGSTROM_TO_BOHR, y*ANGSTROM_TO_BOHR, z*ANGSTROM_TO_BOHR, element, nLine+1) );
    }
    inFile.close();

    // set atomic data is available
    isAtom = true;
}

void Atoms::setCell(double cell[])
{
    for(int d=0; d<9; d++)
        this->cell[d] = cell[d];

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

std::vector<int> Atoms::getListOfIndexForElement(const std::string &element)
{
    std::vector<int> listOfindexForElement;
    for (Atom &atom: atoms) {
        if( atom.getElement() == element )
            listOfindexForElement.push_back( atom.getIndex() );
    }

    if ( listOfindexForElement.size() == 0 )
        throw std::runtime_error("Cannot find the element in list of atoms");

    return listOfindexForElement;
}