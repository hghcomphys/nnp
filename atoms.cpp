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

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */
Atoms::Atoms(): isAtom(false), isCell(false) {}

Atoms::~Atoms() { atoms.clear(); }

void Atoms::readXYZ(std::string filename)
{
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "ERROR: Unable to open file " << filename << "\n";
        exit(1);   // call system to stop
    }
    // TODO: need improvement
    // read number fo atoms
    int nAtoms;
    inFile >> nAtoms;
    // skip the second line
    std::string line;
    inFile >> line;
    // read atomic names and coordinates
    for( int nLine=0; nLine<nAtoms; nLine++)
    {
        std::getline(inFile, line);
        std::stringstream ss(line);
        double x, y, z;
        std::string element;
        ss >> element >> x >> y >> z;
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
    double xij = atom_i.x - atom_j.x;
    double yij = atom_i.y - atom_j.y;
    double zij = atom_i.z - atom_j.z;
    applyPBC(xij, yij, zij);
    return  sqrt( xij*xij + yij*yij + zij*zij );
}