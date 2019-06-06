//
// Atomic Configuration
//

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include "atoms.h"

#define ANGSTROM_TO_BOHR 1.88973

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */
Atom::Atom(double x, double y, double z, std::string element, int index): x(x), y(y), z(z), element(element), index(index) {}

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */
Atoms::Atoms(): is_atom(false), is_cell(false) {}

Atoms::~Atoms() {
    atoms.clear();
}

void Atoms::read_xyz(std::string filename)
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
    for( int nline=0; nline<nAtoms; nline++)
    {
        std::getline(inFile, line);
        std::stringstream ss(line);
        double x, y, z;
        std::string element;
        ss >> element >> x >> y >> z;
        atoms.push_back( Atom(x*ANGSTROM_TO_BOHR, y*ANGSTROM_TO_BOHR, z*ANGSTROM_TO_BOHR, element, nline+1) );
    }
    inFile.close();

    // set atomic data is available
    is_atom = true;
}

void Atoms::set_cell(double cell[]) {
    for(int d=0; d<9; d++)
        this->cell[d] = cell[d];

    // set cell data is available
    is_cell = true;
}

void Atoms::apply_pbc(double &dx, double &dy, double &dz)
{
    // TODO: extend it to non-orthogonal box
    const double lx = cell[1];
    const double ly = cell[4];
    const double lz = cell[8];

    if ( is_cell ) {
        if ( dx > lx*0.5 ) dx -= lx;
        else if ( dx < -lx*0.5 ) dx += lx;  

        if ( dy > ly*0.5 ) dy -= ly;
        else if ( dy < -ly*0.5 ) dy += ly;

        if ( dz > lz*0.5 ) dz -= lz;
        else if ( dz < -lz*0.5 ) dz += lz;   
    }
}

double Atoms::distance(Atom &atom_i, Atom &atom_j) {
    // TODO: apply pbc
    double xij = atom_i.x - atom_j.x;
    double yij = atom_i.y - atom_j.y;
    double zij = atom_i.z - atom_j.z;
    apply_pbc(xij, yij, zij);
    return  sqrt( xij*xij + yij*yij + zij*zij );
}