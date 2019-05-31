//
// Created by hossein on 5/30/19.
//

#include "atoms.h"
#include <fstream>
#include <iostream>
#include <sstream>

/* ----------------------------------------------------------------------
   setup for Atom
------------------------------------------------------------------------- */

Atom::Atom(double x, double y, double z, std::string element): x(x), y(y), z(z), element(element) {}

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */

AtomicConfiguration::AtomicConfiguration(): is_atom(false), is_cell(false) {}

void AtomicConfiguration::read_xyz(std::string filename)
{
    std::ifstream inFile(filename);
    if (!inFile) {
        std::cerr << "Unable to open " << filename;
        exit(1);   // call system to stop
    }
    // TODO: need improvement
    /* read number fo atoms */
    int nAtoms;
    inFile >> nAtoms;
    /* skip the second line */
    std::string line;
    inFile >> line;
    /* read atomic names and coordinates */
    for( int nline=0; nline<nAtoms; nline++)
    {
        std::getline(inFile, line);
        std::stringstream ss(line);

        double x, y, z;
        std::string element;
        ss >> element >> x >> y >> z;
        atoms.push_back( Atom(x, y, z, element) );
    }
    inFile.close();

    /* set atomic data is available */
    is_atom = true;
}

void AtomicConfiguration::set_cell(double cell[])
{
    for(int d=0; d<9; d++)
        this->cell[d] = cell[d];

    /* set cell data is available */
    is_cell = true;
}
