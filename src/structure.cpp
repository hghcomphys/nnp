//
// Atomic Configuration
//

#include "structure.h"
#include "logger.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

/* ----------------------------------------------------------------------
   setup for Atoms
------------------------------------------------------------------------- */
AtomicStructure::AtomicStructure(): isAtom(false), isCell(false), isTableOfDistances(false),
    numberOfAtoms(0), tableOfDistances(NULL), totalCharge(0), totalEnergy(0) {}

AtomicStructure::~AtomicStructure() 
{ 
    // free allocated memory for list of atoms (array of pointers)
    if ( isAtom ) {
        for (int i=0; i<numberOfAtoms; i++)
            delete atoms[i];
        delete[] atoms;
    }
    // free memory for table of distances (matrix of pointets)
    if ( isTableOfDistances ) {
        for (int i=0; i<numberOfAtoms; i++)
            delete[] tableOfDistances[i];
        delete[] tableOfDistances;
    }    
}

void AtomicStructure::writeFileFormatRunner(const std::string& filename)
{
    // open input structure file
    std::ofstream outFile(filename);
    if (!outFile) 
        throw std::runtime_error( (Log(ERROR) << "Unable to open file " << filename).toString() );

    outFile << "begin\n";
    // cell
    if ( isCell )
        for (int d=0; d<9; d+=3)
            outFile << "lattice " << cell[d] << " " << cell[d+1] << " " << cell[d+2] << "\n";
    else
        Log(WARN) << "No cell data available to wrote into " << filename;
    
    // atoms
    if ( isAtom )
    {
        for (int i=0; i<numberOfAtoms; i+=1) 
        {
            Atom *atom = atoms[i];
            outFile << "atom " 
                << atom->x 
                << " " << atom->y << " " << atom->z << " " << atom->element
                << " " << atom->charge << " " << atom->energy 
                << " " << atom->fx << " " << atom->fy << " " << atom->fz 
                << "\n";
        }
    }    
    else
        Log(ERROR) << "No atom data available to wrote into " << filename;
    
    // total energy and charge
    outFile << "energy " << totalEnergy << "\n";
    outFile << "charge " << totalCharge << "\n";
    outFile << "end\n";
    
    outFile.close();

    // log out info regarding number of atoms/elements/cell
    Log(INFO) << "Write structure into " << filename; 
}

void AtomicStructure::readFileFormatRuNNer(const std::string& filename)
{
    // open input structure file
    std::ifstream inFile(filename);
    if (!inFile) 
        throw std::runtime_error( (Log(ERROR) << "Unable to open file " << filename).toString() );

    // variables parsing input file
    std::string line, keyword;

    // pre-reading the structure file (e.g. number of atoms)
    while ( std::getline(inFile, line) ) 
    {
        std::stringstream ss(line);
        ss >> keyword;

        if (keyword == "atom")
        {
            // increase number of atom
            numberOfAtoms++;

            // read element
            double ddummy;
            std::string element;
            ss >> ddummy >> ddummy >> ddummy >> element;
            
            // increase number of atom for element
            numberOfAtomsForElement[element] += 1;
        }

        else if (keyword == "end")
            break; // read only the first frame
    }
    
    // allocate memory for list of atoms
    atoms = new Atom*[numberOfAtoms];

    // allocated memory for each element
    for (auto &it: numberOfAtomsForElement)
    {
        atomsForElement[it.first] = new Atom*[it.second];
        it.second = 0; // reset to zero, will be used as index
    }

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
            // check cell info is correct
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
            // read atom data
            double position[3], force[3];
            double charge, energy;
            std::string element;
            ss >> position[0] >> position[1] >> position[2] >> element >> charge 
                >> energy >> force[0] >> force[1] >> force[2];

            // create atom
            Atom *atom = new Atom(atomIndex, element.c_str(), position, force, charge, energy);

            // add atom to the list of atoms
            atoms[atomIndex++] = atom;

            // add atom to list of atoms for element
            atomsForElement[element][numberOfAtomsForElement[element]++] = atom;
        }
        else if (keyword == "energy") 
        {
            ss >> totalEnergy;
        }
        else if (keyword == "charge") 
        {
            ss >> totalCharge;
        }
        else if (keyword == "end")
        {
            break; // read only the first frame
        }
    }
    inFile.close();

    // set atomic and cell data are available
    isAtom = isCell = true;

    // log out info regarding number of atoms/elements/cell
    Log(INFO) << "Read " << filename << " (" << numberOfAtoms << " atoms)";
    for (auto & each: numberOfAtomsForElement)
        Log(INFO) << "Element " << each.first << ": " << each.second;
    Log(INFO) << "Cell (PBC)";
}

void AtomicStructure::readFileFormatRuNNer() 
{ 
    const char * filename = "input.data";
    readFileFormatRuNNer(filename); 
    Log(WARN) << "Read " << filename << " (as default RuNNer structure file)";
}

void AtomicStructure::setCell(const double cell[9])
{
    for(int d=0; d<9; d++)
        this->cell[d] = cell[d];
    isCell = true;
    Log(INFO) << "Set cell sizes";
}

void AtomicStructure::applyPBC(double &dx, double &dy, double &dz)
{
    // TODO: extend it to non-orthogonal box
    const double lx = cell[0];
    const double ly = cell[4];
    const double lz = cell[8];

    if ( isCell ) {
        // x-direction
        if ( dx > lx*0.5 ) dx -= lx;
        else if ( dx < -lx*0.5 ) dx += lx;  
        // y-direction
        if ( dy > ly*0.5 ) dy -= ly;
        else if ( dy < -ly*0.5 ) dy += ly;
        // z-direction
        if ( dz > lz*0.5 ) dz -= lz;
        else if ( dz < -lz*0.5 ) dz += lz;   
    }
}

double AtomicStructure::distance(Atom *atom_i, Atom *atom_j)
{
    double xij = atom_i->x - atom_j->x;
    double yij = atom_i->y - atom_j->y;
    double zij = atom_i->z - atom_j->z;
    applyPBC(xij, yij, zij);
    return  sqrt( xij*xij + yij*yij + zij*zij );
}

double AtomicStructure::distance(Atom *atom_i, Atom *atom_j, double drij[3])
{
    // TODO: do not repeat yourself
    double xij = atom_i->x - atom_j->x;
    double yij = atom_i->y - atom_j->y;
    double zij = atom_i->z - atom_j->z;
    applyPBC(xij, yij, zij);
    drij[0] = xij; drij[1] = yij; drij[2] = zij; 
    return  sqrt( xij*xij + yij*yij + zij*zij );
}

void AtomicStructure::calculateTableOfDistances(double globalCutOffRadius)
{
    // check atomic data is available
    if ( !isAtom )
        throw std::runtime_error( (Log(ERROR) << "No atomic structure available for table of distances").toString() );

    // check PBC is applied
    if ( !isCell )
        throw std::runtime_error( (Log(WARN) << "No PBC is applied for table of distances").toString() );

    // TODO: optimize memory usage for matrix (e.g. neighbor profiling)
    // allocate memory for table of distances
    tableOfDistances = new Distance*[numberOfAtoms];
    for (int i = 0; i < numberOfAtoms; ++i)
        tableOfDistances[i] = new Distance[numberOfAtoms];
    
    // loop over atom i and j
    for (int i=0; i<numberOfAtoms; i++) 
    {   
        for (int j=0; j<i; j++) 
        {
            // calculate distance between atoms i and j from atomic structure
            double drij[3];
            double rij = distance(atoms[i], atoms[j], drij);

            // skip calculation if it is outside the global cutoff radius 
            // if ( rij > globalCutOffRadius ) continue;

            // fill table of distances
            tableOfDistances[i][j].set(rij, drij);
            tableOfDistances[j][i].set(rij, drij, -1.0);
        }
    }

    // set table of distance is claculated
    isTableOfDistances = true;
}