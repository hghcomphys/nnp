//
// Atomic Descriptor
//

#ifndef NNP_DESCRIPTOR_H
#define NNP_DESCRIPTOR_H

#include "acsf.h"
#include "atoms.h"

//TODO: Descriptor class as template
class Descriptor {
public:
    std::vector<ACSF *> descriptors; /* factory method */
    Descriptor();
    ~Descriptor();
    void add(ACSF *descriptor);
    void calculate(AtomicConfiguration &configuration);
};


#endif //NNP_DESCRIPTOR_H
