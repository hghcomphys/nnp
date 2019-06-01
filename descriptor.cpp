//
// Created by hossein on 6/1/19.
//

#include "descriptor.h"

/* ----------------------------------------------------------------------
   setup for Descriptor
------------------------------------------------------------------------- */

Descriptor::Descriptor() {}

Descriptor::~Descriptor() {
    /* free the allocated memory*/
    for (auto *descriptor: descriptors)
        delete descriptor;
    descriptors.clear();
}

void Descriptor::add(ACSF *descriptor) {
    descriptors.push_back(descriptor);
}

void Descriptor::calculate(AtomicConfiguration &configuration) {}