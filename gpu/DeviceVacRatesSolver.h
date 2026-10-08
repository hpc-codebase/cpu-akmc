#ifndef DEVICEVACRATESSOVLER_H
#define DEVICEVACRATESSOVLER_H

#include "../src/lattice/lattice_types.h"
#include "../src/type_define.h"

#define NN_TOTAL 126
#define RE_SCALE_SIZE 1.2
#define STREAM_SIZE 4

struct dev_Vacancy {
    // LatticeTypes type;
    _type_lattice_id id;
    double rates[8];
};

struct dev_nnLattice {
    LatticeTypes type;
    _type_lattice_id id;
};

struct dev_event {
    _type_lattice_id from_id;
    _type_lattice_id to_id;
    LatticeTypes to_type;
};

#endif /* DEVICEVACRATESSOVLER_H */