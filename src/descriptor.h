#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

namespace descriptors {

struct NoExternalField {
    static const int numScalars = 0;
    static const int numSpecies = 0;
    static const int forceBeginsAt = 0;
    static const int sizeOfForce   = 0;
};

struct NoExternalFieldBase {
    typedef NoExternalField ExternalField;
};

struct CellField
{
    static const int base = 2;
};

struct CellBase{
    typedef CellField Cell;
};

struct CellDescriptor:public CellBase,public NoExternalFieldBase{
};

};

#endif //DESCRIPTOR_H