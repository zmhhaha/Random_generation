#ifndef DESCRIPTOR_H
#define DESCRIPTOR_H

namespace descriptors {

struct NoExternalField {
    static const int numScalars = 1;
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