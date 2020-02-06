#ifndef CELL_H
#define CELL_H

#include<array>
#include"descriptor.h"
using std::array;

/// Helper class: allocation of memory for external fields in a cell
template<typename Descriptor>
class ExternalFieldArray {
public:
    double* get(int index) {
        return data+index;
    }
    double const* get(int index) const {
        return data+index;
    }
private:
    double data[Descriptor::ExternalField::numScalars];
};

template<class Descriptor>
class Cell {
public:
    /// Additional per-cell scalars for external fields
    typedef ExternalFieldArray<Descriptor> External;
public:
    /// Default constructor.
    Cell(){
        iniPara();
        iniExternal();
    }
public:
    /// Read-write access to distribution functions.
    double& operator[](int iPara) {
        return f[iPara];
    }
    /// Read-only access to distribution functions.
    double const& operator[](int iPara) const {
        return f[iPara];
    }
    /// Another way to get direct access to the f's, as in operator[]
    array<double,Descriptor::Cell::base>& getRawParameter() {
        return f;
    }
    /// Another way to get direct, const access to the f's, as in operator[]
    array<double,Descriptor::Cell::base> const& getRawParameter() const {
        return f;
    }

    Cell<Descriptor>& attributeF(Cell<Descriptor> const& rhs) {
        f = rhs.getRawParameter();
        return *this;
    };

    Cell<Descriptor>& attributeValues(Cell<Descriptor> const& rhs) {
        attributeF(rhs);
        for (int iExt=0; iExt < Descriptor::ExternalField::numScalars; ++iExt) {
            *external.get(iExt) = *rhs.external.get(iExt);
        }
        return *this;
    };

    /// Get a pointer to an external field
    double* getExternal(int offset) {
        return external.get(offset);
    }

    /// Get a const pointer to an external field
    double const* getExternal(int offset) const {
        return external.get(offset);
    }
private:
    void iniPara(){
        for (int iData = 0; iData<Descriptor::Cell::base; iData++){
            f[iData]=0;
        }
    }
    void iniExternal(){
        for (int iData=0; iData<Descriptor::ExternalField::numScalars; ++iData) {
            *external.get(iData) = 0;
        }
    }
private:
    array<double,Descriptor::Cell::base>        f;  ///< parameter
    External                             external;  ///< external scalars
};

#endif  // CELL_H