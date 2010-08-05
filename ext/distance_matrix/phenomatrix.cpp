#include "phenomatrix.h"

#include "phenomatrix_base.cpp"

// id = predict species, given_id = source species
Phenomatrix::Phenomatrix(uint id, uint given_id, size_t min_genes = 2)
: PhenomatrixBase(id, id == given_id, min_genes),
  given_id_(given_id)
{
    // this is the part of the construction that is different between the
    // base and derived class.
    if (id != given_id) inherit_construct(min_genes);
}

Phenomatrix::Phenomatrix(const Phenomatrix& rhs)
: PhenomatrixBase(rhs), given_id_(rhs.given_id_)
{ }

// Copy and mask
Phenomatrix::Phenomatrix(const Phenomatrix& rhs, const id_set& mask_rows)
: PhenomatrixBase(rhs, mask_rows), given_id_(rhs.given_id_)
{ }