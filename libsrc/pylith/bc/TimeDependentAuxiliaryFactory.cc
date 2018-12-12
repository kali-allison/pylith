// -*- C++ -*-
//
// ----------------------------------------------------------------------
//
// Brad T. Aagaard, U.S. Geological Survey
// Charles A. Williams, GNS Science
// Matthew G. Knepley, University of Chicago
//
// This code was developed as part of the Computational Infrastructure
// for Geodynamics (http://geodynamics.org).
//
// Copyright (c) 2010-2017 University of California, Davis
//
// See COPYING for license information.
//
// ----------------------------------------------------------------------
//

#include <portinfo>

#include "TimeDependentAuxiliaryFactory.hh" // implementation of object methods

#include "pylith/topology/Field.hh" // HOLDSA AuxiliaryField
#include "pylith/topology/FieldQuery.hh" // USES FieldQuery
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "spatialdata/spatialdb/TimeHistory.hh" // USES TimeHistory
#include "spatialdata/units/Nondimensional.hh" // USES Nondimensional

#include "pylith/utils/journals.hh" // USES PYLITH_JOURNAL*

#include <cassert>

// ---------------------------------------------------------------------------------------------------------------------
namespace pylith {
    namespace bc {
        class _TimeDependentAuxiliaryFactory {
public:

            ///< Names of field components in XYZ coordinate system.
            static const char* componentsXYZ[3];

            ///< Names of field components in 2-D tangential/normal coordinate system.
            static const char* componentsTN[2];

            ///< Names of field components in 3-D tangential/normal coordinate system.
            static const char* componentsTTN[3];

        };

        // _TimeDependentAuxiliaryFactory

        const char* pylith::bc::_TimeDependentAuxiliaryFactory::componentsXYZ[3] = { "_x", "_y", "_z" };
        const char* pylith::bc::_TimeDependentAuxiliaryFactory::componentsTN[2] = { "_tangential", "_normal" };
        const char* pylith::bc::_TimeDependentAuxiliaryFactory::componentsTTN[3] = { "_tangential_1", "_tangential_2", "_normal" };
    } // bc
} // pylith

// ---------------------------------------------------------------------------------------------------------------------
// Default constructor.
pylith::bc::TimeDependentAuxiliaryFactory::TimeDependentAuxiliaryFactory(const ReferenceEnum reference) :
    _auxComponents(reference) { // constructor
    GenericComponent::setName("timedependentauxiliaryfactory");
} // constructor


// ---------------------------------------------------------------------------------------------------------------------
// Destructor.
pylith::bc::TimeDependentAuxiliaryFactory::~TimeDependentAuxiliaryFactory(void) {}


// ---------------------------------------------------------------------------------------------------------------------
// Set names of vector components in auxiliary subfield.
void
pylith::bc::TimeDependentAuxiliaryFactory::_setVectorFieldComponentNames(pylith::topology::FieldBase::Description* description) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("_setVectorFieldComponentNames(description)");

    assert(description);

    const char** componentNames = NULL;
    if (XYZ == _auxComponents) {
        componentNames = _TimeDependentAuxiliaryFactory::componentsXYZ;
    } else if (TANGENTIAL_NORMAL == _auxComponents) {
        if (2 == _spaceDim) {
            componentNames = _TimeDependentAuxiliaryFactory::componentsTN;
        } else if (3 == _spaceDim) {
            componentNames = _TimeDependentAuxiliaryFactory::componentsTTN;
        } // if/else
    } // if/else
    if (!componentNames) {
        PYLITH_JOURNAL_ERROR("Unknown case for auxiliary component reference ("<<_auxComponents<<") and spatial dimension ("<<_spaceDim<<").");
        throw std::logic_error("Unknown case for auxiliary component reference and spatial dimension.");
    } // if

    assert(size_t(_spaceDim) == description->numComponents);
    for (int i = 0; i < _spaceDim; ++i) {
        description->componentNames[i] = description->label + std::string(componentNames[i]);
    } // for

    PYLITH_METHOD_END;
} // _setVectorFieldComponentNames


// ---------------------------------------------------------------------------------------------------------------------
// Add initial amplitude field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::addInitialAmplitude(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addInitialAmplitude(void)");

    const char* fieldName = "initial_amplitude";

    assert(_defaultDescription);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);

    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.validator = NULL;
    switch (subfieldDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        const size_t numComponents = 1;
        assert(numComponents == subfieldDescription.numComponents);
        assert(numComponents == subfieldDescription.componentNames.size());
        subfieldDescription.componentNames[0] = fieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        _setVectorFieldComponentNames(&subfieldDescription);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_ERROR("Unknown vector field case.");
        throw std::logic_error("Unknown vector field case in TimeDependentAuxiliaryFactory::initialAmplitude().");
    } // switch

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addInitialAmplitude


// ---------------------------------------------------------------------------------------------------------------------
// Add rate amplitude field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::addRateAmplitude(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRateAmplitude(void)");

    const char* fieldName = "rate_amplitude";

    assert(_defaultDescription);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.validator = NULL;
    switch (subfieldDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        const size_t numComponents = 1;
        assert(numComponents == subfieldDescription.numComponents);
        assert(numComponents == subfieldDescription.componentNames.size());
        subfieldDescription.componentNames[0] = fieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        _setVectorFieldComponentNames(&subfieldDescription);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_ERROR("Unknown vector field case.");
        throw std::logic_error("Unknown vector field case in TimeDependentAuxiliaryFactory::addRateAmplitude().");
    } // switch

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addRateAmplitude


// ---------------------------------------------------------------------------------------------------------------------
// Add rate start time field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::addRateStartTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addRateStartTime(void)");

    const char* fieldName = "rate_start_time";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = fieldName;
    subfieldDescription.scale = _normalizer->timeScale();
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addRateStartTime


// ---------------------------------------------------------------------------------------------------------------------
// Add time history amplitude field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::addTimeHistoryAmplitude(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryAmplitude(void)");

    const char* fieldName = "time_history_amplitude";

    assert(_defaultDescription);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);

    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    switch (subfieldDescription.vectorFieldType) {
    case pylith::topology::FieldBase::SCALAR: {
        const size_t numComponents = 1;
        assert(numComponents == subfieldDescription.numComponents);
        assert(numComponents == subfieldDescription.componentNames.size());
        subfieldDescription.componentNames[0] = fieldName;
        break;
    } // SCALAR
    case pylith::topology::FieldBase::VECTOR: {
        _setVectorFieldComponentNames(&subfieldDescription);
        break;
    } // VECTOR
    default:
        PYLITH_JOURNAL_ERROR("Unknown vector field case.");
        throw std::logic_error("Unknown vector field case in TimeDependentAuxiliaryFactory::addTimeHistoryAmplitude().");
    } // switch

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addTimeHistoryAmplitude


// ---------------------------------------------------------------------------------------------------------------------
// Add time history start time field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::addTimeHistoryStartTime(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryStartTime(void)");

    const char* fieldName = "time_history_start_time";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = fieldName;
    subfieldDescription.scale = _normalizer->timeScale();
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, pylith::topology::FieldQuery::dbQueryGeneric);

    PYLITH_METHOD_END;
} // addTimeHistoryStartTime


// ---------------------------------------------------------------------------------------------------------------------
// Add time history value field to auxiliary fields.
void
pylith::bc::TimeDependentAuxiliaryFactory::addTimeHistoryValue(void) {
    PYLITH_METHOD_BEGIN;
    PYLITH_JOURNAL_DEBUG("addTimeHistoryValue(void)");

    const char* fieldName = "time_history_value";

    assert(_defaultDescription);
    assert(_normalizer);
    pylith::topology::FieldBase::Description subfieldDescription(*_defaultDescription);
    subfieldDescription.label = fieldName;
    subfieldDescription.alias = fieldName;
    subfieldDescription.vectorFieldType = pylith::topology::FieldBase::SCALAR;
    subfieldDescription.numComponents = 1;
    subfieldDescription.componentNames.resize(1);
    subfieldDescription.componentNames[0] = fieldName;
    subfieldDescription.scale = _normalizer->timeScale();
    subfieldDescription.validator = NULL;

    _field->subfieldAdd(subfieldDescription, getSubfieldDiscretization(fieldName));
    _setSubfieldQueryFn(fieldName, NULL); // populated by integrator or constraint at begining of time step.

    PYLITH_METHOD_END;
} // addTimeHistoryValue


// ---------------------------------------------------------------------------------------------------------------------
void
pylith::bc::TimeDependentAuxiliaryFactory::updateAuxiliaryField(pylith::topology::Field* auxiliaryField,
                                                                const PylithReal t,
                                                                const PylithReal timeScale,
                                                                spatialdata::spatialdb::TimeHistory* const dbTimeHistory) {
    assert(auxiliaryField);
    assert(dbTimeHistory);

    PetscErrorCode err = 0;

    PetscSection auxiliaryFieldSection = auxiliaryField->localSection();assert(auxiliaryFieldSection);
    PetscInt pStart = 0, pEnd = 0;
    err = PetscSectionGetChart(auxiliaryFieldSection, &pStart, &pEnd);PYLITH_CHECK_ERROR(err);
    pylith::topology::VecVisitorMesh auxiliaryFieldVisitor(*auxiliaryField);
    PetscScalar* auxiliaryFieldArray = auxiliaryFieldVisitor.localArray();assert(auxiliaryFieldArray);

    // Compute offset of time history subfields in auxiliary field.
    const PetscInt i_startTime = auxiliaryField->subfieldInfo("time_history_start_time").index;
    const PetscInt i_value = auxiliaryField->subfieldInfo("time_history_value").index;

    // Loop over all points in section.
    for (PetscInt p = pStart; p < pEnd; ++p) {
        // Skip points without values in section.
        if (!auxiliaryFieldVisitor.sectionDof(p)) {continue;}

        // Get starting time and compute relative time for point.
        const PetscInt offStartTime = auxiliaryFieldVisitor.sectionSubfieldOffset(i_startTime, p);
        const PylithScalar tStart = auxiliaryFieldArray[offStartTime];
        const PylithScalar tRel = t - tStart;

        // Query time history for value (normalized amplitude).
        PylithScalar value = 0.0;
        if (tRel >= 0.0) {
            PylithScalar tDim = tRel * timeScale;
            const int err = dbTimeHistory->query(&value, tDim);
            if (err) {
                std::ostringstream msg;
                msg << "Error querying for time '" << tDim << "' in time history database '" << dbTimeHistory->label() << "'.";
                throw std::runtime_error(msg.str());
            } // if
        } // if

        // Update value (normalized amplitude) in auxiliary field.
        const PetscInt offValue = auxiliaryFieldVisitor.sectionSubfieldOffset(i_value, p);
        auxiliaryFieldArray[offValue] = value;
    } // for

} // updateAuilixaryField


// End of file
