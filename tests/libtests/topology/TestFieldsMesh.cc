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

#include "TestFieldsMesh.hh" // Implementation of class methods

#include "pylith/topology/Mesh.hh" // USES Mesh
#include "pylith/topology/Field.hh" // USES Field
#include "pylith/topology/Fields.hh" // USES Fields
#include "pylith/topology/Stratum.hh" // USES Stratum
#include "pylith/topology/VisitorMesh.hh" // USES VecVisitorMesh

#include "pylith/meshio/MeshIOAscii.hh" // USES MeshIOAscii

// ----------------------------------------------------------------------
CPPUNIT_TEST_SUITE_REGISTRATION(pylith::topology::TestFieldsMesh);

// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsMesh::setUp(void) {
    PYLITH_METHOD_BEGIN;

    _mesh = new Mesh;
    meshio::MeshIOAscii importer;
    importer.filename("data/tri3.mesh");
    importer.read(_mesh);

    PYLITH_METHOD_END;
} // setUp


// ----------------------------------------------------------------------
void
pylith::topology::TestFieldsMesh::tearDown(void) {
    PYLITH_METHOD_BEGIN;

    delete _mesh;_mesh = NULL;

    PYLITH_METHOD_END;
} // tearDown


// ----------------------------------------------------------------------
// Test constructor.
void
pylith::topology::TestFieldsMesh::testConstructor(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    Fields fields(*_mesh);

    PYLITH_METHOD_END;
} // testConstructor


// ----------------------------------------------------------------------
// Test add().
void
pylith::topology::TestFieldsMesh::testAdd(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    Fields fields(*_mesh);

    const char* label = "field";
    fields.add(label, "displacement");
    const size_t size = 1;
    CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

    CPPUNIT_ASSERT_THROW(fields.add(label, "displacement"), std::runtime_error);

    PYLITH_METHOD_END;
} // testAdd


// ----------------------------------------------------------------------
// Test del().
void
pylith::topology::TestFieldsMesh::testDelete(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    Fields fields(*_mesh);

    const char* keyA = "field A";
    const char* labelA = "displacement";
    fields.add(keyA, labelA);

    const char* keyB = "field B";
    const char* labelB = "velocity";
    fields.add(keyB, labelB);

    size_t size = 2;
    CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
    fields.del(keyA);
    size = 1;
    CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());
    const Field& field = fields.get(keyB);
    CPPUNIT_ASSERT_EQUAL(std::string(labelB), std::string(field.getLabel()));

    CPPUNIT_ASSERT_THROW(fields.del(keyA), std::runtime_error);

    PYLITH_METHOD_END;
} // testDelete


// ----------------------------------------------------------------------
// Test get().
void
pylith::topology::TestFieldsMesh::testGet(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    Fields fields(*_mesh);

    const char* key = "field";
    const char* label = "velocity";
    fields.add(key, label);
    const Field& field = fields.get(key);
    CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(field.getLabel()));

    CPPUNIT_ASSERT_THROW(fields.get("abc"), std::runtime_error);

    const Fields& fieldsB = fields;
    CPPUNIT_ASSERT_THROW(fieldsB.get("abc"), std::runtime_error);

    PYLITH_METHOD_END;
} // testGet


// ----------------------------------------------------------------------
// Test get() const.
void
pylith::topology::TestFieldsMesh::testGetConst(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    Fields fields(*_mesh);

    const char* key = "field";
    const char* label = "velocity";
    fields.add(key, label);

    const Fields* fieldsPtr = &fields;
    CPPUNIT_ASSERT(fieldsPtr);
    const Field& field = fieldsPtr->get(key);
    CPPUNIT_ASSERT_EQUAL(std::string(label), std::string(field.getLabel()));

    PYLITH_METHOD_END;
} // testGetConst


// ----------------------------------------------------------------------
// Test hasField().
void
pylith::topology::TestFieldsMesh::testHasField(void) {
    PYLITH_METHOD_BEGIN;

    CPPUNIT_ASSERT(_mesh);
    Fields fields(*_mesh);

    fields.add("field A", "velocity");

    CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
    CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field B"));
    CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

    fields.add("field B", "displacement");

    CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field A"));
    CPPUNIT_ASSERT_EQUAL(true, fields.hasField("field B"));
    CPPUNIT_ASSERT_EQUAL(false, fields.hasField("field C"));

    PYLITH_METHOD_END;
} // testHasField


// ----------------------------------------------------------------------
// Test copyLayout().
void
pylith::topology::TestFieldsMesh::testCopyLayout(void) {
    PYLITH_METHOD_BEGIN;

    const int fiberDim = 3;

    CPPUNIT_ASSERT(_mesh);
    Fields fields(*_mesh);

    const char* labelA = "field A";
    fields.add(labelA, "displacement");
    Field& fieldA = fields.get(labelA);
    const int dim = _mesh->dimension();
    const int numComponents = 3;
    const char* components[3] = {"x", "y", "z"};
    const int basisOrder = 1;
    const int quadOrder = 1;
    const double scale = 1.2;
    const pylith::topology::FieldBase::CellBasis cellBasis = _mesh->isSimplex() ?
                                                             pylith::topology::FieldBase::SIMPLEX_BASIS :
                                                             pylith::topology::FieldBase::TENSOR_BASIS;
    fieldA.subfieldAdd("displacement", "displacement", Field::SCALAR, components, numComponents, scale, basisOrder, dim,
                       quadOrder, cellBasis, true, Field::POLYNOMIAL_SPACE);
    fieldA.subfieldsSetup();
    fieldA.createDiscretization();
    fieldA.allocate();

    const char* labelB = "field B";
    fields.add(labelB, "velocity");
    fields.copyLayout(labelA);

    const size_t size = 2;
    CPPUNIT_ASSERT_EQUAL(size, fields._fields.size());

    PetscDM dmMesh = _mesh->dmMesh();CPPUNIT_ASSERT(dmMesh);
    Stratum depthStratum(dmMesh, Stratum::DEPTH, 0);
    const PetscInt vStart = depthStratum.begin();
    const PetscInt vEnd = depthStratum.end();

    const Field& field = fields.get(labelB);
    VecVisitorMesh fieldVisitor(field);
    for (PetscInt v = vStart; v < vEnd; ++v) {
        CPPUNIT_ASSERT_EQUAL(fiberDim, fieldVisitor.sectionDof(v));
    } // for

    CPPUNIT_ASSERT_THROW(fields.copyLayout("zyx"), std::runtime_error);

    PYLITH_METHOD_END;
} // testCopyLayout


// End of file
