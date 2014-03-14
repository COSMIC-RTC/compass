#include "Python.h"
#include "structmember.h"
#include <yoga_host_obj.h>
#include <yoga_obj.h>
#include <sstream>
#include <iomanip>

yoga_context *context_handle;

typedef struct {
  PyObject_HEAD
  PyObject *yoga_object;
  int type;
  int device;
} YogaObj;

static void
YogaObj_dealloc(YogaObj* self) {
  delete (yObjD*) self->yoga_object;
  //Py_XDECREF(self->yoga_object);
  self->ob_type->tp_free((PyObject*) self);
}

static PyObject *
YogaObj_new(PyTypeObject *type, PyObject *args, PyObject *kwds) {
  YogaObj *self;

  self = (YogaObj *) type->tp_alloc(type, 0);
  if (self != NULL) {
    double data_input[5] = { 1, 2, 3, 4, 5 };
    long dims_data[2] = { 1, 5 };
    self->yoga_object = (PyObject *) new yObjD(context_handle,
        (long *) dims_data, (double*) data_input, 0);
    self->type = 0;
    self->device = 0;
  }

  return (PyObject *) self;
}

static int
YogaObj_init(YogaObj *self, PyObject *args, PyObject *kwds) {
  PyObject *yoga_object = NULL;

  //    static char *kwlist[] = {"first", "last", "number", NULL};
  //
  //     if (! PyArg_ParseTupleAndKeywords(args, kwds, "|OOi", kwlist,
  //                                       &first, &last,
  //                                       &self->number))
  //         return -1;
  //
  //     if (first) {
  //         tmp = self->first;
  //         Py_INCREF(first);
  //         self->first = first;
  //         Py_XDECREF(tmp);
  //     }
  //
  //     if (last) {
  //         tmp = self->last;
  //         Py_INCREF(last);
  //         self->last = last;
  //         Py_XDECREF(tmp);
  //     }

  return 0;
}

static PyMemberDef YogaObj_members[] = { { "yoga_object", T_OBJECT_EX,
    offsetof(YogaObj, yoga_object), 0,
    "yoga_object name"},
  { "type", T_INT, offsetof(YogaObj, type), 0,
    "type name"},
  { "device", T_INT, offsetof(YogaObj, device), 0,
    "device number"},
  { NULL} /* Sentinel */
};

static PyObject *
YogaObj_name(YogaObj* self) {
  //static PyObject *format = NULL;
  //PyObject *args, *result = NULL;
  PyObject *result = NULL;

  // if (format == NULL) {
  //     format = PyString_FromString("%s %s");
  //     if (format == NULL)
  //         return NULL;
  // }

  // if (self->first == NULL) {
  //     PyErr_SetString(PyExc_AttributeError, "first");
  //     return NULL;
  // }

  // if (self->last == NULL) {
  //     PyErr_SetString(PyExc_AttributeError, "last");
  //     return NULL;
  // }

  // args = Py_BuildValue("OO", self->first, self->last);
  // if (args == NULL)
  //     return NULL;

  // result = PyString_Format(format, args);
  // Py_DECREF(args);

  return result;
}

static PyMethodDef YogaObj_methods[] = { { "name", (PyCFunction) YogaObj_name,
    METH_NOARGS, "Return nothing" }, { NULL } /* Sentinel */
};

static PyTypeObject YogaObjType = {
  PyObject_HEAD_INIT(NULL)
  0, /*ob_size*/
  "YogaObj.YogaObj", /*tp_name*/
  sizeof(YogaObj), /*tp_basicsize*/
  0, /*tp_itemsize*/
  (destructor)YogaObj_dealloc, /*tp_dealloc*/
  0, /*tp_print*/
  0, /*tp_getattr*/
  0, /*tp_setattr*/
  0, /*tp_compare*/
  0, /*tp_repr*/
  0, /*tp_as_number*/
  0, /*tp_as_sequence*/
  0, /*tp_as_mapping*/
  0, /*tp_hash */
  0, /*tp_call*/
  0, /*tp_str*/
  0, /*tp_getattro*/
  0, /*tp_setattro*/
  0, /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  "YogaObj objects", /* tp_doc */
  0, /* tp_traverse */
  0, /* tp_clear */
  0, /* tp_richcompare */
  0, /* tp_weaklistoffset */
  0, /* tp_iter */
  0, /* tp_iternext */
  YogaObj_methods, /* tp_methods */
  YogaObj_members, /* tp_members */
  0, /* tp_getset */
  0, /* tp_base */
  0, /* tp_dict */
  0, /* tp_descr_get */
  0, /* tp_descr_set */
  0, /* tp_dictoffset */
  (initproc)YogaObj_init, /* tp_init */
  0, /* tp_alloc */
  YogaObj_new, /* tp_new */
};

static PyMethodDef module_methods[] = { { NULL } /* Sentinel */
};

#ifndef PyMODINIT_FUNC	/* declarations for DLL import/export */
#define PyMODINIT_FUNC void
#endif
PyMODINIT_FUNC
initYoga(void) {
  PyObject* m;

  if (PyType_Ready(&YogaObjType) < 0)
    return;

  m = Py_InitModule3("Yoga", module_methods,
      "module that creates Yoga extension type.");

  if (m == NULL)
    return;

  context_handle = new yoga_context();
  Py_INCREF (&YogaObjType);
  PyModule_AddObject(m, "YogaObj", (PyObject *) &YogaObjType);
}

