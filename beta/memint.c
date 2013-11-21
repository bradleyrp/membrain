#include <Python.h>
//#include <iostream>
#include <stdio.h>
#include "bipyramid4.c"
//using namespace std;

//vertex* ver

static PyObject*
thru(PyObject *self, PyObject *args)
{
	int nver_in;
	int nver;
	if (!PyArg_ParseTuple(args, "i", &nver_in)) return NULL;
	//ver = (vertex*)malloc(size*sizeof(vertex))
	nver = nver_in;
	printf("hello, world %d",nver);

 	return Py_BuildValue("");	  
}

static PyMethodDef memint_methods_py[] = {
	{ "thru", thru, METH_VARARGS, "Code development tool."},
	{ NULL, NULL, 0, NULL}
};	

PyMODINIT_FUNC initmemint(void)
{
    (void) Py_InitModule("memint", memint_methods_py);
}

