#include <Python.h>
#include <iostream>
#include "bipyramid4.c"
using namespace std;

//---Global object list definition
vertex* ver;
triangle* tri;

//---Initialize all vertices
static PyObject* 
init(PyObject* self, PyObject* args) {
	int arg1,arg2;
	if (!PyArg_ParseTuple(args, "ii", &arg1, &arg2))
		return NULL;
	ver = new vertex[arg1];
	tri = new triangle[arg2];
	nver = arg1;
	ntr = arg2;
	return Py_BuildValue("");
}

//---Add a single vertex
static PyObject* 
addvertex(PyObject* self, PyObject* args) {
	PyObject* arg2;
	PyObject* arg3;
	PyObject* arg4;
	double* points;
	int vernum,arg1,i;
	int noneighbors,notris;
	if (!PyArg_ParseTuple(args, "iOOO", &arg1, &arg2, &arg3, &arg4))
		return NULL;
	arg2 = PySequence_Fast(arg2, "expected a sequence");
	arg3 = PySequence_Fast(arg3, "expected a sequence");
	arg4 = PySequence_Fast(arg4, "expected a sequence");
 	points = new double[3];
    for(i=0; i < 3; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(arg2, i);
        fitem = PyNumber_Float(item);
        points[i] = PyFloat_AS_DOUBLE(fitem);
        Py_DECREF(fitem);
    }    
    Py_DECREF(arg2);
	vernum = arg1;
    //---Unclear why vcoord is 4-dimensional and starts at 1?
    (ver+vernum)->vcoord[1]=points[0];
	(ver+vernum)->vcoord[2]=points[1];
	(ver+vernum)->vcoord[3]=points[2];
	noneighbors = PySequence_Fast_GET_SIZE(arg3);
	(ver+vernum)->nonei=noneighbors-1;
	for(i=0; i < noneighbors; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(arg3, i);
        fitem = PyNumber_Int(item);
        (ver+vernum)->vneipt[i] = PyLong_AsLong(fitem);
        Py_DECREF(fitem);
    }
    Py_DECREF(arg3);
   	notris = PySequence_Fast_GET_SIZE(arg4);
	for(i=0; i < notris; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(arg4, i);
        fitem = PyNumber_Int(item);
        (ver+vernum)->vneitr[i+1] = PyLong_AsLong(fitem);
        Py_DECREF(fitem);
    }
    Py_DECREF(arg4);
    return Py_BuildValue("");	  
}

//---Add a single simplex
static PyObject* 
addtriangle(PyObject* self, PyObject* args) {
	PyObject* arg2;
	int* triangle_indices;
	int trinum,arg1,i;
	if (!PyArg_ParseTuple(args, "iO", &arg1, &arg2))
		return NULL;
	arg2 = PySequence_Fast(arg2, "expected a sequence");
 	triangle_indices = new int[3];
    for(i=0; i < 3; i++) {
        PyObject *fitem;
        PyObject *item = PySequence_Fast_GET_ITEM(arg2, i);
        fitem = PyNumber_Int(item);
        triangle_indices[i] = PyLong_AsLong(fitem);
        Py_DECREF(fitem);
    }    
    Py_DECREF(arg2);
	trinum = arg1;
	//---Unclear why vert is 4-dimensional and starts at 1?
    (tri+trinum)->vert[1]=triangle_indices[0];
	(tri+trinum)->vert[2]=triangle_indices[1];
	(tri+trinum)->vert[3]=triangle_indices[2];
    return Py_BuildValue("");	  
}

//---Generic diagnostic
static PyObject*
check(PyObject *self, PyObject *args)
{
	int i;
	if (!PyArg_ParseTuple(args,""))
		return NULL;
	cout << "Starting area calculation!" << endl;
	for(i=0;i<ntr;i++) {
		areacal(i,tri,ver);
	}
	cout << "Done!" << endl;
	cout << "Starting curvature calculation!" << endl;
	for(i=0;i<nver;i++) {
		CURVCALC(i,ver,tri);
	}
	cout << "Done!" << endl;
	return Py_BuildValue("");	
}

//---Look at the mesh (for debugging purposes).
static PyObject*
getprop(PyObject *self, PyObject *args, PyObject *keywds)
{
	const char* part;
	int index;
	const char* subpart;
	static char *kwlist[] = {"part","index","subpart",NULL};
	if (!PyArg_ParseTupleAndKeywords(args, keywds,"|sis", kwlist, &part, &index, &subpart)) 
		return NULL;
	cout << "Getting mesh property ..." << endl;
	if (strcmp(part,"vertex")==0) {
		cout << "\tvertex = " << index << endl;
		if (strcmp(subpart,"coords")==0) {
			cout << "\t" << (ver+index)->vcoord[0] << "\t" << (ver+index)->vcoord[1] << "\t" << (ver+index)->vcoord[2] << "\t" << endl;
		} 
		else if (strcmp(subpart,"mcur")==0) {
			cout << "\t" << (ver+index)->mcur << endl;				
		}	
		else if (strcmp(subpart,"pdir")==0) {
			cout << "\t" << (ver+index)->cur1 << " " << (ver+index)->cur2 << endl;						
		}	
	}
	else if (strcmp(part,"triangle")==0) {
		cout << "\ttriangle = " << index << endl;
		if (strcmp(subpart,"vert")==0) {
			cout << "\t" << (tri+index)->vert[1] << "\t" << (tri+index)->vert[2] << "\t" << (tri+index)->vert[3] << "\t" << endl;
		}			
	}
	return Py_BuildValue("");	
}

static PyMethodDef memint_methods_py[] = {
	{ "check", check, METH_VARARGS, "Code development tool."},
	{ "init", init, METH_VARARGS, "Initialize the mesh by providing its size."},
	{ "addvertex", addvertex, METH_VARARGS, "Add a vertex to the mesh."},
	{ "addtriangle", addtriangle, METH_VARARGS, "Add a triangle to the mesh."},
	{ "getprop", (PyCFunction)getprop, METH_VARARGS | METH_KEYWORDS, "Access the mesh (to debug)."},
	{ NULL, NULL, 0, NULL}
};	

PyMODINIT_FUNC initmemint(void)
{
    (void) Py_InitModule("memint", memint_methods_py);
}

