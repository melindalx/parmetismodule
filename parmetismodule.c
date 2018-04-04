#include "Python.h"
#include "numpy/arrayobject.h"
#include "mpi4py/mpi4py.h"
// #include "mpi.h"
#include "parmetis.h"

double *cdoublearray(PyArrayObject *pyarr)
{
	if (NULL == pyarr)
		return NULL;

	if (NPY_DOUBLE != pyarr->descr->type_num)
	{
		printf("Array's type is not double!\n");
	}

	double *pyarrdata = (double *)pyarr->data;

	return pyarrdata;
}

int64_t *cintarray(PyArrayObject *pyarr)
{
	if (NULL == pyarr)
		return NULL;

	if (NPY_INT64 != pyarr->descr->type_num)
	{
		printf("Array's type is not integer!\n");
	}

	int64_t *pyarrdata = (int64_t *)pyarr->data;

	return pyarrdata;
}

static PyObject *
parmetis_PyParMETIS_V3_PartKway(PyObject *self, PyObject *args, PyObject *keywds)
{
	PyArrayObject *pyvtxdist = NULL;
	PyArrayObject *pyxadj = NULL;
	PyArrayObject *pyadjncy = NULL;
	PyArrayObject *pyvwgt = NULL;
	PyArrayObject *pyadjwgt = NULL;
	PyArrayObject *pytpwgts = NULL;
	PyArrayObject *pyubvec = NULL;
	PyArrayObject *pyoptions = NULL;
	idx_t wgtflag = 0;
	idx_t numflag = 0;
	idx_t ncon = 1;
	idx_t nparts = 0; // necessary

	idx_t *vtxdist = NULL;
	idx_t *xadj = NULL;
	idx_t *adjncy = NULL;
	idx_t *vwgt = NULL;
	idx_t *adjwgt = NULL;
	real_t *tpwgts = NULL;
	real_t *ubvec = NULL;
	idx_t dfoptions[] = {0, 0, 0};
	idx_t *options = dfoptions;

	PyObject *pycomm = NULL;
	MPI_Comm *pcomm = NULL;

	static char *kwlist[] = {
		"vtxdist", "xadj", "adjncy", "nparts", "comm", "vwgt", "adjwgt",
		"wgtflag", "numflag", "ncon", "tpwgts", "ubvec", "options", NULL};
	if (!PyArg_ParseTupleAndKeywords(
		args, keywds, "O!O!O!iO|O!O!iiiO!O!O!", kwlist,
		&PyArray_Type, &pyvtxdist, &PyArray_Type, &pyxadj,
		&PyArray_Type, &pyadjncy, &nparts, &pycomm,
		&PyArray_Type, &pyvwgt, &PyArray_Type, &pyadjwgt,
		&wgtflag, &numflag, &ncon,
		&PyArray_Type, &pytpwgts, &PyArray_Type, &pyubvec,
		&PyArray_Type, &pyoptions))
		return NULL;

	pcomm = PyMPIComm_Get(pycomm);
	if (NULL == pcomm)
		return NULL;

	vtxdist = cintarray(pyvtxdist);
	xadj = cintarray(pyxadj);
	adjncy = cintarray(pyadjncy);
	if (NULL != pyvwgt)
		vwgt = cintarray(pyvwgt);
	if (NULL != pyadjwgt)
		adjwgt = cintarray(pyadjwgt);
	if (NULL != pytpwgts)
	{
		tpwgts = cdoublearray(pytpwgts);
	} else {
		tpwgts = (double *)malloc(sizeof(double) * ncon * nparts);
		for (int itp = 0; itp < ncon * nparts; ++itp)
		{
			tpwgts[itp] = 1.0 / (double)nparts;
		}
	}
	if (NULL != pyubvec)
	{
		ubvec = cdoublearray(pyubvec);
	} else {
		ubvec = (double *)malloc(sizeof(double) * ncon);
		for (int iub = 0; iub < ncon; ++iub)
		{
			ubvec[iub] = 1.05;
		}
	}
	if (NULL != pyoptions)
	{
		options = cintarray(pyoptions);
	}

	int nvertex = pyxadj->dimensions[0] - 1;
    idx_t *part = (int64_t *)malloc(sizeof(int64_t) * nvertex);
    idx_t edgecut;

	int err = ParMETIS_V3_PartKway(
		vtxdist, xadj, adjncy, vwgt, adjwgt, &wgtflag,
		&numflag, &ncon, &nparts, tpwgts, ubvec,
		options, &edgecut, part, pcomm);

	PyArrayObject *pypart =
		(PyArrayObject *)PyArray_FromDims(1, &nvertex, NPY_INT64);
	uint64_t *pypartdata = (uint64_t *)pypart->data;
	for (int ipart = 0; ipart < nvertex; ++ipart)
	{
		pypartdata[ipart] = part[ipart];
	}

	free(part);
	if (NULL == pytpwgts) free(tpwgts);
	if (NULL == pyubvec) free(ubvec);

	return Py_BuildValue("iiO", err, edgecut, pypart);
}

// The module's methods table and initialize function
static PyMethodDef ParmetisMethods[] = {
	{"PyParMETIS_V3_PartKway",
		(PyCFunction) parmetis_PyParMETIS_V3_PartKway,
		METH_VARARGS | METH_KEYWORDS,
		"Partation the graph parallelly."},
	{NULL, NULL, 0, NULL}
};

void initparmetis(void)
{
	PyObject *m;

	m = Py_InitModule("parmetis", ParmetisMethods);
	if (m == NULL)
		return;

	import_array();
	import_mpi4py();
}
