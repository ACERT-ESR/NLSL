#include <Python.h>
#include <string.h>

static PyObject *pylog_module = NULL;
static PyObject *pylog_logger = NULL;

static int ensure_pylog_module(void) {
    if (pylog_module != NULL) {
        return 0;
    }
    pylog_module = PyImport_ImportModule("nlsl.logging");
    if (pylog_module == NULL) {
        PyErr_Print();
        return -1;
    }
    return 0;
}

static PyObject *get_logger(void) {
    if (pylog_logger != NULL) {
        Py_INCREF(pylog_logger);
        return pylog_logger;
    }
    PyObject *logging_mod = PyImport_ImportModule("logging");
    if (logging_mod == NULL) {
        PyErr_Print();
        return NULL;
    }
    PyObject *logger = PyObject_CallMethod(logging_mod, "getLogger", "s", "nlsl");
    Py_DECREF(logging_mod);
    if (logger == NULL) {
        PyErr_Print();
        return NULL;
    }
    pylog_logger = logger;
    Py_INCREF(pylog_logger);
    return logger;
}

static int call_helper(const char *name, PyObject *arg) {
    if (ensure_pylog_module() != 0) {
        return -1;
    }
    PyObject *callable = PyObject_GetAttrString(pylog_module, name);
    if (callable == NULL) {
        PyErr_Print();
        return -1;
    }
    PyObject *result;
    if (arg != NULL) {
        result = PyObject_CallFunctionObjArgs(callable, arg, NULL);
    } else {
        result = PyObject_CallFunctionObjArgs(callable, NULL);
    }
    Py_DECREF(callable);
    if (result == NULL) {
        PyErr_Print();
        return -1;
    }
    Py_DECREF(result);
    return 0;
}

int nlsl_pylog_open(const char *filename) {
    PyGILState_STATE gil = PyGILState_Ensure();
    int status = -1;
    PyObject *py_filename = PyUnicode_DecodeUTF8(filename, (Py_ssize_t)strlen(filename), "strict");
    if (py_filename == NULL) {
        PyErr_Print();
        PyGILState_Release(gil);
        return -1;
    }
    status = call_helper("open_log", py_filename);
    Py_DECREF(py_filename);
    if (status == 0) {
        PyObject *logger = get_logger();
        if (logger != NULL) {
            PyObject *setlevel = PyObject_CallMethod(logger, "setLevel", "s", "DEBUG");
            if (setlevel == NULL) {
                PyErr_Print();
                status = -1;
            } else {
                Py_DECREF(setlevel);
            }
            Py_DECREF(logger);
        } else {
            status = -1;
        }
    }
    PyGILState_Release(gil);
    return status;
}

int nlsl_pylog_debug(const char *message) {
    PyGILState_STATE gil = PyGILState_Ensure();
    int status = -1;
    PyObject *logger = get_logger();
    if (logger == NULL) {
        PyGILState_Release(gil);
        return -1;
    }
    PyObject *py_message = PyUnicode_DecodeUTF8(message, (Py_ssize_t)strlen(message), "strict");
    if (py_message == NULL) {
        PyErr_Print();
        Py_DECREF(logger);
        PyGILState_Release(gil);
        return -1;
    }
    PyObject *result = PyObject_CallMethod(logger, "debug", "O", py_message);
    Py_DECREF(py_message);
    Py_DECREF(logger);
    if (result == NULL) {
        PyErr_Print();
        PyGILState_Release(gil);
        return -1;
    }
    Py_DECREF(result);
    status = 0;
    PyGILState_Release(gil);
    return status;
}

int nlsl_pylog_close(void) {
    PyGILState_STATE gil = PyGILState_Ensure();
    int status = call_helper("close_log", NULL);
    PyGILState_Release(gil);
    return status;
}
