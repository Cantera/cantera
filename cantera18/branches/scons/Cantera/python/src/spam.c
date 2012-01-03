/*
 * A trivial Python extension module introduced to force distutils to
 * generate installation packages requiring a specific Python version.
 */

#include <Python.h>

PyMODINIT_FUNC
init_spam(void)
{
  Py_InitModule3("_spam", NULL,
           "Spam, spam, spam, eggs, and spam.");
}
