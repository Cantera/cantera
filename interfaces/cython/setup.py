from setuptools import setup, Extension

extension = Extension("cantera._cantera", sources=[])

setup(ext_modules=[extension])
