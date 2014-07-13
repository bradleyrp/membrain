from distutils.core import setup, Extension
import numpy

setup (	
	name = "memint",
	version = "0.1",
	include_dirs = [numpy.get_include()],
	ext_modules = [Extension('memint',
		language='c++',
		libraries=['m'],
		sources=['memint.cpp'], 
		extra_compile_args=['-O3','-fPIC','-lm'],
		extra_link_args=['-fPIC','-shared'])]
	)
