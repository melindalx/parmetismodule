from distutils.core import setup, Extension
import numpy as np
import mpi4py
import os

# os.environ['CC'] = 'mpicc'

module1 = Extension(
            'parmetis',
			sources = ['parmetismodule.c'],
			include_dirs = [np.get_include(),
                            mpi4py.get_include(),
                            '/usr/include/mpich'],
			libraries = ['parmetis', 'metis', 'mpich'],
			library_dirs = ['/usr/lib/x86_64-linux-gnu'],
			extra_compile_args = ['-Wno-cpp'])

setup (name = 'parmetis',
		version = '1.0',
		description = 'Python interface of ParMETIS library.',
		author = 'Xue Li',
        author_email = 'xli26@nd.edu',
		ext_modules = [module1])
