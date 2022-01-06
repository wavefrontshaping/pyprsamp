try:
    from setuptools import setup
    from setuptools import Extension
except ImportError:
    from distutils.core import setup
    from distutils.extension import Extension


#from distutils.core import setup
#from distutils.extension import Extension
from shutil import move
import sys

python_version = sys.version_info[0]
if not python_version in [2,3]:
    raise "Must use Python 2.X or 3.X."


#os.system('SET VS90COMNTOOLS=%VS150COMNTOOLS%')


#module1 = Extension('calib',
#                    #define_macros = [('MAJOR_VERSION', '1'),
#                    #                 ('MINOR_VERSION', '0')],
#                    include_dirs = ['.','C:\Anaconda\include','C:\Program Files (x86)\MPICH2\include','/usr/lib/openmpi/include', '/home/spopoff/anaconda2/include'],
#                    libraries = ['mpi','m'],#,'libpython2.7'],
#                    library_dirs = ['C:\Anaconda','C:\Program Files (x86)\MPICH2\lib', '/usr/lib/openmpi/lib', '/home/spopoff/anaconda2/lib'],
#                    sources = ['prSAMP_Calibration.c','cgamp.c','cpr.c','common.c','bessel_I0.c','bessel_I1.c','cgb.c'])


Cmodule = Extension('prsamp',
                    define_macros = [('PYTHON_VERSION', sys.version_info[0])],
#                    define_macros = [('PYTHON'+str(python_version),'')],
                    #include_dirs = ['.',"C:\Program Files\Anaconda2\include",'C:\Anaconda\include','/home/spopoff/anaconda2/include'],
                    #libraries = ['m'],#,'python27'],
                    #library_dirs = ['C:\Anaconda',"C:\Program Files\Anaconda2\libs",'/home/spopoff/anaconda2/lib'],
                    sources = ['./src/C/prSAMP_Calibration.c','./src/C/cgamp.c','./src/C/cpr.c',
                               './src/C/common.c','./src/C/bessel_I0.c','./src/C/bessel_I1.c','./src/C/cgb.c'])

setup (name = 'pyprsamp',
       version = '1.0.3',
       description = 'This is a demo package',
       author = 'S.M.P',
       #author_email = '',
       #url = 'https://docs.python.org/extending/building',
       long_description = '''
       Prsamp module for Python.
       ''',
       py_modules = ['./src/pyprsamp'],
       ext_modules = [Cmodule])

#move("./build/lib.win32-2.7/pyprsamp.pyd","./pyprsamp.pyd")
#move("./build/lib.linux-x86_64-2.7/prsamp.so","./prsamp.so")
#move("./build/lib.linux-x86_64-3.5/prsamp.cpython-35m-x86_64-linux-gnu.so","./prsamp.so")
#move("./build/lib.linux-x86_64-3.5/pyprsamp.py","./pyprsamp.py")
