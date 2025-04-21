from setuptools import setup 
  
setup( 
    name='aqua', 
    version='0.1', 
    description='A spectral method Navier-Stokes solver', 
    author='Kento Kaneko', 
    author_email='kento.kaneko@epfl.ch', 
    packages=['aqua'], 
    install_requires=[ 
        'torch', 
        'numpy', 
        'matplotlib'
    ], 
) 