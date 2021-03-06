
from setuptools import setup

setup(name='hgt_detection',
      version='1.0',
      description="Detects horizontal gene transfers from assemblies of bacterial communities.",
      author='Eric Ulrich',
      url='https://github.com/nahanoo/hgt',
      packages=['hgt'],
      install_requires=['pandas',
                        'pysam',
                        'Bio',
                        'dna_features_viewer'],
      entry_points={
          'console_scripts': [
              'detect_hgts = hgt.main:main'
          ]
      }
      )
