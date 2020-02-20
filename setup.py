from setuptools import setup

from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.rst'),'r', encoding='utf-8') as f:
    long_description = f.read()

setup(name='straph',
      version='0.1.1',
      description='Straph is a python package to modelize, analyze and visualize Stream Graphs.',
      author='Leo Rannou',
      author_email='leo.rannou@gmail.com',
      long_description=long_description,
      long_description_content_type="text/x-rst",
      license='Apache 2.0',
      url='https://github.com/StraphX/Straph',
      packages=['straph',
		        'straph.attributed',
                'straph.components',
                'straph.condensation',
                'straph.EulerTourForest',
                'straph.generators',
                'straph.parser',
                'straph.paths',
                'straph.utils'
		],
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Intended Audience :: Science/Research",
          "Programming Language :: Python :: 3",
          "License :: OSI Approved :: Apache Software License",
          "Operating System :: OS Independent",
      ],
      python_requires='>=3.6')
