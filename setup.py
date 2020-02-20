from setuptools import setup

# read the contents of your README file
from os import path
this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.rst'),'r', encoding='utf-8') as f:
    long_description = f.read()

setup(name='straph',
      version='0.1',
      description='Straph is a python package to modelize, analyze and visualize Stream Graphs.',
      author='Leo Rannou',
      author_email='leo.rannou@gmail.com',
      long_description=long_description,
      long_description_content_type="text/x-rst",
      #license='GNUGPLv3 License',
      url='https://github.com/EolOurnan/Straph',
      packages=['straph',
                'straph.analysis',
                'straph.components',
                'straph.dag',
                'straph.EulerTourForest',
                'straph.generators',
                'straph.interval_tree',
                'straph.ML_algorithms',
                'straph.parser',
                'straph.paths',
                'straph.streaming_algorithms',
                'straph.utils',
                'straph.wrack'],
      classifiers=[
          "Development Status :: 3 - Alpha",
          "Intended Audience :: Science/Research",
          "Programming Language :: Python :: 3",
          #"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
          "Operating System :: OS Independent",
      ],
      python_requires='>=3.6')