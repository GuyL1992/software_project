from setuptools import setup, Extension


setup(name='mykmeanssp',
      version='1.0',
      author= 'Guy & Yair',
      description='Implemantation kmeans algorothm By Yair & Guy ',
      ext_modules=[Extension('mykmeanssp', sources=['kmeans.c'])])
