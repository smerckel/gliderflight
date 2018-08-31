#!/usr/bin/env python

import codecs
import os
from setuptools import find_packages, setup


# semantic versioning
VERSION = '0.1.0'


here = os.path.abspath(os.path.dirname(__file__))

with codecs.open(os.path.join(here, 'README.rst'), encoding='utf8') as fh:
    long_description = fh.read()

with open('requirements.txt') as fobj:
    install_requires = [line.strip() for line in fobj]

#scripts = ['']

setup(
    name='gliderflight',
    version=VERSION,
    description='Flight model for Slocum ocean gliders',
    long_description=long_description,
    url='https://github.com/',
    author='Lucas Merckelbach',
    author_email='lucas.merckelbach@hzg.de',
    license='MIT',
    py_modules=['gliderflight'],
    #packages=find_packages(where='.', exclude=['tests', 'docs'])
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Environment :: Console',
        'Natural Language :: English',
        'Operating System :: MacOS',
        'Operating System :: POSIX',
        'Programming Language :: Python :: 3.6',
    ],
    keywords=['ocean gliders', 'glider flight', 'oceanography'],
    install_requires=install_requires,
    include_package_data=True,
    #scripts=scripts,
)
# how to provide scripts etct
#      py_modules = ['fastachar'],
#      entry_points = {'console_scripts':[],#['fastachar_gui = fastachar_gui:main
#                      'gui_scripts':['fastachar = fastachar:main']
#                      },
