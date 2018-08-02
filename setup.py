#!/usr/bin/python3

import shutil
import sys
import sysconfig
import os

cmakeCompiler = None
buildDirectory = "build_python"
ninja_available = False

if sys.version_info.major < 3:
	print("ERROR: NetworKit requires Python 3.")
	sys.exit(1)

if "NETWORKIT_OVERRIDE_CXX" in os.environ:
	cmakeCompiler = os.environ["NETWORKIT_OVERRIDE_CXX"]

if shutil.which("cmake") is None:
	print("ERROR: NetworKit compilation requires cmake.")
	sys.exit(1)

ninja_available = shutil.which("ninja") is not None
if not (ninja_available or shutil.which("make")):
	print("ERROR: NetworKit compilation requires Make or Ninja.")
	sys.exit(1)
try:
	from setuptools import setup # to ensure setuptools is installed
except ImportError:
	print("ERROR: Setuptools is required to install networkit python module.\nInstall via pip3 install setuptools.")
	sys.exit(1)

import os
import subprocess #calling cmake, make and cython

################################################
# get the optional arguments for the compilation
################################################
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("-j", "--jobs", dest="jobs", help="specify number of jobs")
(options,args) = parser.parse_known_args()
if options.jobs is not None:
	jobs = options.jobs
if "NETWORKIT_PARALLEL_JOBS" in os.environ:
	jobs = int(os.environ["NETWORKIT_PARALLEL_JOBS"])
else:
	import multiprocessing
	jobs = multiprocessing.cpu_count()
# make sure sys.argv is correct for setuptools
sys.argv = [__file__] + args

################################################
# compiler identification
################################################

candidates = ["g++", "g++-8", "g++-7", "g++-6.1", "g++-6", "g++-5.3", "g++-5.2", "g++-5.1", "g++-5", "g++-4.9", "g++-4.8", "clang++", "clang++-3.8", "clang++-3.7"]

def determineCompiler(candidates, std, flags):
	#TODO: proper c++11 test?
	sample = open("sample.cpp", "w")
	sample.write("""
	#include <omp.h>
	#include <iostream>
	void helloWorld() {
		std::cout << "Hello world" << std::endl;
	}
	int main (int argc, char *argv[]) {
		helloWorld();
		int nthreads, tid;
		#pragma omp parallel private(nthreads, tid)
		{
			tid = omp_get_thread_num();
			std::cout << \"Hello World from thread = \" << tid << std::endl;
			if (tid == 0) {
				nthreads = omp_get_num_threads();
				std::cout << \"Number of threads = \" << nthreads << std::endl;
			}
		}
	}""")
	sample.close()
	for compiler in candidates:
		cmd = [compiler,"-o","test_build","-std={}".format(std)]
		cmd.extend(flags)
		cmd.append("sample.cpp")
		try:
			if subprocess.call(cmd,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0:
				os.remove("sample.cpp")
				os.remove("test_build")
				return compiler
		except:
			pass
	return ""

# only check for a compiler if none is specified
if cmakeCompiler is None:
	cmakeCompiler = determineCompiler(candidates, "c++11", ["-fopenmp"])
	if cmakeCompiler is None:
		print("ERROR: No suitable compiler found. Install any of these: ",candidates)
		exit(1)

################################################
# functions for cythonizing and building networkit
################################################

def cythonizeFile(filepath):
	cpp_file = filepath.replace("pyx","cpp")

	cython_available = shutil.which("cython") is not None
	if not cython_available:
		if not os.path.isfile(cpp_file):
			print("ERROR: Neither cython nor _NetworKit.cpp is provided. Build cancelled", flush=True)
			exit(1)

		else:
			print("Cython not available, but _NetworKit.cpp provided. Continue build without cythonizing", flush=True)

	elif os.path.isfile(cpp_file) and os.path.getmtime(filepath) < os.path.getmtime(cpp_file):
		print("Cython available; skip as _NetworKit.cpp was create after last modification of _NetworKit.pyx", flush=True)

	else:
		print("Cythonizing _NetworKit.pyx...", flush=True)
		if not os.path.isfile(filepath):
			print("_NetworKit.pyx is not available. Build cancelled.")
			exit(1)
		comp_cmd = ["cython","-3","--cplus","-t",filepath]
		if not subprocess.call(comp_cmd) == 0:
			print("cython returned an error, exiting setup.py")
			exit(1)
		print("_NetworKit.pyx cythonized", flush=True)

def buildNetworKit(install_prefix, withTests = False):
	# Cythonize file
	cythonizeFile("networkit/_NetworKit.pyx")
	try:
		os.makedirs(buildDirectory)
	except FileExistsError:
		pass
	# Build cmake call
	abs_prefix = os.path.join(os.getcwd(), install_prefix)
	comp_cmd = ["cmake","-DCMAKE_BUILD_TYPE=Release"]
	comp_cmd.append("-DCMAKE_INSTALL_PREFIX="+abs_prefix)
	comp_cmd.append("-DCMAKE_CXX_COMPILER="+cmakeCompiler)
	from sysconfig import get_paths, get_config_var
	comp_cmd.append("-DNETWORKIT_PYTHON="+get_paths()['include']) #provide python.h files
	comp_cmd.append("-DNETWORKIT_PYTHON_SOABI="+get_config_var('SOABI')) #provide lib env specification
	if ninja_available:
		comp_cmd.append("-GNinja")
	comp_cmd.append(os.getcwd()) #call CMakeLists.txt from networkit root
	# Run cmake
	print("initializing NetworKit compilation with: '{0}'".format(" ".join(comp_cmd)), flush=True)
	if not subprocess.call(comp_cmd, cwd=buildDirectory) == 0:
		print("cmake returned an error, exiting setup.py")
		exit(1)
	build_cmd = []
	if ninja_available:
		build_cmd = ["ninja", "install", "-j"+str(jobs)]
	else:
		build_cmd = ["make", "install", "-j"+str(jobs)]
	print("Build with: '{0}'".format(" ".join(build_cmd)), flush=True)
	if not subprocess.call(build_cmd, cwd=buildDirectory) == 0:
		print("Build tool returned an error, exiting setup.py")
		exit(1)

################################################
# custom build commands to integrate with setuptools
################################################

# Unfortunately, there is not simple way to invoke external commands from setuptools.
# Thus, we replace the 'build_ext' stage of setuptools completely.
# This looks a bit messy as we need to make sure we satisfy the (undocumented) internal APIs
# of setuptools (and distutils, which setuptools is based on).

# Our build_ext command compiles NetworKit using cmake and installs it into a directory
# supplied by distutils/setuptools. distutils/setuptools picks up all .so files from that
# directory when creating a Python package so that we do not have to do any
# additional installation.

from setuptools import Command
from setuptools import Extension

class build_ext(Command):
	user_options = [ ]

	def initialize_options(self):
		self.build_lib = None # Output directory for libraries
		self.build_temp = None # Temporary directory

		self.extensions = None
		self.package = None

	def finalize_options(self):
		self.set_undefined_options('build',
				('build_lib', 'build_lib'),
				('build_temp', 'build_temp'))

		self.extensions = self.distribution.ext_modules
		self.package = self.distribution.ext_package

	def get_source_files(self):
		sources = [ ]
		for subdir, dirs, files in os.walk('networkit'):
			for filename in files:
				sources.append(os.path.join(subdir, filename))
		return sources

	# Returns the full Python module name of an extension.
	# Prepends the ext_package setup() option to an extension (see distutils).
	def get_ext_fullname(self, extname):
		if self.package is not None:
			return self.package + '.' + extname
		return extname

	# Returns the file name of the DSO implementing a module (see distutils).
	def get_ext_filename(self, fullname):
		return fullname + '.' + sysconfig.get_config_var('SOABI') + '.so'

	def run(self):
		# A generic build_ext command for cmake would iterate over all self.extensions.
		# However, we know that we only want to build a single extension, i.e. NetworKit.
		buildNetworKit(self.build_lib)

################################################
# initialize python setup
################################################
from setuptools import find_packages # in addition to setup
import version
setup(
	name				= version.name,
	version				= version.version,
	author				= version.author,
	author_email		= version.author_email,
	url					= version.url,
	download_url		= version.download_url,
	description			= version.description,
	long_description	= version.long_description,
	license				= version.license,
	packages			= find_packages(),
	package_data		= {'networkit.profiling': ['html/*','latex/*','description/*']},
	keywords			= version.keywords,
	platforms			= version.platforms,
	classifiers			= version.classifiers,
	cmdclass			= {'build_ext': build_ext},
	ext_modules			= [Extension('_NetworKit', sources=[])],
	test_suite			= 'nose.collector',
	install_requires	= version.install_requires,
	zip_safe			= False) # see https://cython.readthedocs.io/en/latest/src/reference/compilation.html

################################################
# check for tkinter installation
# cant be handled by install_requires because
# tkinter is not part of pip
################################################
try:
	import tkinter #python3 only
	del tkinter #python3 only
except:
	print("WARNING: _tkinter is necessary for networkit.profiling module if not used via jupyter.\nInstall tkinter if needed")
