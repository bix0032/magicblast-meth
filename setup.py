import ez_setup
import sys
ez_setup.use_setuptools()
from setuptools import setup

# from mpld3
def get_version(path):
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(path) as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
               if isinstance(node, ast.Assign)
               and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")

install_requires = ['toolshed']
if sys.version_info[:2] < (2, 7):
    install_requires.extend(["argparse", "ordereddict"])

setup(name='mbmeth',
      version=get_version("mbmeth.py"),
      description="align BS-Seq reads with Magic-BLAST",
      py_modules=['mbmeth'],
      author="Shaojun Xie",
      author_email="xie186@purdue.edu",
      license="MIT",
      install_requires=install_requires,
      long_description=open('README.md').read(),
      classifiers=[
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'Programming Language :: Python :: 2',
      'Programming Language :: Python :: 3'
      ],
      scripts=['mbmeth.py']
)
