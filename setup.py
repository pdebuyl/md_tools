from distutils.core import setup
from Cython.Build import cythonize

setup(name="md_tools",
      version="0.1.0",
      description="md_tools - A few tools for Molecular Dynamics analysis",
      author="Pierre de Buyl",
      author_email="pdebuyl at ulb dot ac dot be",
      maintainer="Pierre de Buyl",
      maintainer_email="pdebuyl at ulb dot ac dot be",
      license="BSD", url="http://pdebuyl.be/",
      packages=["md_tools"],
      ext_modules = cythonize(["md_tools/_rdf.pyx", "md_tools/_cell_based_rdf.pyx"]),
)
