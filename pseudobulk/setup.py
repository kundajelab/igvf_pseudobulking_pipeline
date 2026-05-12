from setuptools import setup, Extension
from Cython.Build import cythonize

extensions = [
    Extension(
        name="pseudobulk.tss",
        sources=["src/pseudobulk/tss.pyx"],
    )
]

setup(
    ext_modules=cythonize(
        extensions,
        compiler_directives={
            "language_level": "3",
            "freethreading_compatible": True,  # Crucial flag for Python 3.14t
        },
    )
)
