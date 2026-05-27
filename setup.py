from setuptools import setup, find_packages

setup(
    name="hybrid_quantum_solver",
    packages=find_packages(where="."),
    package_dir={"": "."},
    zip_safe=False,
)