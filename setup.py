from setuptools import (
    setup,
    find_namespace_packages,
)


setup(
    name="idaes-examples",
    version="2.0.dev0",
    package_dir={"": "pkg"},
    packages=find_namespace_packages(where="pkg"),
)
