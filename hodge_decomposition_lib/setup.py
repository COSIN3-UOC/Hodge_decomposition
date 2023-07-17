from setuptools import find_packages, setup

setup(
    name='hodge_decomposition_lib',
    packages=find_packages(include=['hodge_dec']),
    version='0.1.0',
    description='My first Python library',
    author='Robert Benassai Dalmau',
    license='MIT',
    install_requires=['networkx'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)

