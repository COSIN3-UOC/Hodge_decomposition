from setuptools import find_packages, setup

setup(
    name='hodge_decomposition_lib',
    packages=find_packages(include=['hodge_dec']),
    version='0.1.1',
    description='This library performs the Hodge Decomposition of a graph in \
    networkx. The graph must be Directed and weighted with an arbitrary flow \
        called edge_visits',
    author='Robert Benassai Dalmau',
    license='MIT',
    install_requires=['networkx'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests',
)

