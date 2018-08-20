from setuptools import setup, find_packages

setup(
    name='pyabc',
    # Remenber change version in ababe.cmdline.runabalib.exec_from_cmdline()
    version='0.1.0',
    description='A grid structure machine learning tool used\
                 in material science.',
    author='Jason Yu',
    author_email='morty.yu@yahoo.com',
    license='MIT',
    packages=find_packages(),
    include_package_data=True,
    package_data={},
    install_requires=["nose2", "numpy", "spglib==1.9.9.18",
                      "click>=6", "colorama", "ase==3.16"],
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "Intended Audience :: Science/Research",
        "License :: MIT License",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics"],
    scripts=[],
    entry_points={
        'console_scripts': [
            'rexpand = pyabc.cmd.rexpand:cli',
        ],
    }
)
