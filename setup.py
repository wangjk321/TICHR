from setuptools import setup, find_packages

setup(
    name='tichr',
    version='0.1.7',
    author='Jiankang Wang',
    author_email='wangjk321@gmail.com',
    description='TICHR: A computational tool designed to investigate transcriptional regulation',
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=[
    "numpy",
    "pandas",
    "scipy", 
    "matplotlib",
    "seaborn",
    "pyBigWig",
    "scikit-learn",
    "hic-straw",
    "statsmodels",
    "pyranges",
    "tqdm",
    "rpy2",
    # "joblib",
    "adjustText",
    "umap-learn",
    # "multiprocessing",
    # "collections",
    # "concurrent",
    # "functools",
],
    entry_points={
        'console_scripts': [
            'tichr=tichr.__main__:main' 
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',  
        'Operating System :: OS Independent',
    ],
    include_package_data=True,
)

