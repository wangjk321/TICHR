from setuptools import setup, find_packages

setup(
    name='tichr',
    version='0.1.0',
    author='Jiankang Wang',
    author_email='xxx',
    description='TICHR: A computational tool designed to investigate transcriptional regulation',
    # long_description=open('README.md', encoding='utf-8').read(),
    # long_description_content_type='text/markdown',
    # url='https://github.com/wangjk321/tichr', 
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
    # "rpy2",
    # "joblib",
    "adjustText",
    "umap-learn",
    # "multiprocessing",
    # "collections",
    # "concurrent",
    # "functools",
    #这几个好像是python自带的
],
    entry_points={
        'console_scripts': [
            'tichr = __main__:main'  # 指定命令行入口：motsi.py 中需有 main() 函数
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',  
        'Operating System :: OS Independent',
    ],
    include_package_data=True,
)

