import setuptools
with open('README.md', 'r') as fh:
    long_description = fh.read()
    
setuptools.setup(
    name = 'syntiny',
    version = '0.02',
    author = 'Vedika Rajavat', 
    author_email = 'vedika.rajavat.21@ucl.ac.uk', 
    long_description = long_description,
    
    packages = setuptools.find_packages(),
    
    install_requires = ['pandas', 'numpy', 'scipy',
                        'pingouin', 'random', 'matplotlib',
                        'seaborn'],
    
    license = "MIT",
    
    python_requires='>3.5.2',
    
    classifiers=[
        "Programming Language :: Python",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
)