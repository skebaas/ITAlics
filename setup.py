from setuptools import setup, find_packages

setup(
    name='my_package',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'nipype',
        'nibabel',
        'numpy'
    ],
    # Ubuntu specific installations
    extras_require={
        'Ubuntu': ['fsl','freesurfer','matlab','ANTs']
    }
)
