import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(name='extractome',
                 packages=['extractome'],
                 version='0.0.1',
                 description='',
                 long_description=long_description,
                 long_description_content_type="text/markdown",
                 license='MIT',
                 author='Jim Robinson',
                 url='https://github.com/igvteam/extractome',
                 keywords=['igv', 'bioinformatics', 'genomics', 'visualization', 'variant' ],
                 classifiers=[
                     'Programming Language :: Python :: 3',
                     'Development Status :: 4 - Beta ',
                     'Intended Audience :: Science/Research',
                     'Intended Audience :: Developers',
                     'License :: OSI Approved :: MIT License',
                     'Topic :: Scientific/Engineering :: Bio-Informatics '
                 ],
                 install_requires=[
                     'pysam', 'intervalTree', 'requests'
                 ],
                 entry_points={
                     'console_scripts': [
                         'extractome=extractome.extract:main',
                     ],
                 }
                 )
