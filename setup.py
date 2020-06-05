from setuptools import setup
import time

setup(
    name='pulsegen',
    version=time.strftime('%Y%m%d'),
    author='Rohit Navarathna, Markus Jerger',
    author_email='r.navarathna@uq.edu.au',
    description='A python toolbox for generating pulses for AWGs.',
    #long_description=''
    #license='',
    keywords='pulses, pulse generation, AWG',
    url='http://sqd.equs.org/',
    packages=['pulsegen'], 
    zip_safe=False,
    python_requires='''
        >=3.0
    ''', 
    install_requires='''
        uqtools
    '''
)