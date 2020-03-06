from setuptools import setup,find_packages

setup(
    name='Scissor',
    version='1.0',
    packages=find_packages(),
    url='https://github.com/jiadong324/Scissor',
    license='',
    author='Jiadong Lin, Songbo Wang',
    author_email='jiadong324@gmail.com',
    description='A complex genome rearrangement simulator',
    requires=['python (>= 3.6)'],
    entry_points={'console_scripts': ['Scissor=Scissor:main']}
)
