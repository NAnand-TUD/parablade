from setuptools import setup

setup(
   name='parablade',
   version='1.0',
   description='Open-source Python library for the parametrization of turbomachinery blades design using gradient-based optimization algorithms',
   author='',
   author_email='',
   packages=['parablade','parablade.common'],
   install_requires=['numpy','matplotlib','scipy','slackclient','pytecplot'],
)
