from setuptools import find_packages, setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(name='dcg-libgen',
      version='0.0.1',
      py_modules=['dcg_libgen'],
      description=long_description,
      url='GitHub: https://github.com/charlesxu90/dcg-libgen',
      author='Charles Xu',
      author_email='Charlexu90@gmail.com',
      license='MIT',
      packages=find_packages(),
      install_requires=[
            'logomaker',
            'numpy',
            'pandas',
            'loguru',
            'seaborn',
            'matplotlib',
            ],
      )