from setuptools import setup

setup(name='Elementally',
      version='0.0.1',
      packages=['elementally'],
      entry_points={
          'console_scripts': [
              'elementally = elementally.__main__:main'
          ]
      },
      )
